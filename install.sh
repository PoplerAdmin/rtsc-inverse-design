#!/usr/bin/env bash
# One-time setup of everything rtsc-targeted needs, under $RTT_HOME (default ~/.rtt).
# Linux x86-64 or WSL2 (Ubuntu). Re-running is safe: finished steps are skipped.
#
#   ./install.sh                     # full install (QE, Python env, MatterGen, checkpoints, SG15)
#   ./install.sh --with-reference    # also the Alex-MP reference set (lane D + DFT hull; ~5 GB)
#   ./install.sh --no-mattergen      # CPU-only machine: skip MatterGen and its checkpoints
#
# Nothing is installed system-wide and nothing needs root. Downloads: micromamba, conda-forge
# packages, PyPI packages, the MatterGen source (GitHub) and checkpoints (Hugging Face), the SG15
# pseudopotentials (quantum-simulation.org). Every download's sha256 is written to
# $RTT_HOME/install_manifest.txt.
#
# Flow (diagram and reasons: docs/RUNNER.md §5):
#   preflight -> micromamba -> QE env | Python env -> SG15 -> MatterGen env + checkpoints
#   -> MatterSim weights -> [Alex-MP reference] -> env.sh -> rtt init -> rtt doctor
# Design choices in short:
#   * one private prefix ($RTT_HOME): no root, nothing system-wide, `rm -rf ~/.rtt` uninstalls;
#   * versions pinned in tools/versions.env: every measured cost and calibration in docs/design
#     was made with exactly these versions, so a different version means re-measuring;
#   * MatterGen gets its own Python: its pins (Python 3.10, torch 2.2) clash with MatterSim;
#   * all downloads happen here, never during a campaign, so a run is reproducible from the
#     manifest.
set -euo pipefail

REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# shellcheck source=tools/versions.env
source "$REPO/tools/versions.env"
RTT_HOME="${RTT_HOME:-$HOME/.rtt}"
WITH_REF=0
WITH_MG=1
for a in "$@"; do
  case "$a" in
    --with-reference) WITH_REF=1 ;;
    --no-mattergen) WITH_MG=0 ;;
    -h|--help) sed -n '2,/^set -euo/p' "$0" | sed '$d'; exit 0 ;;
    *) echo "unknown option $a"; exit 2 ;;
  esac
done

say() { printf '\n\033[1m== %s\033[0m\n' "$*"; }
die() { printf '\033[31merror:\033[0m %s\n' "$*" >&2; exit 1; }
MANIFEST="$RTT_HOME/install_manifest.txt"
record() { printf '%s\t%s\t%s\n' "$(date -u +%FT%TZ)" "$1" "$2" >> "$MANIFEST"; }

# ------------------------------------------------------------------ preflight
say "preflight"
[ "$(uname -s)" = "Linux" ] || die "Linux or WSL2 required (on Windows: install WSL2 Ubuntu and run this inside it)"
[ "$(uname -m)" = "x86_64" ] || die "x86-64 required"
for t in curl git tar gzip; do command -v "$t" >/dev/null || die "missing tool: $t (e.g. sudo apt install $t)"; done
mkdir -p "$RTT_HOME"/{bin,envs,pseudo,runs,scratch,hf}
free_gb=$(df -Pk "$RTT_HOME" | awk 'NR==2{print int($4/1024/1024)}')
echo "RTT_HOME=$RTT_HOME  free disk: ${free_gb} GB"
[ "$free_gb" -ge 40 ] || die "need at least 40 GB free under $RTT_HOME (DFPT scratch needs much more; 200 GB recommended)"
case "$REPO" in /mnt/[a-z]/*) echo "note: the repository is on a Windows drive; runs and scratch go to $RTT_HOME (Linux file system) for speed and SQLite safety" ;; esac
touch "$MANIFEST"

# ------------------------------------------------------------------ micromamba
# a single static binary: gives conda-forge packages (QE with MPI) without a system conda,
# and needs no bzip2 (absent on stock Ubuntu, unlike the .tar.bz2 distribution)
MM="$RTT_HOME/bin/micromamba"
export MAMBA_ROOT_PREFIX="$RTT_HOME/mamba"
if [ ! -x "$MM" ]; then
  say "micromamba"
  curl -fsSL "https://github.com/mamba-org/micromamba-releases/releases/download/$MICROMAMBA_VERSION/micromamba-linux-64" \
    -o "$MM.tmp"
  record micromamba "$MICROMAMBA_VERSION $(sha256sum "$MM.tmp" | cut -d' ' -f1)"
  chmod +x "$MM.tmp" && mv "$MM.tmp" "$MM"
fi

# ------------------------------------------------------------------ Quantum ESPRESSO
# the reference machine ran conda-forge qe 7.6 built against OpenMPI; the launcher in
# src/rtt/dft/workflow.py uses OpenMPI's `mpirun --bind-to none`, so the MPI flavour matters
QE="$RTT_HOME/envs/qe"
if [ ! -x "$QE/bin/ph.x" ]; then
  say "Quantum ESPRESSO $QE_VERSION (conda-forge, OpenMPI)"
  "$MM" create -y -p "$QE" -c conda-forge "qe=$QE_VERSION=*openmpi*" openmpi \
    || die "could not install qe=$QE_VERSION from conda-forge; edit QE_VERSION in tools/versions.env"
fi
record qe "$("$MM" list -p "$QE" '^qe$' --json 2>/dev/null | tr -d '\n ' | grep -o '"version":"[^"]*"' | head -1)"

# ------------------------------------------------------------------ rtt Python environment
PY="$RTT_HOME/envs/py"
if [ ! -x "$PY/bin/python" ]; then
  say "Python $PYTHON_VERSION environment"
  "$MM" create -y -p "$PY" -c conda-forge "python=$PYTHON_VERSION" pip
fi
say "rtt and Python dependencies"
"$PY/bin/python" -m pip install --upgrade pip >/dev/null
"$PY/bin/python" -m pip install -e "$REPO[ml,tui,test]" "mattersim==$MATTERSIM_VERSION"
record mattersim "$MATTERSIM_VERSION"
record torch "$("$PY/bin/python" -c 'import torch;print(torch.__version__)')"

# ------------------------------------------------------------------ SG15 pseudopotentials
PSEUDO="$RTT_HOME/pseudo/sg15"
if [ ! -f "$PSEUDO/H_ONCV_PBE-1.2.upf" ]; then
  say "SG15 ONCV pseudopotentials"
  mkdir -p "$PSEUDO"
  tgz="$RTT_HOME/pseudo/$(basename "$SG15_URL")"
  curl -fsSL "$SG15_URL" -o "$tgz"
  sum=$(sha256sum "$tgz" | cut -d' ' -f1)
  if [ -n "$SG15_SHA256" ] && [ "$sum" != "$SG15_SHA256" ]; then die "SG15 tarball sha256 $sum != pinned $SG15_SHA256"; fi
  [ -n "$SG15_SHA256" ] || echo "SG15 sha256 (not pinned yet): $sum"
  record sg15 "$sum"
  tar -xzf "$tgz" -C "$PSEUDO" --strip-components=0
  # flatten if the tarball has a top directory
  find "$PSEUDO" -mindepth 2 -name '*.upf' -exec mv -n {} "$PSEUDO"/ \;
fi
[ -f "$PSEUDO/H_ONCV_PBE-1.2.upf" ] || die "SG15 files not found in $PSEUDO (expected <El>_ONCV_PBE-1.2.upf)"

# ------------------------------------------------------------------ MatterGen
MGDIR="$RTT_HOME/mattergen"
MGENV="$RTT_HOME/envs/mattergen"
export HF_HOME="$RTT_HOME/hf"
if [ "$WITH_MG" = 1 ]; then
  if [ ! -d "$MGDIR/.git" ]; then
    say "MatterGen source ($MATTERGEN_COMMIT)"
    GIT_LFS_SKIP_SMUDGE=1 git clone -q "$MATTERGEN_REPO" "$MGDIR"
  fi
  git -C "$MGDIR" -c advice.detachedHead=false checkout -q "$MATTERGEN_COMMIT"
  record mattergen_commit "$(git -C "$MGDIR" rev-parse HEAD)"
  if [ ! -x "$MGENV/bin/mattergen-generate" ]; then
    # MatterGen's pyproject carries uv settings (torch_cluster wheel index, build isolation) that
    # plain pip ignores, so the environment is built with uv, as MatterGen's README does
    say "MatterGen environment (Python $MATTERGEN_PYTHON, uv $UV_VERSION)"
    "$PY/bin/python" -m pip install -q "uv==$UV_VERSION"
    rm -rf "$MGENV"
    UV_PYTHON_INSTALL_DIR="$RTT_HOME/uv-python" "$PY/bin/uv" venv -q --python "$MATTERGEN_PYTHON" "$MGENV"
    (cd "$MGDIR" && GIT_LFS_SKIP_SMUDGE=1 VIRTUAL_ENV="$MGENV" UV_PYTHON_INSTALL_DIR="$RTT_HOME/uv-python" \
      "$PY/bin/uv" pip install -e .)
  fi
  say "MatterGen checkpoints ($MATTERGEN_HF_REPO)"
  "$MGENV/bin/python" - "$MATTERGEN_HF_REPO" $MATTERGEN_CHECKPOINTS <<'EOF'
import hashlib, sys
from pathlib import Path
from huggingface_hub import snapshot_download
repo, names = sys.argv[1], sys.argv[2:]
for n in names:
    root = snapshot_download(repo, allow_patterns=[f"checkpoints/{n}/*", f"checkpoints/{n}/**"])
    ck = Path(root) / "checkpoints" / n / "checkpoints" / "last.ckpt"
    h = hashlib.sha256(ck.read_bytes()).hexdigest() if ck.exists() else "MISSING"
    print(f"checkpoint {n} {h}", flush=True)
EOF
  "$MGENV/bin/python" -c "import huggingface_hub" >/dev/null
fi

# ------------------------------------------------------------------ MatterSim weights (first-use download, done now)
say "MatterSim model weights"
TORCH_DISABLE_NATIVE_JIT=1 "$PY/bin/python" - "$MATTERSIM_MODEL" <<'EOF' || echo "warning: MatterSim weights not fetched (retried on first use)"
import sys
from mattersim.forcefield.potential import Potential
Potential.from_checkpoint(load_path=sys.argv[1], device="cpu")
print("ok")
EOF

# ------------------------------------------------------------------ Alex-MP reference (optional)
REF="$RTT_HOME/reference/reference_MP2020correction.lmdb"
if [ "$WITH_REF" = 1 ] && [ ! -s "$REF" ]; then
  [ "$WITH_MG" = 1 ] || die "--with-reference needs MatterGen (its reference preset)"
  say "Alex-MP reference set (0.9 GB download, 3.9 GB unpacked)"
  mkdir -p "$(dirname "$REF")"
  gz=$("$MGENV/bin/python" -c "from mattergen.evaluation.reference.presets import get_reference_dataset_path as g; print(g('reference_MP2020correction.gz'))")
  record alexmp_reference "$(sha256sum "$gz" | cut -d' ' -f1)"
  zcat "$gz" > "$REF.tmp" && mv "$REF.tmp" "$REF"
fi

# ------------------------------------------------------------------ environment file
say "writing $RTT_HOME/env.sh"
cat > "$RTT_HOME/env.sh" <<EOF
# rtsc-targeted environment (written by install.sh; source it before running rtt)
export RTT_HOME="$RTT_HOME"
export RTT_REPO="$REPO"
export PATH="$PY/bin:\$PATH"
export RTT_QE_BIN="$QE/bin"
export RTT_PSEUDO_DIR="$PSEUDO"
export RTT_MATTERGEN_VENV="$MGENV"
export HF_HOME="$RTT_HOME/hf"
export RTT_HF_HUB="$RTT_HOME/hf/hub"
export RTT_ALEXMP_LMDB="$REF"
export RTT_FRONT_CACHE="$RTT_HOME/scratch/front_cache"
export RTT_PRUNE_MANIFEST_DIR="$RTT_HOME/prune_manifests"
export TORCH_DISABLE_NATIVE_JIT=1
$( [ "$(id -u)" = 0 ] && echo 'export OMPI_ALLOW_RUN_AS_ROOT=1 OMPI_ALLOW_RUN_AS_ROOT_CONFIRM=1' )
EOF
# shellcheck disable=SC1091
source "$RTT_HOME/env.sh"

say "hardware and machine.toml"
rtt init
say "installation check"
rtt doctor || true
cat <<EOF

Installed under $RTT_HOME. Start the dashboard with:

    ./start.sh

EOF
