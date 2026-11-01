"""``rtt doctor``: checks that every external piece the pipeline needs is installed.

Each check returns (name, ok, detail). ``required`` checks must pass for ``rtt campaign start``;
optional ones only switch lanes off (no GPU -> no MatterGen; no Alex-MP reference -> no lane D
and no DFT-reference hull).
"""
from __future__ import annotations

import importlib
import shutil
import subprocess
from pathlib import Path

from .home import machine_path, rtt_home


def _qe(machine) -> tuple[bool, str]:
    from ..dft.profiles import qe_bin
    b = qe_bin()
    missing = [x for x in ("pw.x", "ph.x", "q2r.x", "matdyn.x") if not (b / x).exists()]
    if missing:
        return False, f"{b}: missing {', '.join(missing)}"
    import tempfile
    try:
        with tempfile.TemporaryDirectory() as td:      # pw.x writes CRASH / input_tmp.in into its cwd
            out = subprocess.run([str(b / "pw.x"), "-h"], capture_output=True, text=True, timeout=30,
                                 input="", cwd=td).stdout + ""
    except (OSError, subprocess.SubprocessError):
        out = ""
    ver = next((ln.strip().split(" starts")[0] for ln in out.splitlines() if "v." in ln and "PWSCF" in ln), "")
    mpi = (b / "mpirun").exists()
    return True, f"{b} {ver}".strip() + ("" if mpi else " (no mpirun next to pw.x: single-rank only)")


def _pseudos(machine) -> tuple[bool, str]:
    from ..dft.profiles import pseudo_dir
    d = pseudo_dir()
    n = len(list(d.glob("*.upf"))) + len(list(d.glob("*.UPF")))
    return n >= 60, f"{d}: {n} UPF files"


def _import(mod: str) -> tuple[bool, str]:
    try:
        m = importlib.import_module(mod)
        return True, getattr(m, "__version__", "ok")
    except Exception as ex:  # noqa: BLE001
        return False, repr(ex)[:120]


def _mattergen(machine) -> tuple[bool, str]:
    from ..source.external import path_of
    v = Path(path_of("mattergen_venv", machine))
    exe = v / "bin" / "mattergen-generate"
    return exe.exists(), str(exe)


#: shipped checkpoints the lanes G1-G5 use (src/rtt/source/lanes.py); install.sh fetches these
LANE_CHECKPOINTS = ("mattergen_base", "chemical_system_energy_above_hull", "chemical_system", "space_group",
                    "dft_band_gap", "dft_mag_density", "ml_bulk_modulus")


def _checkpoints(machine) -> tuple[bool, str]:
    from .front import have_checkpoint
    have = [n for n in LANE_CHECKPOINTS if have_checkpoint(n, machine)]
    miss = sorted(set(LANE_CHECKPOINTS) - set(have))
    return not miss, f"{len(have)}/{len(LANE_CHECKPOINTS)} installed" + (f"; missing {miss}" if miss else "")


def _reference(machine) -> tuple[bool, str]:
    from ..source.external import path_of
    p = Path(path_of("alexmp_lmdb", machine))
    return p.exists(), str(p)


def _gpu(_m) -> tuple[bool, str]:
    if not shutil.which("nvidia-smi"):
        return False, "nvidia-smi not found"
    try:
        import torch
        return bool(torch.cuda.is_available()), f"torch {torch.__version__}, cuda={torch.cuda.is_available()}"
    except Exception as ex:  # noqa: BLE001
        return False, repr(ex)[:120]


def run_checks() -> list[tuple[str, bool, bool, str]]:
    """[(name, ok, required, detail)]"""
    from ..config import load_machine
    mp = machine_path()
    machine = load_machine(mp) if mp.exists() else None
    out = [("machine.toml", mp.exists(), True, str(mp) if mp.exists() else f"{mp} missing: run `rtt init`")]
    checks = [
        ("Quantum ESPRESSO", True, lambda: _qe(machine)),
        ("SG15 pseudopotentials", True, lambda: _pseudos(machine)),
        ("python-flint", True, lambda: _import("flint")),
        ("pymatgen", True, lambda: _import("pymatgen.core")),
        ("textual (dashboard)", True, lambda: _import("textual")),
        ("MatterSim", False, lambda: _import("mattersim")),
        ("GPU (CUDA)", False, lambda: _gpu(machine)),
        ("MatterGen venv", False, lambda: _mattergen(machine)),
        ("MatterGen checkpoints", False, lambda: _checkpoints(machine)),
        ("Alex-MP reference (lane D, DFT hull)", False, lambda: _reference(machine)),
    ]
    for name, req, fn in checks:
        try:
            ok, detail = fn()
        except Exception as ex:  # noqa: BLE001
            ok, detail = False, repr(ex)[:120]
        out.append((name, ok, req, detail))
    return out


def report() -> int:
    rows = run_checks()
    print(f"RTT_HOME = {rtt_home()}")
    for name, ok, req, detail in rows:
        tag = "ok  " if ok else ("FAIL" if req else "off ")
        print(f"[{tag}] {name:38s} {detail}")
    bad = [r for r in rows if r[2] and not r[1]]
    print("\nready" if not bad else f"\nnot ready: {len(bad)} required check(s) failed")
    return 1 if bad else 0
