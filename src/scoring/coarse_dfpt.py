#!/usr/bin/env python3
"""Tier-2 coarse DFPT chain (generalisation of coarse_dfpt_probe.sh, which was hard-coded for H3S).

Runs the same QE 7.6 chain as the rtt.core DFPT workflow
    scf_fit (la2F=.true., dense k_fit) -> scf (k_scf) -> ph.x (electron_phonon='interpolated',
    ldisp, q nq^3) -> q2r.x (la2F) -> matdyn.x (la2F, dos)
for ONE structure and ONE set of coarse settings, and writes into <workdir>:
    *.in / *.out / *.err for every step, a2F.dos1..N, lambda, dyn*, ifc.fc, freq.dos
    coarse_run.json : settings, per-step wall time and return code, SHA-256 of every output file,
                      lowest phonon frequency per irreducible q (dyn files) and on the matdyn mesh.

Structure input (one of):
  --pw-template FILE   a pw.x input (e.g. a converged run's scf_fit.in); ATOMIC_SPECIES,
                       CELL_PARAMETERS, ATOMIC_POSITIONS and ecutwfc/ecutrho are taken from it.
  --structure FILE     any ASE-readable structure (CIF, POSCAR, extxyz); pseudopotentials are
                       '<El>_ONCV_PBE-1.2.upf' (SG15) unless --pseudo-pattern is given; --ecutwfc
                       is then required.
Everything numerical is informational (tier 2): coarse spectra never drive a proven status.

Compute policy (project CONTEXT.md): nice 19, OMP_NUM_THREADS=1, --ranks MPI ranks (default 2).
Long runs: setsid nohup python coarse_dfpt.py ... > log 2>&1 < /dev/null &

Example (Al, reference cutoffs, k 12/24, q 3^3, 6 ranks):
  python coarse_dfpt.py --pw-template ref/scf_fit.in --workdir ~/w/al_q3 \
      --kscf 12 --kfit 24 --nq 3 --sigma0 0.005 --nsigma 10 --ranks 6
"""
from __future__ import annotations

import argparse
import hashlib
import json
import os
import re
import subprocess
import sys
import time
from pathlib import Path

QE_BIN = Path(os.environ.get("QE_BIN", "~/micromamba/envs/qe/bin")).expanduser()
PSEUDO_DIR = Path(os.environ.get("PSEUDO_DIR", "~/pseudo/sg15")).expanduser()

# --------------------------------------------------------------------------- structure
_BLOCKS = ("ATOMIC_SPECIES", "CELL_PARAMETERS", "ATOMIC_POSITIONS", "K_POINTS")


def _cards(text: str) -> dict[str, list[str]]:
    """Split the card section of a pw.x input into {card_name: [header, lines...]}."""
    out: dict[str, list[str]] = {}
    cur = None
    for line in text.splitlines():
        head = line.strip().split()[0].upper() if line.strip() else ""
        if head in _BLOCKS:
            cur = head
            out[cur] = [line.strip()]
            continue
        if cur and line.strip():
            out[cur].append(line.rstrip())
    return out


def _nml_value(text: str, key: str) -> str | None:
    m = re.search(rf"\b{key}\s*=\s*([^,\n]+)", text, re.I)
    return m.group(1).strip() if m else None


def structure_from_template(path: Path) -> dict:
    text = path.read_text()
    cards = _cards(text)
    species = [ln.split() for ln in cards["ATOMIC_SPECIES"][1:]]
    nat = len(cards["ATOMIC_POSITIONS"]) - 1
    return {
        "species_block": "\n".join(cards["ATOMIC_SPECIES"]),
        "cell_block": "\n".join(cards["CELL_PARAMETERS"]),
        "pos_block": "\n".join(cards["ATOMIC_POSITIONS"]),
        "masses": [float(s[1]) for s in species],
        "nat": nat, "ntyp": len(species),
        "ecutwfc": float(_nml_value(text, "ecutwfc").replace("d", "e")),
        "ecutrho": float(_nml_value(text, "ecutrho").replace("d", "e")),
        "source": str(path),
    }


def structure_from_file(path: Path, pattern: str) -> dict:
    from ase.data import atomic_masses, atomic_numbers  # noqa: WPS433 (optional dependency)
    from ase.io import read
    at = read(str(path))
    els = list(dict.fromkeys(at.get_chemical_symbols()))
    masses = [float(atomic_masses[atomic_numbers[e]]) for e in els]
    sp = ["ATOMIC_SPECIES"] + [f" {e}  {m:.8f} {pattern.format(el=e)}" for e, m in zip(els, masses)]
    cell = ["CELL_PARAMETERS angstrom"] + ["  " + "  ".join(f"{x:.12f}" for x in v) for v in at.cell]
    pos = ["ATOMIC_POSITIONS crystal"] + [
        f" {s}  " + "  ".join(f"{x:.12f}" for x in f)
        for s, f in zip(at.get_chemical_symbols(), at.get_scaled_positions())]
    return {"species_block": "\n".join(sp), "cell_block": "\n".join(cell),
            "pos_block": "\n".join(pos), "masses": masses, "nat": len(at), "ntyp": len(els),
            "ecutwfc": None, "ecutrho": None, "source": str(path)}


# --------------------------------------------------------------------------- inputs
def pw_input(st: dict, ecut: float, ecutrho: float, k: int, la2f: bool, degauss: float) -> str:
    extra = "    la2F = .true.\n" if la2f else ""
    return f""" &control
    calculation = 'scf'
    prefix = 'pwscf'
    outdir = './tmp/'
    pseudo_dir = '{PSEUDO_DIR}/'
 /
 &system
    ibrav = 0
    nat = {st['nat']}
    ntyp = {st['ntyp']}
    ecutwfc = {ecut}
    ecutrho = {ecutrho}
    occupations = 'smearing'
    smearing = 'mv'
    degauss = {degauss}
{extra} /
 &electrons
    conv_thr = 1.0d-10
    mixing_beta = 0.5
 /
{st['species_block']}
{st['cell_block']}
{st['pos_block']}
K_POINTS automatic
 {k} {k} {k} 0 0 0
"""


def ph_input(st: dict, nq: int, sigma0: float, nsigma: int, tr2: str) -> str:
    am = "\n".join(f"    amass({i + 1}) = {m}" for i, m in enumerate(st["masses"]))
    return f"""e-ph coarse tier-2
 &inputph
    prefix = 'pwscf'
    outdir = './tmp/'
    tr2_ph = {tr2}
    fildyn = 'dyn'
    fildvscf = 'dvscf'
    electron_phonon = 'interpolated'
    el_ph_sigma = {sigma0}
    el_ph_nsigma = {nsigma}
    trans = .true.
    ldisp = .true.
    nq1 = {nq}
    nq2 = {nq}
    nq3 = {nq}
{am}
 /
"""


def q2r_input(nsigma: int) -> str:
    return (" &input\n    fildyn = 'dyn'\n    flfrc = 'ifc.fc'\n    zasr = 'crystal'\n"
            f"    la2F = .true.\n    el_ph_nsigma = {nsigma}\n /\n")


def matdyn_input(st: dict, nk: int, nsigma: int) -> str:
    am = "\n".join(f"    amass({i + 1}) = {m}" for i, m in enumerate(st["masses"]))
    return f""" &input
    asr = 'crystal'
    flfrc = 'ifc.fc'
    flfrq = 'freq.dos'
    dos = .true.
    la2F = .true.
    fldos = 'phonon.dos'
    nk1 = {nk}
    nk2 = {nk}
    nk3 = {nk}
    el_ph_nsigma = {nsigma}
    ndos = 200
{am}
 /
"""


# --------------------------------------------------------------------------- run + record
def sha256(p: Path) -> str:
    h = hashlib.sha256()
    with p.open("rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def run_step(wd: Path, name: str, cmd: list[str], timeout: float | None, rec: dict) -> bool:
    t0 = time.time()
    with (wd / f"{name}.out").open("w") as fo, (wd / f"{name}.err").open("w") as fe:
        try:
            rc = subprocess.run(cmd, cwd=wd, stdout=fo, stderr=fe, timeout=timeout).returncode
            to = False
        except subprocess.TimeoutExpired:
            rc, to = -9, True
    rec["steps"][name] = {"cmd": " ".join(cmd).replace(str(Path.home()), "~"),
                          "wall_s": round(time.time() - t0, 3), "rc": rc, "timed_out": to}
    (wd / "coarse_run.json").write_text(json.dumps(rec, indent=1))
    return rc == 0 and not to


_FREQ = re.compile(r"freq\s*\(\s*\d+\)\s*=\s*(\S+)\s*\[THz\]\s*=\s*(\S+)\s*\[cm-1\]")


def min_freqs(wd: Path) -> dict:
    """Lowest frequency (cm^-1, negative = imaginary) per dyn file (irreducible q) and on the
    matdyn interpolation mesh (freq.dos)."""
    per_q = {}
    for f in sorted(wd.glob("dyn[0-9]*"), key=lambda p: int(p.name[3:])):
        if f.name == "dyn0":
            continue
        vals = [float(m.group(2)) for m in _FREQ.finditer(f.read_text(errors="replace"))]
        if vals:
            per_q[f.name] = min(vals)
    mesh_min = None
    fd = wd / "freq.dos"
    if fd.exists():
        nums = []
        lines = fd.read_text(errors="replace").splitlines()
        for i, ln in enumerate(lines[1:]):
            parts = ln.split()
            if i % 2 == 1:      # every second line holds the frequencies (cm^-1)
                try:
                    nums.extend(float(x) for x in parts)
                except ValueError:
                    pass
        mesh_min = min(nums) if nums else None
    return {"dyn_min_cm1": per_q, "matdyn_mesh_min_cm1": mesh_min}


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    src = ap.add_mutually_exclusive_group(required=True)
    src.add_argument("--pw-template", type=Path)
    src.add_argument("--structure", type=Path)
    ap.add_argument("--pseudo-pattern", default="{el}_ONCV_PBE-1.2.upf")
    ap.add_argument("--workdir", type=Path, required=True)
    ap.add_argument("--ecutwfc", type=float, help="default: template value")
    ap.add_argument("--ecutrho", type=float, help="default: template value, else 4*ecutwfc")
    ap.add_argument("--kscf", type=int, required=True)
    ap.add_argument("--kfit", type=int, required=True)
    ap.add_argument("--nq", type=int, required=True)
    ap.add_argument("--sigma0", type=float, default=0.01, help="el_ph_sigma (Ry)")
    ap.add_argument("--nsigma", type=int, default=10)
    ap.add_argument("--degauss", type=float, default=0.02)
    ap.add_argument("--tr2", default="1.0d-16")
    ap.add_argument("--matdyn-nk", type=int, default=32)
    ap.add_argument("--ranks", type=int, default=2)
    ap.add_argument("--pools", type=int, help="pw.x/ph.x -nk (default = ranks)")
    ap.add_argument("--ph-timeout", type=float, default=None, help="seconds")
    a = ap.parse_args(argv)

    if a.kfit % a.nq:
        ap.error(f"k_fit {a.kfit} must be a multiple of nq {a.nq} (interpolated e-ph needs k+q on the k_fit grid)")
    if a.kfit % a.kscf:
        ap.error(f"k_fit {a.kfit} must be a multiple of k_scf {a.kscf} (QE silently gives garbage otherwise; doc 15)")
    st = structure_from_template(a.pw_template) if a.pw_template else structure_from_file(a.structure, a.pseudo_pattern)
    ecut = a.ecutwfc or st["ecutwfc"]
    if ecut is None:
        ap.error("--ecutwfc is required with --structure")
    ecutrho = a.ecutrho or (st["ecutrho"] if not a.ecutwfc and st["ecutrho"] else 4 * ecut)
    wd = a.workdir.expanduser()
    wd.mkdir(parents=True, exist_ok=True)
    (wd / "scf_fit.in").write_text(pw_input(st, ecut, ecutrho, a.kfit, True, a.degauss))
    (wd / "scf.in").write_text(pw_input(st, ecut, ecutrho, a.kscf, False, a.degauss))
    (wd / "ph.in").write_text(ph_input(st, a.nq, a.sigma0, a.nsigma, a.tr2))
    (wd / "q2r.in").write_text(q2r_input(a.nsigma))
    (wd / "matdyn.in").write_text(matdyn_input(st, a.matdyn_nk, a.nsigma))

    env_ok = os.environ.get("OMP_NUM_THREADS") == "1"
    rec = {"tool": "src/scoring/coarse_dfpt.py", "structure_source": Path(st["source"]).name,
           "settings": {"ecutwfc": ecut, "ecutrho": ecutrho, "kscf": a.kscf, "kfit": a.kfit,
                        "nq": a.nq, "sigma0_Ry": a.sigma0, "nsigma": a.nsigma,
                        "degauss_Ry": a.degauss, "tr2_ph": a.tr2, "matdyn_nk": a.matdyn_nk,
                        "ranks": a.ranks, "omp_threads_1": env_ok},
           "steps": {}}
    np_ = a.ranks
    nk = str(a.pools or a.ranks)
    mpi = ["nice", "-n", "19", str(QE_BIN / "mpirun"), "-np", str(np_), "--bind-to", "none"]
    ser = ["nice", "-n", "19"]
    ok = (run_step(wd, "scf_fit", mpi + [str(QE_BIN / "pw.x"), "-nk", nk, "-in", "scf_fit.in"], None, rec)
          and run_step(wd, "scf", mpi + [str(QE_BIN / "pw.x"), "-nk", nk, "-in", "scf.in"], None, rec)
          and run_step(wd, "ph", mpi + [str(QE_BIN / "ph.x"), "-nk", nk, "-in", "ph.in"], a.ph_timeout, rec)
          and run_step(wd, "q2r", ser + [str(QE_BIN / "q2r.x"), "-in", "q2r.in"], None, rec)
          and run_step(wd, "matdyn", ser + [str(QE_BIN / "matdyn.x"), "-in", "matdyn.in"], None, rec))
    rec["ok"] = ok
    rec["wall_total_s"] = round(sum(s["wall_s"] for s in rec["steps"].values()), 3)
    rec["phonon_min"] = min_freqs(wd)
    outs = sorted([p for p in wd.iterdir() if p.is_file() and p.name != "coarse_run.json"
                   and (p.name.startswith(("a2F.dos", "dyn")) or p.suffix in (".in", ".out", ".fc")
                        or p.name in ("lambda", "freq.dos", "phonon.dos", "matdyn.modes"))])
    rec["sha256"] = {p.name: sha256(p) for p in outs}
    (wd / "coarse_run.json").write_text(json.dumps(rec, indent=1))
    print(json.dumps({k: rec[k] for k in ("ok", "wall_total_s", "steps")}, indent=1))
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
