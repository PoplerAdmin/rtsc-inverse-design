#!/usr/bin/env python3
"""QE 7.6 positive-control driver for Fm-3m Mg2IrH6 (K2PtCl6 type, 9-atom fcc primitive cell).

Reference: Dolui et al., PRL 132, 166001 (2024), arXiv:2310.07562 (PBE a = 6.655 A, harmonic
DFPT stable, Tc ~ 160-175 K).

Structure: ibrav = 2 (QE fcc vectors), Ir 4a (0,0,0), Mg 8c +-(1/4,1/4,1/4), H 24e (x,0,0) and
permutations; positions in alat units. Starting x_H = 0.260649 and a = 6.633 A come from the
Alex-MP DFT reference entry msft-2023062809350933252138247157-GGA (src/ml/check_mg2irh6.py input).
Settings follow the rtt.core DFPT workflow (rtt/core/dfpt/qe_inputs.py): non-spin-polarised,
Marzari-Vanderbilt smearing 0.02 Ry, conv_thr 1e-10, mixing_beta 0.5, SG15 ONCV PBE-1.2 NC pseudos.

Subcommands (all run in the given workdir; nice 19, OMP_NUM_THREADS=1, --ranks MPI ranks):
  ecut   --ecuts 40 50 ... --k 6        scf at each ecutwfc (ecutrho = 4 ecutwfc) -> ecut.json
  kconv  --ecut E --ks 6 8 10 12        scf at each k grid -> kconv.json
  relax  --ecut E --k K                 vc-relax (cell_dofree='ibrav', press 0) -> relax.json
  phon   --ecut E --k K --a A --x X --nq N   scf + ph.x (phonons only, ldisp nq^3) -> phon.json
  pwtemplate --ecut E --k K --a A --x X --out F   write an scf input usable by
                                       src/scoring/coarse_dfpt.py --pw-template
"""
from __future__ import annotations

import argparse
import json
import os
import re
import subprocess
import time
from pathlib import Path

QE_BIN = Path(os.environ.get("QE_BIN", "~/micromamba/envs/qe/bin")).expanduser()
PSEUDO_DIR = Path(os.environ.get("PSEUDO_DIR", "~/pseudo/sg15")).expanduser()
BOHR_A = 0.529177210903
MASS = {"Mg": 24.305, "Ir": 192.217, "H": 1.00794}
A0, X0 = 6.6326, 0.260649


def positions(x: float) -> str:
    rows = [("Ir", 0, 0, 0), ("Mg", .25, .25, .25), ("Mg", -.25, -.25, -.25),
            ("H", x, 0, 0), ("H", -x, 0, 0), ("H", 0, x, 0), ("H", 0, -x, 0),
            ("H", 0, 0, x), ("H", 0, 0, -x)]
    return "ATOMIC_POSITIONS alat\n" + "\n".join(
        f" {s:<3s} {a:16.10f} {b:16.10f} {c:16.10f}" for s, a, b, c in rows) + "\n"


def pw_input(calc: str, a_A: float, x: float, ecut: float, k: int, degauss: float = 0.02,
             la2f: bool = False, extra_ctrl: str = "") -> str:
    species = "ATOMIC_SPECIES\n" + "\n".join(
        f" {e:<3s} {m:.5f} {e}_ONCV_PBE-1.2.upf" for e, m in MASS.items()) + "\n"
    relax = ""
    ctrl = ""
    if calc == "vc-relax":
        ctrl = "    etot_conv_thr = 1.0d-6\n    forc_conv_thr = 1.0d-4\n    nstep = 100\n"
        relax = (" &ions\n    ion_dynamics = 'bfgs'\n /\n &cell\n    cell_dynamics = 'bfgs'\n"
                 "    press = 0.0\n    press_conv_thr = 0.2\n    cell_dofree = 'ibrav'\n /\n")
    return f""" &control
    calculation = '{calc}'
    prefix = 'pwscf'
    outdir = './tmp/'
    pseudo_dir = '{PSEUDO_DIR}/'
    tprnfor = .true.
    tstress = .true.
{ctrl}{extra_ctrl} /
 &system
    ibrav = 2
    celldm(1) = {a_A / BOHR_A:.10f}
    nat = 9
    ntyp = 3
    ecutwfc = {ecut}
    ecutrho = {4 * ecut}
    occupations = 'smearing'
    smearing = 'mv'
    degauss = {degauss}
{"    la2F = .true." + chr(10) if la2f else ""} /
 &electrons
    conv_thr = 1.0d-10
    mixing_beta = 0.5
 /
{relax}{species}{positions(x)}K_POINTS automatic
 {k} {k} {k} 0 0 0
"""


def ph_input(nq: int, tr2: str = "1.0d-16") -> str:
    am = "\n".join(f"    amass({i + 1}) = {m}" for i, m in enumerate(MASS.values()))
    return f"""phonons Mg2IrH6
 &inputph
    prefix = 'pwscf'
    outdir = './tmp/'
    tr2_ph = {tr2}
    fildyn = 'dyn'
    ldisp = .true.
    nq1 = {nq}
    nq2 = {nq}
    nq3 = {nq}
{am}
 /
"""


def run(wd: Path, name: str, code: str, ranks: int, text: str) -> dict:
    (wd / f"{name}.in").write_text(text)
    cmd = ["nice", "-n", "19", str(QE_BIN / "mpirun"), "-np", str(ranks), "--bind-to", "none",
           str(QE_BIN / code), "-nk", str(ranks), "-in", f"{name}.in"]
    env = dict(os.environ, OMP_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1", MKL_NUM_THREADS="1")
    t0 = time.time()
    with (wd / f"{name}.out").open("w") as fo, (wd / f"{name}.err").open("w") as fe:
        rc = subprocess.run(cmd, cwd=wd, stdout=fo, stderr=fe, env=env).returncode
    return {"name": name, "rc": rc, "wall_s": round(time.time() - t0, 2)}


def parse_pw(out: str) -> dict:
    r: dict = {}
    e = re.findall(r"^!\s+total energy\s+=\s+(\S+) Ry", out, re.M)
    if e:
        r["etot_Ry"] = float(e[-1])
    p = re.findall(r"P=\s*(\S+)", out)
    if p:
        r["P_kbar"] = float(p[-1])
    f = re.findall(r"Total force =\s+(\S+)", out)
    if f:
        r["total_force_Ry_bohr"] = float(f[-1])
    fh = re.findall(r"atom\s+4 type\s+3\s+force =\s+(\S+)\s+(\S+)\s+(\S+)", out)
    if fh:
        r["F_H_x_Ry_bohr"] = float(fh[-1][0])
    ef = re.findall(r"the Fermi energy is\s+(\S+) ev", out)
    if ef:
        r["EF_eV"] = float(ef[-1])
    nk = re.search(r"number of k points=\s+(\d+)", out)
    if nk:
        r["nks"] = int(nk.group(1))
    wt = re.findall(r"PWSCF\s+:.*CPU\s+(.*)WALL", out)
    if wt:
        r["pw_wall"] = wt[-1].strip()
    r["job_done"] = "JOB DONE" in out
    return r


def parse_relax(out: str) -> dict:
    r = parse_pw(out)
    a = re.findall(r"CELL_PARAMETERS \(alat=\s*(\S+)\)\s*\n\s*(\S+)\s+(\S+)\s+(\S+)", out)
    if a:
        alat, v1, v2, v3 = (float(t) for t in a[-1])
        # ibrav=2 first vector = alat_new/2 * (-1,0,1) expressed in units of the ORIGINAL alat
        r["a_A"] = 2 * abs(v1) * alat * BOHR_A
    pos = re.findall(r"^H\s+(\S+)\s+0\.0+\s+0\.0+\s*$", out, re.M)
    if pos:
        # final positions block is in alat (original alat) units -> rescale to new a
        r["x_H_alat_orig"] = float(pos[-1])
    r["final_converged"] = "End final coordinates" in out
    return r


_FREQ = re.compile(r"freq\s*\(\s*(\d+)\)\s*=\s*(\S+)\s*\[THz\]\s*=\s*(\S+)\s*\[cm-1\]")
_Q = re.compile(r"q = \(\s*(\S+)\s+(\S+)\s+(\S+)\s*\)")


def parse_dyn(wd: Path) -> list[dict]:
    res = []
    for f in sorted(wd.glob("dyn[0-9]*"), key=lambda p: int(p.name[3:])):
        if f.name == "dyn0":
            continue
        t = f.read_text(errors="replace")
        blocks = t.split("Diagonalizing the dynamical matrix")[1:]
        for b in blocks:
            q = _Q.search(b)
            fr = [float(m.group(3)) for m in _FREQ.finditer(b)]
            res.append({"file": f.name, "q_2pi_over_a": [float(v) for v in q.groups()] if q else None,
                        "freqs_cm1": fr})
    return res


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("cmd", choices=["ecut", "kconv", "relax", "phon", "pwtemplate"])
    ap.add_argument("--wd", type=Path, required=True)
    ap.add_argument("--ranks", type=int, default=4)
    ap.add_argument("--ecuts", type=float, nargs="*", default=[])
    ap.add_argument("--ks", type=int, nargs="*", default=[])
    ap.add_argument("--ecut", type=float)
    ap.add_argument("--k", type=int, default=6)
    ap.add_argument("--a", type=float, default=A0)
    ap.add_argument("--x", type=float, default=X0)
    ap.add_argument("--nq", type=int, default=2)
    ap.add_argument("--degauss", type=float, default=0.02)
    ap.add_argument("--out", type=Path)
    a = ap.parse_args()
    wd = a.wd.expanduser()
    wd.mkdir(parents=True, exist_ok=True)
    rec: dict = {"cmd": a.cmd, "a_A_in": a.a, "x_in": a.x, "ranks": a.ranks, "rows": []}
    if a.cmd == "ecut":
        for e in a.ecuts:
            r = run(wd, f"scf_e{int(e)}", "pw.x", a.ranks, pw_input("scf", a.a, a.x, e, a.k, a.degauss))
            r.update(parse_pw((wd / f"scf_e{int(e)}.out").read_text()), ecutwfc=e, k=a.k)
            rec["rows"].append(r)
            (wd / "ecut.json").write_text(json.dumps(rec, indent=1))
    elif a.cmd == "kconv":
        for k in a.ks:
            r = run(wd, f"scf_k{k}", "pw.x", a.ranks, pw_input("scf", a.a, a.x, a.ecut, k, a.degauss))
            r.update(parse_pw((wd / f"scf_k{k}.out").read_text()), ecutwfc=a.ecut, k=k)
            rec["rows"].append(r)
            (wd / "kconv.json").write_text(json.dumps(rec, indent=1))
    elif a.cmd == "relax":
        r = run(wd, "relax", "pw.x", a.ranks, pw_input("vc-relax", a.a, a.x, a.ecut, a.k, a.degauss))
        r.update(parse_relax((wd / "relax.out").read_text()), ecutwfc=a.ecut, k=a.k)
        if "a_A" in r and "x_H_alat_orig" in r:
            r["x_H"] = r["x_H_alat_orig"] * a.a / r["a_A"]
        rec["rows"].append(r)
        (wd / "relax.json").write_text(json.dumps(rec, indent=1))
    elif a.cmd == "phon":
        r1 = run(wd, "scf", "pw.x", a.ranks, pw_input("scf", a.a, a.x, a.ecut, a.k, a.degauss))
        r1.update(parse_pw((wd / "scf.out").read_text()))
        rec["rows"].append(r1)
        (wd / "phon.json").write_text(json.dumps(rec, indent=1))
        r2 = run(wd, "ph", "ph.x", a.ranks, ph_input(a.nq))
        rec["rows"].append(r2)
        rec["dyn"] = parse_dyn(wd)
        (wd / "phon.json").write_text(json.dumps(rec, indent=1))
    elif a.cmd == "pwtemplate":
        a.out.write_text(pw_input("scf", a.a, a.x, a.ecut, a.k, a.degauss))
    print(json.dumps(rec, indent=1)[:4000])


if __name__ == "__main__":
    main()
