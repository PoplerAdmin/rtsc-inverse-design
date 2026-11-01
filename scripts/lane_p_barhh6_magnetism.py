"""BaRhH6 (campaign-1 c-336a1451e874) H2 re-check: AFM starts and SOC (R30 F6/F7 item 5; doc 32 §6).

    nice -n 19 ~/.rtt/envs/py/bin/python scripts/lane_p_barhh6_magnetism.py --work DIR [--run] [--ranks 2]

Cell: the campaign-1 S4 (PBE, SG15, 60 Ry) vc-relaxed cell, one formula unit (8 atoms; R-3m at
0.1 Å). Rh⁴⁺ is formally low-spin d5 (17 valence electrons per cell). The campaign-1 S4 spin
check started only from a uniform FM moment and converged to 0 µB.

Runs (same S4 settings: SG15 ONCV PBE-1.2, ecutwfc 60 Ry / 240 Ry, MV smearing 0.02 Ry, k mesh
from the S4 run scaled to the supercell), all in the 16-atom cell doubled along a3 or a1:

  nm_<ax>     nspin 1                      (reference)
  fm_<ax>     nspin 2, +0.8 on both Rh
  afm_<ax>    nspin 2, +0.8 / -0.8 (Rh1/Rh2) -- AFM along the doubled axis
  soc_nm_a3   noncollinear + lspinorb, no starting moment     (FR pseudos for Rh, H; Ba has
  soc_afm_a3  noncollinear + lspinorb, AFM start along z       no FR SG15 file: scalar-rel.)

Doubling along a3 makes the two shortest Rh-Rh bonds (a1+a3, a2+a3; 4.63 Å) antiparallel;
along a1 makes a1+a3 and a1-a2 antiparallel. The Rh sublattice is (near-)triangular, so no
collinear AFM is unfrustrated; these are the two simplest bipartitions. No DFT+U.
"""
from __future__ import annotations

import argparse
import json
import os
import re
import subprocess
import time
from pathlib import Path

LAT = [[4.853399452, -2.518737091, -0.293396463],
       [4.853399452, 2.518737091, -0.293396463],
       [-3.478211497, 0.0, 3.922806986]]
SPECIES = ["H"] * 6 + ["Rh", "Ba"]
FRAC = [[0.6407030683, 0.0207829607, 0.8265851126], [0.1107009512, 0.1107009512, 1.0151916739],
        [0.7585903783, 0.3825071009, 0.2666584957], [0.0207829607, 0.6407030683, 0.8265851126],
        [0.898205541, 0.898205541, 0.3089863518], [0.3825071009, 0.7585903783, 0.2666584957],
        [1.0031324552, 1.0031324552, 0.158932026], [0.5020795445, 0.5020795445, 0.6635627317]]
K_PRIM = (7, 7, 7)
MASS = {"H": 1.008, "Rh": 102.9055, "Ba": 137.327}
PSEUDO = os.path.expanduser("~/pseudo/sg15")
QE = os.path.expanduser("~/micromamba/envs/qe/bin")


def supercell(axis: int):
    lat = [list(v) for v in LAT]
    lat[axis] = [2 * x for x in lat[axis]]
    sp, fc = [], []
    for n in (0, 1):
        for s, f in zip(SPECIES, FRAC):
            g = list(f)
            g[axis] = (g[axis] + n) / 2.0
            lab = s if s != "Rh" else ("Rh1" if n == 0 else "Rh2")
            sp.append(lab)
            fc.append(g)
    k = list(K_PRIM)
    k[axis] = max(1, (k[axis] + 1) // 2)
    return lat, sp, fc, k


def pw_input(name: str, axis: int, mode: str) -> str:
    lat, sp, fc, k = supercell(axis)
    soc = mode.startswith("soc")
    spin = {"nm": "", "fm": "nspin = 2\n", "afm": "nspin = 2\n"}.get(mode.split("_")[-1], "")
    if soc:
        spin = "noncolin = .true.\n    lspinorb = .true.\n"
    mag = {"fm": (0.8, 0.8), "afm": (0.8, -0.8), "nm": None}[mode.split("_")[-1]]
    types = ["H", "Rh1", "Rh2", "Ba"]
    sys_extra = ""
    if mag is not None:
        sys_extra += f"    starting_magnetization(2) = {mag[0]}\n    starting_magnetization(3) = {mag[1]}\n"
        if soc:
            sys_extra += "    angle1(2) = 0.0\n    angle1(3) = 0.0\n    angle2(2) = 0.0\n    angle2(3) = 0.0\n"
    def pp(t):
        el = re.sub(r"\d", "", t)
        if soc and el in ("H", "Rh"):
            return f"{el}_ONCV_PBE_FR-1.0.upf"
        return f"{el}_ONCV_PBE-1.2.upf"
    s = f"""&control
    calculation = 'scf'
    prefix = '{name}'
    outdir = './tmp_{name}/'
    pseudo_dir = '{PSEUDO}/'
    tprnfor = .true.
    tstress = .true.
/
&system
    ibrav = 0
    nat = {len(sp)}
    ntyp = 4
    ecutwfc = 60.0
    ecutrho = 240.0
    occupations = 'smearing'
    smearing = 'mv'
    degauss = 0.02
    {spin}{sys_extra}/
&electrons
    conv_thr = 1.0d-9
    mixing_beta = 0.3
    electron_maxstep = 200
/
ATOMIC_SPECIES
"""
    for t in types:
        el = re.sub(r"\d", "", t)
        s += f" {t} {MASS[el]} {pp(t)}\n"
    s += "CELL_PARAMETERS angstrom\n" + "".join(f" {v[0]:.9f} {v[1]:.9f} {v[2]:.9f}\n" for v in lat)
    s += "ATOMIC_POSITIONS crystal\n" + "".join(f" {a} {f[0]:.10f} {f[1]:.10f} {f[2]:.10f}\n" for a, f in zip(sp, fc))
    s += f"K_POINTS automatic\n {k[0]} {k[1]} {k[2]} 0 0 0\n"
    return s


RUNS = [("nm_a3", 2, "nm"), ("fm_a3", 2, "fm"), ("afm_a3", 2, "afm"),
        ("afm_a1", 0, "afm"), ("nm_a1", 0, "nm"),
        ("soc_nm_a3", 2, "soc_nm"), ("soc_afm_a3", 2, "soc_afm")]


def parse(out: Path) -> dict:
    t = out.read_text(errors="replace") if out.exists() else ""
    e = re.findall(r"^!\s+total energy\s+=\s+(-?\d+\.\d+) Ry", t, re.M)
    tm = re.findall(r"total magnetization\s+=\s+(.+?) Bohr", t)
    am = re.findall(r"absolute magnetization\s+=\s+(-?\d+\.\d+) Bohr", t)
    site = re.findall(r"atom\s+(\d+) \(R=.*?\)\s+charge=\s*(-?\d+\.\d+)\s+magn=\s*(-?\d+\.\d+)", t)
    if not site:
        site = re.findall(r"atom:\s+(\d+)\s+charge:\s+(-?\d+\.\d+)\s+magn:\s+(-?\d+\.\d+)", t)
    return {"E_Ry": float(e[-1]) if e else None, "total_mag": tm[-1].strip() if tm else None,
            "abs_mag": float(am[-1]) if am else None, "converged": "convergence has been achieved" in t,
            "site_magn_last": site[-16:] if site else None,
            "wall": (re.findall(r"PWSCF\s+:.*?CPU\s+(.*?) WALL", t) or [None])[-1]}


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--work", required=True)
    ap.add_argument("--run", action="store_true")
    ap.add_argument("--ranks", type=int, default=2)
    ap.add_argument("--only", nargs="*")
    a = ap.parse_args()
    wd = Path(a.work).expanduser()
    wd.mkdir(parents=True, exist_ok=True)
    env = dict(os.environ, OMP_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1", MKL_NUM_THREADS="1")
    res = {}
    for name, ax, mode in RUNS:
        if a.only and name not in a.only:
            continue
        (wd / f"{name}.in").write_text(pw_input(name, ax, mode))
        out = wd / f"{name}.out"
        if a.run and not parse(out)["converged"]:
            t0 = time.time()
            with open(out, "w") as fo:
                subprocess.run(["nice", "-n", "19", f"{QE}/mpirun", "--bind-to", "none", "-np", str(a.ranks), f"{QE}/pw.x", "-nk",
                                str(a.ranks), "-in", f"{name}.in"], cwd=wd, stdout=fo, stderr=subprocess.STDOUT,
                               env=env)
            print(name, "wall", round(time.time() - t0), flush=True)
            subprocess.run(["rm", "-rf", str(wd / f"tmp_{name}")])
        res[name] = parse(out)
        print(name, json.dumps(res[name])[:400], flush=True)
    nat_fu = 2   # formula units per 16-atom cell
    for ax in ("a3", "a1"):
        nm = res.get(f"nm_{ax}", {}).get("E_Ry")
        for m in ("fm", "afm"):
            e = res.get(f"{m}_{ax}", {}).get("E_Ry")
            if nm is not None and e is not None:
                res[f"dE_{m}_{ax}_meV_per_fu"] = round((e - nm) * 13605.693 / nat_fu, 3)
    if res.get("soc_nm_a3", {}).get("E_Ry") is not None and res.get("soc_afm_a3", {}).get("E_Ry") is not None:
        res["dE_soc_afm_vs_soc_nm_meV_per_fu"] = round(
            (res["soc_afm_a3"]["E_Ry"] - res["soc_nm_a3"]["E_Ry"]) * 13605.693 / nat_fu, 3)
    (wd / "barhh6_magnetism.json").write_text(json.dumps(res, indent=1))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
