#!/usr/bin/env python3
"""Calibrate coarse (tier-2) DFPT alpha^2F against a full reference: A, lambda, omega_log, m.

For every spectrum directory (matdyn a2F.dos1..N) this computes, per smearing sigma_i,
    A      = int alpha^2F dw                 (K)    -- QE bin rule: sum_{E>0} a2F_n * Delta
    lambda = 2 int alpha^2F / w dw
    w_log  = exp[(2/lambda) int alpha^2F ln(w) / w dw]   (K)
    m      = lambda <w^2> = 2 int w alpha^2F dw         (K^2)
using the signed a2F_total at E > 0 (the same quantities rtt.core.dfpt tags excluded
smearings with), and applies rtt.core.dfpt's mode-positive smearing admission (standalone copy below).

Usage:
  python coarse_calibration.py --out-sigma per_sigma.csv --out-run per_run.csv \
      Al,ref,<dir> Al,P1_q3,<dir> ... H3S,ref,<dir> H3S,c2,<dir>
Every spec "MATERIAL,TAG,DIR": TAG 'ref' marks the reference of that material.

All numbers are informational (tier 2); nothing here issues a status.
"""
from __future__ import annotations

import csv
import json
import math
import re
import sys
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

import numpy as np

# QE's own constant (rtt/core/dfpt/units.py: RY_TO_K = RYDBERG_SI / K_BOLTZMANN_SI)
RY_TO_K = 157887.512402

# ---------------------------------------------------------------------------------------------
# Standalone copy (trimmed, logic unchanged) of the certified engine's QE a2F handling, so this
# script runs without the package:
#   src/rtt/core/dfpt/parsers.py  (QEA2FDos, parse_qe_a2f_dos, mode_negative_fractions)
#   src/rtt/core/dfpt/workflow.py (admit_smearings, MODE_NEGATIVE_TOL_REL)
# ---------------------------------------------------------------------------------------------
MODE_NEGATIVE_TOL_REL = 1e-6
_FOOTER = re.compile(r"lambda\s*=\s*(\S+)\s+Delta\s*=\s*(\S+)", re.I)


def _ffloat(s: str) -> float:
    return float(s.replace("D", "E").replace("d", "e"))


def _numeric_row(s: str) -> list[float] | None:
    try:
        return [_ffloat(x) for x in s.split()]
    except ValueError:
        return None


@dataclass
class QEA2FDos:
    path: str
    isig: int | None
    E_Ry: np.ndarray
    a2f_total: np.ndarray
    a2f_modes: np.ndarray
    lambda_reported: float | None
    delta_Ry: float | None
    header: list[str] = field(default_factory=list)


def parse_qe_a2f_dos(path: str | Path) -> QEA2FDos:
    p = Path(path)
    m = re.search(r"a2F\.dos(\d+)$", p.name)
    rows, header = [], []
    lam = delta = None
    for line in p.read_text(errors="replace").splitlines():
        s = line.strip()
        if not s:
            continue
        mf = _FOOTER.search(s)
        if mf:
            if lam is not None:
                raise ValueError(f"{p}: more than one footer")
            lam, delta = _ffloat(mf.group(1)), _ffloat(mf.group(2))
            continue
        if lam is not None:
            raise ValueError(f"{p}: content after the footer")
        if s.startswith("#"):
            header.append(s)
            continue
        r = _numeric_row(s)
        if r is None:
            header.append(s)
            continue
        rows.append(r)
    if lam is None:
        raise ValueError(f"{p}: no 'lambda = ... Delta =' footer (truncated file)")
    if not any("Rydberg" in h for h in header):
        raise ValueError(f"{p}: not a QE matdyn a2F.dos file")
    a = np.asarray(rows, dtype=float)
    return QEA2FDos(str(p), int(m.group(1)) if m else None, a[:, 0], a[:, 1], a[:, 2:], lam, delta, header)


def mode_negative_fractions(tab: QEA2FDos) -> list[float]:
    pos = tab.E_Ry > 0
    A = float(np.sum(np.clip(tab.a2f_total[pos], 0.0, None)))
    modes = tab.a2f_modes[pos]
    return [float(-np.sum(np.clip(modes[:, k], None, 0.0))) / max(A, 1e-300)
            for k in range(modes.shape[1])]


def admit_smearings(tabs: dict[int, QEA2FDos], tol_rel: float = MODE_NEGATIVE_TOL_REL
                    ) -> tuple[list[int], dict[str, Any], str | None]:
    """Rule 'mode_positive' (rtt.core.dfpt): smearing i passes iff every mode's alpha^2F_nu is
    non-negative at omega > 0 (within tol_rel of the total mass). Admitted = contiguous block from
    the smallest passing smearing upwards; a failure above a passing smearing is an error."""
    order = sorted(tabs)
    neg = {i: mode_negative_fractions(tabs[i]) for i in order}
    passing = [i for i in order if max(neg[i], default=0.0) <= tol_rel]
    report: dict[str, Any] = {
        "rule": "mode_positive", "tol_rel": tol_rel,
        "max_mode_negative_fraction": {str(i): max(neg[i], default=0.0) for i in order}}
    if not passing:
        report.update(admitted=[], excluded=order)
        return [], report, "no smearing has non-negative per-mode alpha^2F"
    block = [i for i in order if i >= passing[0]]
    holes = [i for i in block if i not in passing]
    report.update(admitted=block, excluded=[i for i in order if i < passing[0]])
    if holes:
        return [], report, f"smearings {holes} fail per-mode positivity above passing {passing[0]}"
    return block, report, None
# ------------------------------------------------------------------------- end of copied code


def moments(tab: QEA2FDos) -> dict[str, float]:
    pos = tab.E_Ry > 0
    w = tab.E_Ry[pos] * RY_TO_K
    a = tab.a2f_total[pos]
    d = tab.delta_Ry * RY_TO_K
    A = float(np.sum(a) * d)
    lam = float(2 * np.sum(a / w) * d)
    wlog = float(math.exp(2 / lam * np.sum(a / w * np.log(w)) * d)) if lam > 0 else float("nan")
    m = float(2 * np.sum(a * w) * d)
    return {"A_K": A, "lambda": lam, "wlog_K": wlog, "m_K2": m, "lambda_qe_footer": tab.lambda_reported}


def sigmas_from(d: Path, n: int) -> list[float]:
    """Smearing values (Ry) from ph.in (el_ph_sigma * i)."""
    txt = (d / "ph.in").read_text() if (d / "ph.in").exists() else ""
    m = re.search(r"el_ph_sigma\s*=\s*([0-9.eEdD+-]+)", txt)
    s0 = _ffloat(m.group(1)) if m else 0.01
    return [round(s0 * i, 6) for i in range(1, n + 1)]


def run_meta(d: Path) -> dict[str, Any]:
    meta: dict[str, Any] = {}
    j = d / "coarse_run.json"
    if j.exists():
        r = json.loads(j.read_text())
        s = r["settings"]
        meta.update(ecutwfc=s["ecutwfc"], kscf=s["kscf"], kfit=s["kfit"], nq=s["nq"], ranks=s["ranks"],
                    ph_wall_s=r["steps"].get("ph", {}).get("wall_s"), wall_total_s=r.get("wall_total_s"))
        pm = r.get("phonon_min", {})
        dq = pm.get("dyn_min_cm1", {})
        meta["min_freq_dyn_nonGamma_cm1"] = min((v for k, v in dq.items() if k != "dyn1"), default=None)
        meta["min_freq_mesh_cm1"] = pm.get("matdyn_mesh_min_cm1")
        return meta
    # older runs: parse inputs (+ timings.txt from coarse_dfpt_probe.sh if present)
    fit = (d / "scf_fit.in").read_text() if (d / "scf_fit.in").exists() else ""
    ph = (d / "ph.in").read_text() if (d / "ph.in").exists() else ""
    mk = re.search(r"K_POINTS.*\n\s*(\d+)", fit)
    mk2 = re.search(r"K_POINTS.*\n\s*(\d+)", (d / "scf.in").read_text()) if (d / "scf.in").exists() else None
    meta.update(ecutwfc=float(_ffloat(re.search(r"ecutwfc\s*=\s*([0-9.]+)", fit).group(1))) if fit else None,
                kfit=int(mk.group(1)) if mk else None, kscf=int(mk2.group(1)) if mk2 else None,
                nq=int(re.search(r"nq1\s*=\s*(\d+)", ph).group(1)) if ph else None)
    t = d / "timings.txt"
    if t.exists():
        tt = t.read_text()
        mp = re.search(r"^ph (\S+) s", tt, re.M)
        meta["ph_wall_s"] = float(mp.group(1)) if mp else None
        mr = re.search(r"ranks=(\d+)", tt)
        meta["ranks"] = int(mr.group(1)) if mr else None
    elif (d / "ph.out").exists():
        mw = re.findall(r"PHONON\s*:.*CPU\s+(.*)WALL", (d / "ph.out").read_text(errors="replace"))
        meta["ph_wall_str"] = mw[-1].strip() if mw else None
    vals = []
    for f in d.glob("dyn[0-9]*"):
        if f.name in ("dyn0", "dyn1"):
            continue
        vals += [float(x) for x in re.findall(r"freq\s*\(\s*\d+\)\s*=\s*\S+\s*\[THz\]\s*=\s*(\S+)\s*\[cm-1\]",
                                               f.read_text(errors="replace"))]
    meta["min_freq_dyn_nonGamma_cm1"] = min(vals) if vals else None
    fd = d / "freq.dos"
    if fd.exists():      # matdyn flfrq: header, then (q line, frequency line) pairs, cm^-1
        lines = fd.read_text(errors="replace").splitlines()[1:]
        nums = [float(x) for ln in lines[1::2] for x in ln.split()]
        meta["min_freq_mesh_cm1"] = min(nums) if nums else None
    return meta


def load(d: Path) -> tuple[dict[int, QEA2FDos], list[int], dict, str | None]:
    files = sorted(d.glob("a2F.dos*"), key=lambda p: int(re.findall(r"\d+$", p.name)[0]))
    tabs = {int(re.findall(r"\d+$", f.name)[0]): parse_qe_a2f_dos(f) for f in files}
    adm, rep, err = admit_smearings(tabs)
    return tabs, adm, rep, err


def main(argv: list[str]) -> int:
    out_sigma = out_run = None
    specs = []
    it = iter(argv)
    for a in it:
        if a == "--out-sigma":
            out_sigma = next(it)
        elif a == "--out-run":
            out_run = next(it)
        else:
            specs.append(a.split(",", 2))
    data: dict[tuple[str, str], dict] = {}
    for mat, tag, d in specs:
        p = Path(d).expanduser()
        tabs, adm, rep, err = load(p)
        sig = sigmas_from(p, len(tabs))
        data[(mat, tag)] = {"dir": p, "tabs": tabs, "adm": adm, "err": err, "rep": rep, "sig": sig,
                            "mom": {i: moments(t) for i, t in tabs.items()}, "meta": run_meta(p)}

    rows_s, rows_r = [], []
    for (mat, tag), c in data.items():
        ref = data.get((mat, "ref"))
        for i in sorted(c["tabs"]):
            r = {"material": mat, "run": tag, "isig": i, "sigma_Ry": c["sig"][i - 1],
                 "admitted": i in c["adm"],
                 "max_mode_neg_frac": c["rep"]["max_mode_negative_fraction"][str(i)]}
            r.update({k: c["mom"][i][k] for k in ("A_K", "lambda", "wlog_K", "m_K2", "lambda_qe_footer")})
            if ref is not None and tag != "ref" and i in ref["mom"]:
                if abs(ref["sig"][i - 1] - c["sig"][i - 1]) > 1e-9:
                    raise SystemExit(f"{mat}/{tag}: sigma ladder differs from reference")
                r["ref_admitted"] = i in ref["adm"]
                for k in ("A_K", "lambda", "wlog_K", "m_K2"):
                    r[f"ratio_{k}"] = c["mom"][i][k] / ref["mom"][i][k]
            rows_s.append(r)
        # per-run summary over sigma admitted in both coarse and reference
        if ref is None or tag == "ref":
            continue
        both = [i for i in c["adm"] if i in ref["adm"]]
        rr = {"material": mat, "run": tag, **{k: c["meta"].get(k) for k in
              ("ecutwfc", "kscf", "kfit", "nq", "ranks", "ph_wall_s", "min_freq_dyn_nonGamma_cm1",
               "min_freq_mesh_cm1")},
              "coarse_admitted": " ".join(map(str, c["adm"])) or "none",
              "ref_admitted": " ".join(map(str, ref["adm"])),
              "n_common": len(both), "admission_error": c["err"] or ""}
        for k in ("A_K", "lambda", "wlog_K", "m_K2"):
            rat = [c["mom"][i][k] / ref["mom"][i][k] for i in both]
            rr[f"{k}_ratio_min"] = min(rat) if rat else None
            rr[f"{k}_ratio_max"] = max(rat) if rat else None
        common = [i for i in c["mom"] if i in ref["mom"]]
        rat_all = [c["mom"][i]["A_K"] / ref["mom"][i]["A_K"] for i in common]
        rr["A_K_ratio_min_allsig"] = min(rat_all)
        rr["A_K_ratio_max_allsig"] = max(rat_all)
        # deferral-relevant: max over the coarse run's own admitted ladder vs reference A_nom and
        # max over the reference's admitted ladder
        if c["adm"]:
            a_c_max = max(c["mom"][i]["A_K"] for i in c["adm"])
            a_ref_nom = ref["mom"][ref["adm"][0]]["A_K"]
            a_ref_max = max(ref["mom"][i]["A_K"] for i in ref["adm"])
            rr["A_coarse_maxadm_K"] = a_c_max
            rr["A_ref_nom_K"] = a_ref_nom
            rr["A_ref_maxadm_K"] = a_ref_max
            rr["defer_ratio_vs_ref_max"] = a_c_max / a_ref_max
        rows_r.append(rr)

    def dump(rows: list[dict], path: str | None) -> None:
        keys: list[str] = []
        for r in rows:
            keys += [k for k in r if k not in keys]
        fh = open(path, "w", newline="") if path else sys.stdout
        w = csv.DictWriter(fh, fieldnames=keys)
        w.writeheader()
        for r in rows:
            w.writerow({k: (f"{v:.6g}" if isinstance(v, float) else v) for k, v in r.items()})
        if path:
            fh.close()

    dump(rows_s, out_sigma)
    dump(rows_r, out_run)
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
