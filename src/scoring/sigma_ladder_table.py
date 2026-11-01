"""Per-smearing A, lambda, omega_log and mode-positive admission for QE a2F.dosN directories.

Used by doc 16 (coarse-screen fix) to compare coarse S5 ladders with the P1 control runs and the
doc 13 references on a common footing. Moments use the QE bin rule on the signed a2F_total at
E > 0 with the bin width from the file footer (same as ``rtt.dft.coarse`` for excluded rows and
doc 13 §2); admission is the core ``admit_smearings`` (rule mode_positive, tol 1e-6).

Usage (WSL venv, from the repo root):
    PYTHONPATH=src python src/scoring/sigma_ladder_table.py <label>=<dir> [<label>=<dir> ...]
Prints one CSV row per (run, sigma): label, i, sigma_Ry, admitted, max_neg_frac, A_K, lambda,
omega_log_K.
"""
from __future__ import annotations

import csv
import math
import sys
from pathlib import Path

import numpy as np

from rtt.core.dfpt.parsers import (mode_negative_fractions, parse_qe_a2f_dos,
                                           parse_qe_lambda_file)
from rtt.core.dfpt.units import to_kelvin_factor
from rtt.core.dfpt.workflow import admit_smearings


def ladder(d: Path) -> list[dict]:
    tabs = {}
    for i in range(1, 31):
        p = d / f"a2F.dos{i}"
        if p.is_file():
            tabs[i] = parse_qe_a2f_dos(p)
    if not tabs:
        return []
    adm, rep, err = admit_smearings(tabs)
    sig = {}
    if (d / "lambda").is_file():
        for k, e in enumerate(parse_qe_lambda_file(d / "lambda"), start=1):
            sig[k] = e.broadening_Ry
    f = to_kelvin_factor("Ry")
    rows = []
    for i, t in sorted(tabs.items()):
        pos = t.E_Ry > 0
        a = t.a2f_total[pos]
        E = t.E_Ry[pos]
        A = float(np.sum(a) * t.delta_Ry * f)
        lam = float(2.0 * np.sum(a / E) * t.delta_Ry)
        wlog = (float(math.exp(2.0 / lam * np.sum(a * np.log(E * f) / E) * t.delta_Ry))
                if lam > 0 else float("nan"))
        rows.append({"i": i, "sigma_Ry": sig.get(i), "admitted": (i in adm) and err is None,
                     "max_neg_frac": max(mode_negative_fractions(t), default=0.0),
                     "A_K": round(A, 1), "lambda": round(lam, 3), "omega_log_K": round(wlog, 1),
                     "admission_error": err is not None})
    return rows


def main(argv: list[str]) -> None:
    w = csv.writer(sys.stdout)
    w.writerow(["label", "i", "sigma_Ry", "admitted", "max_neg_frac", "A_K", "lambda",
                "omega_log_K", "admission_error"])
    for arg in argv:
        label, _, path = arg.partition("=")
        for r in ladder(Path(path).expanduser()):
            w.writerow([label, r["i"], r["sigma_Ry"], int(r["admitted"]), f"{r['max_neg_frac']:.2e}",
                        r["A_K"], r["lambda"], r["omega_log_K"], int(r["admission_error"])])


if __name__ == "__main__":
    main(sys.argv[1:])
