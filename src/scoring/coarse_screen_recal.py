"""Doc 16: recalibration of the S5 coarse screen (T2-SCREEN) on doc 13's references + Mg2IrH6.

For each coarse (q 3^3, q 4^3) pair this computes, from the a2F.dosN files:
  * the doc 13 statistic  A_hi,old = rtt.dft.coarse.coarse_envelope   (admitted rows only)
  * the doc 16 statistic  A_hi,new = rtt.dft.coarse.coarse_envelope_all (all rows, sigma->0)
and divides both by the reference's nominal A and by the reference's A_hi lower bound
(doc 13 §2 comparison 2(b): A_nom,ref (1 + d_sig,ref)(1 + 0.10)/(1 - 0.02)).

Reference numbers for Al, Nb, Pb, H3S are doc 13's (data/13_coarse_envelope.csv). For Mg2IrH6
the reference is the P1 control (doc 15): 80 Ry, q 3^3, k 6^3 / k_fit 36^3, sigma 0.005..0.05 Ry
(all rows; its ladder fails mode-positive admission at tol 1e-6 on every sigma, see doc 16 §2).

Usage (WSL, repo root):
  PYTHONPATH=src python src/scoring/coarse_screen_recal.py > docs/design/data/16_coarse_screen_recal.csv
"""
from __future__ import annotations

import csv
import os
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from sigma_ladder_table import ladder  # noqa: E402

from rtt.dft.coarse import (CoarseMetrics, SigmaRow, a_top, coarse_envelope,  # noqa: E402
                            coarse_envelope_all)

EPS_SYS = 0.10
HOME = Path(os.path.expanduser("~"))
K1 = HOME / "rtsc_targeted_work/K1-calibrate/runs"
T4 = HOME / "rtsc_targeted_work/T4-surrogates"
E1 = HOME / "rtsc_targeted_work/E1-campaign/work/c-9a8fa100a45b/S5"
P1 = Path(__file__).resolve().parents[2] / "runs/controls/mg2irh6"

# doc 13 data/13_coarse_envelope.csv: (A_ref_nom, A_ref_hi_lower_bound)
REF = {"Al": (67.53, 76.31), "Nb": (121.8, 180.2), "Pb": (36.14, 42.45), "H3S": (1648.0, 2482.0)}

PAIRS = [("Al", "P1", K1 / "Al_P1_q3", K1 / "Al_P1_q4"), ("Al", "P2", K1 / "Al_P2_q3", K1 / "Al_P2_q4"),
         ("Nb", "P1", K1 / "Nb_P1_q3", K1 / "Nb_P1_q4"), ("Nb", "P2", K1 / "Nb_P2_q3", K1 / "Nb_P2_q4"),
         ("Pb", "P1", K1 / "Pb_P1_q3", K1 / "Pb_P1_q4"), ("Pb", "P2", K1 / "Pb_P2_q3", K1 / "Pb_P2_q4"),
         ("H3S", "T4c2/c3", T4 / "h3s_c2", T4 / "h3s_c3"),
         ("Mg2IrH6", "E1 60Ry", E1 / "q3x3x3", E1 / "q4x4x4")]
EXTRA = os.environ.get("RECAL_EXTRA")   # "label=q3dir,q4dir" for the 80 Ry re-run (doc 16 §6)


def metrics_from_dir(d: Path, q: int) -> CoarseMetrics:
    rows = [SigmaRow(index=r["i"], sigma_Ry=r["sigma_Ry"], admitted=r["admitted"], A_K=r["A_K"],
                     lambda_=r["lambda"], omega_log_K=r["omega_log_K"]) for r in ladder(d)]
    adm_ok = any(r.admitted for r in rows)
    return CoarseMetrics(rows=rows, plateau_sigma_min_Ry=None, min_freq_matdyn_cm1=None,
                         q=[q, q, q], admission_ok=adm_ok)


def mg2irh6_ref() -> tuple[float, float, str]:
    best = None
    for run in ("eph_q3_k6_kf36", "eph_q3_k8_kf24"):
        rows = ladder(P1 / run)
        a = [r["A_K"] for r in rows]
        nom, lo = a[0], min(a)
        hi = nom * (1 + (nom - lo) / nom) * (1 + EPS_SYS) / (1 - 0.02)
        if best is None or nom > best[0]:
            best = (nom, hi, run)
    return best


def main() -> None:
    refs = dict(REF)
    m_nom, m_hi, m_run = mg2irh6_ref()
    refs["Mg2IrH6"] = (m_nom, m_hi)
    pairs = list(PAIRS)
    if EXTRA:
        lab, _, ds = EXTRA.partition("=")
        a, b = ds.split(",")
        pairs.append(("Mg2IrH6", lab, Path(a), Path(b)))
    w = csv.writer(sys.stdout)
    w.writerow(["material", "profile", "A_ref_nom_K", "A_ref_hi_lb_K", "A_hi_old_K", "old_over_nom",
                "old_over_hi_lb", "A_top_q3_K", "A_top_q4_K", "A_nom_new_K", "d_sig", "d_q", "eps_w",
                "A_hi_new_K", "new_over_nom", "new_over_hi_lb", "admitted_q3", "admitted_q4"])
    for mat, prof, d3, d4 in pairs:
        c, f = metrics_from_dir(d3, 3), metrics_from_dir(d4, 4)
        nom, hilb = refs[mat]
        old = coarse_envelope(f, c, EPS_SYS)
        new = coarse_envelope_all(f, c, EPS_SYS)
        ao = old["A_hi_K"] if old["ok"] else None
        an = new["A_hi_K"] if new["ok"] else None

        def r(x, y):
            return "" if x is None else f"{x / y:.3f}"
        w.writerow([mat, prof, f"{nom:.1f}", f"{hilb:.1f}", "" if ao is None else f"{ao:.1f}",
                    r(ao, nom), r(ao, hilb), f"{a_top(c):.1f}", f"{a_top(f):.1f}",
                    f"{new.get('A_nom_K', float('nan')):.1f}", f"{new.get('d_sig_A', float('nan')):.3f}",
                    f"{new.get('d_q_A', float('nan')):.3f}", f"{new.get('eps_omega', float('nan')):.3f}",
                    "" if an is None else f"{an:.1f}", r(an, nom), r(an, hilb),
                    "-".join(str(x.index) for x in c.rows if x.admitted) or "none",
                    "-".join(str(x.index) for x in f.rows if x.admitted) or "none"])
    print(f"# Mg2IrH6 reference: {m_run} A(0.005 Ry) = {m_nom:.1f} K, A_hi lower bound {m_hi:.1f} K",
          file=sys.stderr)


if __name__ == "__main__":
    main()
