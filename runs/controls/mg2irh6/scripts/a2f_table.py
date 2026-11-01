#!/usr/bin/env python3
"""alpha^2F summary for QE matdyn a2F.dosN files, using src/scoring/a2f_summary.py's formulas
(allen_dynes, RY_TO_K) but integrating on the FULL uniform matdyn grid.

Why not a2f_summary.summarise directly: it first drops every point with a2F <= 0 and then takes
np.gradient of the remaining (non-uniform) frequencies. For a spectrum with gaps (zero a2F
between the Ir/Mg, H-bend and H-stretch bands, as in Mg2IrH6) the points at the gap edges
get dw ~ half the gap width, so lambda, A and Tc are overestimated (here lambda 3.68 vs 2.42
at sigma = 0.005 Ry; QE's own footer gives 2.4225). Here negative a2F values are clipped to 0
(their weight is reported) and dw is the uniform grid step.

Usage: a2f_table.py <dir> [--json out.json]
"""
import glob
import json
import math
import os
import re
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[4] / "src" / "scoring"))
from a2f_summary import RY_TO_K, allen_dynes  # noqa: E402


def load(path):
    w, a, lam_qe = [], [], None
    for line in open(path):
        s = line.strip()
        if s.startswith("lambda"):
            lam_qe = float(s.split()[2])
            continue
        if not s or s.startswith("#"):
            continue
        p = s.split()
        w.append(float(p[0]))
        a.append(float(p[1]))
    return np.array(w) * RY_TO_K, np.array(a), lam_qe


def summ(w, a):
    dw = w[1] - w[0]
    neg = -a[a < 0].sum() / a[a > 0].sum()
    a = np.clip(a, 0, None)
    lam = 2 * np.sum(a / w) * dw
    A = np.sum(a) * dw
    wlog = math.exp(2 / lam * np.sum(a / w * np.log(w)) * dw)
    w2 = math.sqrt(2 / lam * np.sum(a * w) * dw)
    return dict(lam=lam, wlog_K=wlog, w2_K=w2, A_K=A, A_over_4p33148_K=A / 4.33148,
                TcAD_mu010_K=allen_dynes(lam, wlog, w2, 0.10), TcAD_mu013_K=allen_dynes(lam, wlog, w2, 0.13),
                neg_weight_frac=neg)


def main():
    d = sys.argv[1]
    files = sorted(glob.glob(os.path.join(d, "a2F.dos*")), key=lambda p: int(re.findall(r"\d+$", p)[0]))
    rows = []
    print("isig sigma_Ry  lam(QE)  lam    wlog_K  w2_K    A_K     A/4.33  TcAD(0.10) TcAD(0.13) negfrac")
    for f in files:
        i = int(re.findall(r"\d+$", f)[0])
        w, a, lq = load(f)
        r = summ(w, a)
        r.update(isig=i, lambda_qe_footer=lq)
        rows.append(r)
        print(f"{i:4d} {0.005 * i:.3f}   {lq:6.3f}  {r['lam']:6.3f} {r['wlog_K']:7.1f} {r['w2_K']:7.1f} "
              f"{r['A_K']:7.1f} {r['A_over_4p33148_K']:6.1f} {r['TcAD_mu010_K']:9.1f} {r['TcAD_mu013_K']:9.1f} {r['neg_weight_frac']:.1e}")
    if "--json" in sys.argv:
        json.dump(rows, open(sys.argv[sys.argv.index("--json") + 1], "w"), indent=1)


if __name__ == "__main__":
    main()
