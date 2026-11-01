"""Extra extremal computations for 01_eliashberg_targets.md (all at T = 300 K, absolute cutoff).

1. soft-mode (omega -> 0) kernel lam(k) = B/k^2: critical B_c(mu*, omega_c), which gives the
   constant c2 in  Tc <= c2 * sqrt(lambda <omega^2>)  (c2 = 1/(2 pi sqrt(B_c))).
2. Two-atom extremal problem at fixed lambda = L (extreme points of {A fixed, lambda fixed}
   have <= 2 atoms): min A, min lambda<omega^2>, and the omega_log of the best spectra, over a
   grid of (frequency ratio r, high-atom lambda share f); the scale is solved so that Tc = 300 K.
3. Demonstration that Tc/omega_log is unbounded at fixed lambda (two atoms, growing ratio r).
Outputs CSVs in <outdir>.
"""
import csv
import math
import os
import sys

import numpy as np
from scipy.linalg import eigh
from scipy.optimize import brentq

sys.path.insert(0, os.path.dirname(__file__))
from eliashberg import Spectrum, rho, tc  # noqa

OUT = sys.argv[1] if len(sys.argv) > 1 else "."
TT = 300.0
MUS = [0.0, 0.08, 0.10, 0.13, 0.16, 0.20]
WC = 20000.0


def rho_soft(B, T, mu, wc, N=400):
    """Top eigenvalue for lam(k) = B/k^2 (k >= 1). Tail ~ B/k^2 converges; N large."""
    k = np.arange(2 * N + 1, dtype=float)
    L = np.zeros_like(k)
    L[1:] = B / k[1:] ** 2
    n = np.arange(N)
    M = L[np.abs(n[:, None] - n[None, :])] + L[n[:, None] + n[None, :] + 1]
    np.fill_diagonal(M, L[2 * n + 1])
    S = np.concatenate([[0.0], np.cumsum(L[1:N])])
    M[n, n] -= (2 * n + 1) + 2 * S[:N]
    if mu > 0:
        Nc = int(np.sum((2 * n + 1) * math.pi * T <= wc))
        M[:Nc, :Nc] -= 2 * mu
    return float(eigh(M, subset_by_index=[N - 1, N - 1], eigvals_only=True)[0])


def write(path, rows):
    with open(path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        w.writeheader()
        for r in rows:
            w.writerow({k: (f"{v:.6g}" if isinstance(v, float) else v) for k, v in r.items()})


def two_atom(L, r, f, s):
    om = np.array([s, s * r])
    lam = np.array([(1 - f) * L, f * L])
    return Spectrum(om, lam * om / 2)


def scale_for_300(L, r, f, mu, wc):
    g = lambda s: rho(two_atom(L, r, f, s), TT, mu, wc)
    lo, hi = 10.0, 200.0
    while g(hi) < 0:
        hi *= 2
        if hi > 1e6:
            return None
    if g(lo) >= 0:
        return None
    return brentq(g, lo, hi, xtol=1e-6 * hi)


def main():
    # 1. soft-mode constant
    soft = []
    for mu in MUS:
        for N in (200, 400):
            Bc = brentq(lambda B: rho_soft(B, TT, mu, WC, N), 0.1, 5.0, xtol=1e-10)
            soft.append(dict(mu=mu, omega_c=WC, N=N, B_c=Bc, c2=1 / (2 * math.pi * math.sqrt(Bc)),
                             lamw2_req_300=(TT / (1 / (2 * math.pi * math.sqrt(Bc)))) ** 2))
            print(soft[-1], flush=True)
    write(os.path.join(OUT, "soft_mode_constant.csv"), soft)

    # 2. two-atom extremal at fixed lambda
    rs = [1.0, 1.25, 1.5, 2.0, 3.0, 4.0, 6.0, 8.0, 12.0, 20.0]
    fs = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9]
    two = []
    for mu in [0.0, 0.10, 0.13, 0.16]:
        for L in [1.0, 1.5, 2.0, 2.5, 3.0, 4.0]:
            best_A, best_m2, best_wlog = None, None, None
            for r in rs:
                for f in (fs if r > 1 else [0.5]):
                    s = scale_for_300(L, r, f, mu, WC)
                    if s is None:
                        continue
                    sp = two_atom(L, r, f, s)
                    rec = dict(r=r, f=f, A=sp.A, lamw2=sp.lam_omega2, wlog=sp.omega_log,
                               w2=sp.omega2, w1=sp.omega1, wmax=sp.omega_max)
                    if best_A is None or rec["A"] < best_A["A"]:
                        best_A = rec
                    if best_m2 is None or rec["lamw2"] < best_m2["lamw2"]:
                        best_m2 = rec
                    if best_wlog is None or rec["wlog"] < best_wlog["wlog"]:
                        best_wlog = rec
            row = dict(mu=mu, omega_c=WC, lam=L,
                       minA=best_A["A"], minA_r=best_A["r"], minA_f=best_A["f"],
                       minA_wlog=best_A["wlog"], minA_w2=best_A["w2"],
                       min_lamw2=best_m2["lamw2"], min_lamw2_r=best_m2["r"], min_lamw2_f=best_m2["f"],
                       min_w2=math.sqrt(best_m2["lamw2"] / L),
                       einstein_wlog=None, min_wlog_on_grid=best_wlog["wlog"],
                       min_wlog_r=best_wlog["r"], min_wlog_f=best_wlog["f"])
            se = scale_for_300(L, 1.0, 0.5, mu, WC)
            row["einstein_wlog"] = se
            two.append(row)
            print(row, flush=True)
    write(os.path.join(OUT, "two_atom_fixed_lambda.csv"), two)

    # 3. unbounded Tc/omega_log at fixed lambda
    unb = []
    for r in [1, 2, 4, 8, 16, 32, 64, 128]:
        L, f = 2.0, 0.5
        sp = two_atom(L, float(r), f, 100.0)
        T = tc(sp, 0.0)
        unb.append(dict(lam=L, f_hi=f, r=r, Tc_over_wlog=T / sp.omega_log,
                        Tc_over_A=T / sp.A, Tc_over_sqrt_lamw2=T / math.sqrt(sp.lam_omega2)))
        print(unb[-1], flush=True)
    write(os.path.join(OUT, "unbounded_wlog_demo.csv"), unb)
    print("done", flush=True)


if __name__ == "__main__":
    main()
