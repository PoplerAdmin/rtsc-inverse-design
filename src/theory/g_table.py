"""Tables for the lambda-robust moment rule (review R01, F3) and the T-supremum step (F12).

1. mu* = 0, scale-free: f(lam) = Tc/Omega for the Einstein spectrum, on a fine lambda grid.
   g(lam) := lam * Omega_E(lam)^2 with Omega_E = Tc/f(lam), Tc = 300 K
   (minimal lambda<omega^2> at fixed lambda, Sec. 3.3). Checks numerically that
   f is increasing and that h(lam) = f/sqrt(lam) (equivalently g) is monotone.
2. mu* > 0, absolute cutoff omega_c: Omega_E(lam; T) from rho(T) = 0 and
   g_T(lam) = lam * Omega_E^2, for T in [300, 420] K including both sides of every N_c(T) jump;
   reports min_{T >= 300} g_T(lam) * (300/T)^0 -- i.e. the requirement for "some T >= 300 K
   unstable" is lambda<omega^2> >= min_{T>=300} g_T(lam).  Same for the soft-mode limit (N2).
Output: <outdir>/g_table.csv, g_table_mu.csv, tsup_check.txt
"""
import csv
import math
import os
import sys

import numpy as np
from scipy.optimize import brentq

sys.path.insert(0, os.path.dirname(__file__))
from eliashberg import einstein, rho, tc  # noqa
from bounds_extra import rho_soft  # noqa

OUT = sys.argv[1] if len(sys.argv) > 1 else "."
TT = 300.0

# ---------------- 1. mu* = 0
lams = [0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0, 2.25, 2.5, 2.75,
        3.0, 3.25, 3.5, 3.75, 4.0, 4.5, 5.0, 6.0, 7.0, 8.0, 10.0, 15.0, 20.0, 50.0, 100.0]
rows = []
prev_f = prev_h = None
mono_f = mono_h = True
for lam in lams:
    f = tc(einstein(lam, 1.0), xtol_rel=1e-10)
    h = f / math.sqrt(lam)
    if prev_f is not None:
        mono_f &= f > prev_f
        mono_h &= h > prev_h
    prev_f, prev_h = f, h
    Om = TT / f
    rows.append(dict(lam=lam, Tc_over_Omega=f, Tc_over_Omega_sqrtlam=h, Omega_E_300=Om,
                     g_300=lam * Om * Om))
# lam -> infinity: soft-mode limit
Bc = brentq(lambda B: rho_soft(B, TT, 0.0, math.inf, 400), 0.1, 5.0, xtol=1e-12)
c2 = 1 / (2 * math.pi * math.sqrt(Bc))
rows.append(dict(lam="inf", Tc_over_Omega="", Tc_over_Omega_sqrtlam=c2, Omega_E_300="",
                 g_300=(TT / c2) ** 2))
with open(os.path.join(OUT, "g_table.csv"), "w", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
    w.writeheader()
    for r in rows:
        w.writerow({k: (f"{v:.6g}" if isinstance(v, float) else v) for k, v in r.items()})
print(f"mu*=0: f increasing on grid: {mono_f}; f/sqrt(lam) increasing on grid: {mono_h}; c2={c2:.9f}",
      flush=True)
# finer monotonicity check of h on a dense log grid
dense = np.exp(np.linspace(np.log(0.8), np.log(300), 400))
hs = [tc(einstein(l, 1.0), xtol_rel=1e-11) / math.sqrt(l) for l in dense]
dh = np.diff(hs)
print(f"dense grid (400 pts, lam 0.8..300): min diff of f/sqrt(lam) = {dh.min():.3e} "
      f"(all positive: {bool((dh > 0).all())}); max h = {max(hs):.7f} < c2 = {c2:.7f}", flush=True)

# ---------------- 2. mu* > 0, absolute cutoff, T-supremum
def omega_E(lam, T, mu, wc):
    g = lambda Om: rho(einstein(lam, Om), T, mu, wc)
    lo, hi = 50.0, 4000.0
    while g(hi) < 0:
        hi *= 1.5
    return brentq(g, lo, hi, xtol=1e-6)


def Ts_for(wc):
    Ts = set(np.linspace(300, 420, 13))
    for j in range(0, 200):
        Tj = wc / ((2 * j + 1) * math.pi)   # N_c jumps here
        if 300 <= Tj <= 420:
            Ts.update([Tj * (1 - 1e-9), Tj * (1 + 1e-9)])
    return sorted(Ts)


txt = open(os.path.join(OUT, "tsup_check.txt"), "w")
mrows = []
for mu in [0.10, 0.13]:
    for wc in [5000.0, 10000.0, 20000.0]:
        Ts = Ts_for(wc)
        for lam in [1.5, 2.0, 3.0, 4.0]:
            vals = [(T, lam * omega_E(lam, T, mu, wc) ** 2) for T in Ts]
            Tm, gm = min(vals, key=lambda x: x[1])
            g300 = vals[0][1]
            mrows.append(dict(mu=mu, omega_c=wc, lam=lam, g_300=g300, min_over_T=gm, argmin_T=Tm))
            print(f"mu*={mu} wc={wc} lam={lam}: g(300)={g300:.4e}  min_T>=300 g={gm:.4e} at T={Tm:.2f}",
                  file=txt, flush=True)
        # soft-mode (N2) requirement (T/c2(T))^2 = 4 pi^2 T^2 B_c(T)
        sv = []
        for T in Ts:
            B = brentq(lambda B: rho_soft(B, T, mu, wc, 300), 0.1, 5.0, xtol=1e-10)
            sv.append((T, 4 * math.pi ** 2 * T * T * B))
        Tm, vm = min(sv, key=lambda x: x[1])
        print(f"mu*={mu} wc={wc} soft-mode N2: req(300)={sv[0][1]:.4e}  min_T>=300={vm:.4e} at T={Tm:.2f}",
              file=txt, flush=True)
        mrows.append(dict(mu=mu, omega_c=wc, lam="inf", g_300=sv[0][1], min_over_T=vm, argmin_T=Tm))
txt.close()
with open(os.path.join(OUT, "g_table_mu.csv"), "w", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=list(mrows[0].keys()))
    w.writeheader()
    for r in mrows:
        w.writerow({k: (f"{v:.6g}" if isinstance(v, float) else v) for k, v in r.items()})
print("done", flush=True)
