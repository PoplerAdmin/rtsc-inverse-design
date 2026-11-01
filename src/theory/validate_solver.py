"""Validation of src/theory/eliashberg.py.

1. mu* = 0 Einstein optimum: min over w of the critical a = A/Tc must reproduce the certified
   constant a* in [4.3314814, 4.3314815] (wB certificate; rtt/core/certify/PROOFS.md §6) at w* ~ 7.5749.
2. Large-lambda asymptote: Tc / sqrt(lambda <omega^2>) -> 0.1827 (Allen & Dynes 1975).
3. Allen-Dynes comparison for Pb-like and Al-like model spectra and an Einstein lambda scan.
4. Truncation convergence (trunc_mult 30 vs 60).
Run:  nice -n 19 python validate_solver.py > ../../docs/theory/data/validation.txt
"""
import math
import sys
import os

import numpy as np
from scipy.optimize import minimize_scalar

sys.path.insert(0, os.path.dirname(__file__))
from eliashberg import (Spectrum, einstein, gaussians, tc, allen_dynes, rho)  # noqa: E402


def a_crit(w):
    """Einstein, mu*=0: A/Tc at Omega/Tc = w, from Tc of a spectrum with Omega = w (T units)."""
    # Tc of Einstein(lam, Omega) is Omega * F(lam); find lam such that Omega/Tc = w.
    from scipy.optimize import brentq
    g = lambda lam: 1.0 / tc(einstein(lam, 1.0), xtol_rel=1e-11) - w
    lam = brentq(g, 0.3, 20, xtol=1e-12)
    return lam * w / 2.0, lam


print("== 1. Einstein mu*=0 optimum (certified a* in [4.3314814, 4.3314815]) ==")
res = minimize_scalar(lambda w: a_crit(w)[0], bracket=(6.5, 7.5, 8.5), tol=1e-8)
a, lam = a_crit(res.x)
print(f"w* = {res.x:.5f}  a* = {a:.8f}  lambda* = {lam:.5f}  1/a* = {1/a:.8f}")

print("\n== 2. Large-lambda asymptote, Einstein, mu*=0: Tc/(Omega sqrt(lam)) ==")
for lam in [1, 2, 3, 5, 10, 30, 100, 1000]:
    t = tc(einstein(lam, 1.0))
    print(f"lam={lam:6g}  Tc/Omega={t:.5f}  Tc/sqrt(lam<w2>)={t/math.sqrt(lam):.5f}")

print("\n== 3. Allen-Dynes comparison (mu* at omega_c = 10 omega_max; AD with corrections) ==")
models = {
    "Pb-like": gaussians([51.0, 99.0], [7.0, 6.0], [0.72, 0.28], 1.55),
    "Al-like": gaussians([170.0, 410.0], [50.0, 25.0], [0.55, 0.45], 0.43),
}
for name, sp in models.items():
    m = sp.moments()
    for mu in [0.0, 0.10, 0.13]:
        for cm in [5, 10]:
            T = tc(sp, mu, cm * sp.omega_max)
            ad = allen_dynes(sp.lam, mu, sp.omega_log, sp.omega2)
            print(f"{name}: lam={m['lam']:.3f} wlog={m['omega_log']:.1f}K w2={m['omega2']:.1f}K "
                  f"mu*={mu:.2f} wc={cm}wmax  Tc_ME={T:.3f}K  Tc_AD={ad:.3f}K  ratio={T/ad if ad else float('nan'):.3f}")
print("\nEinstein scan, Tc/Omega (ME, mu* at 10 Omega) vs AD:")
for mu in [0.0, 0.1, 0.15]:
    for lam in [0.5, 1.0, 1.5, 2.0, 3.0, 5.0]:
        sp = einstein(lam, 1000.0)
        T = tc(sp, mu, 10 * sp.omega_max)
        ad = allen_dynes(lam, mu, 1000.0, 1000.0)
        print(f"mu*={mu:.2f} lam={lam:4.1f}  Tc_ME={T:8.3f}K  Tc_AD={ad:8.3f}K  ratio={T/ad if ad else float('nan'):.3f}")

print("\n== 4. Truncation convergence ==")
for name, sp in [("Einstein lam=2", einstein(2.0, 1500.0)),
                 ("Pb-like", models["Pb-like"])]:
    for mu in [0.0, 0.15]:
        t1 = tc(sp, mu, 10 * sp.omega_max, trunc_mult=30)
        t2 = tc(sp, mu, 10 * sp.omega_max, trunc_mult=60)
        print(f"{name} mu*={mu}: Tc(30)={t1:.6f}  Tc(60)={t2:.6f}  rel={abs(t1-t2)/t2:.1e}")
