"""R01: is "Tc >= 300 K  =>  A >= 300 * a_E(300 K)" valid when Tc is defined as sup{T unstable}
with an ABSOLUTE block cutoff omega_c?  N_c(T) drops by one at T_j = omega_c / ((2j+1) pi), so the
Coulomb block shrinks discontinuously as T rises.  Compute A_min(T) = T * a_E(T; mu*, omega_c) for
T in [300, 420] K (Einstein family is the extremal one by the doc's own reduction) and report
min_{T >= 300} A_min(T) vs A_min(300).  Also tests a concrete spectrum for re-entrance.
Usage: python reentrance_check.py <src/theory dir> <outfile>
"""
import math
import sys

import numpy as np
from scipy.optimize import brentq, minimize_scalar

sys.path.insert(0, sys.argv[1])
from eliashberg import einstein, rho  # noqa

OUT = open(sys.argv[2], "w")


def log(*a):
    print(*a, file=OUT, flush=True)


def lam_c(T, Om, mu, wc):
    f = lambda lam: rho(einstein(lam, Om), T, mu, wc)
    hi = 2.0
    while f(hi) < 0:
        hi *= 1.6
    return brentq(f, 0.05, hi, xtol=1e-8)


def Amin(T, mu, wc):
    g = lambda w: lam_c(T, T * w, mu, wc) * T * w / 2
    r = minimize_scalar(g, bounds=(4.0, 12.0), method="bounded", options=dict(xatol=1e-5))
    return r.fun, r.x


for mu in (0.10, 0.13):
    for wc in (5000.0, 10000.0, 20000.0):
        A300, w300 = Amin(300.0, mu, wc)
        steps = [wc / ((2 * j + 1) * math.pi) for j in range(0, 40)]
        steps = [t for t in steps if 300 < t < 420]
        Ts = sorted(set(list(np.linspace(300, 420, 25)) + [t * (1 + 1e-6) for t in steps]
                        + [t * (1 - 1e-6) for t in steps]))
        best = min((Amin(T, mu, wc)[0], T) for T in Ts)
        log(f"mu*={mu} wc={wc}: A_min(300)={A300:.2f}  min_(T in [300,420]) A_min(T)={best[0]:.2f} "
            f"at T={best[1]:.2f}  N_c steps in range at T={[round(t,2) for t in steps]}")
        if best[0] < A300 - 1e-6:
            # concrete reentrant spectrum: Einstein optimal at T*, amplitude between best and A300
            Ts_ = best[1]
            _, w = Amin(Ts_, mu, wc)
            Om = Ts_ * w
            A = 0.5 * (best[0] + A300)
            lam = 2 * A / Om
            sp = einstein(lam, Om)
            log(f"   example Einstein Omega={Om:.1f} K lam={lam:.5f} A={A:.2f}: "
                f"rho(300)={rho(sp, 300.0, mu, wc):+.3e}  rho({Ts_:.2f})={rho(sp, Ts_, mu, wc):+.3e}")
            scan = [(T, rho(sp, T, mu, wc)) for T in np.linspace(280, Ts_ + 5, 60)]
            uns = [round(T, 1) for T, r in scan if r >= 0]
            log(f"   unstable T on grid 280..{Ts_+5:.0f}: {uns}")
log("done")
