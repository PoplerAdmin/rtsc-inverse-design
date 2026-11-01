"""R01: try to break N1mu, N2, N3 with continuous optimisation over K-atom spectra (K = 2..4).

Uses the reviewed solver's rho() (verified against indep_solver.py in spot_checks.py) for speed.
For a fixed *shape* the overall amplitude t is scanned upward on a geometric grid and the first
instability at T = 300 K is bracketed and refined (brentq) -> minimal amplitude for that shape.

  N3 : atoms at s*exp(u_i), lambda shares p_i, total lambda = LAM fixed; amplitude = s.
       objective: omega_2 = sqrt(lambda<w^2>/lambda).  Claimed bound: Omega_E(LAM, mu*).
  N2 : atoms at absolute Omega_i, shares p_i of c = int omega a2F; amplitude = c.
       objective: lambda<w^2> = 2c. Claimed bound: 2.6955e6 (mu*=0), 3.634e6 (mu*=0.13, wc=2e4).
  N1mu: atoms at absolute Omega_i, shares p_i of A; amplitude = A.
       objective: A.  Claimed bound: 1811.33 (mu*=0.13, wc=2e4), 1299.44 (mu*=0).
Usage: python adversarial_opt.py <src/theory dir> <outfile>
"""
import math
import sys
import time

import numpy as np
from scipy.optimize import brentq, differential_evolution

sys.path.insert(0, sys.argv[1])
from eliashberg import Spectrum, rho  # noqa

OUT = open(sys.argv[2], "w")
T0 = 300.0
TM = 12.0  # truncation multiple for the search; best points re-checked with 40


def log(*a):
    print(*a, file=OUT, flush=True)


def first_instability(make, lo, hi, mu, wc, tm=TM, nsteps=40):
    grid = np.geomspace(lo, hi, nsteps)
    prev = None
    for t in grid:
        if rho(make(t), T0, mu, wc, trunc_mult=tm) >= 0:
            if prev is None:
                return lo
            return brentq(lambda x: rho(make(x), T0, mu, wc, trunc_mult=tm), prev, t,
                          xtol=1e-7 * t)
        prev = t
    return math.inf


def shares(v):
    e = np.exp(v - v.max())
    return e / e.sum()


def n3_objective(x, K, LAM, mu, wc, tm=TM):
    u = np.concatenate([[0.0], x[:K - 1]])
    p = shares(np.concatenate([[0.0], x[K - 1:]]))
    lam_i = LAM * p

    def make(s):
        om = s * np.exp(u)
        return Spectrum(om, lam_i * om / 2)
    s = first_instability(make, 50.0 / np.exp(u).max(), 2e4 / np.exp(u).max() * 4, mu, wc, tm)
    if not math.isfinite(s):
        return 1e9
    return make(s).omega2


def n2_objective(x, K, mu, wc, tm=TM):
    om = np.exp(x[:K])
    p = shares(np.concatenate([[0.0], x[K:]]))

    def make(c):
        return Spectrum(om, c * p / om)
    c = first_instability(make, 2e5, 2e7, mu, wc, tm)
    return 2 * c if math.isfinite(c) else 1e12


def n1_objective(x, K, mu, wc, tm=TM):
    om = np.exp(x[:K])
    p = shares(np.concatenate([[0.0], x[K:]]))

    def make(A):
        return Spectrum(om, A * p)
    A = first_instability(make, 800.0, 5000.0, mu, wc, tm)
    return A if math.isfinite(A) else 1e9


def run(name, fun, bounds, args, claimed, recheck):
    t0 = time.time()
    res = differential_evolution(fun, bounds, args=args, seed=1, popsize=8, maxiter=40,
                                 tol=1e-6, polish=True, workers=1)
    val = recheck(res.x)
    log(f"{name}: best={res.fun:.6g} recheck(tm=40)={val:.6g} claimed_bound={claimed:.6g} "
        f"ratio={val/claimed:.5f} x={np.round(res.x, 4).tolist()} ({time.time()-t0:.0f}s)")


if __name__ == "__main__":
    INF = math.inf
    # N3 at lambda = 3 and 2, mu* = 0 and 0.10 (wc = 20000)
    for LAM, mu, wc, claimed in [(3.0, 0.0, INF, 1075.33), (2.0, 0.0, INF, 1416.82),
                                 (3.0, 0.10, 20000.0, 1273.31)]:
        for K in (2, 3, 4):
            b = [(-4.0, 4.0)] * (K - 1) + [(-5.0, 5.0)] * (K - 1)
            run(f"N3 lam={LAM} mu*={mu} K={K}", n3_objective, b, (K, LAM, mu, wc), claimed,
                lambda x, K=K, LAM=LAM, mu=mu, wc=wc: n3_objective(x, K, LAM, mu, wc, 40.0))
    # N2
    for mu, wc, claimed in [(0.0, INF, 2.6955e6), (0.13, 20000.0, 3.634e6)]:
        for K in (2, 3):
            b = [(math.log(20.0), math.log(15000.0))] * K + [(-5.0, 5.0)] * (K - 1)
            run(f"N2 mu*={mu} K={K}", n2_objective, b, (K, mu, wc), claimed,
                lambda x, K=K, mu=mu, wc=wc: n2_objective(x, K, mu, wc, 40.0))
    # N1mu
    for mu, wc, claimed in [(0.13, 20000.0, 1811.33), (0.13, 5000.0, 1781.61)]:
        for K in (2, 3):
            b = [(math.log(200.0), math.log(15000.0))] * K + [(-5.0, 5.0)] * (K - 1)
            run(f"N1mu mu*={mu} wc={wc} K={K}", n1_objective, b, (K, mu, wc), claimed,
                lambda x, K=K, mu=mu, wc=wc: n1_objective(x, K, mu, wc, 40.0))
    log("done")
