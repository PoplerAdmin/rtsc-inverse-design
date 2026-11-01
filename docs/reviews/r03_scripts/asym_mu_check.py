"""R03: does mu* >= 0 lower the top eigenvalue also in the ASYMMETRIC (Allen-Mitrovic, cutoff on the
summed index only) convention? The NM-G theorem claims 'every mu* >= 0, every cutoff' via (F1),
which is a quadratic-form statement for the symmetric block. Float check with the independent
solver: top_eig(mu, asym) - top_eig(0) over random spectra, mu*, omega_c, T."""
import math
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "r01_scripts"))
from indep_solver import top_eig  # noqa: E402

rng = np.random.default_rng(3)
worst = -9
for it in range(150):
    k = rng.integers(1, 4)
    om = np.exp(rng.uniform(np.log(50), np.log(5000), k))
    lam = rng.uniform(0.3, 5.0) * rng.dirichlet(np.ones(k))
    w = lam * om / 2
    T = float(np.exp(rng.uniform(np.log(30), np.log(600))))
    mu = float(rng.uniform(0.02, 0.4))
    wc = float(np.exp(rng.uniform(np.log(300), np.log(50000))))
    e0 = top_eig(om, w, T, 0.0, Nmult=20.0)
    ea = top_eig(om, w, T, mu, wc, "asym", Nmult=20.0)
    es = top_eig(om, w, T, mu, wc, "sym", Nmult=20.0)
    worst = max(worst, ea - e0)
    if ea > e0 + 1e-12 or es > e0 + 1e-12:
        print("RAISE it=%d e0=%.6f asym=%.6f sym=%.6f" % (it, e0, ea, es), flush=True)
print("max (asym - mu0) top_ev difference over 150 cases: %.3e" % worst)
