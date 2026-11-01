"""R03: Matsubara-truncation convergence of the independent solver at large lambda, where the
adversarial margins (top_ev - 1 ~ -1e-6) are small. Einstein and soft+hard spectra at
m = frac * g_lo(lambda_k) (T = 300 K), lambda_k in {100, inf}; Nmult controls the truncation."""
import json
import math
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "r01_scripts"))
from indep_solver import top_eig  # noqa: E402

CERT = os.path.join(HERE, "..", "..", "theory", "certificates")
m0 = {k: float(json.load(open(os.path.join(CERT, "nm_%s.json" % k)))["m0_over_T2"]) for k in ("100.0", "inf")}
T = 300.0
for key, lam, frac in (("100.0", 100.0, 0.9999), ("100.0", 100.0, 0.99999), ("100.0", 60.0, 0.99999),
                       ("inf", 1000.0, 0.99999), ("inf", 1e4, 0.99999), ("inf", 1e4, 1.0001)):
    m = frac * m0[key] * T ** 2
    Om = math.sqrt(m / lam)
    for Nwant in (80, 300, 700):
        Nm = (Nwant - 2) * 2 * math.pi * T / Om
        ev = top_eig(np.array([Om]), np.array([lam * Om / 2]), T, Nmult=Nm)
        print("lam_k=%s Einstein lam=%g m=%.5f g_lo  Nmult=%g (N=%d)  top_ev-1 = %+.3e"
              % (key, lam, frac, Nm, max(80, int(Nm * Om / (2 * math.pi * T)) + 2), ev - 1), flush=True)
# two-atom at lambda = 100: soft (80 at Om1) + hard (20 at 20*Om1), m = 0.99999 g_lo(100)
mt = 0.99999 * m0["100.0"] * T ** 2
for r in (5.0, 20.0, 60.0):
    lams = np.array([80.0, 20.0])
    oms = np.array([1.0, r])
    s = math.sqrt(mt / float(np.sum(lams * oms ** 2)))
    oms = oms * s
    ev = top_eig(oms, lams * oms / 2, T, Nmult=600.0)
    print("lam_k=100 2atom 80@%.1fK + 20@%.1fK m=0.99999 g_lo top_ev-1 = %+.3e" % (oms[0], oms[1], ev - 1), flush=True)
