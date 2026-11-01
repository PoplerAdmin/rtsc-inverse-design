"""R01 spot checks of tabulated values with the independent solver (indep_solver.py)."""
import math
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from indep_solver import tc_indep, einstein  # noqa

print("1. Einstein mu*=0 optimum: Omega=7.57486, lam=1.14365 -> expect Tc=1, A/Tc=4.33148")
om, w = einstein(1.14365, 7.57486)
T = tc_indep(om, w)
print(f"   Tc={T:.8f}  A/Tc={w.sum()/T:.7f}")

print("2. Table 2: mu*=0.13 wc=20000, Einstein lam=1.69787 Omega=2133.65 -> expect 300 K")
om, w = einstein(1.69787, 2133.65)
print(f"   sym  Tc={tc_indep(om, w, 0.13, 20000.0, 'sym'):.3f}   asym Tc={tc_indep(om, w, 0.13, 20000.0, 'asym'):.3f}")
print("   Table 2: mu*=0.10 wc=5000, lam=1.47562 Omega=2265.46")
om, w = einstein(1.47562, 2265.46)
print(f"   sym  Tc={tc_indep(om, w, 0.10, 5000.0, 'sym'):.3f}   asym Tc={tc_indep(om, w, 0.10, 5000.0, 'asym'):.3f}")

print("3. Table 4: Einstein lam=3 Omega=1075.33 mu*=0 -> 300; lam=2 Omega=1736.63 mu*=0.10 wc=20000 -> 300")
om, w = einstein(3.0, 1075.33)
print(f"   Tc={tc_indep(om, w):.3f}")
om, w = einstein(2.0, 1736.63)
print(f"   sym Tc={tc_indep(om, w, 0.10, 20000.0):.3f}  asym Tc={tc_indep(om, w, 0.10, 20000.0, 'asym'):.3f}")

print("4. Table 6: two atoms lam=2 split equally r=128 mu*=0 -> Tc/wlog=1.297")
s = 100.0
om = np.array([s, 128 * s]); lam = np.array([1.0, 1.0]); w = lam * om / 2
T = tc_indep(om, w)
g = 2 * w / om
wlog = math.exp((g * np.log(om)).sum() / g.sum())
print(f"   Tc/wlog={T/wlog:.4f}  Tc/A={T/w.sum():.4f}  Tc/sqrt(lamw2)={T/math.sqrt(2*(w*om).sum()):.4f}")
for mu in (0.10, 0.13):
    T = tc_indep(om, w, mu, 10 * om.max())
    print(f"   mu*={mu} (wc=10 wmax): Tc/wlog={T/wlog:.4f}")

print("5. Einstein large lambda: Tc/(Omega sqrt(lam)) -> C_inf = 0.1827262")
for lam in (10.0, 100.0, 1000.0):
    om, w = einstein(lam, 1.0)
    T = tc_indep(om, w, Nmult=40.0)
    print(f"   lam={lam:g}: {T/math.sqrt(lam):.6f}")

print("6. mu* cutoff convention sensitivity (sym = src/theory convention, asym = theta on m only)")
for lam in (0.43, 1.0, 1.7, 3.0):
    om, w = einstein(lam, 1000.0)
    for mu in (0.10, 0.13):
        a = tc_indep(om, w, mu, 10000.0, "sym")
        b = tc_indep(om, w, mu, 10000.0, "asym")
        print(f"   Einstein lam={lam} mu*={mu} wc=10 Omega: sym={a:.3f} asym={b:.3f} asym/sym={b/a:.4f}")
