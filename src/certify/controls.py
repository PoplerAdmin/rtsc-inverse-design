"""Sanity controls for the certificate code (not part of the proof).
(1) cells must FAIL when m0 is above the Einstein value, (2) at lambda = inf above g(inf),
(3) at finite lambda without the Lagrange shift; and succeed just below.
(4) Arb even block == independent float build_M of src/theory/eliashberg.py.
(5) random 1-4 atom spectra with lambda <= 3 and m <= m0(3): float rho < 0.
usage: python src/certify/controls.py"""
import math
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(HERE, "..", "theory"))
from flint import arb  # noqa: E402
from families import MFamily  # noqa: E402
from cellcert import certify_cell  # noqa: E402
from arbcore import cumsums, even_block  # noqa: E402
import run_nm  # noqa: E402
from eliashberg import einstein, rho, build_M, Spectrum  # noqa: E402

eta = arb("1e-7")
lam = 3.0
wE = run_nm.float_wE(lam)
mE = lam * wE ** 2
ws_grid = np.exp(np.linspace(np.log(0.05 * math.sqrt(mE / lam)), np.log(mE), 1500))
b, _ = run_nm.propose_beta(lam, 0.9999 * mE, ws_grid)
b = math.floor(b * 2 ** 30) / 2 ** 30
wc = math.floor(wE * 2 ** 12) / 2 ** 12
res = []
fam = MFamily(run_nm.dec_floor(1.0001 * mE), "3.0", repr(b))
res.append(("(1) m0=1.0001 g(3), cell at w_E", certify_cell(fam, eta, wc, wc + 2 ** -12, 80, 400)[0], False))
fam = MFamily(run_nm.dec_floor(0.9999 * mE), "3.0", repr(b))
res.append(("(1b) m0=0.9999 g(3), same cell", certify_cell(fam, eta, wc, wc + 2 ** -12, 80, 400)[0], True))
minf = 1 / run_nm.float_c2() ** 2
fam = MFamily(run_nm.dec_floor(1.0001 * minf), None, "0")
res.append(("(2) inf, m0=1.0001 g(inf), [0,2^-6]", certify_cell(fam, eta, 0.0, 2 ** -6, 120, 600)[0], False))
fam = MFamily(run_nm.dec_floor(0.9999 * minf), None, "0")
res.append(("(2b) inf, m0=0.9999 g(inf), [0,2^-6]", certify_cell(fam, eta, 0.0, 2 ** -6, 120, 600)[0], True))
fam = MFamily(run_nm.dec_floor(0.9999 * mE), None, "0")
res.append(("(3) m0=0.9999 g(3), beta=0, [0,2^-6]", certify_cell(fam, eta, 0.0, 2 ** -6, 80, 400)[0], False))
allok = True
for name, got, exp in res:
    print("%-42s got %-5s expected %-5s %s" % (name, got, exp, "ok" if got == exp else "MISMATCH"))
    allok &= got == exp
m0 = 0.9999 * mE
for w in [0.7, 3.3, 9.0]:
    N = 60
    f = MFamily(repr(m0), None, "0")
    lv = f.lvals(arb(w), 2 * N)
    A = even_block(lv, cumsums(lv, N), N)
    Af = np.array([[float(A[i, j].mid()) for j in range(N)] for i in range(N)])
    d = np.abs(Af - build_M(einstein(m0 / w ** 2, w), 1.0, N=N)).max()
    print("(4) w=%g max|arb - build_M| = %.2e" % (w, d))
    allok &= d < 1e-10
rng = np.random.default_rng(1)
worst = -9.0
for it in range(400):
    k = rng.integers(1, 5)
    om = np.exp(rng.uniform(np.log(0.05), np.log(60), k))
    lw = rng.dirichlet(np.ones(k)) * lam * rng.uniform(0.7, 1.0)
    sp = Spectrum(om, lw * om / 2)
    if sp.lam_omega2 > m0:
        sp = sp.scaled(lam_factor=m0 / sp.lam_omega2)
    worst = max(worst, rho(sp, 1.0))
print("(5) max float rho, 400 random spectra (lambda<=3, m<=m0): %.4f" % worst)
allok &= worst < 0
print("CONTROLS %s" % ("OK" if allok else "FAILED"))
