"""R01 review: independent linearized isotropic ME Tc solver (no code shared with src/theory).

Differences from src/theory/eliashberg.py (deliberate, for an independent check):
  * full Matsubara axis n in [-N, N-1] (no even/odd reduction, no lambda(0) cancellation done
    by hand: Z contains lambda(0) and the kernel contains lambda(0) on the diagonal);
  * Tc from the non-symmetric eigenproblem  Delta = K diag(1/|w~|) Delta  (largest eigenvalue 1);
  * lambda(k) evaluated in a different (vectorised over atoms, k as float) way;
  * mu* cutoff convention selectable:
       'sym'  : -mu* theta(wc-|w_n|) theta(wc-|w_m|)  (same as src/theory, block form)
       'asym' : -mu* theta(wc-|w_m|)                    (cutoff on the summed index only,
                                                          the textbook Allen-Mitrovic form)
Units k_B = hbar = 1 (kelvin).
"""
import math

import numpy as np
from scipy.optimize import brentq


def lam_of_k(om, w, T, kmax):
    k = np.arange(kmax + 1, dtype=float)
    nu2 = (2.0 * math.pi * T * k) ** 2
    return np.array([np.sum(2.0 * w * om / (om ** 2 + x)) for x in nu2])


def top_eig(om, w, T, mu=0.0, wc=math.inf, conv="sym", Nmult=40.0):
    om = np.asarray(om, float)
    w = np.asarray(w, float)
    top = Nmult * om.max()
    if math.isfinite(wc):
        top = max(top, 2.0 * wc)
    N = max(80, int(top / (2 * math.pi * T)) + 2)
    n = np.arange(-N, N)
    L = lam_of_k(om, w, T, 2 * N + 2)
    sgn = np.where(n >= 0, 1.0, -1.0)
    # |w~_n| / (pi T) = |2n+1| + lambda(0) + 2 sum_{k=1}^{|n|'} lambda(k), exact (no truncation)
    nn = np.where(n >= 0, n, -n - 1)
    S = np.concatenate([[0.0], np.cumsum(L[1:])])
    wt = np.abs(2 * n + 1) + L[0] + 2 * S[nn]
    D = np.abs(n[:, None] - n[None, :])
    K = L[D].copy()
    if mu > 0:
        th = (np.abs(2 * n + 1) * math.pi * T <= wc).astype(float)
        if conv == "sym":
            K -= mu * th[:, None] * th[None, :]
        else:
            K -= mu * th[None, :] * np.ones_like(th)[:, None]
    # symmetrise by similarity for 'sym' is possible but we use the general eigensolver
    Mx = K / wt[None, :]
    ev = np.linalg.eigvals(Mx)
    return float(np.max(ev.real))


def tc_indep(om, w, mu=0.0, wc=math.inf, conv="sym", T_hi=None, Nmult=40.0, xtol=1e-9):
    A = float(np.sum(w))
    if T_hi is None:
        T_hi = A / 4.33 * 1.01
    f = lambda T: top_eig(om, w, T, mu, wc, conv, Nmult) - 1.0
    T2 = T_hi
    T1 = 0.9 * T2
    while f(T1) < 0:
        T2 = T1
        T1 *= 0.9
        if T1 < 1e-4 * T_hi:
            return 0.0
    return brentq(f, T1, T2, xtol=xtol * T1)


def einstein(lam, Om):
    return np.array([Om]), np.array([lam * Om / 2.0])
