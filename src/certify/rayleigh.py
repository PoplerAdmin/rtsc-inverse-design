"""Instability (upper-side) certificates: a dyadic test vector x with Q+(x) > 0, evaluated in arb.
Float proposes x, arb evaluates Q+(x) for the exact parameters.
Q+(x) uses only finitely many exact matrix entries (x has finite support), so Q+(x) > 0 proves
rho+ > 0, i.e. the spectrum is unstable at T (Tc_model >= T)."""
import numpy as np
from flint import arb, arb_mat
from arbcore import TWO_PI, cumsums, even_block


def float_block(lk, N, mu=0.0, Nc=0):
    """float even block from l(k) values lk[0..2N] (lk[0] unused)."""
    n = np.arange(N)
    M = lk[np.abs(n[:, None] - n[None, :])] + lk[n[:, None] + n[None, :] + 1]
    S = np.concatenate([[0.0], np.cumsum(lk[1:N])])
    M[n, n] = lk[2 * n + 1] - (2 * n + 1) - 2 * S
    if Nc:
        M[:Nc, :Nc] -= 2 * mu
    return M


def certify_unstable(l_arb, l_float, N=200, mu=None, Nc=0):
    """l_arb(k) -> arb ball of l(k) (exact parameters), l_float(k) -> float.
    Returns (ok, Q/|x|^2 lower bound)."""
    lk = np.array([0.0] + [l_float(k) for k in range(1, 2 * N + 1)])
    M = float_block(lk, N, 0.0 if mu is None else float(mu), Nc)
    ev, V = np.linalg.eigh(M)
    x = V[:, -1]
    x = np.round(x * 2 ** 40) / 2 ** 40
    lv = [None] + [l_arb(k) for k in range(1, 2 * N + 1)]
    S = cumsums(lv, N)
    A = even_block(lv, S, N)
    if Nc:
        mub = arb(mu)
        for i in range(Nc):
            for j in range(Nc):
                A[i, j] = A[i, j] - 2 * mub
    X = arb_mat([[arb(float(t))] for t in x])
    Q = (X.transpose() * A * X)[0, 0]
    nx2 = sum(arb(float(t)) ** 2 for t in x)
    r = Q / nx2
    return bool(Q > 0), float(r.lower()), float(ev[-1])


def einstein_unstable(lam, w, N=200):
    """Einstein spectrum with coupling lam (decimal string) at w = Omega/T (dyadic float), mu* = 0."""
    lb = arb(lam)
    wb = arb(w)
    la = lambda k: lb * wb * wb / (wb * wb + (TWO_PI * k) ** 2)
    lf = lambda k: float(lam) * w * w / (w * w + (2 * np.pi * k) ** 2)
    return certify_unstable(la, lf, N)
