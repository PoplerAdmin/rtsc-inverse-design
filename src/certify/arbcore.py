"""Rigorous (ball-arithmetic) primitives for the even-sector Eliashberg certificates.

These are the ball-arithmetic building blocks of the wE1 (universal a0) certificate
(`ball_interval`, `cumsums`, `even_block`, `rigorous_negdef`); the kernel-specific functions live
in `families.py`. Lemma numbers refer to the wE1 cell method as stated in
docs/theory/03_certificates.md.

Even sector (k_B = hbar = 1, T = 1 units, lambda(0) terms cancelled):
    M_nm = l(|n-m|) + l(n+m+1)                      (n != m)
    M_nn = l(2n+1) - (2n+1) - 2 sum_{k=1}^n l(k)
Every quantity that enters a verdict is an arb ball; floating point only *proposes*
(a shift s and a Cholesky factor), and the proposal is verified in arb (Lemma 6).
"""
import numpy as np
from flint import arb, arb_mat, ctx

ctx.prec = 128
PI = arb.pi()
TWO_PI = 2 * PI


def ball_interval(w1, w2):
    """Ball containing the closed interval [w1, w2] (w1, w2 exact floats)."""
    return arb(w1).union(arb(w2))


def cumsums(lv, nmax):
    """S[n] = sum_{k=1}^n lv[k], n = 0..nmax."""
    S = [arb(0)]
    for k in range(1, nmax + 1):
        S.append(S[-1] + lv[k])
    return S


def even_block(lv, S, N, deriv=False):
    """N x N leading block of M+ (deriv=True: of dM+/dw, where lv holds l', S its cumsums)."""
    rows = []
    for n in range(N):
        row = []
        for m in range(N):
            if n == m:
                v = lv[2 * n + 1] - 2 * S[n]
                if not deriv:
                    v = v - (2 * n + 1)
            else:
                v = lv[abs(n - m)] + lv[n + m + 1]
            row.append(v)
        rows.append(row)
    return arb_mat(rows)


def rigorous_negdef(S_ball):
    """Prove every symmetric matrix in the ball matrix S_ball is positive definite.
    Returns (ok, certified lower bound on lambda_min of the scaled matrix, float ev_min).
    Lemma 6: exact power-of-two diagonal scaling; float s > 0 and L with S - sI ~ L L^T;
    E := S - sI - L L^T in arb; if ||E||_inf < s then S >= (s - ||E||) I > 0."""
    N = S_ball.nrows()
    Sm = np.array([[float(S_ball[i, j].mid()) for j in range(N)] for i in range(N)])
    dg = np.diag(Sm)
    if np.any(dg <= 0):
        return False, None, float(np.min(dg))
    e = np.round(0.5 * np.log2(dg)).astype(int)
    sc = [arb(2.0 ** (-int(t))) for t in e]
    S_ball = arb_mat([[S_ball[i, j] * sc[i] * sc[j] for j in range(N)] for i in range(N)])
    Sm = np.array([[float(S_ball[i, j].mid()) for j in range(N)] for i in range(N)])
    Sm = 0.5 * (Sm + Sm.T)
    ev = np.linalg.eigvalsh(Sm)
    if ev[0] <= 0:
        return False, None, float(ev[0])
    s = 0.5 * ev[0]
    try:
        Lf = np.linalg.cholesky(Sm - s * np.eye(N))
    except np.linalg.LinAlgError:
        return False, None, float(ev[0])
    Lb = arb_mat([[arb(float(Lf[i, j])) for j in range(N)] for i in range(N)])
    sb = arb(s)
    E = S_ball - Lb * Lb.transpose()
    for i in range(N):
        E[i, i] = E[i, i] - sb
    worst = arb(0)
    for i in range(N):
        r = arb(0)
        for j in range(N):
            r += arb(E[i, j].abs_upper())
        worst = worst.max(r)
    marg = sb - worst
    ok = bool(marg > 0)
    return ok, (float(marg.lower()) if ok else None), float(ev[0])
