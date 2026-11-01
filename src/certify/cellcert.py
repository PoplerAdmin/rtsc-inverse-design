"""Certify one w-cell [w1, w2]:  x^T (M+(w) + s(w) I + C) x <= -eta |x|^2  for every w in the cell and
every finitely supported x, where M+(w) is the even-sector matrix of a kernel family
(families.py), s(w) a concave diagonal shift (or 0) and C a constant negative-semidefinite
block (Coulomb, or None).

This is the wE1 cell method, generalised: the kernel is supplied by a family object, and a
concave diagonal shift s(w) and a constant block C are allowed. Proof (wE1 Lemmas 2-6, with the
shift/Coulomb additions) is in docs/theory/03_certificates.md, Sec. 1.
"""
import time
from flint import arb, arb_mat
from arbcore import ball_interval, cumsums, even_block, rigorous_negdef


def certify_cell(fam, eta, w1, w2, N=120, K=600):
    t0 = time.time()
    wc = (w1 + w2) / 2.0
    h = (w2 - w1) / 2.0
    assert wc - h == w1 and wc + h == w2, "cell endpoints must be dyadic"
    wcb = arb(wc)
    hb = arb(h)
    wb = ball_interval(w1, w2)

    # --- Taylor data of the N x N block at the centre (Lemma 5')
    lc = fam.lvals(wcb, 2 * N)
    Sc = cumsums(lc, N)
    dlc = fam.dlvals(wcb, 2 * N)
    dSc = cumsums(dlc, N)
    A0 = even_block(lc, Sc, N)
    A1 = even_block(dlc, dSc, N, deriv=True)
    s0 = fam.shift(wcb)
    s1 = fam.dshift(wcb)
    for n in range(N):
        A0[n, n] = A0[n, n] + s0
        A1[n, n] = A1[n, n] + s1
    C = fam.const_block(N)
    if C is not None:
        A0 = A0 + C
    mu = fam.d2l_abs_upper(wb, 2 * N)
    Smu = cumsums(mu, N)
    half_h2 = hb * hb / 2
    rem = []
    for n in range(N):
        u = mu[2 * n + 1] + 2 * Smu[n]
        for m in range(N):
            if m != n:
                u += mu[abs(n - m)] + mu[n + m + 1]
        rem.append(arb((half_h2 * u).upper()))

    # --- tail and coupling, valid for every w in the cell (Lemmas 1-5)
    kmax = 2 * N + K + 2
    lb = fam.lvals(wb, kmax)
    Sb = cumsums(lb, kmax)
    Lb = fam.L_total(wb)

    def Tu(j):
        if j < len(lb):
            lj = lb[j]
        else:
            lj = fam.lvals(wb, j)[j]
        best = fam.T_integral(wb, j, lj)
        if Lb is not None and j - 1 < len(Sb):
            best = arb(min(best.upper(), (Lb - Sb[j - 1]).upper()))
        return best

    ssup = fam.shift_sup(w1, w2)
    gam = []
    for m in range(N, N + K):
        g = arb((2 * m + 1 - 2 * Tu(m + 1) - ssup - eta).lower())
        if not g > 0:
            return False, {"why": "gamma<=0", "time": time.time() - t0}
        gam.append(g)
    TNK1 = Tu(N + K + 1)
    gmin = arb((2 * (N + K) + 1 - 2 * TNK1 - ssup - eta).lower())
    if not gmin > 0:
        return False, {"why": "gmin<=0", "time": time.time() - t0}
    Brows, Bgrows = [], []
    for n in range(N):
        r, rg = [], []
        for j in range(K):
            m = N + j
            v = lb[m - n] + lb[n + m + 1]
            r.append(v)
            rg.append(v / gam[j])
        Brows.append(r)
        Bgrows.append(rg)
    G1 = arb_mat(Bgrows) * arb_mat(Brows).transpose()
    cmax = Tu(K + 1) + TNK1
    dfac = cmax / gmin
    d = [arb((dfac * (Tu(N + K - n) + Tu(N + K + n + 1))).upper()) for n in range(N)]

    info = {"N": N, "K": K, "w1": w1, "w2": w2, "rem0": float(rem[0]), "d0": float(d[0])}
    worst = None
    for sgn in (+1, -1):
        P = A0 + (hb * sgn) * A1 + G1
        for n in range(N):
            P[n, n] = P[n, n] + rem[n] + d[n] + eta
        ok, marg, ev0 = rigorous_negdef(-P)
        info["ev_min_float_%+d" % sgn] = ev0
        if not ok:
            info["time"] = time.time() - t0
            return False, info
        worst = marg if worst is None else min(worst, marg)
    info["margin"] = worst
    info["time"] = time.time() - t0
    return True, info
