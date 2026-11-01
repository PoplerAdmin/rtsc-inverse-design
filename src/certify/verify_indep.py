"""INDEPENDENT second checker for the NM-G and N1mu certificates (review R03, finding 2).

Written from the statements in docs/theory/03_certificates.md Sec. 1-4, without importing any of
arbcore / cellcert / families / cover / run_nm / run_n1mu / verify. Only python-flint (arb,
arb_mat), numpy (float midpoints only, never in a verdict) and the JSON files are used.

Differences from the generator's checker (verify.py), on purpose:
  * positive definiteness is proved by an exact-interval LDL^T factorisation in Arb
    (every pivot D_j must have a positive lower bound), not by float Cholesky + residual;
  * truncation N = 110, explicit coupling columns K = 700 (generator: N = 120, K = 600);
  * tail sums T(j) are bounded only by the integral comparison (no closed form for L);
  * the far coupling (m >= N+K) is bounded by Cauchy-Schwarz with column sums <= 2N l(K+1);
  * cells that fail are bisected (up to 8 levels) - a refinement of a stored cell certifies the
    same statement; the number of refinements is reported;
  * end lemmas are re-derived: M-family small w uses L <= m0/24, large w L <= m0/(4w);
    A-family small w uses Lang(u) <= u/3 (no coth evaluation), large w the weighted Schur test.

Statement checked for each NM file (T = 1 units, b_k = 2 pi k, l(k) = m0/(w^2 + b_k^2)):
   (C)   Q+_w(x) + beta (lam - m0/w^2)|x|^2 <= -eta |x|^2   for all w > 0, all finite x,
   (C0)  beta*lam < 1 - eta,  beta >= 0  (beta = 0, no lam, for lambda = inf),
   Q+_w(x) = sum_{n!=m}[l(|n-m|)+l(n+m+1)] x_n x_m + sum_n [l(2n+1) - 2 S_n - (2n+1)] x_n^2.
For each N1mu segment (l(k) = 2 a w/(w^2+b_k^2), Coulomb -2 mu* (sum_{n<j} x_n)^2):
   Q+_w(x) - 2 mu* (sum_{n<j} x_n)^2 <= -eta |x|^2 for all w > 0,
plus the segment bookkeeping (a_j T_lo(j) >= A_floor; every j with T_lo(j) < A_floor/a0 present,
including j = 0).

usage: python src/certify/verify_indep.py [cert_dir] [--only nm|n1mu] [--files a,b]
"""
import glob
from decimal import Decimal
from fractions import Fraction
import json
import math
import os
import sys
import time

import numpy as np
from flint import arb, arb_mat, ctx

ctx.prec = 128
PI = arb.pi()
NN, KK = 110, 700
ETA = arb("1e-9")
MAXDEPTH = 8
A0_UNIV = arb("4.33148145")


# ---------------------------------------------------------------- kernels
class Kern:
    """kind 'M': l = m0/(w^2+b^2), shift beta(lam - m0/w^2);  kind 'A': l = 2aw/(w^2+b^2),
    Coulomb mu on the first j indices."""

    def __init__(self, kind, p, lam=None, beta=None, mu=None, j=0):
        self.kind, self.p = kind, arb(p)
        self.lam = arb(lam) if lam is not None else None
        self.beta = arb(beta) if beta is not None else arb(0)
        self.mu = arb(mu) if mu is not None else arb(0)
        self.j = j

    def l(self, k, w):
        b2 = (2 * PI * k) ** 2
        return self.p / (w * w + b2) if self.kind == "M" else 2 * self.p * w / (w * w + b2)

    def dl(self, k, w):
        b2 = (2 * PI * k) ** 2
        q = w * w + b2
        return -2 * self.p * w / (q * q) if self.kind == "M" else 2 * self.p * (b2 - w * w) / (q * q)

    def d2l_sup(self, k, wb):
        b2 = (2 * PI * k) ** 2
        q = wb * wb + b2
        if self.kind == "M":
            v = self.p * (6 * wb * wb - 2 * b2) / (q * q * q)
        else:
            v = -4 * self.p * wb * (3 * b2 - wb * wb) / (q * q * q)
        return arb(v.abs_upper())

    def tail(self, j, wb):
        """upper bound of T(j) = sum_{k>=j} l(k) over the ball wb (l decreasing in k):
        T(j) <= l(j) + int_j^inf l(t) dt."""
        lj = self.l(j, wb)
        if self.kind == "M":
            u = (lj + self.p / (4 * PI * PI * j)).upper()
            if wb.lower() > 0:
                u = min(u, (lj + self.p * (wb / (2 * PI * j)).atan() / (2 * PI * wb)).upper())
        else:
            u = (lj + (self.p / PI) * (wb / (2 * PI * j)).atan()).upper()
        return arb(u)

    def has_shift(self):
        return self.kind == "M" and self.lam is not None

    def s(self, w):
        return self.beta * (self.lam - self.p / (w * w)) if self.has_shift() else arb(0)

    def ds(self, w):
        return 2 * self.beta * self.p / (w * w * w) if self.has_shift() else arb(0)


def form_block(lv, N, deriv=False):
    """even-sector N x N matrix from lv[k] (k = 1..2N); deriv: drop the bare -(2n+1)."""
    S = [arb(0)]
    for k in range(1, N + 1):
        S.append(S[-1] + lv[k])
    rows = []
    for n in range(N):
        r = []
        for m in range(N):
            if n == m:
                v = lv[2 * n + 1] - 2 * S[n] - (0 if deriv else (2 * n + 1))
            else:
                v = lv[abs(n - m)] + lv[n + m + 1]
            r.append(v)
        rows.append(r)
    return rows


def ldl_posdef(S, N):
    """exact-interval LDL^T of the symmetric ball matrix S (list of lists). True iff every pivot
    has positive lower bound (then every symmetric matrix in the ball is positive definite)."""
    L = [[None] * N for _ in range(N)]
    D = [None] * N
    for jj in range(N):
        acc = S[jj][jj]
        Lj = L[jj]
        for k in range(jj):
            acc = acc - Lj[k] * Lj[k] * D[k]
        if not acc > 0:
            return False
        D[jj] = acc
        for i in range(jj + 1, N):
            acc2 = S[i][jj]
            Li = L[i]
            for k in range(jj):
                acc2 = acc2 - Li[k] * Lj[k] * D[k]
            Li[jj] = acc2 / D[jj]
    return True


def cell_ok(K, w1, w2, N=NN, KC=KK):
    wc, h = (w1 + w2) / 2.0, (w2 - w1) / 2.0
    assert wc - h == w1 and wc + h == w2
    wcb, hb = arb(wc), arb(h)
    wb = arb(w1).union(arb(w2))
    # Taylor data at the centre
    l0 = [None] + [K.l(k, wcb) for k in range(1, 2 * N + 1)]
    l1 = [None] + [K.dl(k, wcb) for k in range(1, 2 * N + 1)]
    l2 = [None] + [K.d2l_sup(k, wb) for k in range(1, 2 * N + 1)]
    P0 = form_block(l0, N)
    P1 = form_block(l1, N, deriv=True)
    R = form_block(l2, N, deriv=True)       # entries are sums of sup|l''|: entrywise bound
    s0, s1 = K.s(wcb), K.ds(wcb)
    ssup = arb(K.s(arb(w2)).upper()) if K.has_shift() else arb(0)   # s increasing in w
    # tail and coupling
    gam = []
    for m in range(N, N + KC):
        g = (2 * m + 1) - 2 * K.tail(m + 1, wb) - ssup - ETA
        if not g > 0:
            return False
        gam.append(arb(g.lower()))
    gmin = (2 * (N + KC) + 1) - 2 * K.tail(N + KC + 1, wb) - ssup - ETA
    if not gmin > 0:
        return False
    gmin = arb(gmin.lower())
    lt = [None] + [K.l(k, wb) for k in range(1, 2 * N + KC + 2)]
    B = arb_mat([[lt[N + c - n] + lt[n + N + c + 1] for c in range(KC)] for n in range(N)])
    Bg = arb_mat([[B[n, c] / gam[c] for c in range(KC)] for n in range(N)])
    G1 = Bg * B.transpose()
    col_far = 2 * N * K.l(KC + 1, wb)          # column sums for m >= N+KC (l decreasing)
    far = [arb((col_far / gmin * (K.tail(N + KC - n, wb) + K.tail(N + KC + n + 1, wb))).upper())
           for n in range(N)]
    rem = []
    for n in range(N):
        r = arb(0)
        for m in range(N):
            r += arb(R[n][m].abs_upper())
        rem.append(arb((hb * hb / 2 * r).upper()))
    for sg in (1, -1):
        S = [[None] * N for _ in range(N)]
        for n in range(N):
            for m in range(n + 1):
                v = P0[n][m] + sg * hb * P1[n][m] + G1[n, m]
                if K.kind == "A" and n < K.j and m < K.j:
                    v = v - 2 * K.mu
                if n == m:
                    v = v + s0 + sg * hb * s1 + rem[n] + far[n] + ETA
                S[n][m] = -v
                S[m][n] = -v
        if not ldl_posdef(S, N):
            return False
    return True


def check_cells(K, cells, lo, hi):
    cells = sorted(tuple(c) for c in cells)
    if not (cells[0][0] == lo and cells[-1][1] == hi and
            all(cells[i][1] == cells[i + 1][0] for i in range(len(cells) - 1))):
        return False, "cells do not tile [%r, %r]" % (lo, hi)
    refined = 0
    bad = []
    for (w1, w2) in cells:
        todo = [(w1, w2, 0)]
        while todo:
            a, b, dpt = todo.pop()
            if cell_ok(K, a, b):
                continue
            if dpt >= MAXDEPTH:
                bad.append((a, b))
                break
            mid = (a + b) / 2
            refined += 1
            todo += [(a, mid, dpt + 1), (mid, b, dpt + 1)]
    if bad:
        return False, "%d cells NOT closed (e.g. %s)" % (len(bad), bad[:2])
    return True, "%d cells closed (%d bisections)" % (len(cells), refined)


# ---------------------------------------------------------------- NM files
def nm_file(d):
    lam_s = d["lambda"]
    m0 = arb(d["m0_over_T2"])
    beta = arb(d["beta"])
    ws, W = d["w_range"]
    msgs, ok = [], True
    if not beta >= 0:
        return False, "beta < 0"
    if lam_s == "inf":
        if not beta == 0 or ws != 0.0:
            return False, "lambda=inf needs beta=0 and cells from w=0"
        lam = arb(0)
        K = Kern("M", d["m0_over_T2"])
    else:
        lam = arb(lam_s)
        K = Kern("M", d["m0_over_T2"], lam=lam_s, beta=d["beta"])
        # (C0): m = 0 endpoint
        if not (beta * lam < 1 - ETA):
            return False, "(C0) fails"
        # small w (0, ws]: row sums 2T(n+1)-(2n+1) <= 2L - 1, L <= m0 sum 1/(2 pi k)^2 = m0/24
        wsb = arb(ws)
        if not (ws > 0 and (m0 / 12 - 1 + beta * lam - beta * m0 / (wsb * wsb) < -ETA)):
            return False, "small-w end fails"
    # large w >= W: L = sum m0/(w^2+b^2) <= m0/(4w)  (proof: L = (m0/4w)(coth(w/2)-2/w), coth u<=1+1/u)
    if not (m0 / (2 * arb(W)) - 1 + beta * lam < -ETA):
        return False, "large-w end fails"
    okc, m = check_cells(K, d["cells"], ws, W)
    # exact certified g_lo
    if not Fraction(float(d["g_lo_K2"])) <= 90000 * Fraction(Decimal(d["m0_over_T2"])):
        return False, "g_lo_K2 > 90000 m0"
    return okc, "ends ok; " + m


# ---------------------------------------------------------------- N1mu files
def a_ends(a, s):
    e = s["ends"]
    ws, W, n0, H = e["w_s"], e["W"], int(e["n0"]), arb(e["H"])
    # small w: Q <= (2L-1)|x|^2, 2L = a Lang(w/2) <= a w/6 (Lang(u) <= u/3), increasing in w
    if not (a * arb(ws) / 6 - 1 < -ETA):
        return None
    eps = 2 * a / arb(W)
    r1 = eps - 1 + 2 * (n0 - 1) * eps + (a + (n0 - 1) * eps) / H
    r2 = a - (2 * n0 + 1) + 2 * (H - 1) * n0 * eps
    if not (r1 < -ETA and r2 < -ETA and H >= 1):
        return None
    return ws, W


def n1mu_file(d):
    A_floor = arb(d["A_floor_K"])
    wc = arb(d["omega_c_K"])
    T_u = A_floor / A0_UNIV
    T0 = arb(300)
    # N_c(300) recomputed: #{n >= 0 : (2n+1) pi 300 <= wc}
    j300 = 0
    while bool((2 * j300 + 1) * PI * T0 <= wc):
        j300 += 1
    if j300 != d["N_c_at_300"]:
        return False, "N_c(300) mismatch"
    segs = {s["j"]: s for s in d["segments"]}
    msgs, ok = [], True
    if not bool(wc / PI >= T_u):
        return False, "segment j=0 not covered"
    for j in range(j300, 0, -1):
        Tlo = T0 if j == j300 else wc / ((2 * j + 1) * PI)
        if j != j300 and not bool(Tlo >= T0):
            return False, "T_lo < 300 for j=%d" % j
        if bool(Tlo >= T_u):
            continue
        s = segs.get(j)
        if s is None or "a_cert" not in s:
            return False, "segment j=%d missing" % j
        a = arb(s["a_cert"])
        if not bool(a * Tlo >= A_floor):
            return False, "a_j T_lo < A_floor at j=%d" % j
        e = a_ends(a, s)
        if e is None:
            return False, "ends fail at j=%d" % j
        K = Kern("A", s["a_cert"], mu=d["mu_star"], j=j)
        okc, m = check_cells(K, s["cells"], e[0], e[1])
        ok &= okc
        msgs.append("j=%d %s" % (j, m))
    return ok, "; ".join(msgs)


def main():
    argv = sys.argv[1:]
    args = [a for i, a in enumerate(argv)
            if not a.startswith("--") and not (i > 0 and argv[i - 1] in ("--only", "--files"))]
    root = args[0] if args else os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                             "..", "..", "docs", "theory", "certificates")
    only = None
    files = None
    for i, a in enumerate(sys.argv):
        if a == "--only":
            only = sys.argv[i + 1]
        if a == "--files":
            files = sys.argv[i + 1].split(",")
    allok = True
    nfiles = 0
    for p in sorted(glob.glob(os.path.join(root, "*.json"))):
        if files and os.path.basename(p) not in files:
            continue
        d = json.load(open(p))
        t0 = time.time()
        if d.get("rule") == "NM-G" and only in (None, "nm"):
            ok, msg = nm_file(d)
        elif d.get("rule") == "N1mu" and only in (None, "n1mu"):
            ok, msg = n1mu_file(d)
        else:
            continue
        allok &= ok
        nfiles += 1
        print("%s %s  %s  (%.0fs)" % ("PASS" if ok else "FAIL", os.path.basename(p), msg,
                                      time.time() - t0), flush=True)
    allok &= nfiles > 0
    print("INDEPENDENT CHECK (%d files): %s" % (nfiles, "ALL VERIFIED" if allok else "NOT VERIFIED"),
          flush=True)
    sys.exit(0 if allok else 1)


if __name__ == "__main__":
    main()
