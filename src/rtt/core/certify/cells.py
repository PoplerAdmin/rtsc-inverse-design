"""R-CELLS: rigorous stability on temperature cells [T1, T2] (even sector, mu* = 0).

Certificate of one cell (all numbers exact doubles / ints):
    {"T1", "T2"           : cell endpoints (exact doubles, stored as float.hex),
     "variant"            : "mhat" (ARCH sec. 4 / Lean CellOK, default) or "cform" (PROOFS.md C2),
     "N"                  : head length,  "u": head entries u_0..u_{N-1} (> 0),
     "C"                  : tail coefficient, v_m = C (2m+1)^-2 for m >= N,
     "eta"                : margin (> 0),
     "N_e"                : rows n < N_e checked one by one, rows n >= N_e by Lemma T1,
     "M_v"                : columns m < M_v summed explicitly, m >= M_v by Lemma R1}

check_cell() re-verifies a certificate from scratch in Arb.  Theorem (PROOFS.md, Thm C):
if check_cell(...) returns True then for every spectrum in the envelope e of the nominal
spectrum and every T in [T1, T2], every row of the (variant's) upper envelope satisfies
sum_m Mhat_nm v_m <= -eta v_n, hence (Lean cell_stable) Stable(T) for all mu* >= 0, all cutoffs.
"""
from __future__ import annotations

import math
import time

import numpy as np
import scipy.linalg as sla
from flint import arb, arb_poly

from .arbcore import ArbSpectrum, PI, poly_coeffs, tail_f
from .floatcore import FloatSpectrum, build_matrix, top_eig

VARIANTS = ("mhat", "cform")


# ============================================================================ Arb check
class LamCache:
    """Cache of nominal lambda(k;T) Arb lists per exact T (shared between adjacent cells)."""

    def __init__(self, asp: ArbSpectrum):
        self.asp = asp
        self.store: dict[float, list] = {}

    def get(self, T: float, kmax: int) -> list:
        cur = self.store.get(T)
        if cur is None or len(cur) <= kmax:
            cur = self.asp.lam(T, kmax)
            self.store[T] = cur
        return cur[: kmax + 1]


def _fail(reason, **kw):
    d = {"ok": False, "reason": reason}
    d.update(kw)
    return d


def check_cell(asp: ArbSpectrum, cert: dict, cache: LamCache | None = None) -> dict:
    """Rigorous check of one cell certificate.  Returns {"ok": bool, ...}."""
    t0 = time.time()
    try:
        T1 = float.fromhex(cert["T1"])
        T2 = float.fromhex(cert["T2"])
        variant = cert["variant"]
        N = int(cert["N"])
        N_e = int(cert["N_e"])
        M_v = int(cert["M_v"])
        u = [float(x) for x in cert["u"]]
        C = float(cert["C"])
        eta = float(cert["eta"])
        p = int(cert.get("p", 3))
    except Exception as exc:  # malformed certificate
        return _fail(f"malformed certificate: {exc!r}")
    # ---- structural conditions (exact comparisons of doubles / ints)
    if variant not in VARIANTS:
        return _fail("unknown variant")
    if not (T1 > 0 and T1 <= T2 and math.isfinite(T2)):
        return _fail("need 0 < T1 <= T2 < inf")
    if not (N >= 1 and len(u) == N and N_e >= 2 * N and N_e >= 2 and M_v >= N_e + 2):
        return _fail("need N >= 1, len(u) = N, N_e >= max(2N,2), M_v >= N_e + 2")
    if not all(math.isfinite(x) and x > 0 for x in u):
        return _fail("head entries must be finite and > 0")
    if not (math.isfinite(C) and C > 0 and math.isfinite(eta) and eta > 0):
        return _fail("need C > 0 and eta > 0")
    if p not in (2, 3):
        return _fail("tail exponent p must be 2 or 3")
    if cache is None:
        cache = LamCache(asp)
    one = arb(1)
    e1 = one + asp.e          # 1 + e  (exact)
    e2 = one - asp.e          # 1 - e
    K = N_e + M_v
    lam1 = cache.get(T1, K)
    hi = [e1 * x for x in lam1]                      # hi(k) = (1+e) lambda(k; T1), hi(0) = 0
    hi[0] = arb(0)
    lam2 = cache.get(T2, N_e)
    lo = [e2 * x for x in lam2]                      # lo(k) = (1-e) lambda(k; T2)
    lo[0] = arb(0)
    b = asp.b_const(T1, e1)                          # hi(k) <= lambda-bound b/k^2  (Lemma T0)
    Ca = arb(C)
    eta_a = arb(eta)
    # ---- test vector on 0..M_v-1 (exact head, ball tail)
    v = [arb(x) for x in u] + [Ca * tail_f(m, p) for m in range(N, M_v)]
    Lp = arb_poly(hi)
    P1 = poly_coeffs(Lp * arb_poly(v), N_e)
    P2 = poly_coeffs(Lp * arb_poly(v[::-1]), M_v + N_e)
    # remainder of columns m >= M_v  (Lemma R1, hi decreasing in k):
    #   sum_{m>=M_v} [hi(m-n) + hi(n+m+1)] C f(m) <= 2 C f(M_v) Tsuf(M_v - n),
    #   Tsuf(J) = sum_{k=J}^{K} hi(k) + b/K  >=  sum_{k>=J} hi(k).
    Tsuf = [arb(0)] * (K + 2)
    Tsuf[K + 1] = b / K
    for k in range(K, 0, -1):
        Tsuf[k] = Tsuf[k + 1] + hi[k]
    rem_c = 2 * Ca * tail_f(M_v, p)
    # Lemma R2: hi(j) <= b/j^2 and (2m+1)^-p <= (2(m-n))^-p give, with J = M_v - n,
    #   remainder <= 2 b C 2^-p [J^-(p+2) + J^-(p+1)/(p+1)]
    two_p = arb(2) ** p
    # ---- correction term for the c-form (Lemma C2)
    corr = None
    if variant == "cform":
        nu = [arb(0)] * (N + 1)
        nz = False
        for i in range(1, N + 1):
            d = v[i - 1] + v[i + 1] - 2 * v[i]
            if not bool(d >= 0):
                nu[i] = abs(d)
                nz = True
        if nz:
            # Corr_n = sum_i nu_i W_n(|i-n|),  W_n(d) = sum_{k=d+1}^{n} (k-d) w(k),  w = hi - lo
            #        = (SW1[d] - SW1[n]) - d (SW0[d] - SW0[n])  for d < n, and 0 for d >= n,
            # with suffix sums SW0[d] = sum_{d<k<N_e} w(k), SW1[d] = sum_{d<k<N_e} k w(k).
            w = [arb(0)] + [hi[k] - lo[k] for k in range(1, N_e)]   # k = 1..N_e-1
            SW0 = [arb(0)] * (N_e + 1)
            SW1 = [arb(0)] * (N_e + 1)
            for d in range(N_e - 2, -1, -1):
                SW0[d] = SW0[d + 1] + w[d + 1]
                SW1[d] = SW1[d + 1] + (d + 1) * w[d + 1]
            idx = [i for i in range(1, N + 1) if not nu[i].is_zero()]
            corr = []
            for n in range(N_e):
                s = arb(0)
                for i in idx:
                    d = abs(i - n)
                    if d < n:
                        s += nu[i] * ((SW1[d] - SW1[n]) - d * (SW0[d] - SW0[n]))
                corr.append(s)
    # ---- head rows n < N_e
    S_hi = arb(0)
    S_lo = arb(0)
    worst = None
    worst_n = -1
    for n in range(N_e):
        if n >= 1:
            S_hi += hi[n]
            S_lo += lo[n]
        TH = P1[n] + P2[M_v - 1 - n] + P2[M_v + n]
        J = arb(M_v - n)
        rem2 = 2 * b * Ca / two_p * (one / J ** (p + 2) + one / ((p + 1) * J ** (p + 1)))
        rem = (rem_c * Tsuf[M_v - n]).min(rem2)
        vn = v[n]
        if variant == "mhat":
            # Mhat row (ERRATA F2: diagonal sum over k = 1..n of lo)
            row = TH + rem - 2 * vn * S_lo - (2 * n + 1) * vn
        else:
            row = TH + rem - 2 * vn * S_hi - (2 * n + 1) * vn
            if corr is not None:
                row = row + corr[n]
        q = row + eta_a * vn
        if not bool(q < 0):
            return _fail("row inequality not certified", row=n, value=str(q), variant=variant,
                         seconds=time.time() - t0)
        rel = q / vn
        if worst is None or float(rel.mid()) > float(worst.mid()):
            worst, worst_n = rel, n
    # ---- tail rows n >= N_e (Lemma T1(p)): G(n) nonincreasing, need G(N_e) < 1
    U0 = arb(0)
    for x in u:
        U0 += arb(x)
    n = N_e
    Z = one / (2 * (p - 1) * arb(2 * N) ** (p - 1))   # Z_N = sum_{j>=N} f(j) <= (2N)^(1-p)/(2(p-1))
    r = U0 / Ca
    tn = arb(2 * n + 1)
    tp1 = tn ** (p - 1)
    D = tn - eta_a
    if variant == "mhat":
        dL = arb(0)
        for k in range(1, N_e + 1):
            dL += hi[k] - lo[k]
        dL += b / N_e                                 # sum_{k>N_e} hi(k) <= b/N_e
        D = D - 2 * dL
    if not bool(D > 0):
        return _fail("tail lemma T1: D(N_e) <= 0", variant=variant, seconds=time.time() - t0)
    last = 4 * (r + Z) if p == 3 else (r + Z) * tp1 / arb(n + 1) ** 2
    G = (tn / D) * b * (2 * p * (p + 1) * tp1 / arb(n + 1) ** (p + 1)
                        + 4 * Z * tp1 / arb(n) ** 2
                        + r * tp1 / arb(n - N + 1) ** 2
                        + tp1 / (arb(4 * n + 3) ** p * n)
                        + last)
    if not bool(G < 1):
        return _fail("tail lemma T1 not certified", G=str(G), variant=variant,
                     seconds=time.time() - t0)
    psi = G
    return {"ok": True, "worst_row": worst_n, "worst_rel": float(worst.mid()),
            "G_tail": float(psi.mid()), "seconds": time.time() - t0}


# ============================================================================ float proposal
def _tail_drive(hi: np.ndarray, N: int, M_v: int, p: int) -> np.ndarray:
    """t_n = sum_{N <= m < M_v} [hi(|n-m|) + hi(n+m+1)] f(m), n < N (float)."""
    m = np.arange(N, M_v)
    f = (2 * m + 1.0) ** (-p)
    n = np.arange(N)
    t = np.zeros(N)
    step = max(1, 4_000_000 // max(1, len(m)))
    for i0 in range(0, N, step):
        nn = n[i0:i0 + step, None]
        t[i0:i0 + step] = ((hi[m[None, :] - nn] + hi[nn + m[None, :] + 1]) * f[None, :]).sum(1)
    return t


def _float_rows(hi: np.ndarray, lo: np.ndarray, v: np.ndarray, n0: int, n1: int) -> np.ndarray:
    """float Mhat rows n0 <= n < n1 over the columns m < len(v) (proposal only)."""
    M = len(v)
    m = np.arange(M)
    out = np.zeros(n1 - n0)
    S = np.concatenate([[0.0], np.cumsum(lo[1:n1 + 1])])
    step = max(1, 2_000_000 // max(1, M))
    for i0 in range(n0, n1, step):
        n = np.arange(i0, min(n1, i0 + step))[:, None]
        T = hi[np.abs(n - m[None, :])] * (n != m[None, :]) + hi[n + m[None, :] + 1]
        r = T @ v - (2 * S[n[:, 0]] + 2 * n[:, 0] + 1) * v[n[:, 0]]
        out[i0 - n0:i0 - n0 + len(r)] = r
    return out


def _G_float(b, r, N, n, p, eta, dL):
    tn = 2 * n + 1.0
    tp1 = tn ** (p - 1)
    Z = (2.0 * N) ** (1 - p) / (2 * (p - 1))
    D = tn - eta - 2 * dL
    if D <= 0:
        return math.inf
    last = 4 * (r + Z) if p == 3 else (r + Z) * tp1 / (n + 1) ** 2
    return (tn / D) * b * (2 * p * (p + 1) * tp1 / (n + 1) ** (p + 1) + 4 * Z * tp1 / n ** 2
                           + r * tp1 / (n - N + 1) ** 2 + tp1 / ((4 * n + 3.0) ** p * n) + last)


def propose_cell(fs: FloatSpectrum, T1: float, T2: float, variant: str, N: int,
                 p: int = 3, eta_prime_frac: float = 0.5, eta_frac: float = 0.25,
                 max_Ne_factor: int = 64, G_target: float = 0.9) -> dict | None:
    """Float proposal of a cell certificate (untrusted).  None if the float data say no."""
    e = fs.e
    b = (1 + e) * fs.sum_pom / (2 * math.pi ** 2 * T1 ** 2)
    C = 1.0
    M_v = 6 * N
    Kf = None
    hi = lo = None
    for _it in range(4):
        if Kf is None or Kf < N + M_v + 2:
            Kf = int((N + M_v + 2) * 1.3)
            hi = fs.lam(T1, Kf, 1 + e)
            lo = fs.lam(T2, Kf, 1 - e)
        if _it == 0:
            M = build_matrix(hi, N, lo if variant == "mhat" else None)
            rho_N, _ = top_eig(M, want_vec=False)
            if not rho_N < 0:
                return None
            eta_p = eta_prime_frac * (-rho_N)
            eta = eta_frac * (-rho_N)
            Afac = M + eta_p * np.eye(N)
        t = _tail_drive(hi, N, M_v, p)
        # crude estimate of the columns m >= M_v (proposal only)
        t += 2 * C * (2 * M_v + 1.0) ** (-p) * b / (M_v - np.arange(N) - 1.0)
        try:
            u = sla.solve(Afac, -t, assume_a="sym")
        except Exception:
            return None
        if not np.all(np.isfinite(u)) or not np.all(u > 0):
            return None
        r = float(u.sum()) / C
        # smallest N_e >= 2N with G(N_e) <= G_target  (G is nonincreasing)
        n = 2 * N
        lim = max_Ne_factor * N
        while True:
            dL = 0.0
            if variant == "mhat":
                kk = min(n, Kf)
                dL = float((hi[1:kk + 1] - lo[1:kk + 1]).sum()) + b / kk
            if _G_float(b, r, N, n, p, eta, dL) <= G_target:
                break
            n = int(n * 1.3) + 1
            if n > lim:
                return None
        N_e = n
        # columns m < M_v explicit; M_v such that the Lemma R1 remainder (~ 2 C f(M_v) b/(M_v-n-1))
        # is <= 2% of the row margin: eta*v_n for head rows, the float margin for rows N..N_e-1
        nn = np.arange(N_e)
        M0 = max(M_v, 2 * N_e)
        if Kf < N_e + M0 + 2:
            Kf = int((N_e + M0 + 2) * 1.3)
            hi = fs.lam(T1, Kf, 1 + e)
            lo = fs.lam(T2, Kf, 1 - e)
        vv = np.concatenate([u, C * (2 * np.arange(N, M0) + 1.0) ** (-p)])
        marg = eta * vv[:N_e]
        tr = _float_rows(hi, lo if variant == "mhat" else hi, vv, N, N_e)
        mt = -tr - eta * vv[N:N_e]
        if np.any(mt <= 0):
            return None
        marg[N:] = mt
        M_new = M0
        while True:
            J = M_new - nn
            rem = np.minimum(2 * C * (2 * M_new + 1.0) ** (-p) * b / (J - 1.0),
                             2 * b * C * 2.0 ** (-p) * (J ** (-p - 2.0) + J ** (-p - 1.0) / (p + 1)))
            if np.all(rem <= 0.02 * marg):
                break
            M_new = int(M_new * 1.25)
            if M_new > 200 * N_e + 200000:
                return None
        if M_new <= M_v:
            break
        M_v = M_new
    return {"T1": float(T1).hex(), "T2": float(T2).hex(), "variant": variant, "p": int(p),
            "N": int(N), "u": [float(x) for x in u], "C": C, "eta": float(eta),
            "N_e": int(N_e), "M_v": int(M_v), "float_rho_N": float(rho_N)}


def head_size(fs: FloatSpectrum, T: float, fac: float = 6.0, nmin: int = 48) -> int:
    k0 = fs.om_max / (2 * math.pi * T)
    return max(nmin, int(math.ceil(fac * k0)))


def certify_cell(asp: ArbSpectrum, fs: FloatSpectrum, T1: float, T2: float, variant: str,
                 cache: LamCache, N0: int | None = None, N_cap: int = 2400, log=None):
    """Propose (float) then check (Arb).  Returns (cert, check_result) or (None, last_result)."""
    N = N0 or head_size(fs, T1, fac=3.0)
    last = None
    while N <= N_cap:
        cert = propose_cell(fs, T1, T2, variant, N)
        if cert is None:
            last = {"ok": False, "reason": "float proposal failed", "N": N}
        else:
            res = check_cell(asp, cert, cache)
            if res["ok"]:
                return cert, res
            last = res
            if res.get("reason", "").startswith("row inequality") and res.get("row", 0) < N:
                # a head row failed: increasing N rarely helps; stop
                return None, last
        if log:
            log(f"    cell [{T1:.9g},{T2:.9g}] N={N}: {last.get('reason')}")
        N *= 2
    return None, last
