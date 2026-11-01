"""C-RAYLEIGH: rigorous instability certificate (lower bound on Tc_model).

Certificate: {"T": hex double, "mu_hi": double, "omega_c": double, "N_c": int, "x": [doubles]}.

check_rayleigh() verifies in Arb (PROOFS.md, Thm R; Lean rayleigh_cert, boxMin_le,
sel_mul_eq_min, Qp_anti_mu):
   N_c = #{n >= 0 : (2n+1) pi T <= omega_c}      (certified comparisons),
   x != 0,
   boxMinQ = -sum (2n+1) x_n^2 - 2 mu_hi (sum_{n < min(N_c, N_x)} x_n)^2
             + sum_j min((1-e) p_j F_j, (1+e) p_j F_j)  >= 0,
   F_j = sum_{k>=1} g_T(k, Omega_j) gamma_k(x)   (finite sum, k <= 2 N_x - 1),
   gamma_k = 2 sum_n x_n x_{n+k} + sum_{n+m=k-1} x_n x_m - 2 sum_{n>=k} x_n^2.
Then every spectrum in the envelope, every mu* in [0, mu_hi], cutoff N_c(T): T is unstable,
hence Tc_model >= T.
"""
from __future__ import annotations

import math
import time

import numpy as np
from flint import arb, arb_poly

from .arbcore import ArbSpectrum, PI, TWO_PI, poly_coeffs
from .floatcore import FloatSpectrum, build_matrix, n_c, top_eig, default_N

R_SERIES = 40


def certified_n_c(T: float, omega_c: float) -> int | None:
    """N_c(T) with both defining inequalities certified in Arb; None if undecidable."""
    n0 = n_c(T, omega_c)
    Ta = arb(float(T))
    wc = arb(float(omega_c))
    # all n < n0 satisfy (2n+1) pi T <= omega_c  <=  n = n0-1 does (monotone in n)
    if n0 >= 1 and not bool((2 * (n0 - 1) + 1) * PI * Ta <= wc):
        return None
    # n = n0 violates it (hence all n >= n0 do)
    if not bool((2 * n0 + 1) * PI * Ta > wc):
        return None
    return n0


def _gamma(x: list) -> list:
    """gamma_k for k = 0..2Nx-1 (index 0 unused), Arb."""
    Nx = len(x)
    X = arb_poly(x)
    A1 = poly_coeffs(X * arb_poly(x[::-1]), 2 * Nx)       # index Nx-1+k: sum_n x_n x_{n+k}
    A2 = poly_coeffs(X * X, 2 * Nx)                        # index k-1: sum_{n+m=k-1} x_n x_m
    suf = [arb(0)] * (Nx + 1)                              # suf[k] = sum_{n>=k} x_n^2
    for k in range(Nx - 1, -1, -1):
        suf[k] = suf[k + 1] + x[k] * x[k]
    g = [arb(0)] * (2 * Nx)
    for k in range(1, 2 * Nx):
        val = A2[k - 1]
        if k <= Nx - 1:
            val = val + 2 * A1[Nx - 1 + k] - 2 * suf[k]
        g[k] = val
    return g


def check_rayleigh(asp: ArbSpectrum, cert: dict) -> dict:
    t0 = time.time()
    try:
        T = float.fromhex(cert["T"])
        mu = float(cert["mu_hi"])
        wc = float(cert["omega_c"])
        Nc_claim = int(cert["N_c"])
        xs = [float(v) for v in cert["x"]]
    except Exception as exc:
        return {"ok": False, "reason": f"malformed certificate: {exc!r}"}
    if not (T > 0 and math.isfinite(T) and mu >= 0 and math.isfinite(mu) and wc > 0):
        return {"ok": False, "reason": "need T > 0, mu_hi >= 0, omega_c > 0"}
    if not all(math.isfinite(v) for v in xs) or not any(v != 0 for v in xs):
        return {"ok": False, "reason": "x must be finite and nonzero"}
    Nc = certified_n_c(T, wc)
    if Nc is None or Nc != Nc_claim:
        return {"ok": False, "reason": f"N_c mismatch/undecided (claimed {Nc_claim}, got {Nc})"}
    Nx = len(xs)
    x = [arb(v) for v in xs]
    gam = _gamma(x)
    Kmax = 2 * Nx - 1
    Ta = arb(T)
    twopiT = TWO_PI * Ta
    ks = min(asp.k_series_start(T), Kmax + 1)
    # direct part: F_j^dir = sum_{k < ks} g_j(k) gamma_k
    x2 = [arb(0)] + [(twopiT * k) ** 2 for k in range(1, ks)]
    # series part (k >= ks, 2 pi k T >= 2 Omega_max):
    #   g_j(k) = (2/Om) sum_s (-1)^s t_j^(2s+1) q_k^(s+1),  q_k = (Om/x_k)^2 <= 1/4, t_j = Om_j/Om
    Gs = None
    Eabs = arb(0)
    if ks <= Kmax:
        Om = asp.om_max
        Gs = [arb(0)] * R_SERIES
        for k in range(ks, Kmax + 1):
            q = (Om / (twopiT * k)) ** 2
            if not bool(q <= 0.25):
                return {"ok": False, "reason": "series ratio check failed"}
            pw = q
            gk = gam[k]
            for s in range(R_SERIES):
                Gs[s] += gk * pw
                pw = pw * q
            Eabs += abs(gk) * pw                         # q^(R+1) |gamma_k|
    e1 = 1 + asp.e
    e2 = 1 - asp.e
    total = arb(0)
    for pj, oj, aj, o2j in zip(asp.p, asp.om, asp.a, asp.o2):
        F = arb(0)
        for k in range(1, ks):
            F += gam[k] * (2 * oj / (o2j + x2[k]))
        if Gs is not None:
            t = oj / asp.om_max
            t2 = t * t
            c = t
            ser = arb(0)
            for s in range(R_SERIES):
                ser = ser + c * Gs[s] if s % 2 == 0 else ser - c * Gs[s]
                c = c * t2
            err = c * Eabs                                # t^(2R+1) * sum |gamma| q^(R+1)
            ser = (ser - err).union(ser + err)
            F += (2 / asp.om_max) * ser
        total += (e2 * pj * F).min(e1 * pj * F)        # Lean sel_mul_eq_min (ERRATA F6)
    diag = arb(0)
    for n in range(Nx):
        diag += (2 * n + 1) * x[n] * x[n]
    sx = arb(0)
    for n in range(min(Nc, Nx)):
        sx += x[n]
    Q = total - diag - 2 * arb(mu) * sx * sx
    ok = bool(Q >= 0)
    nx2 = arb(0)
    for v in x:
        nx2 += v * v
    return {"ok": ok, "Q": str(Q), "Q_over_x2": float((Q / nx2).mid()),
            "reason": None if ok else "boxMinQ >= 0 not certified", "seconds": time.time() - t0}


# ---------------------------------------------------------------------------- proposals
def propose_rayleigh(fs: FloatSpectrum, T: float, mu: float, omega_c: float, N: int | None = None,
                     cap: int = 3000) -> dict | None:
    if N is None:
        N = default_N(fs, T, mu, omega_c, cap)
    lam = fs.lam(T, 2 * N + 2, 1 - fs.e)
    Nc = n_c(T, omega_c)
    M = build_matrix(lam, N, mu=mu, Nc=Nc)
    r, v = top_eig(M)
    if r < 0:
        return None
    v = v * np.sign(v[0] if v[0] != 0 else 1.0)
    v = v / np.abs(v).max()
    return {"T": float(T).hex(), "mu_hi": float(mu), "omega_c": float(omega_c), "N_c": int(Nc),
            "x": [float(t) for t in v], "float_rho": float(r)}
