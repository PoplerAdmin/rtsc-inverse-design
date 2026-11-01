"""Independent float TEST ORACLE for the Tc model of SPEC.md / ARCH.md section 4.

Written independently of rtt.core.certify (team CERT): it only uses the model definition.

Model (k_B = hbar = 1, kelvin).  Finite atomic spectrum alpha^2F = sum_j w_j delta(omega - o_j).
    lam(k) = sum_j w_j * 2 o_j / (o_j^2 + (2 pi k T)^2)
Even-sector matrix on n, m = 0..N-1 (truncation of the infinite operator M+_T):
    M_nm = lam(|n-m|) + lam(n+m+1)                                  (n != m)
    M_nn = lam(2n+1) - (2n+1) - 2 sum_{k=1}^{n} lam(k)             (lam(0) cancels exactly)
    M_nm -= 2 mu*   for n, m < N_c(T),  N_c(T) = #{n >= 0 : (2n+1) pi T <= omega_c}
rho+(T) = sup of the Rayleigh quotient = lim_N top eigenvalue of the N x N truncation
(the truncations are principal submatrices, so the top eigenvalue is NONDECREASING in N; a
truncated value is a lower bound on rho+).

Tc_model := sup{T : rho+(T) >= 0}  (the odd sector is always stable, ARCH R-ODD).

This is FLOAT code.  It is a test oracle, never a certificate.
"""
from __future__ import annotations

import math
from dataclasses import dataclass, field

import numpy as np
from scipy.linalg import eigh, hankel, toeplitz
from scipy.optimize import brentq

A0 = 4.33148145            # universal ceiling Tc <= A / A0 (PROOF.md)
N_CAP = 1500               # max Matsubara truncation
TOL_RHO = 1e-10            # convergence target for rho+ near the root


class OracleCapExceeded(RuntimeError):
    pass


# ----------------------------------------------------------------------------------------------
# model primitives
# ----------------------------------------------------------------------------------------------
def n_c(T: float, omega_c: float) -> int:
    """N_c(T) = #{n >= 0 : (2n+1) pi T <= omega_c}."""
    x = omega_c / (math.pi * T)          # need 2n+1 <= x
    if x < 1.0:
        return 0
    n = int(math.floor((x - 1.0) / 2.0)) + 1
    # guard floating edge effects exactly as the definition reads
    while n > 0 and (2 * (n - 1) + 1) * math.pi * T > omega_c:
        n -= 1
    while (2 * n + 1) * math.pi * T <= omega_c:
        n += 1
    return n


def lam_array(omega: np.ndarray, weight: np.ndarray, T: float, kmax: int) -> np.ndarray:
    """lam(k) for k = 0..kmax (lam(0) = 2 sum w/o)."""
    k = np.arange(kmax + 1, dtype=float)
    nu2 = (2.0 * math.pi * T * k) ** 2
    out = np.empty(kmax + 1)
    # chunk over atoms to bound memory for thousands of atoms
    chunk = max(1, 4_000_000 // (kmax + 1))
    acc = np.zeros(kmax + 1)
    for s in range(0, len(omega), chunk):
        o = omega[s:s + chunk, None]
        w = weight[s:s + chunk, None]
        acc += np.sum(w * 2.0 * o / (o * o + nu2[None, :]), axis=0)
    out[:] = acc
    return out


def even_matrix(omega, weight, T, mu, omega_c, N, lam=None) -> np.ndarray:
    if lam is None:
        lam = lam_array(omega, weight, T, 2 * N + 1)
    M = toeplitz(lam[:N]) + hankel(lam[1:N + 1], lam[N:2 * N])
    S = np.concatenate([[0.0], np.cumsum(lam[1:N])])           # S_n = sum_{k=1}^n lam(k)
    n = np.arange(N)
    M[n, n] = lam[2 * n + 1] - (2 * n + 1) - 2.0 * S
    nc = min(n_c(T, omega_c), N)
    if mu != 0.0 and nc > 0:
        M[:nc, :nc] -= 2.0 * mu
    return M


def top_eig(M: np.ndarray) -> float:
    N = M.shape[0]
    return float(eigh(M, eigvals_only=True, subset_by_index=[N - 1, N - 1], driver="evr")[0])


# ----------------------------------------------------------------------------------------------
@dataclass
class Model:
    omega: np.ndarray
    weight: np.ndarray
    mu: float
    omega_c: float
    n_cap: int = N_CAP
    nevals: int = 0
    cap_delta: float = 0.0

    @classmethod
    def from_spectrum(cls, omega_K, weight_K, mu, omega_c=None, scale=1.0, n_cap=N_CAP):
        o = np.asarray(omega_K, float)
        w = np.asarray(weight_K, float) * scale
        keep = w > 0
        o, w = o[keep], w[keep]
        if omega_c is None:
            omega_c = 10.0 * float(np.max(o))
        return cls(o, w, float(mu), float(omega_c), n_cap)

    @property
    def A(self) -> float:
        return math.fsum(self.weight)

    def n_start(self, T: float) -> int:
        nc = n_c(T, self.omega_c) if self.mu != 0.0 else 0
        x = float(np.max(self.omega)) / (2 * math.pi * T)
        return int(max(40, 1.25 * nc + 24, 3.0 * x + 24))

    def rho_N(self, T: float, N: int) -> float:
        self.nevals += 1
        return top_eig(even_matrix(self.omega, self.weight, T, self.mu, self.omega_c, N))

    def rho(self, T: float, tol: float = TOL_RHO, rel: float = 0.0, N0: int | None = None):
        """Converged rho+(T).  Returns (rho, N, delta) where delta = last change on growing N.

        Stops when |change| < tol, or (rel > 0) when |change| < rel*|rho| (sign decided),
        or immediately-safe: rho_N >= 0 with rel > 0 (truncation only lowers rho, so the sign is
        already certain in floats).
        """
        N = max(N0 or 0, self.n_start(T))
        if N > self.n_cap:
            raise OracleCapExceeded(f"N_start={N} > cap at T={T}")
        r = self.rho_N(T, N)
        if rel > 0 and r >= 0:
            return r, N, float("nan")
        d_prev = float("inf")
        while True:
            N2 = int(math.ceil(1.6 * N))
            if N2 > self.n_cap:
                if N < self.n_cap:
                    N2 = self.n_cap
                else:
                    # at the cap without reaching tol: accept, record the residual change
                    if rel > 0 and d_prev > 0.2 * abs(r):
                        raise OracleCapExceeded(f"sign undecided at cap, T={T} (rho={r}, d={d_prev})")
                    self.cap_delta = max(self.cap_delta, d_prev)
                    return r, N, d_prev
            r2 = self.rho_N(T, N2)
            d = abs(r2 - r)
            N, r, d_prev = N2, r2, d
            if d < tol or (rel > 0 and (d < rel * abs(r) or r >= 0)):
                return r, N, d


@dataclass
class TcResult:
    Tc: float | None                  # None if no instability found above T_floor
    err: float                        # estimated absolute error (K)
    status: str                       # "ok" | "below_floor"
    T_floor: float | None = None      # lowest T scanned (all scanned T >= floor stable)
    N_used: int = 0
    crossings: list = field(default_factory=list)   # [(T_unstable_side, T_stable_side)] brackets
    reentrant: bool = False
    at_Nc_jump: bool = False
    rho_at_ceiling: float = float("nan")
    notes: list = field(default_factory=list)
    nevals: int = 0

    def to_dict(self):
        return {k: getattr(self, k, None) for k in ("Tc", "err", "status", "T_floor", "N_used",
                                                    "crossings", "reentrant", "at_Nc_jump",
                                                    "rho_at_ceiling", "notes", "nevals",
                                                    "jumps_checked")}


def find_tc(model: Model, ratio: float = 0.98, below_factor: float = 0.5,
            T_floor: float | None = None, A_ceiling: float | None = None,
            check_jumps: bool = True, max_jumps: int = 700) -> TcResult:
    """Tc = largest T with rho+(T) >= 0, scanning down from A/A0 on a log grid."""
    A = A_ceiling if A_ceiling is not None else model.A
    T_top = A / A0
    if T_floor is None:
        T_floor = T_top * 1e-4
    # floor also limited by the truncation cap (need N >= n_start)
    notes = []
    r_top, _, _ = model.rho(T_top, rel=1e-3)
    res = TcResult(None, float("nan"), "below_floor", rho_at_ceiling=r_top)
    if r_top >= 0:
        notes.append("THEOREM VIOLATION? rho+ >= 0 at the universal ceiling")
    # ---- downward scan
    grid = []
    T = T_top
    prev = (T_top, r_top)
    found = None
    while True:
        T = T * ratio
        if T < T_floor:
            break
        try:
            r, _, _ = model.rho(T, rel=1e-3)
        except OracleCapExceeded as e:
            notes.append(f"scan stopped: {e}")
            T_floor = T / ratio
            break
        grid.append((T, r))
        if r >= 0 and found is None:
            found = (T, prev[0])           # (unstable T_lo, stable T_hi)
            T_stop = T * below_factor
        if found is not None and T < T_stop:
            break
        prev = (T, r)
    res.notes = notes
    # ---- both sides of every N_c(T) jump above the first unstable grid point (or above the
    # floor if none).  rho+ is discontinuous there and may jump UP when T decreases past a
    # jump, so narrow unstable windows just above a jump temperature can hide between grid
    # points.  T_j = omega_c / ((2j+1) pi) belongs to the piece with N_c = j+1.
    pts = {T_top: r_top}
    pts.update(dict(grid))
    g = found[0] if found is not None else T_floor
    res.jumps_checked = 0
    if model.mu != 0.0 and check_jumps:
        wc = model.omega_c
        jmin = max(0, int(math.ceil((wc / (math.pi * T_top) - 1) / 2)) - 1)
        jmax = int(math.floor((wc / (math.pi * g) - 1) / 2)) + 1
        js = list(range(jmin, jmax + 1))
        if len(js) > max_jumps:
            notes.append(f"jump checks limited to the {max_jumps} lowest-T jumps of {len(js)}")
            js = js[-max_jumps:]
        for j in js:
            Tj = wc / ((2 * j + 1) * math.pi)
            t_at = Tj
            for _ in range(8):
                if n_c(t_at, wc) >= j + 1:
                    break
                t_at = math.nextafter(t_at, 0.0)
            t_ab = Tj * (1 + 1e-12)
            for t in (t_at, t_ab):
                if g < t < T_top and t not in pts:
                    try:
                        pts[t] = model.rho(t, rel=1e-3)[0]
                        res.jumps_checked += 1
                    except OracleCapExceeded as e:
                        notes.append(f"jump check skipped at T={t}: {e}")
    P = sorted(pts.items())
    uns = [t for t, r in P if r >= 0]
    if not uns:
        res.T_floor = T_floor
        res.nevals = model.nevals
        return res
    cross = []
    for (Ta, ra), (Tb, rb) in zip(P, P[1:]):          # ascending T
        if (ra < 0) != (rb < 0):
            cross.append((Ta, Tb) if ra >= 0 else (Tb, Ta))   # (unstable side, stable side)
    res.crossings = sorted(cross, key=lambda c: -c[0])
    res.reentrant = len(cross) > 1
    u = max(uns)
    s_ = min(t for t, _ in P if t > u)
    # ---- refine top crossing with fixed converged N
    T_lo, T_hi = u, s_
    model.cap_delta = 0.0
    r_lo, N_lo, _ = model.rho(T_lo)
    r_hi, N_hi, _ = model.rho(T_hi)
    N = max(N_lo, N_hi)
    f = lambda t: model.rho_N(t, N)
    flo, fhi = f(T_lo), f(T_hi)
    if not (flo >= 0 > fhi):
        notes.append(f"bracket sign issue at fixed N={N}: {flo}, {fhi}")
    jump = model.mu != 0 and n_c(T_lo, model.omega_c) != n_c(T_hi, model.omega_c)
    if jump:
        Tc = T_lo
    else:
        Tc = brentq(f, T_lo, T_hi, xtol=1e-14 * T_hi, rtol=1e-15)
    res.at_Nc_jump = bool(jump)
    # truncation error at the root and slope
    r_N = f(Tc)
    N2 = min(int(1.6 * N), model.n_cap)
    if N2 > N:
        d_trunc = abs(model.rho_N(Tc, N2) - r_N)
    else:
        d_trunc = 1e-10
        notes.append("truncation check at cap")
    h = 1e-6 * Tc
    if jump:
        slope = float("nan")
        notes.append("Tc sits on an N_c(T) jump: Tc is the jump temperature, rho+ >= 0 there "
                     f"(rho={flo:.3e}) and < 0 just above (rho={fhi:.3e})")
        err = (T_hi - T_lo) + 4e-14 * Tc
        if flo < 1e-8 or fhi > -1e-8:
            notes.append("WARNING: jump-side rho values close to 0")
    else:
        slope = (f(Tc + h) - f(Tc - h)) / (2 * h)
        d_trunc = max(d_trunc, 3.0 * model.cap_delta)
        err = (d_trunc + 1e-13 * N) / max(abs(slope), 1e-300) + 4e-14 * Tc
    res.Tc = float(Tc)
    res.err = float(max(err, 1e-12 * Tc))
    res.status = "ok"
    res.N_used = N
    if model.cap_delta > 0:
        notes.append(f"truncation capped at N={model.n_cap}; residual change {model.cap_delta:.2e}")
    res.nevals = model.nevals
    res.notes = notes
    res.slope = slope
    return res


# ----------------------------------------------------------------------------------------------
# independent cross-check 1: kernel eigenvalue = 1 formulation (Allen-Dynes "eta" form)
# ----------------------------------------------------------------------------------------------
def eta_minus_one(model: Model, T: float, N: int) -> float:
    """Top eigenvalue of D^{-1/2} K D^{-1/2} minus 1, where K = Coulomb-corrected kernel incl.
    lam(0) and D_n = (2n+1) + lam(0) + 2 sum_{k<=n} lam(k).  Sylvester inertia: same zero set
    as rho+, but a different function of T."""
    lam = lam_array(model.omega, model.weight, T, 2 * N + 1)
    K = toeplitz(lam[:N]) + hankel(lam[1:N + 1], lam[N:2 * N])
    nc = min(n_c(T, model.omega_c), N)
    if model.mu and nc:
        K[:nc, :nc] -= 2 * model.mu
    S = np.concatenate([[0.0], np.cumsum(lam[1:N])])
    D = (2 * np.arange(N) + 1) + lam[0] + 2 * S
    s = 1 / np.sqrt(D)
    return top_eig(K * s[:, None] * s[None, :]) - 1.0


def tc_eta(model: Model, Tc: float, N: int, width: float = 1e-3) -> float:
    g = lambda t: eta_minus_one(model, t, N)
    a, b = Tc * (1 - width), Tc * (1 + width)
    if not (g(a) >= 0 > g(b)):
        return float("nan")
    return brentq(g, a, b, xtol=1e-14 * Tc, rtol=1e-15)


# ----------------------------------------------------------------------------------------------
# independent cross-check 2: brute-force full Matsubara line (n in Z), folded by hand
# ----------------------------------------------------------------------------------------------
def full_line_matrix_float(omega, weight, T, mu, omega_c, L, Ksum=200000):
    """H on n = -L-1..L from the ORIGINAL equation:
    Z_n = |w_n|/(pi T) with wtilde_n = w_n + pi T sum_m lam(n-m) sgn(w_m) summed directly over
    m in [-Ksum-1, Ksum] (lam computed per atom directly, no telescoping identity used);
    H_nm = lam(n-m) - mu* th_n th_m - Z_n delta_nm, th_n = [|w_n| <= omega_c]."""
    ns = np.arange(-L - 1, L + 1)
    kk = np.arange(0, Ksum + L + 3)
    lam = lam_array(np.asarray(omega), np.asarray(weight), T, len(kk) - 1)
    ms = np.arange(-Ksum - 1, Ksum + 1)
    sg = np.where(ms >= 0, 1.0, -1.0)
    Z = np.empty(len(ns))
    for i, n in enumerate(ns):
        s = np.sum(sg * lam[np.abs(n - ms)])
        wt = (2 * n + 1) + s           # wtilde_n / (pi T)
        Z[i] = abs(wt)
    H = lam[np.abs(ns[:, None] - ns[None, :])].astype(float)
    th = (np.abs(2 * ns + 1) * math.pi * T <= omega_c).astype(float)
    H -= mu * np.outer(th, th)
    H[np.arange(len(ns)), np.arange(len(ns))] -= Z
    return ns, H


def fold_even(ns, H, L):
    """P^T H P with P the explicit even embedding phi_n = phi_{-n-1} = x_n / sqrt2."""
    P = np.zeros((len(ns), L + 1))
    for i, n in enumerate(ns):
        j = n if n >= 0 else -n - 1
        P[i, j] = 1 / math.sqrt(2)
    return P.T @ H @ P


def fold_odd(ns, H, L):
    P = np.zeros((len(ns), L + 1))
    for i, n in enumerate(ns):
        j = n if n >= 0 else -n - 1
        P[i, j] = (1 if n >= 0 else -1) / math.sqrt(2)
    return P.T @ H @ P


def mp_even_matrix_fullline(omega, weight, T, mu, omega_c, L, dps=30):
    """High-precision (mpmath) brute force: build full-line H in mp with Z_n from the direct
    sgn-sum over a finite symmetric window m in [-K-1, K] plus the EXACT finite remainder
    sum_{k=K-n+1}^{K+n+1} lam(k) (n >= 0), then fold onto the even sector."""
    import mpmath as mp
    mp.mp.dps = dps
    om = [mp.mpf(float(o)) for o in omega]
    wt = [mp.mpf(float(w)) for w in weight]
    Tm = mp.mpf(float(T))
    twopiT = 2 * mp.pi * Tm
    cache = {}

    def lam(k):
        k = abs(k)
        if k not in cache:
            nu2 = (twopiT * k) ** 2
            cache[k] = mp.fsum(w * 2 * o / (o * o + nu2) for o, w in zip(om, wt))
        return cache[k]

    K = 4 * L + 40
    Zs = {}
    for n in range(0, L + 1):
        s = mp.fsum((1 if m >= 0 else -1) * lam(n - m) for m in range(-K - 1, K + 1))
        s += mp.fsum(lam(k) for k in range(K - n + 1, K + n + 2))
        Zs[n] = abs((2 * n + 1) + s)
    for n in range(-L - 1, 0):
        Zs[n] = Zs[-n - 1]            # |w_n| symmetric; checked numerically in float version
    ns = list(range(-L - 1, L + 1))
    th = {n: (abs(2 * n + 1) * mp.pi * Tm <= mp.mpf(float(omega_c))) for n in ns}
    # fold: M_ij = sum over (n in {i,-i-1}, m in {j,-j-1}) H_nm / 2
    M = mp.matrix(L + 1, L + 1)
    for i in range(L + 1):
        for j in range(L + 1):
            acc = mp.mpf(0)
            for n in (i, -i - 1):
                for m in (j, -j - 1):
                    h = lam(n - m) - (mu if (th[n] and th[m]) else 0)
                    if n == m:
                        h -= Zs[n]
                    acc += h
            M[i, j] = acc / 2
    return M


def mp_inertia_positive(M, shift, dps=30):
    """Number of eigenvalues of symmetric mp matrix M greater than `shift` (LDL^T pivots,
    Sylvester's law of inertia; no pivoting -- fine for generic shifts)."""
    import mpmath as mp
    mp.mp.dps = dps
    n = M.rows
    A = [[M[i, j] for j in range(n)] for i in range(n)]
    for i in range(n):
        A[i][i] -= shift
    pos = 0
    for k in range(n):
        p = A[k][k]
        if p == 0:
            raise ZeroDivisionError("zero pivot")
        if p > 0:
            pos += 1
        rowk = A[k]
        for i in range(k + 1, n):
            f = A[i][k] / p
            if f != 0:
                Ai = A[i]
                for j in range(k + 1, n):
                    Ai[j] -= f * rowk[j]
    return pos


# ----------------------------------------------------------------------------------------------
def allen_dynes(lam: float, wlog: float, w2: float, mu: float) -> float:
    """Allen-Dynes 1975 modified McMillan formula (literature sanity check only)."""
    den = lam - mu * (1 + 0.62 * lam)
    if den <= 0:
        return 0.0
    L1 = 2.46 * (1 + 3.8 * mu)
    L2 = 1.82 * (1 + 6.3 * mu) * (w2 / wlog)
    f1 = (1 + (lam / L1) ** 1.5) ** (1 / 3)
    f2 = 1 + (w2 / wlog - 1) * lam ** 2 / (lam ** 2 + L2 ** 2)
    return f1 * f2 * wlog / 1.2 * math.exp(-1.04 * (1 + lam) / den)
