"""Float (numpy/scipy) PROPOSAL code.  Nothing here is trusted: every number produced here is
either an input to an Arb check (test vectors, cell endpoints, margins) or informational
(float Tc used for validation and for steering searches).
"""
from __future__ import annotations

import math

import numpy as np
import scipy.linalg as sla

TWO_PI = 2 * math.pi


class FloatSpectrum:
    def __init__(self, omega_K, weight_K, rel_envelope=0.0):
        om = np.asarray(omega_K, float)
        w = np.asarray(weight_K, float)
        keep = w > 0
        self.om = om[keep]
        self.p = w[keep]
        self.e = float(rel_envelope)
        self.a = 2 * self.p * self.om
        self.o2 = self.om ** 2
        self.om_max = float(self.om.max())
        self.A = float(self.p.sum())
        self.sum_pom = float((self.p * self.om).sum())
        import hashlib
        self.digest = hashlib.sha256(self.om.tobytes() + self.p.tobytes()
                                     + np.float64(self.e).tobytes()).hexdigest()

    def lam(self, T: float, kmax: int, factor: float = 1.0) -> np.ndarray:
        """lambda(k; T) for k = 0..kmax (index 0 set to 0)."""
        out = np.zeros(kmax + 1)
        k = np.arange(1, kmax + 1, dtype=float)
        x2 = (TWO_PI * T * k) ** 2
        ks = int(min(kmax, math.ceil(4 * self.om_max / (TWO_PI * T)) + 2))
        # direct, chunked
        step = max(1, int(2_000_000 // max(1, len(self.om))))
        for i0 in range(0, ks, step):
            i1 = min(ks, i0 + step)
            out[1 + i0:1 + i1] = (self.a[None, :] / (self.o2[None, :] + x2[i0:i1, None])).sum(1)
        if ks < kmax:
            # series in r = Omega^2/x^2 <= 1/16 (float only)
            xr = x2[ks:]
            m = [float((self.a * self.o2 ** s).sum()) for s in range(12)]
            acc = np.zeros_like(xr)
            pw = 1.0 / xr
            for s in range(12):
                acc += ((-1) ** s) * m[s] * pw
                pw = pw / xr
            out[1 + ks:] = acc
        return out * factor


def n_c(T: float, omega_c: float) -> int:
    """float proposal of N_c(T) = #{n >= 0 : (2n+1) pi T <= omega_c} (re-checked in Arb)."""
    if omega_c < math.pi * T:
        return 0
    return int(math.floor((omega_c / (math.pi * T) - 1) / 2)) + 1


def build_matrix(lam_hi: np.ndarray, N: int, lam_lo: np.ndarray | None = None,
                 mu: float = 0.0, Nc: int = 0) -> np.ndarray:
    """Even-sector matrix (Mhat if lam_lo given; diagonal sum over k = 1..n, ERRATA F2)."""
    n = np.arange(N)
    M = lam_hi[np.abs(n[:, None] - n[None, :])] + lam_hi[n[:, None] + n[None, :] + 1]
    lo = lam_hi if lam_lo is None else lam_lo
    S = np.concatenate([[0.0], np.cumsum(lo[1:N])])          # S[n] = sum_{k=1}^n lo(k)
    d = lam_hi[2 * n + 1] - 2 * S[:N] - (2 * n + 1)
    M[n, n] = d
    if mu and Nc > 0:
        c = min(Nc, N)
        M[:c, :c] -= 2 * mu
    return M


def top_eig(M: np.ndarray, want_vec: bool = True, v0: np.ndarray | None = None):
    N = M.shape[0]
    if N <= 700:
        w, V = sla.eigh(M, subset_by_index=[N - 1, N - 1])
        return float(w[0]), (V[:, 0] if want_vec else None)
    from scipy.sparse.linalg import eigsh
    if v0 is not None and len(v0) != N:
        v0 = None
    w, V = eigsh(M, k=1, which="LA", tol=1e-12, maxiter=50000, v0=v0)
    return float(w[0]), V[:, 0]


def default_N(fs: FloatSpectrum, T: float, mu: float = 0.0, omega_c: float = 0.0,
              cap: int = 3000, fac: float = 10.0) -> int:
    """Float truncation size.  The truncated top eigenvalue is nondecreasing in N (Cauchy
    interlacing), so a smaller N can only make the float Tc smaller."""
    k0 = fs.om_max / (TWO_PI * T)
    N = max(64, int(math.ceil(fac * k0)))
    if mu > 0 and omega_c > 0:
        N = max(N, n_c(T, omega_c) + int(math.ceil(0.3 * fac * k0)) + 16)
    return min(N, cap)


_TC_CACHE: dict = {}


def rho(fs: FloatSpectrum, T: float, factor: float = 1.0, mu: float = 0.0,
        omega_c: float = 0.0, N: int | None = None, cap: int = 3000, state: dict | None = None):
    """float top eigenvalue of the truncated even-sector matrix at T."""
    if N is None:
        N = default_N(fs, T, mu, omega_c, cap)
    lam = fs.lam(T, 2 * N + 2, factor)
    Nc = n_c(T, omega_c) if mu > 0 else 0
    M = build_matrix(lam, N, mu=mu, Nc=Nc)
    v0 = state.get("v") if state is not None else None
    r, v = top_eig(M, want_vec=True, v0=v0)
    if state is not None and v is not None:
        state["v"] = v
    return r


def unstable(fs: FloatSpectrum, T: float, factor: float = 1.0, mu: float = 0.0,
             omega_c: float = 0.0, N: int | None = None, cap: int = 3000, fac: float = 10.0) -> bool:
    """float test rho(T) >= 0  <=>  -M is not positive definite (Cholesky fails)."""
    if N is None:
        N = default_N(fs, T, mu, omega_c, cap, fac)
    lam = fs.lam(T, 2 * N + 2, factor)
    Nc = n_c(T, omega_c) if mu > 0 else 0
    M = build_matrix(lam, N, mu=mu, Nc=Nc)
    try:
        sla.cholesky(-M, lower=True, check_finite=False, overwrite_a=True)
        return False
    except np.linalg.LinAlgError:
        return True


def float_tc(fs: FloatSpectrum, factor: float = 1.0, mu: float = 0.0, omega_c: float = 0.0,
             T_hi: float | None = None, T_min: float = 1e-3, rtol: float = 1e-10,
             N: int | None = None, cap: int = 3000):
    """Largest T with rho(T) >= 0 (float proposal).  Coarse stage (truncation fac = 3, scan by
    factor 0.5, bisection to 1e-3), then fine stage (fac = 10) bisection to rtol.  Because the
    truncated rho is nondecreasing in N, a coarse-unstable T is fine-unstable.  0.0 if nothing
    unstable above T_min.  Cached per FloatSpectrum object."""
    if T_hi is None:
        T_hi = factor * fs.A / 4.33148145 * (1 + 1e-9)
    key = (fs.digest, factor, mu, omega_c, N, cap)   # content hash (id() may be reused)
    hit = _TC_CACHE.get(key)
    if hit is not None and hit[1] >= T_hi * (1 - 1e-12) and (hit[0] < T_hi or hit[0] == 0.0):
        return hit[0]
    fc = lambda T: unstable(fs, T, factor, mu, omega_c, N, cap, 3.0)
    ff = lambda T: unstable(fs, T, factor, mu, omega_c, N, cap, 10.0)
    if ff(T_hi):
        _TC_CACHE[key] = (T_hi, T_hi)
        return T_hi
    hi = T_hi
    lo = hi * 0.5
    while not fc(lo):
        hi = lo
        lo = lo * 0.5
        if lo < T_min:
            # nothing found by the coarse scan; one fine look at T_min
            if ff(T_min):
                lo = T_min
                break
            _TC_CACHE[key] = (0.0, T_hi)
            return 0.0
    while hi / lo - 1 > 1e-3:
        mid = math.sqrt(lo * hi)
        if fc(mid):
            lo = mid
        else:
            hi = mid
    # fine stage: lo is fine-unstable; find a fine-stable upper end
    d = 2e-3
    hi = min(T_hi, lo * (1 + d))
    while ff(hi) and hi < T_hi:
        lo = hi
        d *= 2
        hi = min(T_hi, lo * (1 + d))
    if ff(hi):
        _TC_CACHE[key] = (hi, T_hi)
        return hi
    while hi / lo - 1 > rtol:
        mid = math.sqrt(lo * hi)
        if ff(mid):
            lo = mid
        else:
            hi = mid
    _TC_CACHE[key] = (lo, T_hi)
    return lo
