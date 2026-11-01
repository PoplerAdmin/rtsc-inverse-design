"""Clean isotropic linearized Migdal-Eliashberg Tc solver (imaginary axis, float64).

Model and conventions follow docs/theory/01_eliashberg_targets.md and the certified engine
(rtt/core/certify/PROOFS.md); this file is an independent float implementation. Units: k_B = hbar = 1, energies in kelvin.

Spectrum: a finite positive measure a2F = sum_i w_i delta(omega - omega_i), omega_i > 0.
    A      = sum_i w_i                       (= \\int a2F d omega)
    lambda = 2 sum_i w_i / omega_i
    lam(k) = sum_i 2 w_i omega_i / (omega_i^2 + (2 pi k T)^2),  k >= 1.

Even-frequency sector (the odd sector is never unstable; wMB Theorem 3, B = 1):
    M_nm = lam(|n-m|) [n != m] + lam(n+m+1) - 2 mu* [n < N_c][m < N_c]
           - delta_nm [ (2n+1) + 2 sum_{k=1}^n lam(k) ]
with N_c(T) = #{n >= 0 : (2n+1) pi T <= omega_c} (block Matsubara cutoff for mu*).
rho(T) = top eigenvalue of M; T unstable iff rho(T) >= 0. Tc = largest root of rho.

Truncation: n < N with (2N+1) pi T >= max(trunc_mult * omega_max, 1.5 * omega_c) (and N >= 64);
the result is checked against a doubled N in the validation script.
"""
from __future__ import annotations

import math
from dataclasses import dataclass

import numpy as np
from scipy.linalg import eigh
from scipy.optimize import brentq


@dataclass
class Spectrum:
    omega: np.ndarray   # K, > 0
    w: np.ndarray       # K (a2F d omega mass per atom), >= 0

    @property
    def A(self) -> float:
        return float(self.w.sum())

    @property
    def lam(self) -> float:
        return float(2 * (self.w / self.omega).sum())

    @property
    def omega_log(self) -> float:
        g = 2 * self.w / self.omega
        return float(math.exp((g * np.log(self.omega)).sum() / g.sum()))

    @property
    def omega2(self) -> float:
        """sqrt(<omega^2>), <omega^2> = (2/lambda) \\int omega a2F."""
        g = 2 * self.w / self.omega
        return float(math.sqrt((g * self.omega ** 2).sum() / g.sum()))

    @property
    def omega1(self) -> float:
        """<omega> = (2/lambda) A."""
        return 2 * self.A / self.lam

    @property
    def lam_omega2(self) -> float:
        """lambda <omega^2> = 2 \\int omega a2F."""
        return float(2 * (self.w * self.omega).sum())

    @property
    def omega_max(self) -> float:
        return float(self.omega[self.w > 0].max())

    def scaled(self, lam_factor: float = 1.0, freq_factor: float = 1.0) -> "Spectrum":
        """Multiply a2F by lam_factor and stretch the frequency axis by freq_factor
        (a2F(omega) -> a2F(omega/s): weights w scale by s, so lambda is unchanged)."""
        return Spectrum(self.omega * freq_factor, self.w * lam_factor * freq_factor)

    def moments(self) -> dict:
        return dict(A=self.A, lam=self.lam, omega_log=self.omega_log, omega2=self.omega2,
                    omega1=self.omega1, lam_omega2=self.lam_omega2, omega_max=self.omega_max)


# ---------------------------------------------------------------- shapes
def einstein(lam: float, Omega: float) -> Spectrum:
    return Spectrum(np.array([float(Omega)]), np.array([lam * Omega / 2.0]))


def from_density(omega: np.ndarray, dens: np.ndarray, lam: float) -> Spectrum:
    """Discretize a density a2F(omega) given on a uniform grid; normalise to lambda."""
    d = omega[1] - omega[0]
    w = dens * d
    sp = Spectrum(omega.astype(float), w.astype(float))
    return sp.scaled(lam_factor=lam / sp.lam)


def gaussians(centers, widths, lam_shares, lam: float, n: int = 400,
              wmin_frac: float = 0.02) -> Spectrum:
    """Sum of Gaussian bands (in a2F, truncated to [c - 3s, c + 3s] and omega > 0) with given
    shares of lambda. Each band is normalised separately so lambda_band = share * lam."""
    om_all, w_all = [], []
    for c, s, f in zip(centers, widths, lam_shares):
        lo = max(c - 3 * s, wmin_frac * c)
        om = np.linspace(lo, c + 3 * s, n)
        dens = np.exp(-0.5 * ((om - c) / s) ** 2)
        b = from_density(om, dens, f * lam)
        om_all.append(b.omega)
        w_all.append(b.w)
    return Spectrum(np.concatenate(om_all), np.concatenate(w_all))


# ---------------------------------------------------------------- kernel
def lam_k(sp: Spectrum, T: float, kmax: int) -> np.ndarray:
    k = np.arange(kmax + 1)
    x = (2 * math.pi * T * k) ** 2
    L = (2 * sp.w[None, :] * sp.omega[None, :] / (sp.omega[None, :] ** 2 + x[:, None])).sum(1)
    return L  # L[0] = lambda (unused in M)


def n_trunc(sp: Spectrum, T: float, omega_c: float, trunc_mult: float) -> int:
    top = max(trunc_mult * sp.omega_max, 1.5 * omega_c)
    return max(64, int(math.ceil((top / (math.pi * T) - 1) / 2)) + 1)


def build_M(sp: Spectrum, T: float, mu: float = 0.0, omega_c: float = math.inf,
            N: int | None = None, trunc_mult: float = 30.0) -> np.ndarray:
    if N is None:
        N = n_trunc(sp, T, omega_c if math.isfinite(omega_c) else 0.0, trunc_mult)
    L = lam_k(sp, T, 2 * N)
    n = np.arange(N)
    D = np.abs(n[:, None] - n[None, :])
    M = L[D] + L[n[:, None] + n[None, :] + 1]
    np.fill_diagonal(M, L[2 * n + 1])
    S = np.concatenate([[0.0], np.cumsum(L[1:N])])  # S[n] = sum_{k=1}^n L(k)
    M[n, n] -= (2 * n + 1) + 2 * S[:N]
    if mu > 0:
        if math.isfinite(omega_c):
            Nc = int(np.sum((2 * n + 1) * math.pi * T <= omega_c))
            if (2 * N + 1) * math.pi * T <= omega_c:
                raise ValueError("truncation below cutoff")
        else:
            Nc = N
        M[:Nc, :Nc] -= 2 * mu
    return M


def rho(sp: Spectrum, T: float, mu: float = 0.0, omega_c: float = math.inf,
        N: int | None = None, trunc_mult: float = 30.0, vec: bool = False):
    M = build_M(sp, T, mu, omega_c, N, trunc_mult)
    m = M.shape[0]
    vals, vecs = eigh(M, subset_by_index=[m - 1, m - 1])
    if vec:
        return float(vals[0]), vecs[:, 0]
    return float(vals[0])


def tc(sp: Spectrum, mu: float = 0.0, omega_c: float = math.inf, trunc_mult: float = 30.0,
       T_lo: float | None = None, xtol_rel: float = 1e-7) -> float:
    """Largest T with rho(T) = 0. Scans down from T_hi = A/4.3 (universal bound, so rho < 0
    above it) in factors of 0.8 until rho >= 0, then brentq. Returns 0.0 if none above T_lo."""
    T_hi = sp.A / 4.3314814 * 1.0000001
    if T_lo is None:
        T_lo = 1e-4 * sp.omega_log
    f = lambda T: rho(sp, T, mu, omega_c, None, trunc_mult)
    T2 = T_hi
    T1 = T2 * 0.8
    while f(T1) < 0:
        T2 = T1
        T1 *= 0.8
        if T1 < T_lo:
            return 0.0
    return brentq(f, T1, T2, xtol=xtol_rel * T1, rtol=1e-12)


# ---------------------------------------------------------------- Bergmann-Rainer derivative
def functional_derivative(sp: Spectrum, Tc: float, omegas: np.ndarray, mu: float = 0.0,
                          omega_c: float = math.inf, trunc_mult: float = 30.0) -> np.ndarray:
    """delta Tc / delta a2F(omega) (dimensionless, per unit a2F mass) by Hellmann-Feynman:
    d rho = x^T (dM/d a2F(omega)) x, d rho/dT by central difference; dTc = -d rho/(d rho/dT)."""
    N = n_trunc(sp, Tc, omega_c if math.isfinite(omega_c) else 0.0, trunc_mult)
    r0, x = rho(sp, Tc, mu, omega_c, N, vec=True)
    h = 1e-5 * Tc
    drdT = (rho(sp, Tc + h, mu, omega_c, N) - rho(sp, Tc - h, mu, omega_c, N)) / (2 * h)
    out = []
    for om in omegas:
        unit = Spectrum(np.array([om]), np.array([1.0]))
        Mu = build_M(unit, Tc, 0.0, math.inf, N)
        n = np.arange(N)
        Mu[n, n] += (2 * n + 1)            # remove the bare omega_n term: only the a2F part
        out.append(-(x @ Mu @ x) / drdT)
    return np.array(out)


# ---------------------------------------------------------------- Allen-Dynes
def allen_dynes(lam: float, mu: float, omega_log: float, omega2: float,
                corrections: bool = True) -> float:
    """Allen & Dynes, PRB 12, 905 (1975), Eqs. as reproduced in standard reviews:
    Tc = f1 f2 (omega_log/1.2) exp[-1.04(1+lam)/(lam - mu(1+0.62 lam))],
    f1 = [1 + (lam/L1)^{3/2}]^{1/3}, L1 = 2.46(1+3.8 mu),
    f2 = 1 + (omega2/omega_log - 1) lam^2/(lam^2 + L2^2), L2 = 1.82(1+6.3 mu)(omega2/omega_log)."""
    den = lam - mu * (1 + 0.62 * lam)
    if den <= 0:
        return 0.0
    t = omega_log / 1.2 * math.exp(-1.04 * (1 + lam) / den)
    if not corrections:
        return t
    L1 = 2.46 * (1 + 3.8 * mu)
    L2 = 1.82 * (1 + 6.3 * mu) * (omega2 / omega_log)
    f1 = (1 + (lam / L1) ** 1.5) ** (1 / 3)
    f2 = 1 + (omega2 / omega_log - 1) * lam ** 2 / (lam ** 2 + L2 ** 2)
    return f1 * f2 * t


def mu_rescale(mu: float, omega_c_from: float, omega_c_to: float) -> float:
    """Tolmachev/Morel-Anderson renormalisation: 1/mu*(w2) = 1/mu*(w1) + ln(w1/w2)."""
    return 1.0 / (1.0 / mu + math.log(omega_c_from / omega_c_to))
