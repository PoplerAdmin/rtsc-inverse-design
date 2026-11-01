"""SYNTHETIC alpha^2F spectra for tests/validation (not data for any real compound)."""
from __future__ import annotations

import math

import numpy as np

from rtt.core.spectrum import Spectrum


def multipeak(seed: int, lam_target: float, n_pts: int = 600, om_lo: float = 20.0,
              om_hi: float = 1500.0, n_peaks: int | None = None, e: float = 0.0) -> Spectrum:
    """Random sum of Gaussian peaks on a grid (kelvin), scaled to lambda = lam_target."""
    rng = np.random.default_rng(seed)
    k = n_peaks or int(rng.integers(2, 6))
    w = np.linspace(om_lo * 0.5, om_hi, n_pts)
    a = np.zeros_like(w)
    for _ in range(k):
        c = math.exp(rng.uniform(math.log(om_lo), math.log(om_hi * 0.9)))
        s = c * rng.uniform(0.03, 0.25)
        a += rng.uniform(0.2, 1.0) * np.exp(-0.5 * ((w - c) / s) ** 2)
    a *= (w / w.max()) ** 0.5            # soft acoustic onset
    sp = Spectrum.from_grid(w, a, unit="K")
    lam = sp.lambda_()
    f = lam_target / lam
    return Spectrum(sp.omega_K, tuple(x * f for x in sp.weight_K), e,
                    label=f"synthetic-multipeak(seed={seed},lambda={lam_target})",
                    provenance={"synthetic": True, "seed": seed})


def einstein(Omega: float, A: float, e: float = 0.0) -> Spectrum:
    return Spectrum((float(Omega),), (float(A),), e, label=f"synthetic-einstein(Omega={Omega},A={A})",
                    provenance={"synthetic": True})


def float_tc_reference(spec: Spectrum, mu: float = 0.0, omega_c: float | None = None,
                       N: int = 600, scale: float = 1.0, T_lo: float = 1e-3) -> float:
    """Independent float reference: largest T with top eigenvalue >= 0 of the N x N even-sector
    matrix (block cutoff N_c(T) = #{n: (2n+1) pi T <= omega_c}), bisection in log T."""
    om = np.asarray(spec.omega_K)
    p = np.asarray(spec.weight_K) * scale
    n = np.arange(N)

    def rho(T):
        k = np.arange(2 * N + 2)
        x2 = (2 * np.pi * T * k[1:]) ** 2
        lam = np.zeros(2 * N + 2)
        for i0 in range(0, len(x2), 400):
            lam[1 + i0:1 + i0 + 400] = (2 * p * om / (om ** 2 + x2[i0:i0 + 400, None])).sum(1)
        M = lam[np.abs(n[:, None] - n[None, :])] + lam[n[:, None] + n[None, :] + 1]
        S = np.concatenate([[0], np.cumsum(lam[1:N])])
        M[n, n] = lam[2 * n + 1] - 2 * S - (2 * n + 1)
        if mu > 0:
            Nc = 0 if omega_c < np.pi * T else int((omega_c / (np.pi * T) - 1) // 2) + 1
            M[:Nc, :Nc] -= 2 * mu
        return np.linalg.eigvalsh(M)[-1]

    hi = p.sum() / 4.33148145 * 1.000001
    lo = hi
    while rho(lo) < 0:
        hi = lo
        lo *= 0.7
        if lo < T_lo:
            return 0.0
    for _ in range(60):
        mid = math.sqrt(lo * hi)
        if rho(mid) >= 0:
            lo = mid
        else:
            hi = mid
        if hi / lo - 1 < 1e-10:
            break
    return lo
