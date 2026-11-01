"""Alternatives for the core engine (architecture §3.2).

* ``freq_scaled``: for every admitted spectrum add the Hopfield-scaled copies Ω' = sΩ, p' = p/s
  for s in {1-ε_ω, 1+ε_ω} (doc 01 §5.6: A ∝ 1/ω, λ ∝ 1/ω², m invariant), so the frequency
  uncertainty enters as alternatives (ARCH errata F5).
* ``envelope_e``: e = (1+δ_A^q)(1+ε_sys) - 1, the weight envelope of every alternative.
  If e >= 1 no core rejection is possible (``Spectrum`` requires e < 1): returns None.
"""
from __future__ import annotations

from typing import Optional, Sequence

from ..core.spectrum import Spectrum


def freq_scaled(spectra: Sequence[Spectrum], eps_omega: float) -> list[Spectrum]:
    if not (0 <= eps_omega < 1):
        raise ValueError("eps_omega must lie in [0, 1)")
    out: list[Spectrum] = []
    for s in spectra:
        out.append(s)
        if eps_omega == 0:
            continue
        for tag, f in (("lo", 1.0 - eps_omega), ("hi", 1.0 + eps_omega)):
            prov = dict(s.provenance)
            prov["freq_scaled"] = {"factor": f, "eps_omega": eps_omega, "rule": "Hopfield A~1/omega"}
            out.append(Spectrum(tuple(o * f for o in s.omega_K), tuple(w / f for w in s.weight_K),
                                s.rel_envelope, f"{s.label}|omega*{f:g}", prov))
    return out


def envelope_e(delta_A_q: float, eps_sys: float) -> Optional[float]:
    if delta_A_q < 0 or eps_sys < 0:
        raise ValueError("negative envelope inputs")
    e = (1.0 + delta_A_q) * (1.0 + eps_sys) - 1.0
    return e if e < 1.0 else None


def with_envelope(spectra: Sequence[Spectrum], e: float) -> list[Spectrum]:
    return [s.with_envelope(e) for s in spectra]
