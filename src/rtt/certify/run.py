"""The only entry into the core certify engine (contract R1b, architecture §3.2, §5.2).

``certify`` calls ``core.certify_policy.certify_spectra`` and nothing else, after forcing
``allow_unverified_rejection`` / ``allow_unverified_certification`` to false so that the core
policy and Lean ``decideT`` (which has no override) agree.
"""
from __future__ import annotations

from typing import Any, Sequence

from ..core.certify_policy import certify_spectra, validate_certify_config
from ..core.spectrum import Spectrum
from ..core.status import Verdict

ORIGIN = "rtt.core.certify_policy.certify_spectra"
FORCED_FALSE = ("allow_unverified_rejection", "allow_unverified_certification")
MODES = ("decide", "bracket")


def certify_config(T_star_K: float, mu_hi: float, mode: str = "bracket", omega_c_factor: float = 10.0,
                   omega_c_K: float | None = None, **extra: Any) -> dict:
    """Build a core ``[certify]`` table."""
    cc = {"T_target_K": float(T_star_K), "mu_hi": float(mu_hi), "omega_c_factor": float(omega_c_factor),
          "cutoff_rule": "block", "mode": mode}
    if omega_c_K is not None:
        cc["omega_c_K"] = float(omega_c_K)
    cc.update(extra)
    return cc


def certify(material_id: str, spectra: Sequence[Spectrum], ccfg: dict, mode: str | None = None,
            band_couplings: Sequence[Any] | None = None) -> Verdict:
    if mode is not None:
        if mode not in MODES:
            raise ValueError(f"mode must be one of {MODES}")
        ccfg = {**ccfg, "mode": mode}
    for k in FORCED_FALSE:
        if ccfg.get(k):
            raise ValueError(f"[certify] {k} = true is not allowed in rtt (R1b/R1c)")
    ccfg = {**ccfg, **{k: False for k in FORCED_FALSE}}
    validate_certify_config(ccfg)
    return certify_spectra(material_id, list(spectra), ccfg, band_couplings=band_couplings)
