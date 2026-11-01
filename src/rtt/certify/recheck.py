"""Independent re-verification of a core verdict before it may enter the ledger (R4)."""
from __future__ import annotations

from typing import Any, Sequence

from ..core.certify import CertifyParams, recheck, spectrum_digest  # noqa: F401
from ..core.certify_policy import _certify_params
from ..core.spectrum import Spectrum
from ..core.status import Verdict


class RecheckError(RuntimeError):
    pass


def recheck_or_raise(verdict: Verdict, spectra: Sequence[Spectrum], ccfg: dict | None = None,
                     band_couplings: Sequence[Any] | None = None) -> dict:
    """Run core ``recheck``; with ``ccfg`` the run parameters are bound too."""
    exp = None
    if ccfg is not None:
        exp = _certify_params(dict(ccfg))
    res = recheck(verdict, list(spectra), band_couplings=list(band_couplings or []) or None,
                  expected_params=exp)
    if not res.get("ok"):
        raise RecheckError("; ".join(res.get("errors", [])) or "recheck failed")
    return res


def digests(spectra: Sequence[Spectrum]) -> list[str]:
    return [spectrum_digest(s) for s in spectra]
