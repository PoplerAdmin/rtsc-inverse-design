"""rtt.core.certify: rigorous Eliashberg engine (team CERT).

Public API
----------
classify(spectra, params, material_id, band_couplings=None) -> Verdict
                                                        the ONLY producer of REJECTED_PROVEN /
                                                        CERTIFIED_CANDIDATE (ARCH S1)
recheck(verdict, spectra, band_couplings=None) -> dict  re-verify a stored certificate from scratch
CertifyParams.rejection_model = "isotropic"|"multiband" ARCH 4c (R-MULTI audited SOUND; opt-in)
rmulti_for / rmulti_evidence                            R-MULTI certificate / re-check
CertifyParams                                           run parameters (recorded in the Verdict)
AltEvidence, classifyAlt, classify_evidence             Lean PipelineLogic mirror
check_univ(spectrum, T_star) -> bool                    R-UNIV (exact rational)
certify_upper / certify_lower                           R-CELLS cover / C-RAYLEIGH search
check_cell / check_rayleigh                             Arb checks of single certificates
Proofs: rtt/core/certify/PROOFS.md.
"""
from .classify import (AltEvidence, CertifyParams, classify, classify_evidence, classifyAlt,
                       evidence_for_alternative, recheck, spectrum_digest)
from .arbcore import ArbSpectrum
from .floatcore import FloatSpectrum, float_tc
from .cells import check_cell, propose_cell
from .rayleigh import check_rayleigh, propose_rayleigh
from .search import certify_lower, certify_upper, t_end_for
from .multiband import rmulti_evidence, rmulti_for
from .classify import combine_multiband


def check_univ(spectrum, T_star: float) -> bool:
    """R-UNIV: (1+e) A < a0 T*, exact rational arithmetic (Lean R_UNIV)."""
    from .arbcore import q_of_float
    from .search import A0_Q
    asp = ArbSpectrum(spectrum.omega_K, spectrum.weight_K, spectrum.rel_envelope)
    return bool(asp.A_hi_q < A0_Q * q_of_float(float(T_star)))


__all__ = ["classify", "recheck", "CertifyParams", "AltEvidence", "classifyAlt",
           "classify_evidence", "evidence_for_alternative", "spectrum_digest", "check_univ",
           "certify_upper", "certify_lower", "check_cell", "check_rayleigh", "propose_cell",
           "propose_rayleigh", "t_end_for", "float_tc", "ArbSpectrum", "FloatSpectrum",
           "rmulti_for", "rmulti_evidence", "combine_multiband"]
