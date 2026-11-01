"""Certify post-processing policy: the only bridge from spectra to certified verdicts.

Contains ConfigError, validate_certify_config, _certify_params, certify_spectra, trust_base,
blocking_tags, _sigma_exclusion_directions, apply_input_policy, annotate_model,
cutoff_description, ensure_cutoff_assumption (architecture §5.2, review R12 #1).
The two marked blocks are frozen (sha256 in PINS.json, checked by tests/test_core_pins.py).
Do not edit them without re-deriving the soundness arguments and updating the pins.
"""
from __future__ import annotations

import copy  # noqa: F401  (used by _merge)
import dataclasses
import math
import re
from typing import Any, Sequence

from .spectrum import Spectrum
from .status import PROVEN_STATUSES, Status, Verdict  # noqa: F401

# ==== BEGIN FROZEN BLOCK L83-L133
class ConfigError(ValueError):
    pass


class AccountingError(AssertionError):
    """S4 violated: some input has no status, several, or an unknown id appeared."""


def _merge(base: dict, over: dict) -> dict:
    out = copy.deepcopy(base)
    for k, v in over.items():
        if isinstance(v, dict) and isinstance(out.get(k), dict):
            out[k] = _merge(out[k], v)
        else:
            out[k] = copy.deepcopy(v)
    return out


def _real(cc: dict[str, Any], key: str) -> float:
    v = cc.get(key)
    if isinstance(v, bool) or not isinstance(v, (int, float)):
        raise ConfigError(f"certify.{key} must be a real number, got {v!r}")
    v = float(v)
    if not math.isfinite(v):
        raise ConfigError(f"certify.{key} must be finite, got {v!r}")
    return v


def validate_certify_config(cc: dict[str, Any]) -> None:
    """[certify] values that change what a verdict claims are checked strictly (audit F4/F9)."""
    if cc.get("cutoff_rule", "block") != "block":
        raise ConfigError("certify.cutoff_rule: only 'block' (ARCH section 4 N_c rule) is defined")
    T = _real(cc, "T_target_K")
    if not T > 0:
        raise ConfigError(f"certify.T_target_K must be > 0, got {T}")
    mu = _real(cc, "mu_hi")
    if not mu >= 0:
        raise ConfigError(f"certify.mu_hi must be >= 0, got {mu}")
    if cc.get("omega_c_K") is not None:
        oc = _real(cc, "omega_c_K")
        if not oc >= math.pi * T:
            raise ConfigError(f"certify.omega_c_K = {oc} K < pi*T* = {math.pi * T:.6g} K: the "
                              "Coulomb block would be empty at T* (N_c(T*) = 0)")
    for k in ("allow_unverified_rejection", "allow_unverified_certification"):
        if k in cc and not isinstance(cc[k], bool):
            raise ConfigError(f"certify.{k} must be true or false")
    if "omega_c_factor" in cc:
        f = _real(cc, "omega_c_factor")
        if not f >= 1.0:
            raise ConfigError(f"certify.omega_c_factor must be >= 1 (omega_c >= max phonon "
                              f"frequency), got {f}")
# ==== END FROZEN BLOCK L83-L133

# ==== BEGIN FROZEN BLOCK L420-L650

# ----------------------------------------------------------------------------- certify
#: [certify] config key -> rtt.core.certify.CertifyParams field (other keys pass through unchanged)
_CERTIFY_KEYMAP = {"T_target_K": "T_star", "omega_c_K": "omega_c"}
_CERTIFY_ORCH_ONLY = {"cutoff_rule",           # validated here ("block" is the only rule)
                      "allow_unverified_rejection", "allow_unverified_certification"}


def _certify_params(ccfg: dict[str, Any]):
    """Build ``rtt.core.certify.CertifyParams`` from the ``[certify]`` table."""
    from .certify import CertifyParams
    names = {f.name for f in dataclasses.fields(CertifyParams)}
    kw = {}
    for k, v in ccfg.items():
        if k in _CERTIFY_ORCH_ONLY:
            continue
        name = _CERTIFY_KEYMAP.get(k, k)
        if name not in names:
            raise ConfigError(f"unknown [certify] key {k!r} (CertifyParams fields: {sorted(names)})")
        kw[name] = v
    return CertifyParams(**kw)


def certify_spectra(material_id: str, spectra: Sequence[Spectrum], ccfg: dict[str, Any],
                    band_couplings: Sequence[Any] | None = None) -> Verdict:
    """The single call into rtt.core.certify.classify (S1).  All alternatives are passed together;
    the S5 conjunction is certify's job.  Afterwards only assumptions and flags are added
    (cutoff rule, discretisation, single-band model, Coulomb-inactive); status and bounds are
    never touched."""
    validate_certify_config(ccfg)
    spectra = list(spectra)
    if not spectra:   # ARCH 4b F1: the S5 conjunction over zero alternatives is vacuous
        raise ValueError("refusing to classify zero alternative spectra (ARCH 4b F1)")
    from .certify import classify
    kw: dict[str, Any] = {"material_id": material_id}
    if band_couplings:
        kw["band_couplings"] = list(band_couplings)
    v = classify(spectra, _certify_params(ccfg), **kw)
    if not isinstance(v, Verdict):
        raise TypeError(f"certify.classify returned {type(v).__name__}, expected Verdict")
    if v.material_id != material_id:
        v = dataclasses.replace(v, material_id=material_id)
    v = ensure_cutoff_assumption(v, ccfg, spectra)
    v = annotate_model(v, ccfg, spectra)
    v = dataclasses.replace(v, assumptions=list(v.assumptions) + [trust_base()])
    return apply_input_policy(v, spectra, ccfg)


#: python-flint version the trusted base was verified with (correspondence audit G5)
VERIFIED_PYTHON_FLINT = "0.9.0"


def trust_base() -> str:
    """'TRUST: python-flint X / FLINT Y' for Verdict.assumptions; warns on an unverified version."""
    import warnings
    try:
        import flint
        pv, fv = getattr(flint, "__version__", "?"), getattr(flint, "__FLINT_VERSION__", "?")
    except ImportError:
        pv, fv = "missing", "?"
    if pv != VERIFIED_PYTHON_FLINT:
        warnings.warn(f"python-flint {pv} installed; the trusted base was verified with "
                      f"{VERIFIED_PYTHON_FLINT}", RuntimeWarning, stacklevel=2)
    return f"TRUST: python-flint {pv} / FLINT {fv}"


_TAG_MASS = re.compile(r"^A-(DROP-NONPOS|NEG-CLIP)\(([-+0-9.eE]+)K\)")


_EXCL_ITEM = re.compile(r"sigma(\d+):lam=([-+0-9.eE]+):A=([-+0-9.eE]+)K")


def _sigma_exclusion_directions(tag: str, spectra: Sequence[Spectrum]) -> tuple[bool, bool]:
    """(blocks REJECTED, blocks CERTIFIED) for an A-SIGMA-EXCLUDED tag.

    The excluded smearings violate alpha^2F_nu >= 0, so they are not admissible alternatives; the
    question is only whether they point outside what the admitted alternatives already cover. Each
    admitted alternative s is treated by the engine as the box (1 +- e_s) w_j, so its lambda and A
    range over [(1-e_s) X_s, (1+e_s) X_s]. An excluded smearing whose lambda or A lies ABOVE the
    union of these ranges blocks REJECTED (it might be harder to reject); BELOW blocks CERTIFIED.
    This is a moment-level comparison (lambda, A), not a pointwise one: inside the range the
    exclusion is recorded but not blocking. A malformed tag blocks both."""
    items = _EXCL_ITEM.findall(tag)
    if not items or not spectra:
        return True, True
    rng = {}
    for name, fn in (("lam", lambda s: s.lambda_()), ("A", lambda s: s.A())):
        rng[name] = (min((1 - s.rel_envelope) * fn(s) for s in spectra),
                     max((1 + s.rel_envelope) * fn(s) for s in spectra))
    above = below = False
    for _, lam, A in items:
        for name, x in (("lam", float(lam)), ("A", float(A))):
            above |= x > rng[name][1]
            below |= x < rng[name][0]
    return above, below


def blocking_tags(spectra: Sequence[Spectrum]) -> tuple[list[str], list[str]]:
    """Input-quality tags (team DFT, Spectrum.provenance['assumptions']) that make one
    direction of proof unsafe.  Returns (tags blocking REJECTED, tags blocking CERTIFIED):
      REJECTED : A-FILE-COMPLETE(unverified), A-FILE-INCOMPLETE..., A-DROP-NONPOS(m K) with m > 0
                 (missing weight lowers A, which is anti-conservative for an upper bound);
      CERTIFIED: A-NEG-CLIP(m K) with m > 0 (clipping raises A), A-FILE-INCOMPLETE...
      either   : A-SIGMA-EXCLUDED(...) (smearings excluded by smearing_admission='mode_positive'):
                 see _sigma_exclusion_directions
    """
    rej, cert = [], []
    seen_excl: set[str] = set()
    for s in spectra:
        for a in s.provenance.get("assumptions", []) or []:
            a = str(a)
            m = _TAG_MASS.match(a)
            mass = float(m.group(2)) if m else 0.0
            tag = f"{a} [{s.label}]"
            if a.startswith("A-FILE-INCOMPLETE"):
                rej.append(tag)
                cert.append(tag)
            elif a.startswith("A-FILE-COMPLETE(unverified)"):
                rej.append(tag)
            elif m and m.group(1) == "DROP-NONPOS" and mass > 0:
                rej.append(tag)
            elif m and m.group(1) == "NEG-CLIP" and mass > 0:
                cert.append(tag)
            elif a.startswith("A-SIGMA-EXCLUDED(") and a not in seen_excl:
                seen_excl.add(a)
                r, c = _sigma_exclusion_directions(a, spectra)
                if r:
                    rej.append(tag)
                if c:
                    cert.append(tag)
    return rej, cert


def apply_input_policy(v: Verdict, spectra: Sequence[Spectrum], ccfg: dict[str, Any]) -> Verdict:
    """Coordinator policy (S1-compatible: it can only DOWNGRADE).  A proven status whose
    direction is blocked by an input-quality tag on any alternative becomes UNDECIDED, unless
    [certify] allow_unverified_rejection / allow_unverified_certification is true (then the
    tag is kept as an assumption).  The certificate is kept unchanged; the untrusted bound of the
    blocked direction is dropped from the bracket."""
    rej, cert = blocking_tags(spectra)
    reasons, assumptions = list(v.reasons), list(v.assumptions)
    status, T_low, T_up = v.status, v.T_low_K, v.T_up_K
    for kind, tags, key in (("rejection", rej, "allow_unverified_rejection"),
                            ("certification", cert, "allow_unverified_certification")):
        hit = (status is Status.REJECTED_PROVEN and kind == "rejection") or \
              (status is Status.CERTIFIED_CANDIDATE and kind == "certification")
        if not (hit and tags):
            continue
        if ccfg.get(key, False):
            assumptions.extend(f"UNVERIFIED-INPUT-ACCEPTED({key}): {t}" for t in tags)
            continue
        reasons.append(f"certify.classify returned {status.value}; {kind} blocked: "
                       + "; ".join(tags))
        if kind == "rejection":
            T_up = None
        else:
            T_low = None
        status = Status.UNDECIDED
    if status is v.status and reasons == v.reasons and assumptions == v.assumptions:
        return v
    out = dataclasses.replace(v, status=status, T_low_K=T_low, T_up_K=T_up, reasons=reasons,
                              assumptions=assumptions)
    assert out.status is v.status or out.status is Status.UNDECIDED   # downgrade only
    assert out.certificate == v.certificate
    return out


def annotate_model(v: Verdict, ccfg: dict[str, Any], spectra: Sequence[Spectrum]) -> Verdict:
    """Add the model assumptions the report must state (audit low items, F9).  Only
    ``assumptions`` and ``flags`` change."""
    assumptions = list(v.assumptions)
    flags = list(v.flags)
    seen = set()
    for s in spectra:
        rule = s.provenance.get("rule") or s.provenance.get("discretisation")
        src = s.provenance.get("source_file") or s.provenance.get("attached_file") or \
            ("inline/synthetic" if s.provenance.get("synthetic") else "inline")
        if rule:
            a = f"DISCRETISATION({str(rule).split(':')[0]}, {src})"
            if a not in seen:
                seen.add(a)
                assumptions.append(a)
    for s in spectra:                       # e.g. team DFT's A-FILE-COMPLETE(unverified)
        for a in s.provenance.get("assumptions", []) or []:
            if a not in assumptions:
                assumptions.append(a)
    has_sb = any(a.startswith("M-ISO-SINGLEBAND") for a in assumptions)
    if v.status is Status.REJECTED_PROVEN and "R-MULTI" not in v.rule_ids and not has_sb:
        assumptions.append("M-ISO-SINGLEBAND (rejection uses the DOS-averaged single-band A)")
    elif v.status in (Status.UNDECIDED, Status.CERTIFIED_CANDIDATE) and not has_sb and \
            (v.T_low_K is not None or v.T_up_K is not None):
        # Lean note MB-N1: the bracket [T_low, T_up] is the single-band (DOS-averaged) one,
        # also when rejection_model = "multiband"
        assumptions.append("M-ISO-SINGLEBAND (the bracket [T_low, T_up] is the single-band "
                           "DOS-averaged one)")
    # F9: is the Coulomb block empty at the certified lower bound?  N_c(T) = 0 iff pi*T > omega_c
    oc = (v.certificate or {}).get("params", {}).get("omega_c")
    if oc is None:
        oc = (float(ccfg["omega_c_K"]) if ccfg.get("omega_c_K") is not None else
              float(ccfg.get("omega_c_factor", 10.0)) * max(s.omega_max() for s in spectra))
    if v.T_low_K is not None and math.pi * float(v.T_low_K) > float(oc) and \
            "COULOMB_INACTIVE" not in flags:
        flags.append("COULOMB_INACTIVE")
        assumptions.append(f"COULOMB_INACTIVE: N_c(T_low={v.T_low_K} K) = 0 for omega_c = {oc} K; "
                           "the lower bound does not depend on mu*")
    out = dataclasses.replace(v, assumptions=assumptions, flags=flags)
    assert (out.status, out.T_low_K, out.T_up_K, out.certificate) == \
        (v.status, v.T_low_K, v.T_up_K, v.certificate)
    return out


def cutoff_description(ccfg: dict[str, Any], spectra: Sequence[Spectrum] = ()) -> str:
    """ARCH 4b F3: the Coulomb block-cutoff rule and omega_c used (config [certify])."""
    rule = ccfg.get("cutoff_rule", "block")
    if ccfg.get("omega_c_K") is not None:
        oc = f"omega_c={float(ccfg['omega_c_K']):g} K (fixed)"
    else:
        f = float(ccfg.get("omega_c_factor", 10.0))
        vals = sorted({round(f * s.omega_max(), 6) for s in spectra})
        oc = f"omega_c={f:g}*max(omega) per spectrum" + (f" = {vals} K" if vals else "")
    return f"CUTOFF({rule}: N_c(T)=#{{n>=0 : (2n+1)*pi*T <= omega_c}}; {oc})"


def ensure_cutoff_assumption(v: Verdict, ccfg: dict[str, Any], spectra: Sequence[Spectrum]) -> Verdict:
    """Make sure Verdict.assumptions names the cutoff rule and omega_c (ARCH 4b F3).  Only the
    assumptions list is extended; status, bounds and certificate are untouched."""
    if any("omega_c" in a for a in v.assumptions):
        return v
    out = dataclasses.replace(v, assumptions=list(v.assumptions) + [cutoff_description(ccfg, spectra)])
    assert (out.status, out.T_low_K, out.T_up_K) == (v.status, v.T_low_K, v.T_up_K)
    return out
# ==== END FROZEN BLOCK L420-L650

