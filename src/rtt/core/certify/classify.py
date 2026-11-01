"""The single decision function (ARCH S1, section 4; ERRATA 4b) and the standalone re-checker.

Control flow mirrors Lean PipelineLogic (Classify.lean):

    structure AltEvidence where univFires, cellsBelowTarget, rayleighAtOrAboveTarget : Bool
    def classifyAlt ev := if ev.univFires then rejectedProven
                          else if ev.cellsBelowTarget then rejectedProven
                          else if ev.rayleighAtOrAboveTarget then certifiedCandidate
                          else undecided
    def classify evs := if evs.isEmpty then undecided
                        else if evs.all (classifyAlt . == rejectedProven) then rejectedProven
                        else if evs.all (classifyAlt . == certifiedCandidate) then certifiedCandidate
                        else undecided

Every boolean of AltEvidence is produced by ``evidence_for_alternative`` from a stored
certificate, using only exact rational comparisons and Arb ball checks (S1).  ``recheck``
calls the very same function, so a stored Verdict can be re-verified from scratch.
"""
from __future__ import annotations

import hashlib
import json
import math
import time
from dataclasses import asdict, dataclass, field

from flint import fmpq

from ..constants import A0_ETA_FRACTION, A0_FRACTION
from ..spectrum import Spectrum
from ..status import Status, Verdict
from .arbcore import ArbSpectrum, q_of_float
from .cells import VARIANTS, LamCache, check_cell
from .floatcore import FloatSpectrum
from .rayleigh import check_rayleigh
from .search import A0_Q, certify_lower, certify_upper, t_end_for

CERT_VERSION = "rtsc-certify-1"

LEAN = {
    "R-UNIV": ["PipelineLogic.R_UNIV", "PipelineLogic.univ_stable", "PipelineLogic.stable_of_A_le",
               "TcR1.stable_of_ge"],
    "R-CELLS": ["PipelineLogic.R_CELLS", "PipelineLogic.cells_stable", "PipelineLogic.cell_stable",
                "PipelineLogic.envelope_form_bound", "PipelineLogic.schur_row_test_inf",
                "PipelineLogic.rows_of_eventually", "PipelineLogic.lam_box", "PipelineLogic.cover"],
    "R-MU": ["PipelineLogic.R_MU", "PipelineLogic.stable_of_even_bound"],
    "R-ODD": ["TcR1.L6 via PipelineLogic.stable_of_even_bound"],
    "C-RAYLEIGH": ["PipelineLogic.rayleigh_cert", "PipelineLogic.boxMin_le",
                   "PipelineLogic.boxMin_attained", "PipelineLogic.sel_mul_eq_min",
                   "PipelineLogic.Qp_anti_mu"],
    "CLASSIFY": ["PipelineLogic.classify_rejected_sound", "PipelineLogic.classify_certified_sound",
                 "PipelineLogic.classify_eq_rejected_iff", "PipelineLogic.classify_eq_certified_iff",
                 "PipelineLogic.classify_cases"],
}


# ============================================================================ Lean mirror
@dataclass(frozen=True)
class AltEvidence:
    univFires: bool
    cellsBelowTarget: bool
    rayleighAtOrAboveTarget: bool


def classifyAlt(ev: AltEvidence) -> Status:
    if ev.univFires:
        return Status.REJECTED_PROVEN
    elif ev.cellsBelowTarget:
        return Status.REJECTED_PROVEN
    elif ev.rayleighAtOrAboveTarget:
        return Status.CERTIFIED_CANDIDATE
    else:
        return Status.UNDECIDED


def classify_evidence(evs: list[AltEvidence]) -> Status:
    if len(evs) == 0:
        return Status.UNDECIDED                                  # ERRATA F1
    elif all(classifyAlt(ev) == Status.REJECTED_PROVEN for ev in evs):
        return Status.REJECTED_PROVEN
    elif all(classifyAlt(ev) == Status.CERTIFIED_CANDIDATE for ev in evs):
        return Status.CERTIFIED_CANDIDATE
    else:
        return Status.UNDECIDED


# ============================================================================ parameters
@dataclass
class CertifyParams:
    T_star: float = 300.0          # target temperature T* (K)
    mu_hi: float = 0.15            # mu* range [0, mu_hi] (CERTIFIED needs all of it)
    omega_c: float | None = None   # block-cutoff frequency (K); None -> factor * max omega
    omega_c_factor: float = 10.0
    row_check: str = "mhat"        # "mhat" (Lean CellOK) or "cform" (PROOFS.md C2; not in Lean)
    rel_tol: float = 1e-4          # target relative gap of T_up/T_low to the float Tc
    mode: str = "bracket"          # "bracket": best brackets; "decide": stop once decided
    N_cap: int = 2400              # max head size for cells
    time_budget_s: float = 120.0   # per spectrum and per search (refinement stops, never the decision)
    verbose: bool = False
    rejection_model: str = "isotropic"   # "isotropic" | "multiband" (ARCH 4c; R-MULTI audited SOUND, opt-in)

    def resolved_omega_c(self, spectra: list[Spectrum]) -> float:
        if self.omega_c is not None:
            return float(self.omega_c)
        if not spectra:
            return float("nan")
        return float(self.omega_c_factor) * max(max(s.omega_K) for s in spectra)


def spectrum_digest(s: Spectrum) -> str:
    payload = json.dumps({"omega_K": [float(x).hex() for x in s.omega_K],
                          "weight_K": [float(x).hex() for x in s.weight_K],
                          "rel_envelope": float(s.rel_envelope).hex()}, sort_keys=True)
    return hashlib.sha256(payload.encode()).hexdigest()


# ============================================================================ evidence
def _cert_key(c: dict) -> str:
    d = {k: v for k, v in c.items() if k not in ("check", "float_rho_N", "float_rho")}
    return hashlib.sha256(json.dumps(d, sort_keys=True).encode()).hexdigest()


def evidence_for_alternative(spec: Spectrum, alt: dict, T_star: float, mu_hi: float,
                             omega_c: float, verified: set | None = None,
                             row_check: str = "mhat") -> tuple[AltEvidence, dict]:
    """Re-derive AltEvidence from a stored per-alternative certificate.  Only exact rational
    comparisons and Arb checks.  Returns (evidence, info) where info has the certified bracket
    of this alternative and a log of every check."""
    info: dict = {"checks": []}
    asp = ArbSpectrum(spec.omega_K, spec.weight_K, spec.rel_envelope)
    if alt.get("spectrum_sha256") != spectrum_digest(spec):
        info["checks"].append("spectrum digest MISMATCH")
        return AltEvidence(False, False, False), info
    Ts_q = q_of_float(T_star)
    # ---- R-UNIV: A_hi < a0 T*   (A_hi = (1+e) A exactly, rational)
    univ = bool(asp.A_hi_q < A0_Q * Ts_q)
    info["checks"].append(f"R-UNIV: A_hi={int(asp.A_hi_q.p) / int(asp.A_hi_q.q):.10g} K, "
                          f"a0*T*={int((A0_Q * Ts_q).p) / int((A0_Q * Ts_q).q):.10g} K -> {univ}")
    # ---- R-CELLS
    cells_ok = False
    T_up_cert = None
    cblk = alt.get("cells")
    if cblk is not None:
        try:
            t_end = float.fromhex(cblk["t_end"])
            T_up = float.fromhex(cblk["T_up"])
            cells = cblk["cells"]
            # D1: A_hi <= a0 * t_end (exact)
            d1 = bool(A0_Q * q_of_float(t_end) >= asp.A_hi_q) and t_end > 0
            info["checks"].append(f"D1 A_hi <= a0*t_end: {d1}")
            if d1:
                T_up_cert = t_end           # R-UNIV form: stable for T >= t_end (univ_stable)
            if len(cells) == 0:
                # Lean: Tiles _ _ [] = False -> no R-CELLS evidence; t_end stays a valid bound
                tiles = False
                all_ok = False
                info["checks"].append("no cells (bound t_end via univ_stable only)")
            else:
                # Tiles T_up t_end cells (exact equality of doubles)
                ends = [(float.fromhex(c["T1"]), float.fromhex(c["T2"])) for c in cells]
                tiles = (len(ends) > 0 and ends[0][0] == T_up and ends[-1][1] == t_end
                         and all(ends[i][1] == ends[i + 1][0] for i in range(len(ends) - 1)))
                info["checks"].append(f"Tiles(T_up, t_end, {len(ends)} cells): {tiles}")
                all_ok = tiles
                if tiles:
                    cache = LamCache(asp)
                    for i, c in enumerate(cells):
                        if c.get("variant") not in VARIANTS or c.get("variant") != row_check:
                            # G8: only the run's declared row check is accepted
                            info["checks"].append(f"cell {i}: variant {c.get('variant')!r} != "
                                                  f"row_check {row_check!r}")
                            all_ok = False
                            break
                        if verified is not None and ("cell", alt.get("spectrum_sha256"), _cert_key(c)) in verified:
                            r = {"ok": True}      # same certificate already passed check_cell
                        else:
                            r = check_cell(asp, c, cache)
                        if not r["ok"]:
                            info["checks"].append(f"cell {i} FAILED: {r.get('reason')}")
                            all_ok = False
                            break
                    info["checks"].append(f"all {len(cells)} cells certified: {all_ok}")
            if d1 and tiles and all_ok:
                T_up_cert = T_up
                cells_ok = True
        except Exception as exc:
            info["checks"].append(f"cells block malformed: {exc!r}")
    cells_below = bool(cells_ok and T_up_cert is not None and T_up_cert < T_star)
    info["T_up"] = T_up_cert
    info["cells_valid"] = cells_ok
    # ---- C-RAYLEIGH
    ray_ok = False
    T_low = None
    rc = alt.get("rayleigh")
    if rc is not None:
        try:
            same = (float(rc["mu_hi"]) == float(mu_hi) and float(rc["omega_c"]) == float(omega_c))
            if not same:
                info["checks"].append("rayleigh: mu_hi/omega_c differ from the run parameters")
            else:
                if verified is not None and ("ray", alt.get("spectrum_sha256"), _cert_key(rc)) in verified:
                    r = {"ok": True, "Q_over_x2": rc.get("check", {}).get("Q_over_x2")}
                else:
                    r = check_rayleigh(asp, rc)
                info["checks"].append(f"C-RAYLEIGH at T={float.fromhex(rc['T']):.10g}: {r['ok']} "
                                      f"(Q/|x|^2={r.get('Q_over_x2')})")
                if r["ok"]:
                    ray_ok = True
                    T_low = float.fromhex(rc["T"])
        except Exception as exc:
            info["checks"].append(f"rayleigh block malformed: {exc!r}")
    ray_above = bool(ray_ok and T_low >= T_star)
    info["T_low"] = T_low
    ev = AltEvidence(univFires=univ, cellsBelowTarget=cells_below, rayleighAtOrAboveTarget=ray_above)
    info["evidence"] = asdict(ev)
    info["classifyAlt"] = classifyAlt(ev).value
    return ev, info


# ============================================================================ classify
def _combine_bracket(infos: list[dict]):
    """ERRATA F7: [min T_low, max T_up]; None if any alternative lacks a bound."""
    lows = [i.get("T_low") for i in infos]
    ups = [i.get("T_up") for i in infos]
    T_low = min(lows) if lows and all(x is not None for x in lows) else None
    T_up = max(ups) if ups and all(x is not None for x in ups) else None
    return T_low, T_up


def _rmulti_run(bcs: list, T_star: float):
    """R-MULTI (multiband.py) on every band-coupling alternative; evidence re-derived from the
    stored blocks by rmulti_evidence (exact + Arb only)."""
    from . import multiband as mb
    blocks, fires, infos = [], [], []
    for bc in bcs:
        blk = mb.rmulti_for(bc, T_star)
        f, inf = mb.rmulti_evidence(blk, T_star)
        blk["evidence"] = inf
        blocks.append(blk)
        fires.append(f)
        infos.append(inf)
    return blocks, fires, infos, mb.LEAN_MB


def _rmulti_recheck(blocks: list, bcs: list, T_star: float):
    from . import multiband as mb
    if len(bcs) != len(blocks):
        return None, None, f"{len(blocks)} certified band couplings but {len(bcs)} given"
    fires, infos = [], []
    for bc, blk in zip(bcs, blocks):
        fresh = mb.bc_payload(bc)
        if mb.bc_digest(fresh) != blk.get("digest"):
            return None, None, "band-coupling data differ from the certified data"
        # referee R2: re-derive with the payload rebuilt from the SUPPLIED data
        blk = dict(blk, payload=fresh)
        f, inf = mb.rmulti_evidence(blk, T_star)
        fires.append(f)
        infos.append(inf)
    return fires, infos, None


def combine_multiband(status_sb: Status, mb_fires: list[bool], rejection_model: str,
                      n_spectra: int | None = None) -> Status:
    """ARCH 4c.  "isotropic": the single-band status (Lean classify) unchanged.
    "multiband": a rejection needs R-MULTI on every band-coupling alternative (non-empty list);
    a single-band rejection alone gives UNDECIDED; CERTIFIED (single-band lower bound) stays
    CERTIFIED only if R-MULTI does not fire; any conflict gives UNDECIDED."""
    if rejection_model == "isotropic":
        return status_sb
    if n_spectra is not None and n_spectra == 0:
        return Status.UNDECIDED                     # LOW-2: proven status needs spectra (F1)
    mb_rej = len(mb_fires) > 0 and all(mb_fires)
    if mb_rej and status_sb != Status.CERTIFIED_CANDIDATE:
        return Status.REJECTED_PROVEN
    if status_sb == Status.CERTIFIED_CANDIDATE and not mb_rej:
        return Status.CERTIFIED_CANDIDATE
    return Status.UNDECIDED


def _prm_dict(params: "CertifyParams", omega_c: float) -> dict:
    return {"T_star": float(params.T_star), "mu_hi": float(params.mu_hi), "omega_c": float(omega_c),
            "row_check": params.row_check, "rejection_model": params.rejection_model}


#: assumption prefixes whose exact text is derived from the run parameters (recheck: equality)
BOUND_PREFIXES = ("TARGET(", "RUN-PARAMS(", "U-MU[", "CUTOFF-BLOCK(")


def _annotations(spectra, prm: dict, status: Status, evs, infos, mb_infos, n_mb: int,
                 T_low, T_up):
    """Deterministic rule_ids / assumptions / flags (and the final bracket) from the run
    parameters and the re-derived evidence.  Used by classify AND recheck (audit MEDIUM-1, LOW-1)."""
    rule_ids: list[str] = []
    flags: list[str] = []
    if any(e.univFires for e in evs):
        rule_ids.append("R-UNIV")
    if any(i.get("cells_valid") for i in infos):
        rule_ids += ["R-CELLS", "R-MU", "R-ODD"]
    if any(i.get("T_low") is not None for i in infos):
        rule_ids.append("C-RAYLEIGH")
    rule_ids.append("CLASSIFY-S5")
    T_star, mu_hi, omega_c = prm["T_star"], prm["mu_hi"], prm["omega_c"]
    assumptions = ["M-ISO-MIGDAL (isotropic Migdal-Eliashberg, linearised gap equation, given a2F)",
                   f"H(a0={float(A0_FRACTION)}, eta={float(A0_ETA_FRACTION)}) [PROOF.md Prop. 6, wE1]",
                   "TRUST: Arb/FLINT (python-flint 0.9.0) ball arithmetic and exact conversions",
                   f"TARGET(T*={T_star!r} K)",
                   f"RUN-PARAMS(T*={T_star!r}, mu_hi={mu_hi!r}, omega_c={omega_c!r}, "
                   f"row_check={prm['row_check']}, rejection_model={prm['rejection_model']})"]
    for spec in spectra:
        assumptions.append(f"E-ENV({spec.rel_envelope}) weights only (ERRATA F5) [{spec.label}]")
    assumptions.append(f"U-MU[0,{mu_hi}]")
    if status == Status.CERTIFIED_CANDIDATE or T_low is not None:
        assumptions.append(f"CUTOFF-BLOCK(N_c(T)=#{{n>=0:(2n+1)piT<=omega_c}}, omega_c={omega_c} K) "
                           "(ERRATA F3: lower bounds / CERTIFIED hold for this cutoff rule)")
    if status == Status.REJECTED_PROVEN:
        assumptions.append("REJECTED holds for every mu*>=0 and every cutoff (Lean D2)")
    if prm["row_check"] == "cform":
        assumptions.append("R-CELLS[cform]: relies on PROOFS.md Thm C2 (hand proof, NOT in Lean)")
    if prm["rejection_model"] == "multiband":
        if status == Status.REJECTED_PROVEN:
            rule_ids = ["R-MULTI", "CLASSIFY-S5"]
            for inf in mb_infos:
                for a in inf.get("assumptions", []):
                    if a not in assumptions:
                        assumptions.append(a)
            if n_mb != len(spectra):
                assumptions.append(f"BAND-DATA-NOT-PER-ALTERNATIVE ({n_mb} band couplings for "
                                   f"{len(spectra)} spectra; every band coupling must fire)")
            T_up = max(inf["T_bound"] for inf in mb_infos)
            T_low = None
            assumptions = [a for a in assumptions if not a.startswith("REJECTED holds")]
            assumptions.append("REJECTED (multiband model) holds for every PSD N_i mu*_ij and every cutoff")
        elif T_low is not None or T_up is not None:
            assumptions.append("M-ISO-SINGLEBAND (bracket / certification from the single-band "
                               "DOS-averaged model)")
    elif status == Status.REJECTED_PROVEN:
        assumptions.append("M-ISO-SINGLEBAND (rejection uses the DOS-averaged single-band A)")
        flags.append("MULTIBAND_RISK")
    return rule_ids, assumptions, flags, T_low, T_up


def classify(spectra: list[Spectrum], params: CertifyParams | None = None,
             material_id: str = "", log=None, band_couplings: list | None = None) -> Verdict:
    """ARCH section 4 classify over alternative spectra (S5 conjunction, ERRATA F1-F7)."""
    params = params or CertifyParams()
    if params.row_check not in VARIANTS:
        raise ValueError(f"row_check must be one of {VARIANTS}")
    if log is None and params.verbose:
        log = print
    t0 = time.time()
    omega_c = params.resolved_omega_c(spectra)
    T_star = float(params.T_star)
    mu_hi = float(params.mu_hi)
    alts = []
    for idx, spec in enumerate(spectra):
        ta = time.time()
        asp = ArbSpectrum(spec.omega_K, spec.weight_K, spec.rel_envelope)
        fs = FloatSpectrum(spec.omega_K, spec.weight_K, spec.rel_envelope)
        alt = {"index": idx, "label": spec.label, "spectrum_sha256": spectrum_digest(spec),
               "rel_envelope": spec.rel_envelope,
               "A_hi": f"{asp.A_hi_q.p}/{asp.A_hi_q.q}"}
        univ = bool(asp.A_hi_q < A0_Q * q_of_float(T_star))
        decide = params.mode == "decide"
        if log:
            log(f"[{material_id}] alt {idx} ({spec.label}): {spec.n_atoms()} atoms, "
                f"A_hi={spec.A_hi():.6g} K, R-UNIV={univ}")
        # ---- R-CELLS (upper bound)
        if decide and univ:
            alt["cells"] = {"t_end": t_end_for(asp).hex(), "T_up": t_end_for(asp).hex(), "cells": []}
        else:
            up = certify_upper(asp, fs, params.row_check, params.rel_tol, T_star=T_star,
                               refine=not decide, N_cap=params.N_cap,
                               time_budget=params.time_budget_s, log=log)
            alt["cells"] = {"t_end": up["t_end"].hex(), "T_up": up["T_up"].hex(),
                            "cells": up["cells"]}
            alt["search_upper"] = {"tc_float_hi": up["tc_float_hi"], "seconds": up["seconds"],
                                   "attempts": up["attempts"]}
        # ---- C-RAYLEIGH (lower bound)
        T_up_f = float.fromhex(alt["cells"]["T_up"])
        if not (decide and (univ or T_up_f < T_star)):
            lo = certify_lower(asp, fs, mu_hi, omega_c, params.rel_tol,
                               T_floor=T_star if decide else None, T_hi=T_up_f,
                               time_budget=params.time_budget_s, log=log)
            alt["rayleigh"] = lo["cert"]
            alt["search_lower"] = {"tc_float_lo": lo["tc_float_lo"], "seconds": lo["seconds"],
                                   "attempts": lo["attempts"]}
        alt["seconds"] = time.time() - ta
        alts.append(alt)
    # ---- evidence re-derived from the stored certificates (the only input to the decision)
    # certificates accepted by the search passed check_cell / check_rayleigh on exactly these
    # data (search.py only keeps certificates whose Arb check returned ok); keyed by content hash
    verified = set()
    for alt in alts:
        for c in alt["cells"]["cells"]:
            verified.add(("cell", alt["spectrum_sha256"], _cert_key(c)))
        if alt.get("rayleigh"):
            verified.add(("ray", alt["spectrum_sha256"], _cert_key(alt["rayleigh"])))
    evs, infos = [], []
    for spec, alt in zip(spectra, alts):
        ev, info = evidence_for_alternative(spec, alt, T_star, mu_hi, omega_c, verified,
                                            row_check=params.row_check)
        alt["evidence"] = info["evidence"]
        alt["classifyAlt"] = info["classifyAlt"]
        alt["checks"] = info["checks"]
        evs.append(ev)
        infos.append(info)
    status_sb = classify_evidence(evs)
    T_low, T_up = _combine_bracket(infos)
    # ---- R-MULTI (ARCH 4c), evidence re-derived from stored blocks
    if params.rejection_model not in ("isotropic", "multiband"):
        raise ValueError("rejection_model must be 'isotropic' or 'multiband'")
    mb_blocks, mb_fires, mb_infos, LEAN_MB = _rmulti_run(band_couplings or [], T_star)
    status = combine_multiband(status_sb, mb_fires, params.rejection_model, n_spectra=len(spectra))
    reasons = []
    if not spectra:
        reasons.append("no alternative spectra (ERRATA F1): UNDECIDED")
    for spec, alt, info in zip(spectra, alts, infos):
        reasons.append(f"alt {alt['index']} ({spec.label or 'unlabelled'}): {alt['classifyAlt']}; "
                       f"bracket [{info['T_low']}, {info['T_up']}] K; evidence {info['evidence']}")
    for k, (blk, inf) in enumerate(zip(mb_blocks, mb_infos)):
        reasons.append(f"band coupling {k} ({blk['payload'].get('label') or blk['payload']['kind']}): "
                       f"R-MULTI kappa_hi={inf.get('kappa_hi')} K, fires={inf.get('fires')} "
                       f"({inf.get('reason') or (blk['route'] or {}).get('route')})")
    if params.rejection_model == "multiband":
        reasons.append(f"rejection_model=multiband: single-band status {status_sb.value}, "
                       f"R-MULTI over {len(mb_blocks)} band-coupling alternative(s) -> {status.value}")
        if not mb_blocks:
            reasons.append("no band data: multiband rejection impossible (ARCH 4c)")
    prm = _prm_dict(params, omega_c)
    rule_ids, assumptions, flags, T_low, T_up = _annotations(
        spectra, prm, status, evs, infos, mb_infos, len(mb_blocks), T_low, T_up)
    certificate = {
        "version": CERT_VERSION,
        "a0": f"{A0_FRACTION.numerator}/{A0_FRACTION.denominator}",
        "params": {"T_star": T_star, "mu_hi": mu_hi, "omega_c": omega_c,
                   "row_check": params.row_check, "rel_tol": params.rel_tol, "mode": params.mode},
        "alternatives": alts,
        "lean": dict(LEAN, **{"R-MULTI": LEAN_MB}),
        "rejection_model": params.rejection_model,
        "multiband": mb_blocks,
        "status_single_band": status_sb.value,
        "seconds": time.time() - t0,
    }
    return Verdict(material_id=material_id, status=status, T_low_K=T_low, T_up_K=T_up,
                   reasons=reasons, rule_ids=rule_ids, assumptions=assumptions, flags=flags,
                   certificate=certificate)


# ============================================================================ recheck
def recheck(verdict: Verdict, spectra: list[Spectrum], band_couplings: list | None = None,
            expected_params: "CertifyParams | None" = None) -> dict:
    """Re-verify a stored Verdict from scratch (Arb + exact arithmetic only; no float search).
    ok == True iff the status and the bracket re-derived from the certificate equal the stored
    ones and every proven claim is backed."""
    t0 = time.time()
    cert = verdict.certificate or {}
    out = {"ok": False, "errors": [], "alternatives": []}
    if cert.get("version") != CERT_VERSION:
        out["errors"].append("unknown certificate version")
        return out
    if cert.get("a0") != f"{A0_FRACTION.numerator}/{A0_FRACTION.denominator}":
        out["errors"].append("a0 mismatch")
        return out
    prm = cert["params"]
    alts = cert.get("alternatives", [])
    if len(alts) != len(spectra):
        out["errors"].append(f"{len(alts)} certified alternatives but {len(spectra)} spectra given")
        return out
    evs, infos = [], []
    for spec, alt in zip(spectra, alts):
        ev, info = evidence_for_alternative(spec, alt, float(prm["T_star"]), float(prm["mu_hi"]),
                                            float(prm["omega_c"]),
                                            row_check=prm.get("row_check", "mhat"))
        evs.append(ev)
        infos.append(info)
        out["alternatives"].append(info)
    status_sb = classify_evidence(evs)
    T_low, T_up = _combine_bracket(infos)
    blocks = cert.get("multiband", []) or []
    given = band_couplings or []
    if blocks or given:
        mb_fires, mb_infos, err = _rmulti_recheck(blocks, given, float(prm["T_star"]))
        if err:
            out["errors"].append(err)
            return out
    else:
        mb_fires, mb_infos = [], []
    model = cert.get("rejection_model", "isotropic")
    status = combine_multiband(status_sb, mb_fires, model, n_spectra=len(spectra))
    # referee R2: every assumption rebuilt from the supplied band data must be in the verdict
    # ---- MEDIUM-1: bind the run parameters
    sprm = {"T_star": float(prm["T_star"]), "mu_hi": float(prm["mu_hi"]),
            "omega_c": float(prm["omega_c"]), "row_check": prm.get("row_check", "mhat"),
            "rejection_model": model}
    if expected_params is not None:
        exp = _prm_dict(expected_params, expected_params.resolved_omega_c(spectra))
        for k in sprm:
            if exp[k] != sprm[k]:
                out["errors"].append(f"run parameter {k}: certificate {sprm[k]!r} != expected {exp[k]!r}")
    else:
        out["warnings"] = ["no expected_params given: run parameters taken from the certificate"]
    # ---- MEDIUM-1 / LOW-1: rebuild rule ids, assumptions and flags deterministically
    _rules, req_assum, req_flags, T_low, T_up = _annotations(
        spectra, sprm, status, evs, infos, mb_infos or [], len(blocks), T_low, T_up)
    have = set(verdict.assumptions or [])
    missing = [a for a in req_assum if a not in have]
    if missing:
        out["errors"].append(f"verdict lacks required assumptions: {missing[:4]}")
    req_set = set(req_assum)
    stray = [a for a in have if a.startswith(BOUND_PREFIXES) and a not in req_set]
    if stray:
        out["errors"].append(f"verdict assumptions not derivable from the run parameters: {stray[:4]}")
    mf = [f for f in req_flags if f not in set(verdict.flags or [])]
    if mf:
        out["errors"].append(f"verdict lacks required flags: {mf}")
    out["status"] = status.value
    out["T_low_K"], out["T_up_K"] = T_low, T_up
    if status != verdict.status:
        out["errors"].append(f"status mismatch: stored {verdict.status.value}, rechecked {status.value}")
    if T_low != verdict.T_low_K or T_up != verdict.T_up_K:
        out["errors"].append(f"bracket mismatch: stored [{verdict.T_low_K}, {verdict.T_up_K}], "
                             f"rechecked [{T_low}, {T_up}]")
    out["ok"] = not out["errors"]
    out["seconds"] = time.time() - t0
    return out
