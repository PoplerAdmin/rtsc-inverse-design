"""Build :class:`~rtt.records.DecideInputs` (Lean ``Inputs``) from stage records
(architecture §6.1; doc 01 §0.2b, §0.6).

Every field is filled from a named record; a field without a source stays ``None`` and
``decide_t`` then returns ``undecided`` (MISSING). The mappers do the sign conventions of §0.6a:
imaginary modes and anharmonic shifts become magnitudes here, so a negative value reaching
``decide_t`` is malformed.

Metric keys read from ``CandidateRecord.stages[...].metrics`` (contract with the S4/S5/S6/S7
producers, B3):

S4:  ``E_F_eV`` (bottom of the shallowest crossing band), ``n_sheets``
S5/S6: ``ladder`` = {"q_meshes": [[n,n,n], ...] (coarse first),
                     "fine":   {"A_K": [...], "lambda": [...], "m_K2": [...]}  admitted σ, finest q,
                     "coarse": {"A_K": [...], "lambda": [...], "m_K2": [...]}  same σ, next q,
                     "omega_log_K_by_q": [coarse, fine]},
       ``omega_max_K``, ``min_freq_K`` (signed; negative = imaginary), ``run_id`` (optional),
       ``mu_star_negative_suspected`` (bool)
S7:  ``anharmonic_check`` (bool), ``anharmonic_shift`` (signed relative shift, max |.| over
     coupled modes), ``alpha2F_recomputed_anharmonic`` (bool),
     ``mb`` = {"calA": [[...]], "N": [...], "dCalA": float, "w": [...], "dW": float,
               "run_id": ..., "tier": "fullDFPT"|...}
"""
from __future__ import annotations

import re
import statistics
from typing import Any, Iterable, Optional, Sequence

import numpy as np

from ..records import (CandidateRecord, DecideInputs, DosData, EpsData, InputTag, IsoData, MatData, Prov,
                       Validity)
from ..units import EV_TO_K

H_LI_BE = frozenset({"H", "Li", "Be"})

#: stage tier -> Lean Prov tier (Lean has no "converged"; an EPW run is fullDFPT-grade)
PROV_TIER = {"fullDFPT": "fullDFPT", "converged": "fullDFPT", "coarseDFPT": "coarseDFPT", "ml": "ml"}


# ------------------------------------------------------------------------------------ tags
_MASS = re.compile(r"^A-(DROP-NONPOS|NEG-CLIP)\(([-+0-9.eE]+)K\)")


def tag_from_string(s: str, core_rej: Iterable[str] = ()) -> InputTag:
    """Map a ``Spectrum.provenance['assumptions']`` string to a Lean ``InputTag``.

    ``core_rej``: first list returned by core ``blocking_tags(spectra)`` (entries carry a
    `` [label]`` suffix); A-SIGMA-EXCLUDED blocks rejection iff the core policy says so.
    An ``A-...`` string the mapper does not recognise becomes ``unknown`` (blocks; RL1 #7b)."""
    s = str(s)
    rej = list(core_rej)
    if s.startswith("A-FILE-INCOMPLETE"):
        return InputTag(s, "fileIncomplete")
    if s.startswith("A-FILE-COMPLETE(unverified)"):
        return InputTag(s, "fileCompleteUnverified")
    m = _MASS.match(s)
    if m:
        return InputTag(s, "dropNonpos" if m.group(1) == "DROP-NONPOS" else "negClip", massK=float(m.group(2)))
    if s.startswith("A-SIGMA-EXCLUDED("):
        blocks = any(r == s or r.startswith(s + " [") for r in rej)
        return InputTag(s, "sigmaExcluded", blocksRejection=blocks)
    if s.startswith("A-FILE-COMPLETE"):
        return InputTag(s, "informational")
    if s.startswith("A-"):
        return InputTag(s, "unknown")
    return InputTag(s, "informational")


def tags_from_spectra(spectra: Sequence[Any]) -> list[InputTag]:
    """Input-quality tags of all alternatives; BlockedRej agrees with core ``blocking_tags``
    (architecture §3.2): every core rejection-blocking tag maps to a blocking InputTag."""
    from ..core.certify_policy import blocking_tags
    rej, _ = blocking_tags(spectra)
    seen, out = set(), []
    for sp in spectra:
        for a in sp.provenance.get("assumptions", []) or []:
            if a not in seen:
                seen.add(a)
                out.append(tag_from_string(a, rej))
    for r in rej:
        base = r.rsplit(" [", 1)[0]
        if base not in seen:          # fail-closed: anything core blocks stays blocking
            seen.add(base)
            out.append(InputTag(base, "unknown"))
    return out


# ------------------------------------------------------------------------------------ envelope
def rel_spread(values: Sequence[float], nom: float) -> float:
    """δ^σ = max_i |X_i - X_nom| / X_nom."""
    return max(abs(v - nom) for v in values) / nom


def envelope_from_ladder(fine: dict, coarse: dict, prov: Prov, *, eps_sys: float = 0.10,
                         eps_omega_min: float = 0.02, omega_log_by_q: Sequence[float] | None = None,
                         anharm_shift_abs: float = 0.0) -> tuple[IsoData, EpsData]:
    """doc 01 §0.2b. ``fine``/``coarse``: {"A_K", "lambda", "m_K2"} over the *admitted* smearings
    (same σ order). X_nom = median over admitted σ on the fine q mesh (A_mid, §3.6);
    δ_X = δ^σ + δ^q; ε_ω = max(ε_ω,min, |shift|, q-spread of ω_log)."""
    vals = {}
    for key in ("A_K", "lambda", "m_K2"):
        xs, cs = [float(v) for v in fine[key]], [float(v) for v in coarse[key]]
        nom = float(statistics.median(xs))
        vals[key] = (nom, rel_spread(xs, nom) + abs(nom - float(statistics.median(cs))) / nom)
    qspread = 0.0
    if omega_log_by_q and len(omega_log_by_q) >= 2:
        qspread = abs(omega_log_by_q[-1] - omega_log_by_q[-2]) / omega_log_by_q[-1]
    eps_omega = max(eps_omega_min, abs(anharm_shift_abs), qspread)
    iso = IsoData(A=vals["A_K"][0], lam=vals["lambda"][0], m=vals["m_K2"][0], dA=vals["A_K"][1],
                  dLam=vals["lambda"][1], dM=vals["m_K2"][1], provA=prov, provLam=prov, provM=prov)
    return iso, EpsData(epsSys=eps_sys, epsOmega=eps_omega, qSpread=qspread, prov=prov)


# ------------------------------------------------------------------------------------ multiband
def perron_vector(K: np.ndarray, iters: int = 500) -> list[float]:
    """Nominal Perron vector of a non-negative matrix (float power iteration), strictly > 0.
    Any v > 0 gives a valid Collatz-Wielandt upper bound; this choice makes it tight."""
    K = np.asarray(K, dtype=float)
    v = np.ones(K.shape[0])
    for _ in range(iters):
        w = K @ v + v * 1e-12
        v = w / np.max(w)
    return [float(max(x, 1e-9)) for x in v]


def mat_from_matrix(calA: Sequence[Sequence[float]], N: Sequence[float], dCalA: float, prov: Prov,
                    v: Sequence[float] | None = None) -> MatData:
    """MatData of doc 01 §0.6b. Lean computes κ_hi = max_i (calA_hi v)_i / v_i exactly; v defaults
    to the nominal Perron vector of calA (same spectrum as sqrt(N_i/N_j) calA_ij by similarity)."""
    A = np.asarray(calA, dtype=float)
    Nn = np.asarray(N, dtype=float)
    if v is None:
        v = list(perron_vector(A))
    return MatData(calA=[[float(a) for a in r] for r in A], n=[float(x) for x in Nn], v=[float(x) for x in v],
                   dCalA=float(dCalA), prov=prov)


# ------------------------------------------------------------------------------------ build
def imag_abs_max(min_freq_K: Optional[float]) -> Optional[float]:
    """Largest |ω| among negative (imaginary) eigenvalues, 0 if none (§0.6a)."""
    if min_freq_K is None:
        return None
    return abs(float(min_freq_K)) if min_freq_K < 0 else 0.0


def build_inputs(rec: CandidateRecord, tier: str, *, stage: str = "S6", tags: Sequence[InputTag] = (),
                 eps_sys: float = 0.10, eps_omega_min: float = 0.02,
                 cfg_record: dict | None = None) -> DecideInputs:
    """Assemble DecideInputs from the S4, ``stage`` (S5/S6/S8) and S7 records of ``rec``."""
    s4, sd, s7 = rec.stages.get("S4"), rec.stages.get(stage), rec.stages.get("S7")
    m4 = s4.metrics if s4 else {}
    md = sd.metrics if sd else {}
    m7 = s7.metrics if s7 else {}
    contains = bool(set(rec.structure.species) & H_LI_BE)
    shift = m7.get("anharmonic_shift")
    shift_abs = abs(float(shift)) if shift is not None else (0.0 if not m7.get("anharmonic_check") else None)
    run = md.get("run_id", sd.tier_run_id if sd else None)
    prov = Prov(PROV_TIER[tier], run) if run is not None else None
    lad = md.get("ladder") or {}
    iso, eps = IsoData(), EpsData()
    admitted = None
    if prov is not None and lad.get("fine") and lad.get("coarse"):
        iso, eps = envelope_from_ladder(lad["fine"], lad["coarse"], prov, eps_sys=eps_sys,
                                        eps_omega_min=eps_omega_min, omega_log_by_q=lad.get("omega_log_K_by_q"),
                                        anharm_shift_abs=shift_abs or 0.0)
        admitted = len(lad["fine"]["A_K"])
    ef = m4.get("E_F_eV")
    val = Validity(EF=None if ef is None else float(ef) * EV_TO_K, omegaMax=md.get("omega_max_K"),
                   imagAbsMaxK=(abs(float(md["max_imag_mode_K"])) if md.get("max_imag_mode_K") is not None
                                else imag_abs_max(md.get("min_freq_K"))), containsHLiBe=contains,
                   anharmonicCheck=bool(m7.get("anharmonic_check", False)), anharmShiftAbs=shift_abs,
                   alpha2FRecomputedAnharmonic=bool(m7.get("alpha2F_recomputed_anharmonic", False)),
                   muStarNegativeSuspected=md.get("mu_star_negative_suspected"),
                   admittedSmearings=admitted, qMeshes=(len(lad["q_meshes"]) if lad.get("q_meshes") else None),
                   nSheets=m4.get("n_sheets"))
    mat = dos = None
    mbd = m7.get("mb") or {}
    if mbd:
        mprov = Prov(PROV_TIER.get(mbd.get("tier"), mbd.get("tier") or "external"), mbd.get("run_id"))
        if mbd.get("calA") is not None:
            mat = mat_from_matrix(mbd["calA"], mbd["N"], mbd["dCalA"], mprov)
        if mbd.get("w") is not None:
            dos = DosData(w=[float(w) for w in mbd["w"]], dW=float(mbd["dW"]), prov=mprov)
    prov_rec = {"tier_run_ids": [e.tier_run_id for k, e in sorted(rec.stages.items())
                                 if k in ("S5", "S6", "S7", "S8") and e.tier_run_id is not None],
                "record_ids": {k: e.tier_run_id for k, e in rec.stages.items()}}
    return DecideInputs(tier=tier, val=val, iso=iso, eps=eps, mat=mat, dos=dos, tags=list(tags),
                        cfg=cfg_record or {}, provenance=prov_rec)
