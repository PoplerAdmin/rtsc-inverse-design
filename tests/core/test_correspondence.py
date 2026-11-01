"""Correspondence tests: Python decision code vs the Lean definitions (docs/CORRESPONDENCE.md).

1. ``classifyAlt`` / ``classify_evidence`` are compared, on every boolean combination and on
   every alternative list of length <= 5 (plus long random lists), with a Python transcription of
   ``PipelineLogic.classifyAlt`` / ``PipelineLogic.classify``.  The Lean source is copied
   verbatim below; a drift test checks that it still occurs verbatim in Classify.lean.
2. The Lean theorems about ``classify`` (classify_cases, classify_eq_rejected_iff,
   classify_eq_certified_iff, classify_singleton, status_exclusive) and the five Lean
   ``example`` lines are checked on the Python implementation.
3. ``combine_multiband`` and the full decision ``combine_multiband(classify_evidence(..))`` are
   compared with a transcription of Lean ``PipelineLogic.combineMultiband`` / ``classifyFull``
   (Multiband.lean, copied verbatim, drift-checked), and with the characterisation theorems
   ``combine_mb_rejected_iff`` / ``combine_mb_certified_iff`` / ``combine_mb_empty_not_rejected``
   / ``combine_mb_zero_spectra`` / ``classifyFull_zero_spectra``, including the ``n_spectra``
   argument (Lean ``nSpectra : Option ℕ``; ``classifyFull`` passes ``some evs.length``).
4. ``_combine_bracket`` (ERRATA F7) and the strictness of each boolean (R-UNIV strict, D1
   non-strict, R-CELLS strict, C-RAYLEIGH non-strict) are checked at exact boundaries.
5. ``check_rayleigh``'s Q is compared with Lean's ``boxMinQ`` evaluated literally from the Lean
   definitions (C, muTerm, Kp/Fp double sums, g, sel), and ``certified_n_c`` with the literal
   cutoff definition.
6. An accepted R-CELLS certificate is checked against Lean's ``Mhat`` / ``rowSum`` literally
   (every truncation N <= 400, every row n < N), i.e. the conclusion of ``CellOK`` on those N.
"""
from __future__ import annotations

import enum
import itertools
import math
import random
import sys
from fractions import Fraction
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from flint import arb  # noqa: E402

from rtt.core.constants import A0_FRACTION  # noqa: E402
from rtt.core.spectrum import Spectrum  # noqa: E402
from rtt.core.status import Status  # noqa: E402
from rtt.core.certify import classify as certify_classify  # noqa: E402
from rtt.core.certify import CertifyParams  # noqa: E402
from rtt.core.certify.classify import (AltEvidence, _combine_bracket, classify_evidence,  # noqa: E402
                                   classifyAlt, combine_multiband, evidence_for_alternative,
                                   spectrum_digest)
from rtt.core.certify.arbcore import ArbSpectrum  # noqa: E402
from rtt.core.certify.floatcore import FloatSpectrum, float_tc  # noqa: E402
from rtt.core.certify.cells import LamCache, certify_cell, check_cell  # noqa: E402
from rtt.core.certify.rayleigh import certified_n_c, check_rayleigh, propose_rayleigh  # noqa: E402
from rtt.core.certify.search import t_end_for  # noqa: E402


def _lean_file(rel):
    """Locate a Lean source: $RTT_LEAN_DIR or the repository's lean/ directory."""
    import os
    cands = [os.environ.get("RTT_LEAN_DIR"), str(Path(__file__).resolve().parents[2] / "lean")]
    for c in cands:
        if c and (Path(c) / rel).exists():
            return Path(c) / rel
    return Path(cands[1]) / rel

LEAN_CLASSIFY = _lean_file("PipelineLogic/Classify.lean")

# ============================================================================ Lean source (verbatim)
LEAN_STATUS_SRC = """\
inductive Status
  | rejectedProven
  | certifiedCandidate
  | undecided
  | excludedAssumption
  | deferred
  | error
  deriving DecidableEq, Repr"""

LEAN_ALTEVIDENCE_SRC = """\
structure AltEvidence where
  univFires : Bool
  cellsBelowTarget : Bool
  rayleighAtOrAboveTarget : Bool
  deriving DecidableEq, Repr"""

LEAN_CLASSIFYALT_SRC = """\
def classifyAlt (ev : AltEvidence) : Status :=
  if ev.univFires then .rejectedProven
  else if ev.cellsBelowTarget then .rejectedProven
  else if ev.rayleighAtOrAboveTarget then .certifiedCandidate
  else .undecided"""

LEAN_CLASSIFY_SRC = """\
def classify (evs : List AltEvidence) : Status :=
  if evs.isEmpty then .undecided
  else if evs.all (fun ev => classifyAlt ev == .rejectedProven) then .rejectedProven
  else if evs.all (fun ev => classifyAlt ev == .certifiedCandidate) then .certifiedCandidate
  else .undecided"""

LEAN_EXAMPLES_SRC = """\
example : classify [] = .undecided := rfl
example : classify [⟨true, false, false⟩, ⟨false, true, true⟩] = .rejectedProven := rfl
example : classify [⟨false, false, true⟩, ⟨false, false, true⟩] = .certifiedCandidate := rfl
example : classify [⟨true, false, false⟩, ⟨false, false, true⟩] = .undecided := rfl
example : classify [⟨false, false, true⟩, ⟨false, false, false⟩] = .undecided := rfl"""


# ============================================================================ transcription
class LStatus(enum.Enum):
    """Lean `inductive Status` (constructor order kept)."""
    rejectedProven = 0
    certifiedCandidate = 1
    undecided = 2
    excludedAssumption = 3
    deferred = 4
    error = 5


# Lean constructor <-> rtt.core.status.Status member (by documented meaning)
L2P = {LStatus.rejectedProven: Status.REJECTED_PROVEN,
       LStatus.certifiedCandidate: Status.CERTIFIED_CANDIDATE,
       LStatus.undecided: Status.UNDECIDED,
       LStatus.excludedAssumption: Status.EXCLUDED_ASSUMPTION,
       LStatus.deferred: Status.DEFERRED,
       LStatus.error: Status.ERROR}


def lean_List_isEmpty(l: tuple) -> bool:
    # Lean core: def List.isEmpty : List α → Bool | [] => true | _ :: _ => false
    match l:
        case ():
            return True
        case _:
            return False


def lean_List_all(l: tuple, p) -> bool:
    # Lean core: def List.all : List α → (α → Bool) → Bool
    #   | [], _ => true | h :: t, p => p h && all t p
    if l == ():
        return True
    return p(l[0]) and lean_List_all(l[1:], p)


def lean_classifyAlt(ev: tuple[bool, bool, bool]) -> LStatus:
    univFires, cellsBelowTarget, rayleighAtOrAboveTarget = ev     # anonymous constructor order
    if univFires:                                                  # if ev.univFires then
        return LStatus.rejectedProven                              #   .rejectedProven
    elif cellsBelowTarget:                                         # else if ev.cellsBelowTarget then
        return LStatus.rejectedProven                              #   .rejectedProven
    elif rayleighAtOrAboveTarget:                                  # else if ev.rayleighAtOrAboveTarget then
        return LStatus.certifiedCandidate                          #   .certifiedCandidate
    else:
        return LStatus.undecided                                   # else .undecided


def lean_classify(evs: tuple) -> LStatus:
    if lean_List_isEmpty(evs):                                     # if evs.isEmpty then .undecided
        return LStatus.undecided
    elif lean_List_all(evs, lambda ev: lean_classifyAlt(ev) == LStatus.rejectedProven):
        return LStatus.rejectedProven
    elif lean_List_all(evs, lambda ev: lean_classifyAlt(ev) == LStatus.certifiedCandidate):
        return LStatus.certifiedCandidate
    else:
        return LStatus.undecided


ALL_EV = list(itertools.product((False, True), repeat=3))       # 8 AltEvidence values


def py(ev: tuple) -> AltEvidence:
    return AltEvidence(univFires=ev[0], cellsBelowTarget=ev[1], rayleighAtOrAboveTarget=ev[2])


def all_lists(max_len: int):
    for n in range(max_len + 1):
        yield from itertools.product(ALL_EV, repeat=n)


# ============================================================================ 1. drift + structure
@pytest.mark.skipif(not LEAN_CLASSIFY.exists(), reason="Lean sources not available")
@pytest.mark.parametrize("src", [LEAN_STATUS_SRC, LEAN_ALTEVIDENCE_SRC, LEAN_CLASSIFYALT_SRC,
                                 LEAN_CLASSIFY_SRC, LEAN_EXAMPLES_SRC])
def test_lean_source_verbatim(src):
    """The transcribed Lean text still occurs verbatim in Classify.lean (no drift)."""
    assert src in LEAN_CLASSIFY.read_text()


def test_status_constructors_bijective():
    assert [s.name for s in LStatus] == ["rejectedProven", "certifiedCandidate", "undecided",
                                         "excludedAssumption", "deferred", "error"]
    assert set(L2P.values()) == set(Status) and len(set(L2P.values())) == len(Status) == 6


def test_altevidence_fields_match_lean():
    import dataclasses
    assert [f.name for f in dataclasses.fields(AltEvidence)] == \
        ["univFires", "cellsBelowTarget", "rayleighAtOrAboveTarget"]


# ============================================================================ 2. control flow
@pytest.mark.parametrize("ev", ALL_EV)
def test_classifyAlt_equals_lean(ev):
    assert classifyAlt(py(ev)) == L2P[lean_classifyAlt(ev)]


def test_classify_equals_lean_exhaustive_len_le_5():
    n = 0
    for evs in all_lists(5):
        assert classify_evidence([py(e) for e in evs]) == L2P[lean_classify(evs)], evs
        n += 1
    assert n == sum(8 ** k for k in range(6)) == 37449


def test_classify_equals_lean_random_long_lists():
    rng = random.Random(20261027)
    for _ in range(3000):
        n = rng.randint(6, 60)
        mode = rng.random()
        if mode < 0.33:      # mostly-rejecting lists with possibly one odd element
            pool = [e for e in ALL_EV if e[0] or e[1]]
        elif mode < 0.66:    # mostly-certifying lists
            pool = [(False, False, True)]
        else:
            pool = ALL_EV
        evs = [rng.choice(pool) for _ in range(n)]
        if rng.random() < 0.5:
            evs[rng.randrange(n)] = rng.choice(ALL_EV)
        evs = tuple(evs)
        assert classify_evidence([py(e) for e in evs]) == L2P[lean_classify(evs)]


def test_lean_examples_hold_in_python():
    T, F = True, False
    assert classify_evidence([]) == Status.UNDECIDED
    assert classify_evidence([py((T, F, F)), py((F, T, T))]) == Status.REJECTED_PROVEN
    assert classify_evidence([py((F, F, T)), py((F, F, T))]) == Status.CERTIFIED_CANDIDATE
    assert classify_evidence([py((T, F, F)), py((F, F, T))]) == Status.UNDECIDED
    assert classify_evidence([py((F, F, T)), py((F, F, F))]) == Status.UNDECIDED


def test_lean_classify_theorems_hold_in_python():
    three = {Status.REJECTED_PROVEN, Status.CERTIFIED_CANDIDATE, Status.UNDECIDED}
    assert len(three) == 3                                                  # status_exclusive
    for evs in all_lists(4):
        p = [py(e) for e in evs]
        s = classify_evidence(p)
        assert s in three                                                   # classify_cases
        assert (s == Status.REJECTED_PROVEN) == \
            (len(p) > 0 and all(classifyAlt(e) == Status.REJECTED_PROVEN for e in p))   # _rejected_iff
        assert (s == Status.CERTIFIED_CANDIDATE) == \
            (len(p) > 0 and all(classifyAlt(e) == Status.CERTIFIED_CANDIDATE for e in p))  # _certified_iff
        if len(p) == 1:
            assert s == classifyAlt(p[0])                                   # classify_singleton
    for ev in ALL_EV:
        assert classifyAlt(py(ev)) in three                                 # classifyAlt_cases


def test_top_level_classify_empty_list_is_undecided():
    """ERRATA F1 at the public entry point (both rejection models, no band data)."""
    for model in ("isotropic", "multiband"):
        v = certify_classify([], CertifyParams(rejection_model=model))
        assert v.status == Status.UNDECIDED
        assert v.T_low_K is None and v.T_up_K is None


# ============================================================================ 3. combine_multiband
LEAN_MULTIBAND = _lean_file("PipelineLogic/Multiband.lean")

LEAN_COMBINE_SRC = """\
def combineMultiband (sb : Status) (mbFires : List Bool) (model : RejectionModel)
    (nSpectra : Option ℕ := none) : Status :=
  match model with
  | .isotropic => sb
  | .multiband =>
    if nSpectra == some 0 then .undecided
    else
    let mbRej : Bool := !mbFires.isEmpty && mbFires.all id
    if mbRej && sb != .certifiedCandidate then .rejectedProven
    else if sb == .certifiedCandidate && !mbRej then .certifiedCandidate
    else .undecided"""

LEAN_CLASSIFYFULL_SRC = """\
def classifyFull (evs : List AltEvidence) (mbFires : List Bool) (model : RejectionModel) :
    Status :=
  combineMultiband (classify evs) mbFires model (some evs.length)"""

LEAN_MB_EXAMPLES_SRC = """\
example : combineMultiband .rejectedProven [] .multiband = .undecided := rfl
example : combineMultiband .undecided [true, true] .multiband = .rejectedProven := rfl
example : combineMultiband .certifiedCandidate [true] .multiband = .undecided := rfl
example : combineMultiband .certifiedCandidate [true, false] .multiband = .certifiedCandidate :=
  rfl
example : combineMultiband .rejectedProven [] .isotropic = .rejectedProven := rfl
example : combineMultiband .undecided [true] .multiband (some 0) = .undecided := rfl
example : combineMultiband .undecided [true] .multiband (some 2) = .rejectedProven := rfl
example : combineMultiband .rejectedProven [true] .isotropic (some 0) = .rejectedProven := rfl"""

P2L = {v: k for k, v in L2P.items()}


def lean_RejectionModel(s: str) -> str:
    # Python string -> Lean RejectionModel: "isotropic" -> .isotropic, any other -> .multiband
    return "isotropic" if s == "isotropic" else "multiband"


def lean_combineMultiband(sb: LStatus, mbFires: tuple, model: str,
                          nSpectra: int | None = None) -> LStatus:
    if model == "isotropic":                                   # | .isotropic => sb
        return sb
    if nSpectra == 0:                                          # if nSpectra == some 0 then .undecided
        return LStatus.undecided
    mbRej = (not lean_List_isEmpty(mbFires)) and lean_List_all(mbFires, lambda f: f)
    if mbRej and sb != LStatus.certifiedCandidate:             # if mbRej && sb != .certifiedCandidate
        return LStatus.rejectedProven
    elif sb == LStatus.certifiedCandidate and not mbRej:       # else if sb == .cert && !mbRej
        return LStatus.certifiedCandidate
    else:
        return LStatus.undecided


def lean_classifyFull(evs: tuple, mbFires: tuple, model: str) -> LStatus:
    # combineMultiband (classify evs) mbFires model (some evs.length)
    return lean_combineMultiband(lean_classify(evs), mbFires, model, len(evs))


N_SPECTRA = (None, 0, 1, 2, 5)


@pytest.mark.skipif(not LEAN_MULTIBAND.exists(), reason="Lean sources not available")
@pytest.mark.parametrize("src", [LEAN_COMBINE_SRC, LEAN_CLASSIFYFULL_SRC, LEAN_MB_EXAMPLES_SRC])
def test_lean_multiband_source_verbatim(src):
    assert src in LEAN_MULTIBAND.read_text()


def test_lean_multiband_examples_hold_in_python():
    S = Status
    assert combine_multiband(S.REJECTED_PROVEN, [], "multiband") == S.UNDECIDED
    assert combine_multiband(S.UNDECIDED, [True, True], "multiband") == S.REJECTED_PROVEN
    assert combine_multiband(S.CERTIFIED_CANDIDATE, [True], "multiband") == S.UNDECIDED
    assert combine_multiband(S.CERTIFIED_CANDIDATE, [True, False], "multiband") == S.CERTIFIED_CANDIDATE
    assert combine_multiband(S.REJECTED_PROVEN, [], "isotropic") == S.REJECTED_PROVEN
    assert combine_multiband(S.UNDECIDED, [True], "multiband", n_spectra=0) == S.UNDECIDED
    assert combine_multiband(S.UNDECIDED, [True], "multiband", n_spectra=2) == S.REJECTED_PROVEN
    assert combine_multiband(S.REJECTED_PROVEN, [True], "isotropic", n_spectra=0) == S.REJECTED_PROVEN


def test_combine_multiband_equals_lean_exhaustive():
    """All 6 input statuses (incl. ones classify never yields) x fire lists of length <= 6 x
    model strings (incl. an unknown one, which Lean maps to .multiband) x n_spectra in
    (None, 0, 1, 2, 5) (Python None <-> Lean `none`, the default argument)."""
    n = 0
    for model in ("isotropic", "multiband", "garbage"):
        for sb in Status:
            for k in range(7):
                for mb in itertools.product((False, True), repeat=k):
                    for ns in N_SPECTRA:
                        got = (combine_multiband(sb, list(mb), model) if ns is None else
                               combine_multiband(sb, list(mb), model, n_spectra=ns))
                        want = L2P[lean_combineMultiband(P2L[sb], mb, lean_RejectionModel(model), ns)]
                        assert got == want, (sb, mb, model, ns)
                        n += 1
    assert n == 3 * 6 * 127 * len(N_SPECTRA)


def test_classifyFull_equals_lean():
    """The composition used by classify(): combine_multiband(classify_evidence(evs), fires, m)."""
    for model in ("isotropic", "multiband"):
        for evs in all_lists(3):
            for k in range(4):
                for mb in itertools.product((False, True), repeat=k):
                    got = combine_multiband(classify_evidence([py(e) for e in evs]), list(mb), model,
                                            n_spectra=len(evs))
                    assert got == L2P[lean_classifyFull(evs, mb, model)]


SB = list(Status)


def mb_lists(max_len: int):
    for n in range(max_len + 1):
        yield from itertools.product((False, True), repeat=n)


def test_combine_multiband_isotropic_is_identity():
    for sb in SB:
        for mb in mb_lists(3):
            for ns in N_SPECTRA:
                assert combine_multiband(sb, list(mb), "isotropic", n_spectra=ns) == sb


def test_zero_spectra_never_rejected():
    """Lean combine_mb_zero_spectra / classifyFull_zero_spectra(_not_rejected)."""
    for sb in SB:
        for mb in mb_lists(4):
            assert combine_multiband(sb, list(mb), "multiband", n_spectra=0) == Status.UNDECIDED
    for model in ("isotropic", "multiband"):
        for mb in mb_lists(4):
            assert combine_multiband(classify_evidence([]), list(mb), model, n_spectra=0) == Status.UNDECIDED


def test_combine_multiband_characterisation():
    """docs/CORRESPONDENCE.md section 6: the statement a Lean mirror `combineMB` must prove."""
    for sb in SB:
        for mb in mb_lists(4):
            for ns in N_SPECTRA:
                mb_rej = len(mb) > 0 and all(mb)
                ok_n = ns != 0                      # Lean: n ≠ some 0
                out = combine_multiband(sb, list(mb), "multiband", n_spectra=ns)
                # (MB1) REJECTED  <=>  n_spectra != 0, nonempty band list, every R-MULTI fires,
                #                      single band not CERTIFIED          (combine_mb_rejected_iff)
                assert (out == Status.REJECTED_PROVEN) == (ok_n and mb_rej and sb != Status.CERTIFIED_CANDIDATE)
                # (MB2) CERTIFIED <=>  n_spectra != 0, single-band CERTIFIED, not all R-MULTI fire
                assert (out == Status.CERTIFIED_CANDIDATE) == (ok_n and sb == Status.CERTIFIED_CANDIDATE and not mb_rej)
                # (MB3) never anything but the three classify statuses
                assert out in (Status.REJECTED_PROVEN, Status.CERTIFIED_CANDIDATE, Status.UNDECIDED)
                # (MB4) a single-band rejection alone never survives without R-MULTI
                if sb == Status.REJECTED_PROVEN and not mb_rej:
                    assert out == Status.UNDECIDED


def test_combine_multiband_soundness_shape():
    """The two implications the soundness proof needs (with classify_* and R_UNIV_MB)."""
    for sb in (Status.REJECTED_PROVEN, Status.CERTIFIED_CANDIDATE, Status.UNDECIDED):
        for mb in mb_lists(4):
            out = combine_multiband(sb, list(mb), "multiband")
            if out == Status.REJECTED_PROVEN:
                assert len(mb) > 0 and all(mb)          # every band alternative has R_UNIV_MB
            if out == Status.CERTIFIED_CANDIDATE:
                assert sb == Status.CERTIFIED_CANDIDATE  # classify_certified_sound applies


# ============================================================================ 4. bracket / strictness
def test_combine_bracket_F7():
    vals = [None, 1.0, 2.5, 7.0]
    for n in range(0, 4):
        for lows in itertools.product(vals, repeat=n):
            for ups in itertools.product(vals, repeat=n):
                infos = [{"T_low": a, "T_up": b} for a, b in zip(lows, ups)]
                lo, up = _combine_bracket(infos)
                if n > 0 and all(x is not None for x in lows):
                    assert lo == min(lows) and all(lo <= x for x in lows)
                else:
                    assert lo is None
                if n > 0 and all(x is not None for x in ups):
                    assert up == max(ups) and all(up >= x for x in ups)
                else:
                    assert up is None


# a0 * 5^8 = 433148145 / 2^8 is a double, so A_hi == a0 T* can be hit exactly
T_EXACT = 5.0 ** 8
A_EXACT = float(A0_FRACTION * Fraction(5 ** 8))
assert Fraction(A_EXACT) == A0_FRACTION * 5 ** 8


def _evidence(spec, alt_extra, T_star, mu_hi=0.1, omega_c=1000.0):
    alt = {"spectrum_sha256": spectrum_digest(spec)}
    alt.update(alt_extra)
    return evidence_for_alternative(spec, alt, T_star, mu_hi, omega_c)


def test_R_UNIV_strict():
    """Lean UnivCert: A_hi < a0 T* (strict)."""
    s = Spectrum((100.0,), (A_EXACT,), 0.0)
    ev, _ = _evidence(s, {}, T_EXACT)
    assert ev.univFires is False                                   # equality: does not fire
    ev, _ = _evidence(s, {}, math.nextafter(T_EXACT, math.inf))
    assert ev.univFires is True
    # the envelope enters exactly: A_hi = (1+e) A with e = 1/4, A = a0 T*/(5/4) (both doubles)
    s2 = Spectrum((100.0,), (A_EXACT * 0.8,), 0.25)
    assert Fraction(A_EXACT * 0.8) * Fraction(5, 4) == Fraction(A_EXACT)
    ev, _ = _evidence(s2, {}, T_EXACT)
    assert ev.univFires is False                                   # (1+e)A == a0 T*: not strict
    s3 = Spectrum((100.0,), (math.nextafter(A_EXACT * 0.8, -math.inf),), 0.25)
    ev, _ = _evidence(s3, {}, T_EXACT)
    assert ev.univFires is True


def test_D1_nonstrict_and_no_cells():
    """Lean CellsCert: A_hi <= a0 t_end (non-strict); Tiles _ _ [] = False (no R-CELLS)."""
    s = Spectrum((100.0,), (A_EXACT,), 0.0)
    asp = ArbSpectrum(s.omega_K, s.weight_K, 0.0)
    assert t_end_for(asp) == T_EXACT
    blk = {"cells": {"t_end": T_EXACT.hex(), "T_up": T_EXACT.hex(), "cells": []}}
    ev, info = _evidence(s, blk, math.nextafter(T_EXACT, math.inf))
    assert info["T_up"] == T_EXACT and info["cells_valid"] is False
    assert ev.cellsBelowTarget is False                            # no cells -> no R-CELLS evidence
    below = math.nextafter(T_EXACT, -math.inf)
    blk = {"cells": {"t_end": below.hex(), "T_up": below.hex(), "cells": []}}
    ev, info = _evidence(s, blk, 1e9)
    assert info["T_up"] is None                                    # A_hi > a0 t_end: no bound


SPEC_SMALL = Spectrum((60.0, 100.0, 140.0), (15.0, 25.0, 10.0), 0.05)


@pytest.fixture(scope="module")
def small_cell():
    s = SPEC_SMALL
    asp = ArbSpectrum(s.omega_K, s.weight_K, s.rel_envelope)
    fs = FloatSpectrum(s.omega_K, s.weight_K, s.rel_envelope)
    te = t_end_for(asp)
    tc = float_tc(fs, factor=1 + s.rel_envelope, T_hi=te, T_min=1e-3, rtol=1e-6)
    T1 = 0.5 * (tc + te)
    cert, res = certify_cell(asp, fs, T1, te, "mhat", LamCache(asp))
    assert cert is not None and res["ok"]
    return s, asp, te, T1, cert


def test_R_CELLS_strict(small_cell):
    """Lean CellsCert: T_up < T* (strict, ERRATA F4); T_up = cells[0].T1, cells end at t_end."""
    s, asp, te, T1, cert = small_cell
    blk = {"cells": {"t_end": te.hex(), "T_up": T1.hex(), "cells": [cert]}}
    ev, info = _evidence(s, blk, T1)
    assert info["cells_valid"] is True and info["T_up"] == T1
    assert ev.cellsBelowTarget is False                            # T_up == T*: not below
    ev, _ = _evidence(s, blk, math.nextafter(T1, math.inf))
    assert ev.cellsBelowTarget is True
    # a gap in the tiling (exact contiguity) kills the evidence
    blk_gap = {"cells": {"t_end": te.hex(), "T_up": math.nextafter(T1, -math.inf).hex(), "cells": [cert]}}
    ev, info = _evidence(s, blk_gap, 1e9)
    assert ev.cellsBelowTarget is False and info["cells_valid"] is False


@pytest.fixture(scope="module")
def small_rayleigh():
    s = SPEC_SMALL
    fs = FloatSpectrum(s.omega_K, s.weight_K, s.rel_envelope)
    asp = ArbSpectrum(s.omega_K, s.weight_K, s.rel_envelope)
    T = 5.847876662508529
    c = propose_rayleigh(fs, T, 0.1, 1400.0)
    assert c is not None and check_rayleigh(asp, c)["ok"]
    return s, asp, T, c


def test_C_RAYLEIGH_nonstrict(small_rayleigh):
    """Lean RayleighCert: T* <= T (non-strict); run mu_hi / omega_c must match."""
    s, asp, T, c = small_rayleigh
    ev, info = _evidence(s, {"rayleigh": c}, T, mu_hi=0.1, omega_c=1400.0)
    assert ev.rayleighAtOrAboveTarget is True and info["T_low"] == T
    ev, _ = _evidence(s, {"rayleigh": c}, math.nextafter(T, math.inf), mu_hi=0.1, omega_c=1400.0)
    assert ev.rayleighAtOrAboveTarget is False
    ev, _ = _evidence(s, {"rayleigh": c}, T, mu_hi=0.11, omega_c=1400.0)   # other mu_hi
    assert ev.rayleighAtOrAboveTarget is False
    ev, _ = _evidence(s, {"rayleigh": c}, T, mu_hi=0.1, omega_c=1300.0)    # other cutoff
    assert ev.rayleighAtOrAboveTarget is False


# ============================================================================ 5. boxMinQ literal
TWO_PI = 2 * arb.pi()


def lean_g(T, k, w):                      # def g (κ T) (k) (w) := 2 * w / (w ^ 2 + (κ * k * T) ^ 2)
    return 2 * w / (w ** 2 + (TWO_PI * k * T) ** 2)


def lean_adist(n, m):                     # def adist (n m) := (n - m) + (m - n)   (ℕ subtraction)
    return max(n - m, 0) + max(m - n, 0)


def lean_Ssum(L, n):                      # def Ssum (L) (n) := ∑ k ∈ range n, L (k + 1)
    s = arb(0)
    for k in range(n):
        s += L(k + 1)
    return s


def lean_Kp(L, n, m):                     # if n = m then L (2n+1) - 2 * Ssum L n else L (adist n m) + L (n+m+1)
    if n == m:
        return L(2 * n + 1) - 2 * lean_Ssum(L, n)
    return L(lean_adist(n, m)) + L(n + m + 1)


def lean_Fp(L, N, x):                     # ∑ n ∈ range N, ∑ m ∈ range N, Kp L n m * x n * x m
    s = arb(0)
    for n in range(N):
        for m in range(N):
            s += lean_Kp(L, n, m) * x[n] * x[m]
    return s


def lean_C(N, x):                         # -∑ n ∈ range N, (2 * n + 1) * x n ^ 2
    s = arb(0)
    for n in range(N):
        s += (2 * n + 1) * x[n] ** 2
    return -s


def lean_muTerm(mu, Nc, N, x):            # 2 * μ * (∑ n ∈ (range N).filter (n < Nc), x n) ^ 2
    s = arb(0)
    for n in range(N):
        if Nc is None or n < Nc:          # Nc = None stands for ⊤
            s += x[n]
    return 2 * mu * s ** 2


def lean_boxMinQ(omegas, weights, e, mu, Nc, T, x):
    """C N x - muTerm μ Nc N x + ∑ j, boxWeight j * atomF (Ω j);
    boxWeight = sel (atomF) ((1-e) p) ((1+e) p), sel F lo hi = if 0 ≤ F then lo else hi.
    When the sign of the ball F is undecided, sel·F = min(lo F, hi F) (Lean sel_mul_eq_min)."""
    N = len(x)
    xa = [arb(v) for v in x]
    Ta = arb(T)
    ea = arb(e)
    tot = lean_C(N, xa) - lean_muTerm(arb(mu), Nc, N, xa)
    for w, p in zip(omegas, weights):
        wa, pa = arb(w), arb(p)
        F = lean_Fp(lambda k: lean_g(Ta, k, wa), N, xa)
        lo, hi = (1 - ea) * pa, (1 + ea) * pa
        if F >= 0:
            tot += lo * F
        elif F < 0:
            tot += hi * F
        else:
            tot += (lo * F).min(hi * F)
    return tot


def _arb_of(s: str) -> arb:
    return arb(s.strip().strip("[]"))


def lean_Nc(T: float, omega_c: float) -> int:
    """#{n ≥ 0 : (2n+1) π T ≤ ω_c}, decided with exact rationals and a 60-digit π enclosure."""
    from fractions import Fraction as Fr
    pi_lo = Fr(31415926535897932384626433832795028841971693993751, 10 ** 49)
    pi_hi = pi_lo + Fr(1, 10 ** 49)
    Tq, wq = Fr(T), Fr(omega_c)
    n = 0
    while True:
        lo_val, hi_val = (2 * n + 1) * pi_lo * Tq, (2 * n + 1) * pi_hi * Tq
        if hi_val <= wq:
            n += 1
            continue
        if lo_val > wq:
            return n
        raise AssertionError("undecidable at this precision")


def test_certified_n_c_matches_definition():
    rng = random.Random(7)
    for _ in range(300):
        T = rng.uniform(0.5, 400.0)
        wc = rng.uniform(1.0, 5000.0)
        got = certified_n_c(T, wc)
        assert got is None or got == lean_Nc(T, wc)


@pytest.mark.parametrize("seed", range(12))
def test_check_rayleigh_Q_equals_lean_boxMinQ(seed):
    """check_rayleigh's Q (gamma_k / moment-series route) encloses Lean's boxMinQ, and its
    decision equals 'boxMinQ >= 0' whenever the latter is decided."""
    rng = random.Random(seed)
    nat = rng.randint(1, 4)
    omegas = [rng.uniform(5.0, 200.0) for _ in range(nat)]
    weights = [rng.choice([0.0, rng.uniform(1.0, 80.0)]) for _ in range(nat)]
    if not any(w > 0 for w in weights):
        weights[0] = 30.0
    e = rng.choice([0.0, 0.05, 0.3])
    s = Spectrum(tuple(omegas), tuple(weights), e)
    asp = ArbSpectrum(s.omega_K, s.weight_K, s.rel_envelope)
    T = rng.choice([rng.uniform(0.5, 5.0), rng.uniform(20.0, 120.0)])   # second: series branch
    mu = rng.choice([0.0, 0.1, 0.3])
    wc = rng.uniform(10.0, 3000.0)
    Nc = certified_n_c(T, wc)
    if Nc is None:
        pytest.skip("cutoff undecidable at this T")
    Nx = rng.randint(1, 14)
    if rng.random() < 0.5:
        x = [rng.uniform(-1, 1) for _ in range(Nx)]
    else:   # near-optimal vector: Q of either sign near the threshold
        fs = FloatSpectrum(s.omega_K, s.weight_K, s.rel_envelope)
        c = propose_rayleigh(fs, T, mu, wc, N=Nx)
        x = c["x"] if c is not None else [1.0] + [0.0] * (Nx - 1)
    cert = {"T": float(T).hex(), "mu_hi": mu, "omega_c": wc, "N_c": Nc, "x": x}
    r = check_rayleigh(asp, cert)
    assert "Q" in r, r
    Qpy = _arb_of(r["Q"])
    Qlean = lean_boxMinQ(s.omega_K, s.weight_K, e, mu, Nc, T, x)
    assert Qpy.overlaps(Qlean), (str(Qpy), str(Qlean))
    scale = sum(v * v for v in x) + 1.0
    assert float(Qpy.rad()) < 1e-12 * scale * 100 and float(Qlean.rad()) < 1e-12 * scale * 100
    if Qlean >= 0 and not Qlean.contains(0):
        assert r["ok"] is True
    if Qlean < 0:
        assert r["ok"] is False


@pytest.mark.parametrize("case", ["fixture", "strong"])
def test_check_rayleigh_accepting_case_equals_lean_boxMinQ(case, small_rayleigh):
    """Accepted certificates (ok = True): Lean's boxMinQ is certainly >= 0 as well."""
    if case == "fixture":
        s, asp, T, c = small_rayleigh
    else:
        s = Spectrum((50.0, 90.0), (400.0, 150.0), 0.1)
        asp = ArbSpectrum(s.omega_K, s.weight_K, s.rel_envelope)
        fs = FloatSpectrum(s.omega_K, s.weight_K, s.rel_envelope)
        T = 25.0                                            # series branch for k >= 2
        c = propose_rayleigh(fs, T, 0.15, 500.0, N=24)
        assert c is not None
    r = check_rayleigh(asp, c)
    assert r["ok"] is True
    Qlean = lean_boxMinQ(s.omega_K, s.weight_K, s.rel_envelope, float(c["mu_hi"]), int(c["N_c"]),
                         float.fromhex(c["T"]), c["x"])
    assert Qlean >= 0
    assert _arb_of(r["Q"]).overlaps(Qlean)


# ============================================================================ 6. Mhat literal
def lean_Mhat(hi, lo, n, m):
    """if n = m then hi (2n+1) - 2 * Ssum lo n - (2n+1) else hi (adist n m) + hi (n+m+1)"""
    if n == m:
        return hi(2 * n + 1) - 2 * lean_Ssum(lo, n) - (2 * n + 1)
    return hi(lean_adist(n, m)) + hi(n + m + 1)


def test_cell_certificate_satisfies_lean_CellOK_rows(small_cell):
    """For the accepted certificate: rowSum (Mhat hi lo) v N n <= -η v n for N in {N_e, 400}
    and every n < N, with hi = (1+e) lam T1, lo = (1-e) lam T2 (lam by direct atom sums).
    Rows n >= N_e exercise the tail lemma T1 (PROOFS.md), rows n < N_e the head check."""
    s, asp, te, T1, cert = small_cell
    assert check_cell(asp, cert)["ok"]
    T2 = float.fromhex(cert["T2"])
    e = arb(s.rel_envelope)
    Nmax = 400
    Kmax = 2 * Nmax + 1

    def lam(T):
        Ta = arb(T)
        out = [arb(0)] * (Kmax + 1)
        for k in range(Kmax + 1):
            acc = arb(0)
            for w, p in zip(s.omega_K, s.weight_K):
                acc += arb(p) * lean_g(Ta, k, arb(w))
            out[k] = acc
        return out
    l1, l2 = lam(T1), lam(T2)
    hi = lambda k: (1 + e) * l1[k]          # noqa: E731  (k = 0 is never used by Mhat)
    lo = lambda k: (1 - e) * l2[k]          # noqa: E731
    N, u, C, p = int(cert["N"]), cert["u"], float(cert["C"]), int(cert.get("p", 3))
    v = [arb(u[m]) if m < N else arb(C) / arb(2 * m + 1) ** p for m in range(Nmax)]
    eta = arb(float(cert["eta"]))
    S_lo = [arb(0)] * (Nmax + 1)
    for n in range(1, Nmax + 1):
        S_lo[n] = S_lo[n - 1] + lo(n)
    for Ntr in sorted({int(cert["N_e"]), Nmax}):
        for n in range(Ntr):
            row = arb(0)
            for m in range(Ntr):
                if m == n:
                    M = hi(2 * n + 1) - 2 * S_lo[n] - (2 * n + 1)
                    assert M.overlaps(lean_Mhat(hi, lo, n, m))
                else:
                    M = hi(lean_adist(n, m)) + hi(n + m + 1)
                row += M * v[m]
            assert row <= -eta * v[n], (Ntr, n, str(row), str(-eta * v[n]))
