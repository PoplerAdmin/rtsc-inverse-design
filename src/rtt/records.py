"""Shared record types (architecture §6.1). Schema id ``rtt.candidate/2``.

This module is the interface between builders: stages exchange these dataclasses, the ledger
stores them as JSON, and ``rtt.decide`` consumes :class:`DecideInputs`.

Conventions
* All dataclasses have ``to_json()`` (a JSON-able dict) and ``from_json(d)``.
* Field names of :class:`Validity`, :class:`Envelope`, :class:`MBInputs` are those of Lean
  ``TargetLogic.Decide`` (camelCase on purpose, so the mirror can be diff-checked).
* A ``None`` in a field that a rule or G-VALID needs makes G-VALID fail (§6.1: missing data
  always yields ``undecided``).
* Per-number provenance (review RL1 #4): the envelope and each multiband input carry their own
  :data:`Tier`. A rejection needs every input of the fired rule to be proof grade.
"""
from __future__ import annotations

import dataclasses
import json
from dataclasses import dataclass, field
from typing import Any, Literal, Optional

SCHEMA_CANDIDATE = "rtt.candidate/2"
SCHEMA_DECIDE_INPUTS = "rtt.decide_inputs/1"

# --------------------------------------------------------------------------- vocabularies
STAGES = ("S0", "S1", "S2", "S3", "S4", "S5", "S6", "S7", "S8", "S9")

# stage outcomes (§3.1)
PASS = "PASS"
SCREENED_OUT = "SCREENED_OUT"
EXCLUDED_ASSUMPTION = "EXCLUDED_ASSUMPTION"
DEFERRED = "DEFERRED"
INVALID = "INVALID"
DUPLICATE = "DUPLICATE"
ERROR = "ERROR"
STAGE_OUTCOMES = (PASS, SCREENED_OUT, EXCLUDED_ASSUMPTION, DEFERRED, INVALID, DUPLICATE, ERROR)

# decide outcomes (Lean TargetLogic.Outcome)
D_REJECTED_PROVEN = "rejectedProven"
D_REJECTED_ASSUMING = "rejectedAssuming"
D_SCREENED_OUT = "screenedOut"
D_UNDECIDED = "undecided"
D_PASS = "pass"
DECIDE_OUTCOMES = (D_REJECTED_PROVEN, D_REJECTED_ASSUMING, D_SCREENED_OUT, D_UNDECIDED, D_PASS)
REMOVAL_DECIDE = frozenset({D_REJECTED_PROVEN, D_REJECTED_ASSUMING, D_SCREENED_OUT})

#: (stage -> outcomes it may emit); §3.1 table. Mirrored by the ledger trigger
#: ``outcomes_allowed`` (table ``stage_outcomes`` in schema.sql). Deviation (20_build_core.md):
#: ERROR is also allowed at S7 and S8 (a compute step can fail there too).
ALLOWED_OUTCOMES: dict[str, frozenset[str]] = {
    "S1": frozenset({PASS, SCREENED_OUT, DEFERRED}),     # DEFERRED: R30 F1(b), doc 31 §7
    "S2": frozenset({PASS, INVALID, DUPLICATE}),
    "S3": frozenset({PASS, SCREENED_OUT, DEFERRED}),
    "S4": frozenset({PASS, EXCLUDED_ASSUMPTION, SCREENED_OUT, DEFERRED, ERROR}),
    "S5": frozenset({PASS, SCREENED_OUT, DEFERRED, ERROR, D_SCREENED_OUT, D_UNDECIDED, D_PASS}),
    "S6": frozenset({D_REJECTED_PROVEN, D_REJECTED_ASSUMING, D_UNDECIDED, D_PASS, ERROR}),
    "S7": frozenset({EXCLUDED_ASSUMPTION, D_REJECTED_PROVEN, D_REJECTED_ASSUMING, D_UNDECIDED,
                     D_PASS, ERROR}),
    "S8": frozenset({D_REJECTED_PROVEN, D_REJECTED_ASSUMING, D_UNDECIDED, D_PASS, ERROR}),
}

Tier = Literal["ml", "coarseDFPT", "fullDFPT", "converged"]
TIERS = ("ml", "coarseDFPT", "fullDFPT", "converged")
PROOF_GRADE = frozenset({"fullDFPT", "converged"})

LANES = ("G1", "G2", "G3", "G4", "G5", "GA", "D", "C", "P")   # P: prototype substitution (doc 32)
COST_CLASSES = ("C1", "C2", "C3")

#: tier_runs.profile -> verdicts.spectrum_tier and -> decide Tier (ledger triggers use the same)
PROFILE_SPECTRUM_TIER = {"coarse": "coarse", "full": "full", "epw": "converged"}
PROFILE_DECIDE_TIER = {"coarse": "coarseDFPT", "full": "fullDFPT", "epw": "converged"}

ORIGIN_CERTIFY = "rtt.core.certify_policy.certify_spectra"
ORIGIN_DECIDE = "rtt.decide.core.decide_t"


def proof_grade(t: Optional[str]) -> bool:
    return t in PROOF_GRADE


# --------------------------------------------------------------------------- helpers
class _J:
    """Mixin: dataclass <-> JSON dict (nested dataclasses handled via ``_nested``)."""

    _nested: dict[str, Any] = {}

    def to_json(self) -> dict:
        out = {}
        for f in dataclasses.fields(self):  # type: ignore[arg-type]
            v = getattr(self, f.name)
            out[f.name] = _to_jsonable(v)
        return out

    @classmethod
    def from_json(cls, d: dict | None):
        if d is None:
            return None
        names = {f.name for f in dataclasses.fields(cls)}  # type: ignore[arg-type]
        kw = {}
        for k, v in d.items():
            if k not in names:
                continue
            sub = cls._nested.get(k)
            if sub is not None and v is not None:
                if isinstance(sub, tuple) and sub[0] == "list":
                    v = [sub[1].from_json(x) for x in v]
                elif isinstance(sub, tuple) and sub[0] == "dict":
                    v = {kk: sub[1].from_json(x) for kk, x in v.items()}
                else:
                    v = sub.from_json(v)
            kw[k] = v
        return cls(**kw)

    def dumps(self) -> str:
        return json.dumps(self.to_json(), sort_keys=True, separators=(",", ":"))


def _to_jsonable(v: Any) -> Any:
    if isinstance(v, _J):
        return v.to_json()
    if isinstance(v, (list, tuple)):
        return [_to_jsonable(x) for x in v]
    if isinstance(v, dict):
        return {k: _to_jsonable(x) for k, x in v.items()}
    if isinstance(v, frozenset):
        return sorted(v)
    return v


# --------------------------------------------------------------------------- structures
@dataclass
class StructureData(_J):
    lattice: list[list[float]]
    species: list[str]
    frac_coords: list[list[float]]
    spacegroup: Optional[int] = None
    n_atoms: Optional[int] = None

    def __post_init__(self):
        if self.n_atoms is None:
            self.n_atoms = len(self.species)
        if len(self.frac_coords) != len(self.species):
            raise ValueError("species and frac_coords differ in length")

    @property
    def formula(self) -> str:
        from collections import Counter
        c = Counter(self.species)
        return "".join(f"{el}{n if n > 1 else ''}" for el, n in sorted(c.items()))


@dataclass
class Origin(_J):
    lane: str
    control: Optional[str] = None
    T_star_surrogate_K: Optional[float] = None
    generator: Optional[str] = None
    checkpoint: Optional[str] = None
    checkpoint_sha256: Optional[str] = None
    conditioning: Optional[dict] = None
    guidance: Optional[float] = None
    batch_id: Optional[str] = None
    seed: Optional[int] = None
    flags: list[str] = field(default_factory=list)
    pressure_override_GPa: Optional[float] = None

    def __post_init__(self):
        if self.lane not in LANES:
            raise ValueError(f"unknown lane {self.lane!r}")


@dataclass
class RawStructure(_J):
    """Output of S1 sources (MatterGen, database, controls) before canonicalisation."""
    structure: StructureData
    origin: Origin
    pressure_GPa: float = 0.0
    source_ref: Optional[str] = None      # file path / LMDB key / MatterGen batch item
    meta: dict = field(default_factory=dict)
    _nested = {"structure": StructureData, "origin": Origin}


@dataclass
class Context(_J):
    targets_sha256: str
    tiers_sha256: str
    code_git_sha: Optional[str] = None
    core_version: str = "core-1"


@dataclass
class ArtifactRef(_J):
    path: str
    sha256: str
    bytes: int
    id: Optional[int] = None
    retained: bool = True


@dataclass
class CostEstimate(_J):
    core_h: float = 0.0
    gpu_h: float = 0.0
    wall_s: float = 0.0
    disk_GB: float = 0.0
    cost_class: Optional[str] = None
    lane: str = "none"


# --------------------------------------------------------------------------- stage outcome
@dataclass
class Outcome(_J):
    """One stage outcome (a row of ``outcomes``). Heuristic removals carry criterion, value,
    threshold, source_id and provisional (R2, R7)."""
    stage: str
    outcome: str
    criterion: Optional[str] = None
    value: Any = None
    threshold: Any = None
    source_id: Optional[str] = None
    provisional: bool = False
    reason: Optional[str] = None          # DEFERRED reason / ERROR text / decide reason
    retry: int = 0
    tier_run_id: Optional[int] = None
    decide: Optional[dict] = None         # DecideOutcome.to_json() at S5-S8
    tags: list[str] = field(default_factory=list)   # e.g. PROVISIONAL_MARGIN, CROSSCHECK_DISAGREE

    def __post_init__(self):
        allowed = ALLOWED_OUTCOMES.get(self.stage)
        if allowed is not None and self.outcome not in allowed:
            raise ValueError(f"stage {self.stage} may not emit {self.outcome!r} (§3.1)")
        if self.outcome in (SCREENED_OUT, EXCLUDED_ASSUMPTION) and not (self.criterion and self.source_id):
            raise ValueError(f"{self.outcome} needs criterion and source_id (R2/R3/R7)")
        if self.outcome == DEFERRED and not (self.reason or self.criterion):
            raise ValueError("DEFERRED needs a reason")

    @property
    def removes(self) -> bool:
        return self.outcome in (SCREENED_OUT, EXCLUDED_ASSUMPTION, INVALID, DUPLICATE,
                                D_REJECTED_PROVEN, D_REJECTED_ASSUMING, D_SCREENED_OUT)

    @property
    def passes(self) -> bool:
        return self.outcome in (PASS, D_PASS, D_UNDECIDED)


@dataclass
class StageResult(_J):
    """What ``Stage.run`` returns; ``Stage.decide`` turns it into an :class:`Outcome`."""
    stage: str
    candidate_id: str
    ok: bool = True
    metrics: dict = field(default_factory=dict)
    tier_run_id: Optional[int] = None
    artifacts: list[ArtifactRef] = field(default_factory=list)
    error: Optional[str] = None
    payload: Any = None                   # in-memory objects (spectra, structures); not stored
    _nested = {"artifacts": ("list", ArtifactRef)}

    def to_json(self) -> dict:
        d = super().to_json()
        d.pop("payload", None)
        return d


# --------------------------------------------------------------------------- decide inputs
# Mirror of Lean TargetLogic.Decide (v2, after review RL1): Prov, IsoData, EpsData, MatData,
# DosData, Validity, InputTag, Inputs. Field names are the Lean ones.
LEAN_TIERS = ("fullDFPT", "coarseDFPT", "ml", "external")


@dataclass
class Prov(_J):
    """Provenance of one input group (doc 01 §0.6c): Lean tier and run / structure id.
    A converged (EPW) run is recorded as ``fullDFPT`` (Lean has no separate converged tier)."""
    tier: str
    run: Any

    def __post_init__(self):
        if self.tier not in LEAN_TIERS:
            raise ValueError(f"Prov.tier must be one of {LEAN_TIERS}, got {self.tier!r}")


@dataclass
class IsoData(_J):
    A: Optional[float] = None
    lam: Optional[float] = None
    m: Optional[float] = None
    dA: Optional[float] = None
    dLam: Optional[float] = None
    dM: Optional[float] = None
    provA: Optional[Prov] = None
    provLam: Optional[Prov] = None
    provM: Optional[Prov] = None
    _nested = {"provA": Prov, "provLam": Prov, "provM": Prov}


@dataclass
class EpsData(_J):
    epsSys: Optional[float] = None
    epsOmega: Optional[float] = None
    qSpread: Optional[float] = None
    prov: Optional[Prov] = None
    _nested = {"prov": Prov}


@dataclass
class MatData(_J):
    """Band/k-resolved coupling-weight matrix calA_ij (nominal), DOS weights N_i, a
    Collatz-Wielandt test vector v > 0, spread dCalA (doc 01 §0.6b)."""
    calA: list[list[float]]
    n: list[float]
    v: list[float]
    dCalA: float
    prov: Prov
    _nested = {"prov": Prov}

    @property
    def B(self) -> int:
        return len(self.calA)


@dataclass
class DosData(_J):
    """Sheet DOS fractions w_i = N_i/N with spread dW."""
    w: list[float]
    dW: float
    prov: Prov
    _nested = {"prov": Prov}


@dataclass
class Validity(_J):
    """G-VALID inputs. ``imagAbsMaxK`` and ``anharmShiftAbs`` are MAGNITUDES (§0.6a); the mapper
    (``rtt.decide.inputs``) takes absolute values, WELL-FORMED rejects negative entries.
    ``None`` = no source record (only ``EF`` may be None in Lean; any other None -> undecided)."""
    EF: Optional[float] = None
    omegaMax: Optional[float] = None
    imagAbsMaxK: Optional[float] = None
    containsHLiBe: Optional[bool] = None
    anharmonicCheck: Optional[bool] = None
    anharmShiftAbs: Optional[float] = None
    alpha2FRecomputedAnharmonic: Optional[bool] = None
    muStarNegativeSuspected: Optional[bool] = None
    admittedSmearings: Optional[int] = None
    qMeshes: Optional[int] = None
    nSheets: Optional[int] = None


@dataclass
class InputTag(_J):
    """Lean ``InputTag``: kind in fileIncomplete | fileCompleteUnverified | dropNonpos | negClip |
    sigmaExcluded | informational | unknown. ``unknown`` blocks rejection (fail-closed)."""
    name: str = ""
    kind: str = "informational"
    massK: Optional[float] = None
    blocksRejection: Optional[bool] = None       # sigmaExcluded only (from core blocking_tags)
    blocksCertification: Optional[bool] = None


@dataclass
class DecideInputs(_J):
    """Lean ``Inputs`` (val, iso, eps, mat, dos, tags) plus ledger metadata: ``tier`` is the
    stage tier of the record (ledger binding, §6.2; decideT itself uses the per-input Prov),
    ``cfg`` and ``provenance`` (tier_run_ids, artifact_ids, record_ids) are recorded only."""
    tier: str
    val: Validity
    iso: IsoData
    eps: EpsData
    mat: Optional[MatData] = None
    dos: Optional[DosData] = None
    tags: list[InputTag] = field(default_factory=list)
    cfg: dict = field(default_factory=dict)
    provenance: dict = field(default_factory=dict)
    schema: str = SCHEMA_DECIDE_INPUTS
    _nested = {"val": Validity, "iso": IsoData, "eps": EpsData, "mat": MatData, "dos": DosData,
               "tags": ("list", InputTag)}

    def __post_init__(self):
        if self.tier not in TIERS:
            raise ValueError(f"unknown tier {self.tier!r}")


@dataclass
class DecideOutcome(_J):
    """Result of ``rtt.decide.core.decide_t`` (Lean ``Outcome`` + recorded data): outcome in
    rejectedProven | screenedOut | undecided | pass; ``rule`` = fired rule (rejectedProven) or the
    gate that failed; ``assumptions`` = Lean ``Outcome.assumptions``; ``would_fire`` = rules whose
    floor holds (on any tier)."""
    outcome: str
    rule: Optional[str] = None
    assumptions: list[str] = field(default_factory=list)
    inputs_used: list[dict] = field(default_factory=list)
    reasons: list[str] = field(default_factory=list)
    would_fire: list[str] = field(default_factory=list)
    core_outcome: Optional[str] = None
    blocked: bool = False

    def __post_init__(self):
        if self.outcome not in DECIDE_OUTCOMES:
            raise ValueError(self.outcome)


# --------------------------------------------------------------------------- candidate record
@dataclass
class StageEntry(_J):
    outcome: str
    tier_run_id: Optional[int] = None
    metrics: dict = field(default_factory=dict)
    profile: Optional[str] = None
    score: Optional[float] = None
    margin: Optional[dict] = None
    decide: Optional[dict] = None
    verdict_ids: list[int] = field(default_factory=list)
    flags: list[str] = field(default_factory=list)
    criterion: Optional[str] = None
    reason: Optional[str] = None


@dataclass
class CandidateRecord(_J):
    candidate_id: str
    structure_sha256: str
    formula: str
    structure: StructureData
    origin: Origin
    pressure_GPa: float = 0.0
    cost_class: Optional[str] = None
    context: Optional[Context] = None
    stages: dict[str, StageEntry] = field(default_factory=dict)
    flags: list[str] = field(default_factory=list)
    current: dict = field(default_factory=dict)
    schema: str = SCHEMA_CANDIDATE
    _nested = {"structure": StructureData, "origin": Origin, "context": Context,
               "stages": ("dict", StageEntry)}

    def passed(self, stage: str) -> bool:
        e = self.stages.get(stage)
        return e is not None and e.outcome in (PASS, D_PASS, D_UNDECIDED)

    def is_control(self) -> bool:
        return self.origin.lane == "C"


def new_candidate_id(structure_sha256: str, lane: str, pressure_GPa: float = 0.0) -> str:
    import hashlib
    h = hashlib.sha256(f"{structure_sha256}|{lane}|{pressure_GPa:.3f}".encode()).hexdigest()
    return "c-" + h[:12]
