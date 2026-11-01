"""Stage interface (architecture §5.3).

Every stage S1..S8 implements :class:`Stage`. ``run`` does the compute (and writes ``tier_runs``
and ``artifacts`` through ``ctx.ledger``); ``decide`` is pure: it maps a stored
:class:`~rtt.records.StageResult` plus thresholds and target spec to an
:class:`~rtt.records.Outcome`. ``rtt redecide --stage S5`` calls only ``decide``.

The allowed outcomes per stage are ``rtt.records.ALLOWED_OUTCOMES`` (enforced again by the
ledger trigger ``outcomes_allowed``).
"""
from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Literal, Optional, Protocol, runtime_checkable

from ..config import LoadedConfig, Targets
from ..records import ALLOWED_OUTCOMES, CandidateRecord, CostEstimate, Outcome, StageResult

StageKind = Literal["heuristic", "physical", "dfpt", "report"]
LaneId = Literal["G", "S", "F", "none"]


@dataclass
class Context:
    """Everything a stage needs at run time."""
    ledger: Any                          # rtt.ledger.store.Ledger
    targets: Targets
    tiers: LoadedConfig
    machine: LoadedConfig
    workdir: Path
    dry_run: bool = False
    launcher: Any = None                 # dft launcher (B3), optional
    extra: dict = field(default_factory=dict)

    def targets_for(self, rec: CandidateRecord) -> Targets:
        """Control lane C uses its surrogate target T*_s (§3.3)."""
        if rec.origin.lane == "C" and rec.origin.T_star_surrogate_K:
            return self.targets.with_T_star(rec.origin.T_star_surrogate_K)
        return self.targets


Thresholds = LoadedConfig    # tiers.toml (Threshold objects via .threshold("S5.margin_M_c"))


@runtime_checkable
class Stage(Protocol):
    id: str
    kind: StageKind
    lane: LaneId
    requires: tuple[str, ...]
    allowed: frozenset[str]

    def estimate_cost(self, rec: CandidateRecord, machine: LoadedConfig) -> CostEstimate: ...

    def run(self, rec: CandidateRecord, ctx: Context) -> StageResult: ...

    def decide(self, res: StageResult, th: Thresholds, tgt: Targets) -> Outcome: ...


class StageBase:
    """Convenience base: fills ``allowed`` from the §3.1 table and checks emitted outcomes."""
    id: str = "S?"
    kind: StageKind = "heuristic"
    lane: LaneId = "none"
    requires: tuple[str, ...] = ()

    @property
    def allowed(self) -> frozenset[str]:
        return ALLOWED_OUTCOMES.get(self.id, frozenset())

    def estimate_cost(self, rec: CandidateRecord, machine: LoadedConfig) -> CostEstimate:
        return CostEstimate(lane=self.lane, cost_class=rec.cost_class)

    def check(self, o: Outcome) -> Outcome:
        if o.stage != self.id or o.outcome not in self.allowed:
            raise ValueError(f"{self.id} may not emit {o.stage}/{o.outcome}")
        return o

    def ready(self, rec: CandidateRecord) -> bool:
        return all(rec.passed(s) for s in self.requires)

    def run(self, rec: CandidateRecord, ctx: Context) -> StageResult:   # batch stages override run_batch
        raise NotImplementedError(f"{self.id} is a batch stage (run_batch) or not implemented yet")

    def decide(self, res: StageResult, th: Thresholds, tgt: Targets) -> Outcome:
        raise NotImplementedError(f"{self.id}.decide")
