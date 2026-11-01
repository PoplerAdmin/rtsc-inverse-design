"""Status enum and Verdict record: the soundness contract (ARCH.md section 3).

Shared by every team.  Only ``rtt.core.certify.classify`` may produce
``Status.REJECTED_PROVEN`` or ``Status.CERTIFIED_CANDIDATE`` (contract S1).
"""
from __future__ import annotations

import json
from dataclasses import asdict, dataclass, field
from enum import Enum
from typing import Any


class Status(Enum):
    REJECTED_PROVEN = "REJECTED_PROVEN"          # Tc_model < T* proved (Lean-backed rule + Arb checks)
    CERTIFIED_CANDIDATE = "CERTIFIED_CANDIDATE"  # Tc_model >= T* proved
    UNDECIDED = "UNDECIDED"                      # certified bracket [T_low, T_up] straddles T*
    EXCLUDED_ASSUMPTION = "EXCLUDED_ASSUMPTION"  # excluded by a physical, non-proved criterion
    DEFERRED = "DEFERRED"                        # not computed (budget / ordering); never dropped
    ERROR = "ERROR"

    def __str__(self) -> str:  # pragma: no cover - cosmetic
        return self.value


#: statuses that claim a theorem about the model; only rtt.core.certify.classify may emit them
PROVEN_STATUSES = frozenset({Status.REJECTED_PROVEN, Status.CERTIFIED_CANDIDATE})


@dataclass
class Verdict:
    material_id: str
    status: Status
    T_low_K: float | None = None          # certified lower bound on Tc_model (None: no bound)
    T_up_K: float | None = None           # certified upper bound on Tc_model (None: no bound)
    reasons: list[str] = field(default_factory=list)       # human-readable
    rule_ids: list[str] = field(default_factory=list)      # "R-UNIV", "R-CELLS", "C-RAYLEIGH", ...
    assumptions: list[str] = field(default_factory=list)   # "M-ISO-MIGDAL", "E-ENV(0.10)", ...
    flags: list[str] = field(default_factory=list)         # out-of-model risks
    certificate: dict = field(default_factory=dict)        # data needed to re-check

    def to_dict(self) -> dict[str, Any]:
        d = asdict(self)
        d["status"] = self.status.value
        return d

    def to_json(self, **kw) -> str:
        return json.dumps(self.to_dict(), **kw)

    @classmethod
    def from_dict(cls, d: dict[str, Any]) -> "Verdict":
        d = dict(d)
        d["status"] = Status(d["status"])
        return cls(**d)

    @classmethod
    def from_json(cls, s: str) -> "Verdict":
        return cls.from_dict(json.loads(s))
