"""Final ranking (architecture §3.6). Lexicographic, never by T_up:

1. decide outcome: pass > undecided > screenedOut > rejectedAssuming > rejectedProven
   (screenedOut at S6 = "a floor holds on lower-tier inputs", Lean decideT v2);
2. certified T_low: lower end of the combined core bracket at mu*_hi = 0.13 (missing = 0);
3. A_mid: median over admitted sigma, finest q mesh, of A_nom (no envelope).
T_up, A_hi and A_lo are carried along as reported ceilings/floors only.
"""
from __future__ import annotations

import json
import statistics
from dataclasses import asdict, dataclass, field
from typing import Any, Iterable, Optional

OUTCOME_ORDER = {"pass": 4, "undecided": 3, "screenedOut": 2, "rejectedAssuming": 1, "rejectedProven": 0}


@dataclass
class RankRow:
    candidate_id: str
    formula: str = ""
    decide_outcome: str = "undecided"
    firing_rule: Optional[str] = None
    T_low_K: Optional[float] = None
    T_up_K: Optional[float] = None
    A_mid_K: Optional[float] = None
    A_lo_K: Optional[float] = None
    verdict_status: Optional[str] = None
    flags: list[str] = field(default_factory=list)
    rank: int = 0

    def key(self):
        return (-OUTCOME_ORDER.get(self.decide_outcome, -1), -(self.T_low_K or 0.0), -(self.A_mid_K or 0.0),
                self.candidate_id)

    def to_json(self) -> dict:
        return asdict(self)


def final_rank(records: Iterable[Any]) -> list[RankRow]:
    """``records``: RankRow objects or dicts with RankRow field names."""
    rows = [r if isinstance(r, RankRow) else RankRow(**{k: v for k, v in dict(r).items()
                                                       if k in RankRow.__dataclass_fields__}) for r in records]
    rows.sort(key=RankRow.key)
    for i, r in enumerate(rows, 1):
        r.rank = i
    return rows


def rows_from_ledger(ledger: Any, mu_hi: float = 0.13) -> list[RankRow]:
    """One row per candidate that reached S6 (current S6/S7/S8 decision; T_low from the
    full/converged core verdicts at ``mu_hi``, min over verdicts of the latest run)."""
    rows = []
    for (cid,) in ledger.qa("SELECT DISTINCT candidate_id FROM decisions WHERE stage IN ('S6','S7','S8')"):
        d = ledger.q1("SELECT outcome, firing_rule FROM decisions WHERE candidate_id=? AND stage IN ('S6','S7','S8')"
                      " ORDER BY id DESC LIMIT 1", (cid,))
        cur = None
        for st in ("S8", "S7", "S6"):
            cur = ledger.current_outcome(cid, st)
            if cur is not None:
                break
        outcome = cur["outcome"] if cur and cur["outcome"] in OUTCOME_ORDER else d[0]
        vs = ledger.qa("SELECT status, T_low_K, T_up_K FROM verdicts WHERE candidate_id=? AND spectrum_tier IN"
                       " ('full','converged') AND abs(mu_hi - ?) < 1e-9 ORDER BY id DESC", (cid, mu_hi))
        T_low = min((v[1] for v in vs if v[1] is not None), default=None)
        T_up = max((v[2] for v in vs if v[2] is not None), default=None)
        A_mid = A_lo = None
        tr = ledger.q1("SELECT metrics_json FROM tier_runs WHERE candidate_id=? AND kind IN ('dfpt','epw') AND"
                       " profile IN ('full','epw') AND state='done' ORDER BY id DESC LIMIT 1", (cid,))
        if tr and tr[0]:
            m = json.loads(tr[0])
            As = (m.get("ladder") or {}).get("fine", {}).get("A_K") or m.get("A_K_by_sigma")
            if As:
                A_mid, A_lo = float(statistics.median(As)), float(min(As))
        rec = ledger.candidate_record(cid)
        rows.append(RankRow(candidate_id=cid, formula=rec.formula, decide_outcome=outcome, firing_rule=d[1],
                            T_low_K=T_low, T_up_K=T_up, A_mid_K=A_mid, A_lo_K=A_lo,
                            verdict_status=vs[0][0] if vs else None))
    return final_rank(rows)
