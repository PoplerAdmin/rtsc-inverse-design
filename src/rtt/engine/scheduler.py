"""Scheduler (architecture §5.1, §8): lanes G/S/F, budgets, rank cuts; ``next_jobs``.

The queue table is derived and rebuildable: :func:`rebuild_queue` recomputes it from current
outcomes. A candidate is queued for the first stage it has not reached whose ``requires`` all
passed; controls (lane C) get top priority (R5/R6); otherwise priority = latest score label
(S3 ``score_s3``, S4 ``score_s4``, S5 ``A_coarse/floor``). Per-day budgets come from
machine.toml (capacity x derate) and the per-stage cost estimate by cost class (§8.1).
"""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any, Optional

from ..records import D_PASS, D_UNDECIDED, PASS
from .retry import backoff_ready, failed_attempts, needs_run

PER_CANDIDATE_STAGES = ("S4", "S5", "S6", "S7")
STAGE_LANE = {"S1": "G", "S2": "G", "S3": "G", "S4": "S", "S5": "S", "S6": "F", "S7": "F", "S8": "F"}
PASSING = {PASS, D_PASS, D_UNDECIDED}
REQUIRES = {"S2": "S1", "S3": "S2", "S4": "S3", "S5": "S4", "S6": "S5", "S7": "S6", "S8": "S7"}
#: §8.1 core-hours per candidate by cost class (upper ends of the table; est.)
COST_CORE_H = {"S4": {"C1": 0.2, "C2": 2.0, "C3": 3.0},
               # doc 31 §6: class fallback from the campaign-1 refit (S5_COST_FIT) at 80 Ry (H):
               # C1 = 5 atoms, C2 = 10 atoms; C3 keeps the doc 12 §8.1 upper end (low symmetry is
               # not covered by the fit, which has only cubic/tetragonal cells)
               "S5": {"C1": 0.062 * 5 ** 3 * (80 / 60) ** 1.5, "C2": 0.062 * 10 ** 3 * (80 / 60) ** 1.5,
                      "C3": 700.0 * 1.5},
               "S6": {"C1": 80.0, "C2": 2000.0, "C3": 30000.0},
               "S7": {"C1": 10.0, "C2": 50.0, "C3": 50.0}}

# ------------------------------------------------------------------ S5 cost model (doc 31 §6)
#: Refit from campaign 1 measured S5 costs (doc 30 §5.2/§5.3; q 3^3 + q 4^3, 60 Ry, core-h):
#: ScH2 1.67, LaH2 1.35 (3 atoms), NaIrH3 3.30, BaPtH3 6.47 (5), Mg2IrH6 29.7 (9); elements
#: Al 0.49, Nb 0.84, Pb 0.96 (1 atom). Model: core_h = max(floor, a * N^3) * (ecut/60)^1.5.
#: a = max over the 5 hydrides of core_h / N^3 (0.050, 0.062, 0.026, 0.052, 0.041; geometric mean
#: 0.044) = 0.062: the upper envelope, as the old table used the upper ends. N^3 fits all five
#: hydrides within x2.4 with no separate 5d factor. floor = 1.0 (1-atom elements 0.5-1.0).
#: (ecut/60)^1.5 = plane-wave count scaling when H forces 80 Ry (doc 15/16): x1.54.
S5_COST_FIT = {"a_core_h": 0.062, "exponent": 3.0, "floor_core_h": 1.0, "ecut_ref_Ry": 60.0,
               "ecut_exponent": 1.5, "default_ecut_floor_Ry": {"H": 80.0},
               "source": "doc30#5.2-5.3 measured S5 (campaign 1); doc31#6 refit"}


def s5_ecut_Ry(elements: Any, tiers: Any = None) -> float:
    """Coarse-profile cutoff for these elements: max(profile ecut, tiers S5.ecut_floor_Ry or the
    default H floor 80 Ry)."""
    base = S5_COST_FIT["ecut_ref_Ry"]
    floors = None
    if tiers is not None:
        try:
            floors = tiers.get("S5.ecut_floor_Ry", None)
            from ..dft.profiles import profiles_from_tiers
            base = float(profiles_from_tiers(tiers)["coarse"].ecutwfc_Ry)
        except Exception:  # noqa: BLE001
            pass
    floors = floors if floors else S5_COST_FIT["default_ecut_floor_Ry"]
    els = set(elements or ())
    vals = [float(v) for k, v in dict(floors).items()
            if k in els and isinstance(v, (int, float)) and not isinstance(v, bool)]
    return max([base, *vals])


def estimate_s5_core_h(n_atoms: int, elements: Any = (), *, ecut_Ry: float | None = None, tiers: Any = None) -> float:
    """S5 (q 3^3 + q 4^3) core-hours from the campaign-1 refit :data:`S5_COST_FIT`."""
    f = S5_COST_FIT
    ec = float(ecut_Ry) if ecut_Ry is not None else s5_ecut_Ry(elements, tiers)
    base = max(f["floor_core_h"], f["a_core_h"] * float(n_atoms) ** f["exponent"])
    return base * (ec / f["ecut_ref_Ry"]) ** f["ecut_exponent"]


def estimate_core_h(ledger: Any, cid: str, stage: str, cost_class: str | None, tiers: Any = None) -> float:
    """Per-candidate estimate: the S5 refit from the cell when the structure is known; the class
    table otherwise (and for the other stages)."""
    if stage == "S5":
        try:
            rec = ledger.candidate_record(cid)
            sp = list(rec.structure.species)
            return estimate_s5_core_h(len(sp), set(sp), tiers=tiers)
        except Exception:  # noqa: BLE001
            pass
    return COST_CORE_H[stage].get(cost_class or "C2", COST_CORE_H[stage]["C2"])


SCORE_LABEL = {"S4": "score_s3", "S5": "score_s4", "S6": "A_coarse_over_floor", "S7": "A_coarse_over_floor"}


@dataclass
class Job:
    candidate_id: str
    stage: str
    lane: str
    priority: float
    est_core_h: float = 0.0
    control: bool = False
    meta: dict = field(default_factory=dict)


def lane_budget_core_h(machine: Any, lane: str) -> float:
    """Effective core-h per day of a CPU lane: mpi_ranks x 24 h x contention derate (§8.2)."""
    ranks = float(machine.get(f"lanes.{lane}.mpi_ranks", 0) or 0)
    derate = float(machine.get("derate.contention", 0.85))
    return ranks * 24.0 * derate


def _latest_label(ledger: Any, cid: str, name: str) -> Optional[float]:
    r = ledger.q1("SELECT value FROM labels WHERE candidate_id=? AND name=? ORDER BY rowid DESC LIMIT 1", (cid, name))
    return None if r is None else float(r[0])


def next_stage(ledger: Any, cid: str, tiers: Any = None) -> Optional[str]:
    """First per-candidate stage to (re)run: not reached yet, or a retryable DEFERRED/ERROR within
    the attempt budget (RC1 F12/F13). Only terminal removals and exhausted retries stop a
    candidate. Back-off timing is checked by :func:`next_jobs`."""
    cur = ledger.current_outcomes(cid)
    for st in PER_CANDIDATE_STAGES:
        req = REQUIRES[st]
        if not (req in cur and cur[req]["outcome"] in PASSING):
            return None
        if st not in cur or needs_run(ledger, cid, st, tiers):
            return st
        if cur[st]["outcome"] not in PASSING:
            return None
    return None


def rebuild_queue(ledger: Any, tiers: Any = None) -> int:
    ledger.queue_clear()
    n = 0
    for cid in ledger.candidate_ids():
        st = next_stage(ledger, cid, tiers)
        if st is None:
            continue
        lane = ledger.q1("SELECT lane FROM candidates WHERE id=?", (cid,))[0]
        pr = 1e9 if lane == "C" else (_latest_label(ledger, cid, SCORE_LABEL[st]) or 0.0)
        pr -= failed_attempts(ledger, cid, st)          # retries after first attempts
        if not backoff_ready(ledger, cid, st, tiers):
            continue
        ledger.queue_put(cid, st, pr)
        n += 1
    return n


def next_jobs(ledger: Any, machine: Any, *, stages: tuple[str, ...] = PER_CANDIDATE_STAGES,
              spent_core_h: dict[str, float] | None = None, max_jobs: int | None = None,
              tiers: Any = None) -> list[Job]:
    """Jobs for today within each lane's remaining budget, highest priority first. A job that
    alone exceeds the remaining budget is skipped (it stays queued: DEFERRED by budget is a queue
    state, not an outcome)."""
    rebuild_queue(ledger, tiers)
    spent = dict(spent_core_h or {})
    left = {ln: lane_budget_core_h(machine, ln) - spent.get(ln, 0.0) for ln in ("S", "F")}
    jobs: list[Job] = []
    for cid, st, pr in ledger.queue_items():
        if st not in stages:
            continue
        cc = ledger.q1("SELECT cost_class, lane FROM candidates WHERE id=?", (cid,))
        cost = estimate_core_h(ledger, cid, st, cc[0], tiers)
        lane = STAGE_LANE[st]
        if cost > left[lane] and cc[1] != "C":
            continue
        left[lane] -= cost
        jobs.append(Job(cid, st, lane, pr, cost, control=cc[1] == "C"))
        if max_jobs and len(jobs) >= max_jobs:
            break
    return jobs
