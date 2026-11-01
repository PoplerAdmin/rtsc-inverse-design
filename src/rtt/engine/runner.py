"""Run one scheduled job (architecture §5.1): ``run_job(job, ctx)``.

Every stage writes its ``tier_runs`` row as 'started' before launching compute (W12), then the
outcome. Errors become ERROR outcomes, never silent drops (R5)."""
from __future__ import annotations

import traceback
from typing import Any

from ..records import ERROR, Outcome
from .retry import failed_attempts
from .scheduler import Job
from .stage import Context
from .stages import STAGES, write_outcome


def run_job(job: Job, ctx: Context) -> Any:
    st = STAGES[job.stage]
    rec = ctx.ledger.candidate_record(job.candidate_id)
    tgt = ctx.targets_for(rec)
    if job.stage == "S6":
        res = st.run(rec, ctx)
        return st.certify_and_decide(rec, ctx, res)
    if job.stage == "S7":
        return st.redecide(rec, ctx)
    try:
        res = st.run(rec, ctx)
        o = st.decide(res, ctx.tiers, tgt)
    except Exception:  # noqa: BLE001
        o = Outcome(job.stage, ERROR, reason=traceback.format_exc(limit=3)[:500])
        res = None
    if o.outcome in (ERROR, "DEFERRED"):
        o.retry = failed_attempts(ctx.ledger, rec.candidate_id, job.stage) + 1      # RC1 F13
    write_outcome(ctx, rec.candidate_id, o, rec=rec)
    if job.stage == "S5" and res is not None and res.ok:
        try:
            st.evidence(rec, ctx, res)
        except Exception:  # noqa: BLE001 - evidence only; recorded in the log, never changes the outcome
            traceback.print_exc()
    return o
