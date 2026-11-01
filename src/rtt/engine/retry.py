"""Re-entry of non-terminal outcomes (review RC1 F12, F13).

DEFERRED and ERROR are not removals (architecture §3.1): a candidate whose current outcome at a
stage is DEFERRED with a retryable reason, or ERROR, is run again at that stage, up to
``[retry] max_attempts`` attempts (tiers.toml, sourced), with exponential back-off between
attempts (``backoff_h * 2**(n-1)`` after the n-th failure; applied by the scheduler, not by
eligibility). Exhausted candidates stay at their last outcome and are listed by
:func:`exhausted` (S9 reports them). Terminal outcomes are only the removals (SCREENED_OUT,
EXCLUDED_ASSUMPTION, INVALID, DUPLICATE, rejectedProven, rejectedAssuming) and DEFERRED with a
non-retryable reason.

Batch runners (S2, S3) and the per-candidate scheduler both use :func:`needs_run`.
"""
from __future__ import annotations

import datetime as _dt
from typing import Any, Optional

#: DEFERRED reasons that are queue states, not failures: re-entered without using the attempt budget
QUEUE_DEFERRED = ("rank_cut", "budget", "awaiting_dependency")
NON_RETRYABLE_DEFERRED = ("uncalibrated_window",)   # deterministic re-run (doc 31 §7); waits for a calibration
DEFAULT_MAX_ATTEMPTS = 3
DEFAULT_BACKOFF_H = 6.0


def _policy(tiers: Any) -> tuple[int, float]:
    if tiers is None:
        return DEFAULT_MAX_ATTEMPTS, DEFAULT_BACKOFF_H
    try:
        n = int(tiers.threshold("retry.max_attempts").value)
        b = float(tiers.threshold("retry.backoff_h").value)
        return n, b
    except Exception:  # noqa: BLE001 - older tiers.toml without [retry]
        return DEFAULT_MAX_ATTEMPTS, DEFAULT_BACKOFF_H


def _reason_text(o: dict) -> str:
    return " ".join(str(o.get(k) or "") for k in ("reason", "criterion", "value"))


def is_retryable(o: Optional[dict]) -> bool:
    if o is None:
        return False
    if o["outcome"] == "ERROR":
        return True
    if o["outcome"] == "DEFERRED":
        txt = _reason_text(o)
        if any(k in txt for k in NON_RETRYABLE_DEFERRED):
            return False
        return True             # DEFERRED is never a removal; unknown reasons fail open (re-run)
    return False


def failed_attempts(ledger: Any, cid: str, stage: str) -> int:
    """Failed attempts (ERROR, or DEFERRED for a failure reason) since the last other outcome."""
    rows = ledger.qa("SELECT outcome, reason, criterion, value FROM outcomes WHERE candidate_id=? AND stage=?"
                     " ORDER BY id DESC", (cid, stage))
    n = 0
    for o, reason, crit, val in rows:
        if o not in ("DEFERRED", "ERROR"):
            break
        txt = f"{reason or ''} {crit or ''} {val or ''}"
        if o == "DEFERRED" and any(k in txt for k in QUEUE_DEFERRED):
            continue                       # queue state, not a failed attempt
        n += 1
    return n


def needs_run(ledger: Any, cid: str, stage: str, tiers: Any = None) -> bool:
    """True if the stage has no outcome yet, or its current outcome is a retryable DEFERRED or
    ERROR within the attempt budget. (Back-off timing is applied by the scheduler.)"""
    cur = ledger.current_outcome(cid, stage)
    if cur is None:
        return True
    if not is_retryable(cur):
        return False
    max_attempts, _ = _policy(tiers)
    return failed_attempts(ledger, cid, stage) < max_attempts


def backoff_ready(ledger: Any, cid: str, stage: str, tiers: Any = None, now: _dt.datetime | None = None) -> bool:
    cur = ledger.current_outcome(cid, stage)
    if cur is None or not is_retryable(cur):
        return True
    _, base = _policy(tiers)
    n = failed_attempts(ledger, cid, stage)
    wait = _dt.timedelta(hours=base * 2 ** max(0, n - 1))
    t = _dt.datetime.fromisoformat(cur["created"])
    now = now or _dt.datetime.now(_dt.timezone.utc)
    return now - t >= wait


def exhausted(ledger: Any, tiers: Any = None) -> list[tuple[str, str, str, int]]:
    """(candidate, stage, outcome, attempts) whose retries are used up (S9 lists them)."""
    out = []
    max_attempts, _ = _policy(tiers)
    for (cid,) in ledger.qa("SELECT id FROM candidates"):
        for st, o in ledger.current_outcomes(cid).items():
            if is_retryable(o):
                n = failed_attempts(ledger, cid, st)
                if n >= max_attempts:
                    out.append((cid, st, o["outcome"], n))
    return out
