"""Daily summary of a ledger (architecture §5.1 ``daily_summary``; R5, R8, §3.5, §8.3).

Accounts for every candidate (current outcome per stage), splits core REJECTED_PROVEN by
spectrum tier (R5), lists provisional outcomes (§3.4: S9 lists them), the compute spent per
stage/profile (actual funnel vs the §8.3 assumptions) and the control recall per stage at the
surrogate floors (R8). Read-only.
"""
from __future__ import annotations

import json
from typing import Any

REMOVING = {"SCREENED_OUT", "EXCLUDED_ASSUMPTION", "INVALID", "DUPLICATE", "rejectedProven",
            "rejectedAssuming", "screenedOut"}


PASSING = {"PASS", "pass", "undecided"}


def control_recall(ledger: Any) -> dict[str, dict[str, int]]:
    """Per stage: lane-C controls that reached it; ``passed`` counts only PASS-like outcomes
    (PASS, pass, undecided). DEFERRED and ERROR are counted separately and are NOT passes
    (review RC1 F14): a stalled control is a recall failure, so R8 can fire."""
    out: dict[str, dict[str, int]] = {}
    for (cid,) in ledger.qa("SELECT id FROM candidates WHERE lane='C'"):
        for stage, o in ledger.current_outcomes(cid).items():
            d = out.setdefault(stage, {"passed": 0, "total": 0, "removed": 0})
            d["total"] += 1
            if o["outcome"] in PASSING:
                d["passed"] += 1
            elif o["outcome"] in REMOVING:
                d["removed"] += 1
            elif o["outcome"] == "ERROR":
                d["error"] = d.get("error", 0) + 1         # keys present only when nonzero
            else:
                d["deferred"] = d.get("deferred", 0) + 1
    return out


def provisional_outcomes(ledger: Any) -> list[dict]:
    rows = ledger.qa("SELECT candidate_id, stage, outcome, criterion, source_id, tags FROM outcomes "
                     "WHERE provisional=1 OR (tags IS NOT NULL AND tags LIKE '%PROVISIONAL%') "
                     "ORDER BY id")
    return [{"candidate_id": r[0], "stage": r[1], "outcome": r[2], "criterion": r[3],
             "source_id": r[4], "tags": json.loads(r[5]) if r[5] else []} for r in rows]


def compute_spent(ledger: Any) -> list[tuple]:
    return ledger.qa("SELECT stage, COALESCE(profile,''), state, COUNT(*), COALESCE(SUM(wall_s),0), "
                     "COALESCE(SUM(core_h),0), COALESCE(SUM(gpu_h),0) FROM tier_runs "
                     "GROUP BY stage, profile, state ORDER BY stage, profile, state")


def daily_summary(ledger: Any) -> str:
    c = ledger.summary_counts()
    lines = ["# rtt daily summary", "", f"Candidates: {c['candidates']}", "",
             "## Current outcome per stage (R5)", "", "| stage | outcome | n |", "|---|---|---|"]
    for stage in sorted(c["outcomes"]):
        for oc, n in sorted(c["outcomes"][stage].items()):
            lines.append(f"| {stage} | {oc} | {n} |")
    vt = c.get("core_rejected_proven_by_tier") or {}
    lines += ["", "Core REJECTED_PROVEN by spectrum tier (evidence, not outcomes): "
              + (", ".join(f"{k}: {v}" for k, v in sorted(vt.items())) or "none")]
    lines += ["", "## Compute spent", "", "| stage | profile | state | runs | wall h | core h | GPU h |",
              "|---|---|---|---|---|---|---|"]
    for st, pr, state, n, wall, ch, gh in compute_spent(ledger):
        lines.append(f"| {st} | {pr} | {state} | {n} | {wall / 3600:.2f} | {ch:.2f} | {gh:.2f} |")
    prov = provisional_outcomes(ledger)
    lines += ["", f"## Outcomes made with provisional thresholds ({len(prov)})", ""]
    for p in prov:
        lines.append(f"- {p['candidate_id']} {p['stage']} {p['outcome']} "
                     f"{p['criterion'] or ''} ({p['source_id'] or ''}) {' '.join(p['tags'])}".rstrip())
    cr = control_recall(ledger)
    lines += ["", "## Control recall at surrogate floors (R8)", "",
              "| stage | controls passed | reached | removed | deferred | error |",
              "|---|---|---|---|---|---|"]
    for stage, d in sorted(cr.items()):
        flag = "  **recall < 1**" if d["passed"] < d["total"] else ""
        lines.append(f"| {stage} | {d['passed']} | {d['total']} | {d['removed']} | "
                     f"{d.get('deferred', 0)} | {d.get('error', 0)}{flag} |")
    if not cr:
        lines.append("| - | 0 | 0 | 0 | 0 | 0 |")
    lines += ["", "At T* = 300 K about zero positives are expected, so a false-negative *rate* at "
              "the real target cannot be estimated; the audits give proxy error distributions and "
              "the controls give recall at the surrogate floors (§3.5)."]
    return "\n".join(lines) + "\n"
