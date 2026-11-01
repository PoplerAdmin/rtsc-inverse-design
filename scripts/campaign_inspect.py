"""Dump a per-candidate summary of a campaign ledger (all stages) as JSON lines, for reports.

    ~/.rtt/envs/py/bin/python scripts/campaign_inspect.py LEDGER [--stage S5] [--cids ...]
"""
from __future__ import annotations

import argparse
import json
import sqlite3
import sys


def metrics(c, tr):
    if tr is None:
        return {}
    r = c.execute("SELECT metrics_json, wall_s, core_h, gpu_h FROM tier_runs WHERE id=?", (tr,)).fetchone()
    if not r:
        return {}
    m = json.loads(r[0] or "{}") or {}
    m["_wall_s"], m["_core_h"], m["_gpu_h"] = r[1], r[2], r[3]
    return m


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("ledger")
    ap.add_argument("--stage", default=None)
    ap.add_argument("--cids", nargs="*")
    a = ap.parse_args()
    c = sqlite3.connect(a.ledger)
    cids = a.cids or [r[0] for r in c.execute("SELECT id FROM candidates ORDER BY rowid")]
    for cid in cids:
        lane, = c.execute("SELECT lane FROM candidates WHERE id=?", (cid,)).fetchone()
        rows = c.execute("SELECT stage, outcome, criterion, value, tags, tier_run_id, reason FROM outcomes "
                         "WHERE candidate_id=? ORDER BY id", (cid,)).fetchall()
        cur = {}
        for st, o, crit, val, tags, tr, reason in rows:
            cur[st] = {"outcome": o, "criterion": crit, "value": val, "tags": tags, "tr": tr,
                       "reason": (reason or "")[:200]}
        if a.stage and a.stage not in cur:
            continue
        out = {"cid": cid, "lane": lane, "stages": {k: {kk: vv for kk, vv in v.items() if vv} for k, v in cur.items()}}
        if a.stage:
            m = metrics(c, cur[a.stage]["tr"])
            m.pop("relaxed_structure", None)
            out["metrics"] = m
        print(json.dumps(out, default=str))
    return 0


if __name__ == "__main__":
    sys.exit(main())
