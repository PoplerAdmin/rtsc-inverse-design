"""Funnel: current outcome per stage x lane (counts), from a campaign ledger."""
import collections, json, sqlite3, sys
c = sqlite3.connect(sys.argv[1])
F = collections.defaultdict(collections.Counter)
for cid, lane in c.execute("SELECT id, lane FROM candidates"):
    for st in ("S1", "S2", "S3", "S4", "S5", "S6"):
        o = c.execute("SELECT outcome, criterion FROM outcomes WHERE candidate_id=? AND stage=? ORDER BY id DESC LIMIT 1",
                      (cid, st)).fetchone()
        if o:
            F[(st, lane)][o[0] + (f"({o[1]})" if o[1] else "")] += 1
for k in sorted(F):
    print(k[0], k[1], dict(F[k]))
