"""Space groups of the S3 (MatterSim-relaxed) cells per lane at symprec 0.1 A (spglib)."""
import collections, json, sqlite3, sys
import spglib
c = sqlite3.connect(sys.argv[1])
S6 = {191, 194, 139, 123, 166, 221, 223, 225, 227, 229}
out = collections.defaultdict(collections.Counter)
for cid, lane in c.execute("SELECT id, lane FROM candidates"):
    o = c.execute("SELECT tier_run_id FROM outcomes WHERE candidate_id=? AND stage='S3' ORDER BY id DESC LIMIT 1", (cid,)).fetchone()
    if not o or not o[0]:
        continue
    m = json.loads(c.execute("SELECT metrics_json FROM tier_runs WHERE id=?", (o[0],)).fetchone()[0] or "{}") or {}
    rs = m.get("relaxed_structure")
    if not rs:
        continue
    nums = [sorted(set(rs["species"])).index(s) + 1 for s in rs["species"]]
    try:
        sg = spglib.get_symmetry_dataset((rs["lattice"], rs["frac_coords"], nums), symprec=0.1).number
    except Exception:
        sg = None
    out[lane]["n"] += 1
    out[lane]["P1"] += sg == 1
    out[lane]["S6_list"] += sg in S6
    out[lane]["cubic"] += (sg or 0) >= 195
    out[lane]["sg>=75"] += (sg or 0) >= 75
print(json.dumps(out, indent=1))
