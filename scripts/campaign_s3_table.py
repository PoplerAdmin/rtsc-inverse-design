"""S3 table of a campaign ledger: one JSON line per S3-computed candidate (score, hull, features)."""
import json, sqlite3, sys
c = sqlite3.connect(sys.argv[1])
for cid, lane, flags_json, sha, cc in c.execute("SELECT id, lane, origin_json, structure_sha256, cost_class FROM candidates"):
    o = c.execute("SELECT outcome, criterion, tier_run_id, tags FROM outcomes WHERE candidate_id=? AND stage='S3' ORDER BY id DESC LIMIT 1", (cid,)).fetchone()
    if not o:
        continue
    m = json.loads(c.execute("SELECT metrics_json FROM tier_runs WHERE id=?", (o[2],)).fetchone()[0] or "{}") or {}
    rs = m.get("relaxed_structure") or {}
    f = m.get("features") or {}
    ph = m.get("phonons") or {}
    st = json.loads(c.execute("SELECT json FROM structures WHERE sha256=?", (sha,)).fetchone()[0])
    org = json.loads(flags_json or "{}")
    print(json.dumps({"cid": cid, "lane": lane, "S3": o[0], "crit": o[1], "tags": o[3], "score": m.get("score_s3"),
                      "ehull": m.get("e_hull_ml_self_eV"), "ehull_dft": m.get("e_hull_ml_dftref_eV"),
                      "species": st.get("species"), "sg": st.get("spacegroup"), "n": len(st.get("species", [])),
                      "H_f": f.get("H_f"), "motif": f.get("motif"), "minf": ph.get("min_freq_cm1"),
                      "o2l": ph.get("omega2_light_K"), "cond": org.get("conditioning"), "flags": org.get("flags"),
                      "conv": (m.get("relax") or {}).get("converged"), "cc": cc}))
