"""S5 table: per candidate, the T2-SCREEN outcome, per-mesh A/lambda/omega_log per admitted sigma, costs."""
import json, sqlite3, sys
c = sqlite3.connect(sys.argv[1])
FLOOR = 433148145 / 1e8 * 300.0
for cid, lane in c.execute("SELECT id, lane FROM candidates"):
    o = c.execute("SELECT outcome, criterion, tags, tier_run_id, reason, targets_sha256 FROM outcomes WHERE candidate_id=? "
                  "AND stage='S5' ORDER BY id DESC LIMIT 1", (cid,)).fetchone()
    if not o:
        continue
    trs = c.execute("SELECT id, wall_s, core_h, params_json, state FROM tier_runs WHERE candidate_id=? AND stage='S5' ORDER BY id",
                    (cid,)).fetchall()
    m = json.loads(c.execute("SELECT metrics_json FROM tier_runs WHERE id=?", (o[3],)).fetchone()[0] or "{}") or {}
    meshes = []
    for b in m.get("by_mesh", []):
        meshes.append({"q": b.get("q"), "admitted": b.get("admitted_sigmas"),
                       "A_K": [round(x, 1) for x in b.get("A_K_by_sigma") or []],
                       "lambda": [round(x, 3) for x in b.get("lambda_by_sigma") or []],
                       "omega_log_K": [round(x) for x in b.get("omega_log_K_by_sigma") or []],
                       "A_mid_K": b.get("A_mid_K"), "min_freq_cm1": b.get("min_freq_cm1")})
    v = c.execute("SELECT status, T_low_K, T_up_K, T_star_K FROM verdicts WHERE candidate_id=? ORDER BY id DESC LIMIT 1", (cid,)).fetchone()
    print(json.dumps({"cid": cid, "lane": lane, "S5": o[0], "crit": o[1], "tags": o[2], "reason": o[4],
                      "A_mid_K": m.get("A_mid_K"), "A_lo_K": m.get("A_lo_K"), "A_coarse_over_1299": (m.get("A_mid_K") or 0) / FLOOR,
                      "omega_max_K": m.get("omega_max_K"), "min_freq_cm1": m.get("min_freq_cm1"), "zpe_meV_atom": m.get("zpe_meV_atom"),
                      "n_sheets": m.get("n_sheets"), "meshes": meshes,
                      "runs": [{"id": t[0], "wall_s": t[1], "core_h": t[2], "params": t[3], "state": t[4]} for t in trs],
                      "coarse_verdict": v}))
