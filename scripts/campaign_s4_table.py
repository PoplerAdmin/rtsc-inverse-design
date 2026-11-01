"""S4 table: current S4 outcome + key metrics per candidate (JSON lines)."""
import json, sqlite3, sys
c = sqlite3.connect(sys.argv[1])
keys = ("gap_eV", "N_EF_per_eV_atom", "moment_per_atom_muB", "n_sheets", "sheet_w_min", "E_F_eV", "H_DOS", "phi",
        "Tc_phi_K", "e_form_eV", "H_f", "scf_pressure_kbar")
for cid, lane, cc in c.execute("SELECT id, lane, cost_class FROM candidates"):
    o = c.execute("SELECT outcome, criterion, value, tags, tier_run_id, reason FROM outcomes WHERE candidate_id=? AND stage='S4' ORDER BY id DESC LIMIT 1", (cid,)).fetchone()
    if not o:
        continue
    tr = c.execute("SELECT metrics_json, wall_s, core_h, state FROM tier_runs WHERE id=?", (o[4],)).fetchone() if o[4] else None
    m = json.loads((tr or [None])[0] or "{}") or {}
    print(json.dumps({"cid": cid, "lane": lane, "cc": cc, "S4": o[0], "crit": o[1], "value": o[2], "tags": o[3],
                      "reason": (o[5] or "")[:160], "formula": m.get("formula"),
                      **{k: m.get(k) for k in keys}, "wall_s": tr and tr[1], "core_h": tr and tr[2]}))
