"""Per-lane generation statistics of a campaign ledger (S1-S3): composition, H fraction,
space groups, size, on-target fraction, novelty, S3 E_hull. Prints JSON."""
import collections, json, sqlite3, statistics, sys
c = sqlite3.connect(sys.argv[1])
TM = {"Rh", "Pd", "Ir", "Pt", "Os", "Au", "Ag"}
DON = {"Li", "Na", "Mg", "Ca", "Sr", "Ba", "Al", "Sc", "Y", "La"}
CRYS = [(2, "triclinic"), (15, "monoclinic"), (74, "orthorhombic"), (142, "tetragonal"), (167, "trigonal"),
        (194, "hexagonal"), (230, "cubic")]
S6 = {191, 194, 139, 123, 166, 221, 223, 225, 227, 229}


def cur(cid, st):
    return c.execute("SELECT outcome, criterion, value, tier_run_id FROM outcomes WHERE candidate_id=? AND stage=? "
                     "ORDER BY id DESC LIMIT 1", (cid, st)).fetchone()


lanes = collections.defaultdict(list)
for cid, lane, oj, sha, cc in c.execute("SELECT id, lane, origin_json, structure_sha256, cost_class FROM candidates"):
    st = json.loads(c.execute("SELECT json FROM structures WHERE sha256=?", (sha,)).fetchone()[0])
    o = json.loads(oj)
    sp = st["species"]
    els = set(sp)
    s1, s2, s3 = cur(cid, "S1"), cur(cid, "S2"), cur(cid, "S3")
    eh = None
    if s3 and s3[3]:
        m = json.loads(c.execute("SELECT metrics_json FROM tier_runs WHERE id=?", (s3[3],)).fetchone()[0] or "{}") or {}
        eh = m.get("e_hull_ml_self_eV")
    nov = None
    if s2 and s2[0] == "PASS" and s2[2]:
        try:
            nov = json.loads(s2[2]).get("novel")
        except Exception:
            pass
    lanes[lane].append({"els": els, "n": len(sp), "Hf": sp.count("H") / len(sp), "sg": st.get("spacegroup"),
                        "flags": o.get("flags", []), "cc": cc, "s1": s1 and (s1[0], s1[1]), "s2": s2 and s2[0],
                        "s3": s3 and (s3[0], s3[1]), "novel": nov, "eh": eh, "nel": len(els)})
out = {}
for lane, R in sorted(lanes.items()):
    sgs = [r["sg"] for r in R if r["sg"]]
    hs = [r["Hf"] for r in R]
    ehs = [r["eh"] for r in R if r["eh"] is not None]
    cs = collections.Counter()
    for r in R:
        for lo, name in CRYS:
            if r["sg"] and r["sg"] <= lo:
                cs[name] += 1
                break
    d = {"n": len(R),
         "with_H": sum(1 for r in R if "H" in r["els"]),
         "donor_TM_H_ternary": sum(1 for r in R if "H" in r["els"] and r["els"] & DON and r["els"] & TM and r["nel"] == 3),
         "off_requested_system": sum(1 for r in R if "OFF_REQUESTED_SYSTEM" in r["flags"]),
         "subsystem_of_requested": sum(1 for r in R if "SUBSYSTEM_OF_REQUESTED" in r["flags"]),
         "H_f_mean": round(statistics.mean(hs), 3) if hs else None,
         "H_f_median": round(statistics.median(hs), 3) if hs else None,
         "H_f_ge_0.5": sum(1 for h in hs if h >= 0.5), "H_f_ge_0.6": sum(1 for h in hs if h >= 0.6),
         "n_atoms_median": statistics.median([r["n"] for r in R]),
         "n_elements": dict(collections.Counter(r["nel"] for r in R)),
         "cost_class": dict(collections.Counter(r["cc"] for r in R)),
         "P1": sum(1 for s in sgs if s == 1), "S6_list_groups": sum(1 for s in sgs if s in S6),
         "crystal_system": dict(cs), "top_sg": collections.Counter(sgs).most_common(6),
         "S1": dict(collections.Counter(str(r["s1"]) for r in R)),
         "S2": dict(collections.Counter(str(r["s2"]) for r in R)),
         "novel_of_S2pass": [sum(1 for r in R if r["novel"] is True), sum(1 for r in R if r["s2"] == "PASS")],
         "S3": dict(collections.Counter(str(r["s3"]) for r in R)),
         "ehull_ml_median": round(statistics.median(ehs), 3) if ehs else None,
         "ehull_ml_le_0.1": sum(1 for e in ehs if e <= 0.1), "n_ehull": len(ehs)}
    out[lane] = d
print(json.dumps(out, indent=1, default=str))
