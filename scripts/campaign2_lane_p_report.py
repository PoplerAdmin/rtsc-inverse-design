"""Lane P report tables (doc 32): S3 table with composition novelty, S4 selection, S4 table, S5 queue.

    nice -n 19 ~/.rtt/envs/py/bin/python scripts/campaign2_lane_p_report.py \
        --ledger ~/rtsc_targeted_work/P2-lane/ledger.sqlite --out runs/campaign2_laneP [--select 10]

Read-only on the ledger. Outputs (repo, small, no paths): ``s3_table.jsonl``, ``s4_selection.json``,
``s4_table.jsonl``, ``s5_queue.json``, ``funnel.json``.

S4 selection (ranking only, R2): lane-P candidates with an S3 PASS or S3 rank-queue outcome,
ordered by the stored S3 score (desc), ties broken by the larger of the two ML hull distances
(MatterSim-self and DFT-reference; asc; R30 F4 ``ML_BELOW_DFT_HULL`` caution), then by
|N - 19|, then by cid. Candidates whose composition already has a campaign-1 S4 result are
skipped (Mg2IrH6, Ba2IrH6) and listed.

S5 cost model (est., doc 32 §5): anchor Mg2IrH6 = 30 core-h at 60 Ry (campaign 1 measured), x
(ecut/60)^1.5, x V/V0 (plane waves), x nat/9 (perturbations), x n_el/43 (bands), x 48/|G| (irreducible
k and q; non-centrosymmetric cells: low end uses min(48, 2|G|), time reversal). V from the S4 DFT
relax when present, else the MatterSim S3 relax.
"""
from __future__ import annotations

import argparse
import collections
import json
import os
import re
import sqlite3
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

ANCHOR = {"core_h_60Ry": 30.0, "V_A3": 6.652064 ** 3 / 4, "nat": 9, "nel": 43, "nsym": 48}
CAMPAIGN1_S4_DONE = {"Mg2IrH6": "c-9a8fa100a45b (control)", "Ba2IrH6": "c-539acaa4a6e7 (lane D)"}


def zval_table(pseudo_dir: Path) -> dict[str, float]:
    out = {}
    for f in pseudo_dir.glob("*_ONCV_PBE-1.2.upf"):
        m = re.search(r'z_valence="\s*([\d.Ee+-]+)"', f.read_text(errors="replace")[:20000])
        if m:
            out[f.name.split("_")[0]] = float(m.group(1))
    return out


def nsym(sd: dict) -> int | None:
    """Point-group order; for a non-centrosymmetric cell the k-mesh reduction also uses time
    reversal, so the effective order is min(48, 2|G|) (returned negative-free as the lower-cost
    bound by :func:`s5_cost_range`)."""
    from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
    from rtt.structure.canonical import to_pymatgen
    try:
        sga = SpacegroupAnalyzer(to_pymatgen(sd), symprec=0.01)
        n = len(sga.get_point_group_operations())
        return n if sga.is_laue() else -n
    except Exception:  # noqa: BLE001
        return None


def s5_cost(nat: int, vol: float, nel: float, ns: int | None, ecut: float) -> float:
    a = ANCHOR
    return round(a["core_h_60Ry"] * (ecut / 60.0) ** 1.5 * (vol / a["V_A3"]) * (nat / a["nat"])
                 * (nel / a["nel"]) * (a["nsym"] / (ns or 1)), 1)


def s5_cost_range(nat: int, vol: float, nel: float, ns: int | None, ecut: float) -> tuple[float, float]:
    """(low, high): high uses 48/|G|; low credits time reversal for non-centrosymmetric cells."""
    if ns is None:
        return (None, None)
    g = abs(ns)
    lo_g = g if ns > 0 else min(48, 2 * g)
    return (s5_cost(nat, vol, nel, lo_g, ecut), s5_cost(nat, vol, nel, g, ecut))


def comp_novel(formula_red: str, chemsys: str, cache: Path) -> bool | None:
    from pymatgen.core import Composition
    f = cache / f"ref_{chemsys}.json"
    if not f.exists():
        return None
    fr = Composition(formula_red).reduced_formula
    rows = json.loads(f.read_text())["rows"]
    return not any(Composition(r["formula"]).reduced_formula == fr for r in rows)


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--ledger", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--select", type=int, default=10)
    ap.add_argument("--work", default="~/rtsc_targeted_work/P2-lane/work", help="S4 workdir (relaxed cells)")
    ap.add_argument("--cache", default="~/rtt_scratch/front_cache/alexmp_ref")
    a = ap.parse_args()
    out = Path(a.out)
    out.mkdir(parents=True, exist_ok=True)
    L = sqlite3.connect(os.path.expanduser(a.ledger))
    zv = zval_table(Path(os.path.expanduser("~/pseudo/sg15")))
    cache = Path(os.path.expanduser(a.cache))

    def cur(cid, st):
        r = L.execute("SELECT outcome, criterion, value, reason, tags FROM outcomes WHERE candidate_id=? AND stage=? "
                      "ORDER BY id DESC LIMIT 1", (cid, st)).fetchone()
        return None if r is None else dict(zip(("outcome", "criterion", "value", "reason", "tags"), r))

    def metrics(cid, st):
        r = L.execute("SELECT metrics_json, core_h, wall_s FROM tier_runs WHERE candidate_id=? AND stage=? AND "
                      "state='done' ORDER BY id DESC LIMIT 1", (cid, st)).fetchone()
        return (json.loads(r[0] or "{}"), r[1], r[2]) if r else ({}, None, None)

    rows = []
    for cid, lane, oj, sj in L.execute("SELECT c.id, c.lane, c.origin_json, s.json FROM candidates c "
                                       "JOIN structures s ON s.sha256=c.structure_sha256"):
        o, sd = json.loads(oj), json.loads(sj)
        cond = o.get("conditioning") or {}
        species = sd["species"]
        cnt = collections.Counter(species)
        from pymatgen.core import Composition
        comp = Composition(dict(cnt))
        chemsys = "-".join(sorted(cnt))
        r = {"cid": cid, "lane": lane, "formula": cond.get("formula") or comp.reduced_formula,
             "prototype": cond.get("prototype"), "substitution": cond.get("substitution"),
             "N": cond.get("N_per_complex"), "parent_N": cond.get("parent_N"),
             "parity": cond.get("valence_parity"), "nat": len(species)}
        for st in ("S1", "S2", "S3", "S4"):
            c = cur(cid, st)
            r[st] = c and c["outcome"]
            if c and st == "S2" and c["value"]:
                v = json.loads(c["value"])
                r["novel_structure"] = v.get("novel")
                r["alexmp_match"] = v.get("alexmp_match")
            if c and st in ("S3", "S4"):
                r[f"{st}_criterion"] = c["criterion"]
                r[f"{st}_tags"] = c["tags"]
        r["novel_composition"] = comp_novel(comp.reduced_formula, chemsys, cache)
        m3, _, _ = metrics(cid, "S3")
        r["e_hull_ml_self"] = m3.get("e_hull_ml_self_eV")
        r["e_hull_ml_dftref"] = m3.get("e_hull_ml_dftref_eV")
        r["score_s3"] = m3.get("score_s3")
        rl = m3.get("relax") or {}
        r["relax_converged"] = rl.get("converged")
        r["vol_ml_A3"] = (rl.get("volume_per_atom_A3") or 0) * len(species) or None
        r["ml_below_dft_hull"] = (r["e_hull_ml_dftref"] is not None and r["e_hull_ml_dftref"] < -0.05)
        if m3.get("relaxed_structure"):
            from rtt.structure.canonical import spacegroup_number
            r["sg_relaxed_0.01"] = spacegroup_number(m3["relaxed_structure"], 0.01)
            r["sg_relaxed_0.1"] = spacegroup_number(m3["relaxed_structure"], 0.1)
        r["nsym"] = nsym(sd)
        r["n_el_sg15"] = sum(zv.get(s, 0) for s in species)
        m4, ch4, w4 = metrics(cid, "S4")
        if m4:
            r["S4_metrics"] = {k: m4.get(k) for k in ("gap_eV", "N_EF_per_eV_atom", "moment_per_atom_muB",
                                                      "n_sheets", "H_DOS", "e_form_eV", "scf_pressure_kbar")}
            r["S4_core_h"], r["S4_wall_s"] = ch4, w4
            rs = Path(os.path.expanduser(a.work)) / cid / "S4" / "relaxed_structure.json"
            if rs.exists():
                import numpy as np
                r["vol_dft_A3"] = abs(float(np.linalg.det(np.array(json.loads(rs.read_text())["lattice_A"]))))
        rows.append(r)

    with open(out / "s3_table.jsonl", "w") as fh:
        for r in sorted(rows, key=lambda r: (r["lane"] != "C", -(r["score_s3"] or 0), r["cid"])):
            fh.write(json.dumps(r, default=str) + "\n")

    # S4 selection
    def hmax(r):
        xs = [x for x in (r["e_hull_ml_self"], r["e_hull_ml_dftref"]) if x is not None]
        return max(xs) if xs else 9.0
    pool = [r for r in rows if r["lane"] == "P" and r["S3"] in ("PASS", "DEFERRED") and
            (r.get("S3_criterion") in (None, "rank_cut"))]
    pool.sort(key=lambda r: (-(r["score_s3"] or 0), hmax(r), abs((r["N"] or 19) - 19), r["cid"]))
    sel, skipped = [], []
    for r in pool:
        if r["formula"] in CAMPAIGN1_S4_DONE:
            skipped.append((r["cid"], r["formula"], CAMPAIGN1_S4_DONE[r["formula"]]))
            continue
        if len(sel) < a.select:
            sel.append(r)
    (out / "s4_selection.json").write_text(json.dumps(
        {"rule": "score_s3 desc, max(e_hull_ml_self, e_hull_ml_dftref) asc, |N-19| asc, cid",
         "selected": [{k: r[k] for k in ("cid", "formula", "prototype", "N", "score_s3", "e_hull_ml_self",
                                          "e_hull_ml_dftref", "novel_structure", "novel_composition", "S3")}
                      for r in sel],
         "skipped_campaign1_done": skipped, "pool_size": len(pool)}, indent=1, default=str))

    # S4 table + S5 queue
    s4 = [r for r in rows if r.get("S4")]
    with open(out / "s4_table.jsonl", "w") as fh:
        for r in s4:
            fh.write(json.dumps({k: r.get(k) for k in ("cid", "lane", "formula", "prototype", "N", "parity", "S4",
                                                       "S4_criterion", "S4_tags", "S4_metrics", "S4_core_h",
                                                       "S4_wall_s", "novel_structure", "novel_composition")},
                                default=str) + "\n")
    q = []
    for r in s4:
        if r["S4"] != "PASS":
            continue
        vol = r.get("vol_dft_A3") or r.get("vol_ml_A3")
        if not vol:
            continue
        m = r.get("S4_metrics") or {}
        q.append({"cid": r["cid"], "lane": r["lane"], "formula": r["formula"], "prototype": r["prototype"],
                  "N": r["N"], "nat": r["nat"], "nsym": r["nsym"], "vol_A3": round(vol, 2),
                  "n_el": r["n_el_sg15"], "N_EF_per_eV_atom": m.get("N_EF_per_eV_atom"), "H_DOS": m.get("H_DOS"),
                  "n_sheets": m.get("n_sheets"), "e_hull_ml_self": r["e_hull_ml_self"],
                  "e_hull_ml_dftref": r["e_hull_ml_dftref"], "novel_composition": r["novel_composition"],
                  "core_h_60Ry_est_range": s5_cost_range(r["nat"], vol, r["n_el_sg15"], r["nsym"], 60.0),
                  "core_h_80Ry_est_range": s5_cost_range(r["nat"], vol, r["n_el_sg15"], r["nsym"], 80.0)})
    q.sort(key=lambda x: (-(x["N_EF_per_eV_atom"] or 0), -(x["H_DOS"] or 0), x["cid"]))
    for i, x in enumerate(q, 1):
        x["rank"] = i
    (out / "s5_queue.json").write_text(json.dumps({"cost_model": __doc__.split("S5 cost model")[1].strip(),
                                                   "anchor": ANCHOR, "queue": q}, indent=1, default=str))
    funnel = {"S1": collections.Counter(r["S1"] for r in rows if r["lane"] == "P"),
              "S2": collections.Counter(r["S2"] for r in rows if r["lane"] == "P"),
              "S3": collections.Counter(f"{r['S3']}({r.get('S3_criterion')})" for r in rows if r["lane"] == "P"),
              "S4": collections.Counter(f"{r['S4']}({r.get('S4_criterion')})" for r in rows if r.get("S4")),
              "novel_structure": collections.Counter(str(r.get("novel_structure")) for r in rows if r["lane"] == "P"),
              "novel_composition": collections.Counter(str(r.get("novel_composition")) for r in rows if r["lane"] == "P"),
              "by_prototype": collections.Counter(r["prototype"] for r in rows if r["lane"] == "P")}
    (out / "funnel.json").write_text(json.dumps(funnel, indent=1))
    print(json.dumps(funnel, indent=1))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
