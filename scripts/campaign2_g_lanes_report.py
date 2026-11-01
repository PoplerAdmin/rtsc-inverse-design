"""Doc 34: per-lane "is this meaningful MatterGen output" statistics across ledgers, and the
ranked S4 queue of the G2/G3/G4-bulk lanes. All ledgers are opened read-only.

    ~/.rtt/envs/py/bin/python scripts/campaign2_g_lanes_report.py \
        --ledger G=~/rtsc_targeted_work/G-lanes/ledger.sqlite \
        --ledger E1=~/rtsc_targeted_work/E1-campaign/ledger.sqlite \
        --ledger P=~/rtsc_targeted_work/P2-lane/ledger.sqlite \
        --timing G=~/rtsc_targeted_work/G-lanes/timing_G.jsonl \
        --timing E1=~/rtsc_targeted_work/E1-campaign/timing_front.jsonl \
        --queue-from G --out runs/campaign2_G

Lane group key: ``lane`` for D/C/P/G1/G5, ``lane:checkpoint`` for G2/G3/G4 (G3 has two
checkpoints; campaign-1 G4 is ``chemical_system``, this campaign's is ``ml_bulk_modulus``).

Columns (all per lane group; "of S2 PASS" unless said otherwise):
* validity: S1 PASS, S1 no_pseudo, S2 INVALID / DUPLICATE;
* novelty vs Alex-MP (S2 structure match), novel composition (reduced formula absent from the
  Alex-MP references of the system), and a **training-set proxy**: the matched Alex-MP entry has
  <= 20 atoms and E_hull,DFT <= 0.1 eV/atom (the alex_mp_20 selection rule of the MatterGen
  paper, recalled, UNVERIFIED; the alex_mp_20 file itself is a git-LFS pointer here);
* H content, elements, the doc 02 S10 palette;
* metallicity prior (rtt.score.metallicity, on the PRIMITIVE cell's counts so parity is exact);
* E_hull,ML-self (S3), and the doc 02 hard conditions checkable at S3:
  H1 by prior (p_metal > 0.5), H7 by ML hull (<= 0.35 eV/atom; no ZPE), H in the cell.
* G2 only: requested space group recovered (raw cell 0.1 A; MatterSim-relaxed cell 0.1 A).
"""
from __future__ import annotations

import argparse
import collections
import json
import sqlite3
import statistics
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

PALETTE_DON = {"Li", "Na", "Mg", "Ca", "Sr", "Ba", "Al", "Sc", "Y", "La"}
PALETTE_TM = {"Rh", "Pd", "Ir", "Pt", "Os", "Au", "Ag"}
LIGHT = {"H", "Li", "Be", "B", "C", "N"}
MAGNETIC_3D = {"Cr", "Mn", "Fe", "Co", "Ni"}


def ro(p: str) -> sqlite3.Connection:
    return sqlite3.connect(f"file:{Path(p).expanduser()}?mode=ro", uri=True)


def cur(c, cid, st):
    r = c.execute("SELECT outcome, criterion, value, tier_run_id, tags FROM outcomes WHERE candidate_id=? AND "
                  "stage=? ORDER BY id DESC LIMIT 1", (cid, st)).fetchone()
    return r


def prim_counts(st: dict) -> dict[str, int]:
    import spglib
    els = sorted(set(st["species"]))
    nums = [els.index(s) + 1 for s in st["species"]]
    try:
        cell = spglib.standardize_cell((st["lattice"], st["frac_coords"], nums), to_primitive=True,
                                       no_idealize=True, symprec=0.1)
        if cell is not None:
            return dict(collections.Counter(els[n - 1] for n in cell[2]))
    except Exception:  # noqa: BLE001
        pass
    return dict(collections.Counter(st["species"]))


def sg(st: dict, prec: float):
    from rtt.structure.canonical import spacegroup_number
    try:
        return spacegroup_number(st, prec)
    except Exception:  # noqa: BLE001
        return None


def reduced(comp: dict) -> str:
    from pymatgen.core import Composition
    return Composition(comp).reduced_formula


_REF: dict[str, dict] = {}


def refrows(cs: str) -> dict:
    """Alex-MP reference rows of a system from the front cache (written by S2); never queried here."""
    if cs in _REF:
        return _REF[cs]
    from rtt.source.external import path_of
    f = path_of("front_cache") / "alexmp_ref" / f"ref_{cs}.json"
    d = {"by_id": {}, "formulas": set()}
    if f.exists():
        from pymatgen.core import Composition
        for r in json.loads(f.read_text())["rows"]:
            d["by_id"][r["entry_id"]] = {"n_atoms": r.get("n_atoms"), "eh": r.get("e_above_hull_dft_eV")}
            try:
                d["formulas"].add(Composition(r["formula"]).reduced_formula)
            except Exception:  # noqa: BLE001
                pass
        d["ok"] = True
    else:
        d["ok"] = False
    _REF[cs] = d
    return d


def rows_of(tag: str, path: str, posthoc: dict | None = None) -> list[dict]:
    from rtt.score.metallicity import metallicity_prior
    posthoc = posthoc or {}
    c = ro(path)
    out = []
    for cid, lane, oj, sha in c.execute("SELECT id, lane, origin_json, structure_sha256 FROM candidates"):
        st = json.loads(c.execute("SELECT json FROM structures WHERE sha256=?", (sha,)).fetchone()[0])
        o = json.loads(oj)
        ck = o.get("checkpoint")
        key = f"{lane}:{ck}" if lane in ("G2", "G3", "G4") else lane
        sp = st["species"]
        els = set(sp)
        pc = prim_counts(st)
        mp = metallicity_prior(pc, cell=True)
        r = {"src": tag, "cid": cid, "lane": lane, "key": key, "cond": o.get("conditioning"), "n": len(sp),
             "n_prim": sum(pc.values()), "els": sorted(els), "nel": len(els), "Hf": sp.count("H") / len(sp),
             "formula": reduced(pc), "p_metal": mp["p_metal"], "prior_rule": mp["rule"],
             "n_val": mp.get("n_valence_e"), "flags": o.get("flags", [])}
        s1, s2, s3, s4 = (cur(c, cid, s) for s in ("S1", "S2", "S3", "S4"))
        r["S1"] = s1 and (s1[0] + (f"({s1[1]})" if s1[1] else ""))
        r["S2"] = s2 and s2[0]
        r["S3"] = s3 and (s3[0] + (f"({s3[1]})" if s3[1] else ""))
        r["S4"] = s4 and (s4[0] + (f"({s4[1]})" if s4[1] else ""))
        r["novel"] = r["match"] = None
        if s2 and s2[0] == "PASS" and s2[2]:
            try:
                v = json.loads(s2[2])
                r["novel"], r["match"] = v.get("novel"), v.get("alexmp_match")
            except Exception:  # noqa: BLE001
                pass
        cs = "-".join(sorted(els))
        rr = refrows(cs)
        r["ref_ok"] = rr["ok"]
        r["novel_comp"] = (r["formula"] not in rr["formulas"]) if rr["ok"] else None
        if r["match"] and rr["ok"]:
            e = rr["by_id"].get(r["match"])
            r["in_train_proxy"] = bool(e and (e["n_atoms"] or 99) <= 20 and e["eh"] is not None and e["eh"] <= 0.1)
        else:
            r["in_train_proxy"] = False if r["novel"] else None
        r["eh"] = r["eh_dft"] = r["score"] = r["sg_rel"] = None
        if s3 and s3[3]:
            m = json.loads(c.execute("SELECT metrics_json FROM tier_runs WHERE id=?", (s3[3],)).fetchone()[0] or "{}") or {}
            r["eh"], r["eh_dft"], r["score"] = m.get("e_hull_ml_self_eV"), m.get("e_hull_ml_dftref_eV"), m.get("score_s3")
            if m.get("relaxed_structure"):
                r["sg_rel"] = sg(m["relaxed_structure"], 0.1)
                r["relax_conv"] = (m.get("relax") or {}).get("converged")
        r["eh_src"] = "S3" if r["eh"] is not None else None
        if r["eh"] is None and (posthoc.get(cid) or {}).get("e_hull_ml_self_eV") is not None:
            r["eh"], r["eh_src"] = posthoc[cid]["e_hull_ml_self_eV"], "posthoc"
        r["sg_raw"] = sg(st, 0.1)
        out.append(r)
    return out


def med(x):
    x = [v for v in x if v is not None]
    return round(statistics.median(x), 3) if x else None


def summarize(R: list[dict]) -> dict:
    n = len(R)
    s2p = [r for r in R if r["S2"] == "PASS"]
    s3 = [r for r in s2p if r["eh"] is not None]
    hyd = [r for r in s2p if "H" in r["els"]]
    d = {
        "n": n,
        "S1_PASS": sum(1 for r in R if (r["S1"] or "").startswith("PASS")),
        "S1_no_pseudo": sum(1 for r in R if "no_pseudo" in (r["S1"] or "")),
        "S1_other_removed": sum(1 for r in R if r["S1"] and not r["S1"].startswith("PASS") and "no_pseudo" not in r["S1"]),
        "S2_PASS": len(s2p),
        "S2_INVALID": sum(1 for r in R if r["S2"] == "INVALID"),
        "S2_DUPLICATE": sum(1 for r in R if r["S2"] == "DUPLICATE"),
        "novel_struct": sum(1 for r in s2p if r["novel"] is True),
        "novel_comp": sum(1 for r in s2p if r["novel_comp"] is True),
        "in_train_proxy": sum(1 for r in s2p if r["in_train_proxy"] is True),
        "with_H": len(hyd),
        "H_f_median": med([r["Hf"] for r in s2p]),
        "H_f_ge_0.5": sum(1 for r in s2p if r["Hf"] >= 0.5),
        "n_elements": dict(sorted(collections.Counter(r["nel"] for r in s2p).items())),
        "n_atoms_median": med([r["n"] for r in s2p]),
        "palette_donor_TM_H": sum(1 for r in s2p if "H" in r["els"] and set(r["els"]) & PALETTE_DON and set(r["els"]) & PALETTE_TM),
        "light_element_cells": sum(1 for r in s2p if set(r["els"]) & LIGHT),
        "magnetic_3d_cells": sum(1 for r in s2p if set(r["els"]) & MAGNETIC_3D),
        "prior_p_ge_0.8": sum(1 for r in s2p if r["p_metal"] >= 0.8),
        "prior_odd_parity": sum(1 for r in s2p if r["prior_rule"] == "odd_valence_electron_count"),
        "prior_closed_shell_le_0.15": sum(1 for r in s2p if r["p_metal"] <= 0.15),
        "prior_neutral_0.5": sum(1 for r in s2p if r["p_metal"] == 0.5),
        "prior_rules": dict(collections.Counter(r["prior_rule"] for r in s2p).most_common(6)),
        "n_ehull": len(s3),
        "ehull_ml_median": med([r["eh"] for r in s3]),
        "ehull_le_0.1": sum(1 for r in s3 if r["eh"] <= 0.1),
        "ehull_le_0.35_H7ml": sum(1 for r in s3 if r["eh"] <= 0.35),
        "H1prior_and_H7ml": sum(1 for r in s3 if r["p_metal"] > 0.5 and r["eh"] <= 0.35),
        "H1prior_H7ml_withH": sum(1 for r in s3 if r["p_metal"] > 0.5 and r["eh"] <= 0.35 and "H" in r["els"]),
        "useful_H_metal_novel_le0.1": sum(1 for r in s3 if r["p_metal"] > 0.5 and r["eh"] <= 0.1 and "H" in r["els"]
                                          and r["novel"] is True),
        "S3": dict(collections.Counter(r["S3"] for r in s2p).most_common(8)),
        "top_elements": collections.Counter(e for r in s2p for e in r["els"]).most_common(12),
    }
    return d


def g2_sg(R):
    out = {}
    for r in R:
        if r["lane"] != "G2" or r["S2"] != "PASS":
            continue
        want = (r["cond"] or {}).get("space_group")
        d = out.setdefault(want, {"n": 0, "raw_match": 0, "relaxed_match": 0, "n_relaxed": 0, "with_H": 0})
        d["n"] += 1
        d["raw_match"] += r["sg_raw"] == want
        d["with_H"] += "H" in r["els"]
        if r["sg_rel"] is not None:
            d["n_relaxed"] += 1
            d["relaxed_match"] += r["sg_rel"] == want
    return dict(sorted(out.items()))


def gpu_h(timing: dict[str, str]) -> dict:
    out = collections.defaultdict(float)
    for tag, p in timing.items():
        p = Path(p).expanduser()
        if not p.exists():
            continue
        for line in p.read_text().splitlines():
            t = json.loads(line)
            if t.get("kind") == "mattergen" and t.get("wall_s") and not t.get("reused"):
                ck = t.get("checkpoint")
                lane = t["lane"]
                key = f"{lane}:{ck}" if lane in ("G2", "G3", "G4") else lane
                out[f"{tag}/{key}"] += t["wall_s"] / 3600
    return {k: round(v, 3) for k, v in out.items()}


def s4_queue(R: list[dict], k: int = 15) -> list[dict]:
    """Ranked S4 queue (ranking only, R2): S3 survivors (PASS or rank queue) of the G2/G3/G4b
    lanes with p_metal > 0.5 and E_hull <= 0.35 eV/atom (ML-self, else the DFT-reference hull).
    Key: (1) H in the cell (doc 02 target class; non-hydrides only after hydrides), (2) no ionic
    anion, (3) no magnetic 3d metal, (4) metallicity prior p_metal desc (odd parity 0.95 >
    open-shell counting 0.85 > ...), (5) E_hull asc, (6) H fraction desc, (7) primitive size."""
    pool = [r for r in R if r["lane"] in ("G2", "G3", "G4") and r["S3"] and
            (r["S3"].startswith("PASS") or r["S3"] == "DEFERRED(rank_cut)")]
    for r in pool:
        r["eh_any"] = r["eh"] if r["eh"] is not None else r["eh_dft"]
        els = set(r["els"])
        # electron-counting caveats on the parity rule: an odd count forces a metal only if the
        # normal state is non-magnetic. Ionic anions (O, F, Cl, Br, I, S, Se, Te) with an odd count
        # usually mean a local moment (e.g. MnF2), and 3d Cr-Ni carry an H2 risk.
        r["ionic_anion"] = bool(els & {"O", "F", "Cl", "Br", "I", "S", "Se", "Te"})
        r["mag3d"] = bool(els & MAGNETIC_3D)
    pool = [r for r in pool if r["p_metal"] > 0.5 and (r["eh_any"] is None or r["eh_any"] <= 0.35)]
    pool.sort(key=lambda r: ("H" not in r["els"], r["ionic_anion"], r["mag3d"], -r["p_metal"],
                             r["eh_any"] if r["eh_any"] is not None else 9, -r["Hf"], r["n_prim"], r["cid"]))
    return [{k2: r[k2] for k2 in ("cid", "key", "formula", "n_prim", "Hf", "p_metal", "prior_rule", "n_val",
                                  "eh", "eh_src", "eh_dft", "ionic_anion", "mag3d", "novel", "novel_comp", "sg_raw",
                                  "sg_rel", "S3", "S4")}
            for r in pool[:k]]


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--ledger", action="append", required=True, help="TAG=path")
    ap.add_argument("--timing", action="append", default=[], help="TAG=path")
    ap.add_argument("--queue-from", default="G")
    ap.add_argument("--posthoc", action="append", default=[], help="TAG=posthoc_hull.json")
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    out = Path(a.out)
    out.mkdir(parents=True, exist_ok=True)
    allrows = []
    ph = {t: json.loads(Path(p).expanduser().read_text()) for t, p in (x.split("=", 1) for x in a.posthoc)
          if Path(p).expanduser().exists()}
    for spec in a.ledger:
        tag, p = spec.split("=", 1)
        allrows += rows_of(tag, p, ph.get(tag))
    groups = collections.defaultdict(list)
    for r in allrows:
        groups[f"{r['src']}/{r['key']}"].append(r)
    summ = {k: summarize(v) for k, v in sorted(groups.items())}
    res = {"lanes": summ,
           "G2_space_group": g2_sg([r for r in allrows if r["src"] == a.queue_from]),
           "gpu_h_mattergen": gpu_h(dict(s.split("=", 1) for s in a.timing)),
           "s4_queue": s4_queue([r for r in allrows if r["src"] == a.queue_from])}
    (out / "lane_stats.json").write_text(json.dumps(res, indent=1, default=str))
    with open(out / "rows.jsonl", "w") as fh:
        for r in allrows:
            if r["src"] == a.queue_from:
                fh.write(json.dumps(r, default=str) + "\n")
    print(json.dumps(res, indent=1, default=str))
    return 0


if __name__ == "__main__":
    sys.exit(main())
