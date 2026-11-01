"""Campaign 2, lane P (doc 32): prototype substitution -> S1 -> S2 -> S3 on its own ledger.

    export OMP_NUM_THREADS=2 OPENBLAS_NUM_THREADS=2 MKL_NUM_THREADS=2 TORCH_DISABLE_NATIVE_JIT=1
    nice -n 19 ~/.rtt/envs/py/bin/python scripts/campaign2_lane_p.py \
        --ledger ~/rtsc_targeted_work/P2-lane/ledger.sqlite --work ~/rtsc_targeted_work/P2-lane \
        [--phase s1|s2|s3|table|all] [--rank-cut 40] [--chunk 48]

* S1: ``rtt.source.prototype.plan_lane_p`` (configs/lane_p.toml) + the Mg2IrH6 control (lane C).
* S2: the standard validity / dedupe / Alex-MP novelty path (``rtt.structure.stage.run_s2``).
* S3: ``rtt.mlff.screen.run_s3`` unchanged, with two run options:
  - ``phonons=False``: MatterSim phonons are never used for H compounds (rank_metrics drops them,
    decide_s3 never reads them; doc 15 shows they fake the Mg2IrH6 T1g instability), so they
    are not computed (GPU time);
  - ``hull_fn``: the MatterSim-only hull of :mod:`rtt.mlff.hull` with an **entry-level** cache
    (Alex-MP entry id -> MatterSim enthalpy, same model, same pressure, same relax settings), so
    a quaternary A-A'-M-H hull relaxes only the reference entries no earlier hull relaxed.
    Same numbers as ``build_self_hull`` (each reference is relaxed independently); written to
    ``<work>/mlhull`` (the shared cache is read, never written).
* table: ``<work>/s3_table_laneP.jsonl`` (one row per candidate: provenance, S2/S3, metrics).
"""
from __future__ import annotations

import argparse
import collections
import json
import os
import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))


def make_cached_hull(work: Path, log):
    from rtt.mlff import hull as H
    from rtt.source.external import path_of, refquery
    seeds: dict[tuple, dict] = {}

    def seed_from(d: Path, key):
        for f in d.glob("*.json"):
            try:
                for e in json.loads(f.read_text())["entries"]:
                    seeds.setdefault(key, {})[e["entry_id"]] = e
            except Exception:  # noqa: BLE001
                pass

    def cached_hull(chemsys, *, pressure_GPa=0.0, model="MatterSim-v1.0.0-5M", device=None, machine=None,
                    relax_fn=None, max_ref_atoms=64):
        cs = "-".join(sorted(set(chemsys.split("-"))))
        key = (model, float(pressure_GPa))
        mine = work / "mlhull" / f"{model}_P{pressure_GPa:g}"
        mine.mkdir(parents=True, exist_ok=True)
        shared = path_of("front_cache", machine) / "mlhull" / f"{model}_P{pressure_GPa:g}"
        for d in (mine, shared):
            if (d / f"{cs}.json").exists():
                return H.SelfHull(**json.loads((d / f"{cs}.json").read_text()))
        if key not in seeds:
            seeds[key] = {}
            for d in (shared, mine):
                if d.is_dir():
                    seed_from(d, key)
            log(f"hull entry cache seeded: {len(seeds[key])} entries")
        cache = seeds[key]
        data = refquery(cs, machine=machine)
        rows = [r for r in data["rows"] if r["n_atoms"] <= max_ref_atoms]
        todo = [r for r in rows if r["entry_id"] not in cache]
        t0 = time.time()
        if todo:
            if relax_fn is None:
                from rtt.mlff.relax import relax_batch as relax_fn  # noqa: N813
            res = relax_fn([r["structure"] for r in todo], model=model, pressure_GPa=pressure_GPa, device=device)
            for r, x in zip(todo, res):
                if x.error or x.enthalpy_per_atom_eV != x.enthalpy_per_atom_eV:
                    cache[r["entry_id"]] = None
                    continue
                cache[r["entry_id"]] = {"entry_id": r["entry_id"], "formula": r["formula"],
                                        "enthalpy_per_atom_eV": x.enthalpy_per_atom_eV,
                                        "converged": bool(x.converged)}
        entries = [cache[r["entry_id"]] for r in rows if cache.get(r["entry_id"])]
        covered = sorted({H._element_of(e["formula"]) for e in entries
                          if e["converged"] and H._is_element(e["formula"])})
        h = H.SelfHull(chemsys=cs, pressure_GPa=pressure_GPa, model=model, entries=entries,
                       n_ref_total=len(data["rows"]), n_ref_skipped_size=len(data["rows"]) - len(rows),
                       n_ref_unconverged=sum(1 for e in entries if not e["converged"]), elements_covered=covered)
        f = mine / f"{cs}.json"
        f.with_suffix(".tmp").write_text(json.dumps({k: v for k, v in h.__dict__.items() if k != "_pd"}))
        os.replace(f.with_suffix(".tmp"), f)
        log(f"hull {cs}: {len(rows)} refs, {len(todo)} newly relaxed in {time.time() - t0:.0f} s")
        return h

    return cached_hull


def table(L, work: Path) -> Path:
    out = work / "s3_table_laneP.jsonl"
    rows = []
    for cid in L.candidate_ids():
        rec = L.candidate_record(cid)
        cond = rec.origin.conditioning or {}
        r = {"cid": cid, "lane": rec.origin.lane, "formula": cond.get("formula") or rec.structure.formula,
             "prototype": cond.get("prototype"), "substitution": cond.get("substitution"),
             "N": cond.get("N_per_complex"), "parent_N": cond.get("parent_N"),
             "parity": cond.get("valence_parity"), "n_atoms": rec.structure.n_atoms, "cc": rec.cost_class}
        for st in ("S1", "S2", "S3", "S4"):
            o = L.current_outcome(cid, st)
            r[st] = None if o is None else o["outcome"]
            if o is not None and st == "S2":
                v = o.get("value")
                v = json.loads(v) if isinstance(v, str) else (v or {})
                r["novel"] = (v or {}).get("novel")
                r["alexmp_match"] = (v or {}).get("alexmp_match")
            if o is not None and st == "S3":
                r["S3_crit"] = o.get("criterion")
                r["S3_reason"] = (o.get("reason") or "")[:120]
        tr = L.q1("SELECT metrics_json FROM tier_runs WHERE candidate_id=? AND stage='S3' AND state='done' "
                  "ORDER BY id DESC LIMIT 1", (cid,))
        if tr and tr[0]:
            m = json.loads(tr[0])
            for k in ("e_hull_ml_self_eV", "e_hull_ml_dftref_eV", "score_s3"):
                r[k] = m.get(k)
            rl = m.get("relax") or {}
            r["relax_converged"] = rl.get("converged")
            if m.get("relaxed_structure"):
                from rtt.structure.canonical import spacegroup_number
                r["sg_relaxed_0.01"] = spacegroup_number(m["relaxed_structure"], 0.01)
                r["sg_relaxed_0.1"] = spacegroup_number(m["relaxed_structure"], 0.1)
                r["relaxed_match_of"] = m.get("relaxed_match_of")
            r["vol_per_atom"] = rl.get("volume_per_atom_A3")
        rows.append(r)
    with open(out, "w") as fh:
        for r in rows:
            fh.write(json.dumps(r, default=str) + "\n")
    return out


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--ledger", required=True)
    ap.add_argument("--work", required=True)
    ap.add_argument("--run-id", default="campaign2-laneP")
    ap.add_argument("--phase", default="all", choices=["s1", "s2", "s3", "rank", "table", "all"])
    ap.add_argument("--rank-cut", type=int, default=40)
    ap.add_argument("--chunk", type=int, default=48)
    ap.add_argument("--reverse", action="store_true", help="S3: compute from the end of the list (2nd worker)")
    ap.add_argument("--max-n", type=int, default=None, help="S3: at most N candidates")
    ap.add_argument("--no-close", action="store_true", help="S3: do not close the rank pool")
    ap.add_argument("--limit", type=int, default=None, help="first N planned candidates only (smoke)")
    a = ap.parse_args()
    work = Path(a.work).expanduser()
    work.mkdir(parents=True, exist_ok=True)
    os.environ["HF_HUB_OFFLINE"] = "1"
    logf = open(work / "laneP.log", "a")

    def log(msg):
        line = f"{time.strftime('%Y-%m-%dT%H:%M:%S')} {msg}"
        print(line, flush=True)
        logf.write(line + "\n")
        logf.flush()

    from rtt.config import load_machine, load_targets, load_tiers
    from rtt.ledger.store import Ledger
    from rtt.source import prototype as PR
    t, ti, m = load_targets(), load_tiers(), load_machine()
    L = Ledger(Path(a.ledger).expanduser(), run_id=a.run_id)
    L.register_configs(t, ti, m)
    cfg = PR.load_lane_p_config()

    if a.phase in ("s1", "all"):
        plan = PR.plan_lane_p(cfg)
        if a.limit:
            plan = plan[:a.limit]
        (work / "laneP_plan.json").write_text(json.dumps(
            {"config_sha256": cfg.sha256, "n": len(plan),
             "by_prototype": collections.Counter(c.prototype for c in plan),
             "by_N": collections.Counter(str(float(c.N)) for c in plan),
             "candidates": [c.provenance() for c in plan]}, indent=1, default=str))
        t0 = time.time()
        cids = PR.run_lane_p(L, targets=t, tiers=ti, cfg=cfg, cands=plan)
        from rtt.source.controls import inject_controls
        cc = inject_controls(L, ("Mg2IrH6",), targets=t, tiers=ti)
        s1 = collections.Counter(L.current_outcome(c, "S1")["outcome"] for c in cids)
        log(f"S1: {len(cids)} lane-P cells ({len(set(cids))} unique ids) {dict(s1)}; control {cc}; "
            f"{time.time() - t0:.0f} s")

    if a.phase in ("s2", "all"):
        from rtt.structure.stage import run_s2
        t0 = time.time()
        s2 = run_s2(L, targets=t, tiers=ti, novelty=True, machine=m)
        log(f"S2: {collections.Counter(o.outcome for o in s2.values())} in {time.time() - t0:.0f} s")

    if a.phase in ("s3", "all"):
        from rtt.engine.retry import needs_run
        from rtt.mlff.screen import rerank_s3, run_s3
        try:   # GPU etiquette (brief: < 16 GB): cap the caching allocator; 0.6 x 24 GB = 14.4 GB
            import torch
            if torch.cuda.is_available():
                torch.cuda.set_per_process_memory_fraction(float(os.environ.get("RTT_GPU_FRACTION", "0.6")))
        except Exception as ex:  # noqa: BLE001
            log(f"GPU memory cap not set: {ex}")
        hull_fn = make_cached_hull(work, log)
        todo = [cid for cid in L.candidate_ids()
                if (L.current_outcome(cid, "S2") or {}).get("outcome") == "PASS" and needs_run(L, cid, "S3", ti)]
        lane = {c: L.q1("SELECT lane FROM candidates WHERE id=?", (c,))[0] for c in todo}
        todo.sort(key=lambda c: (lane[c] != "C", L.candidate_record(c).structure.n_atoms, c))
        if a.reverse:
            todo = todo[::-1]
        if a.max_n:
            todo = todo[:a.max_n]
        log(f"S3 to compute: {len(todo)} (reverse={a.reverse})")
        for i in range(0, len(todo), a.chunk):
            chunk = todo[i:i + a.chunk]
            t0 = time.time()
            res = run_s3(L, chunk, targets=t, tiers=ti, machine=m, device="cuda", workdir=work / "S3",
                         phonons=False, alignn=False, close_pool=False, pool_id="S3:campaign2-laneP",
                         hull_fn=hull_fn)
            c = collections.Counter(o.outcome + (f"({o.criterion})" if o.criterion else "") for o in res.values())
            log(f"S3 chunk {i // a.chunk}: {len(chunk)} in {time.time() - t0:.0f} s: {dict(c)}")
    if a.phase in ("s3", "all", "rank") and not a.no_close:
        from rtt.mlff.screen import rerank_s3
        fin = rerank_s3(L, targets=t, tiers=ti, rank_cut=a.rank_cut, pool_id="S3:campaign2-laneP")
        log(f"S3 rank cut {a.rank_cut}: promoted {sum(1 for o in fin.values() if o.outcome == 'PASS')}")

    if a.phase in ("table", "all", "s3", "rank"):
        log(f"table: {table(L, work)}")
    log(f"counts: {json.dumps(L.summary_counts())[:1500]}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
