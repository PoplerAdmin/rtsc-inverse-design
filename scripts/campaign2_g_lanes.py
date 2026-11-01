"""Campaign 2, lanes G2 / G3 / G4-bulk (doc 34): MatterGen -> S1 -> S2 -> S3 on its own ledger.

    export OMP_NUM_THREADS=2 OPENBLAS_NUM_THREADS=2 MKL_NUM_THREADS=2 TORCH_DISABLE_NATIVE_JIT=1
    nice -n 19 ~/.rtt/envs/py/bin/python scripts/campaign2_g_lanes.py \
        --ledger ~/rtsc_targeted_work/G-lanes/ledger.sqlite --work ~/rtsc_targeted_work/G-lanes \
        [--phase plan|gen|s2|s3|rank|table|all] [--per-job 32] [--batch 16] [--rank-cut 40]

* Lane plan (architecture §3.3, ``rtt.source.lanes``; one LaneSpec per job, each checkpoint
  conditioned on exactly its own property, guidance 2):
  - G2 ``space_group``: one job per group of the doc 02 S6 list ``SPACE_GROUPS_S6`` (10 jobs);
  - G3 ``dft_band_gap = 0`` and ``dft_mag_density = 0`` (``--g3-jobs`` jobs each);
  - G4-bulk ``ml_bulk_modulus = ML_BULK_MODULUS_HIGH_GPA`` (300 GPa; ``--g4b-jobs`` jobs).
  ``--per-job`` structures per job in batches of ``--batch`` (GPU brief: batch <= 16).
  The day's shares (G2 0.15 : G3 0.10 : G4-bulk 0.0375) are not reproduced exactly: G4-bulk is
  given two jobs so that it has a sample at all.
* Mg2IrH6 control (lane C) is injected at S1 and runs every stage.
* S2: ``rtt.structure.stage.run_s2`` (validity, dedupe, Alex-MP novelty).
* S3: ``rtt.mlff.screen.run_s3`` unchanged, as lane P ran it (doc 32 §2): ``phonons=False``
  (MatterSim phonons never decide or rank a hydride; not computed), the entry-level cached
  MatterSim hull of ``campaign2_lane_p.make_cached_hull`` (same numbers as build_self_hull),
  one pool ``S3:campaign2-G`` closed once by ``rerank_s3``. Without ``--hull-nonH`` the ML-self
  hull is built only for H systems (a non-H system reuses a cached hull or gets NO_ML_HULL);
  measured 1-2 min per new system, ~400 mostly distinct systems would not fit the GPU budget.
* table: ``<work>/s3_table_G.jsonl``.
Resumable: MatterGen outputs are reused (same argv), ingest is idempotent, S2/S3 skip finished.
"""
from __future__ import annotations

import argparse
import collections
import hashlib
import json
import os
import random
import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
sys.path.insert(0, str(Path(__file__).resolve().parent))


def build_jobs(day: str, seed: int, per_job: int, batch: int, g3_jobs: int, g4b_jobs: int):
    from rtt.source import lanes as LN
    rng = random.Random(int(hashlib.sha256(f"{day}|{seed}".encode()).hexdigest()[:8], 16))
    jobs = []

    def add(spec, k):
        jobs.append(LN.LaneJob(spec=spec, n=per_job, seed=rng.randrange(1, 2**31 - 1),
                               batch_id=f"{day}-{spec.lane}-{k:03d}"))
    for k, sg in enumerate(LN.SPACE_GROUPS_S6):
        add(LN.LaneSpec("G2", "space_group", {"space_group": int(sg)}, batch_size=batch, source_id="doc02#5d2"), k)
    for k in range(2 * g3_jobs):
        if k % 2 == 0:
            add(LN.LaneSpec("G3", "dft_band_gap", {"dft_band_gap": 0.0}, batch_size=batch,
                            source_id="doc02#5d3:H1"), k)
        else:
            add(LN.LaneSpec("G3", "dft_mag_density", {"dft_mag_density": 0.0}, batch_size=batch,
                            source_id="doc02#5d3:H2"), k)
    for k in range(g4b_jobs):
        add(LN.LaneSpec("G4", "ml_bulk_modulus", {"ml_bulk_modulus": LN.ML_BULK_MODULUS_HIGH_GPA},
                        batch_size=batch, source_id="doc02#5c:S9"), 100 + k)
    return jobs


def table(L, work: Path) -> Path:
    from rtt.structure.canonical import spacegroup_number
    out = work / "s3_table_G.jsonl"
    with open(out, "w") as fh:
        for cid in L.candidate_ids():
            rec = L.candidate_record(cid)
            o = rec.origin
            r = {"cid": cid, "lane": o.lane, "checkpoint": o.checkpoint, "conditioning": o.conditioning,
                 "batch_id": o.batch_id, "formula": rec.structure.formula, "species": rec.structure.species,
                 "n_atoms": rec.structure.n_atoms, "cc": rec.cost_class, "sg_raw_0.01": None}
            try:
                r["sg_raw_0.01"] = spacegroup_number(rec.structure.to_json(), 0.01)
                r["sg_raw_0.1"] = spacegroup_number(rec.structure.to_json(), 0.1)
            except Exception:  # noqa: BLE001
                pass
            for st in ("S1", "S2", "S3", "S4"):
                oc = L.current_outcome(cid, st)
                r[st] = None if oc is None else oc["outcome"]
                if oc is None:
                    continue
                r[st + "_crit"] = oc.get("criterion")
                if st == "S2":
                    v = oc.get("value")
                    try:
                        v = json.loads(v) if isinstance(v, str) else v
                    except ValueError:
                        v = None
                    v = v if isinstance(v, dict) else {}
                    r["novel"] = (v or {}).get("novel")
                    r["alexmp_match"] = (v or {}).get("alexmp_match")
                    r["S2_tags"] = oc.get("tags")
            tr = L.q1("SELECT metrics_json FROM tier_runs WHERE candidate_id=? AND stage='S3' AND state='done' "
                      "ORDER BY id DESC LIMIT 1", (cid,))
            if tr and tr[0]:
                m = json.loads(tr[0])
                for k in ("e_hull_ml_self_eV", "e_hull_ml_dftref_eV", "score_s3", "metallicity_prior"):
                    r[k] = m.get(k)
                f = m.get("features") or {}
                if isinstance(f, dict) and "metallicity_prior" in f and r.get("metallicity_prior") is None:
                    r["metallicity_prior"] = f.get("metallicity_prior")
                rl = m.get("relax") or {}
                r["relax_converged"] = rl.get("converged")
                r["vol_per_atom"] = rl.get("volume_per_atom_A3")
                if m.get("relaxed_structure"):
                    try:
                        r["sg_relaxed_0.1"] = spacegroup_number(m["relaxed_structure"], 0.1)
                    except Exception:  # noqa: BLE001
                        pass
            fh.write(json.dumps(r, default=str) + "\n")
    return out


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--ledger", required=True)
    ap.add_argument("--work", required=True)
    ap.add_argument("--run-id", default="campaign2-G")
    ap.add_argument("--day", default="campaign2-G")
    ap.add_argument("--seed", type=int, default=20261030)
    ap.add_argument("--per-job", type=int, default=32)
    ap.add_argument("--batch", type=int, default=16)
    ap.add_argument("--g3-jobs", type=int, default=3)
    ap.add_argument("--g4b-jobs", type=int, default=2)
    ap.add_argument("--rank-cut", type=int, default=40)
    ap.add_argument("--chunk", type=int, default=48)
    ap.add_argument("--gpu-max-mib", type=int, default=8000)
    ap.add_argument("--phase", default="all", choices=["plan", "gen", "s2", "s3", "rank", "table", "all"])
    ap.add_argument("--max-n", type=int, default=None, help="S3: at most N candidates this call")
    ap.add_argument("--no-close", action="store_true")
    ap.add_argument("--hull-nonH", action="store_true", help="S3: also build ML-self hulls for non-H systems")
    a = ap.parse_args()
    work = Path(a.work).expanduser()
    work.mkdir(parents=True, exist_ok=True)
    os.environ["HF_HUB_OFFLINE"] = "1"          # never download a checkpoint implicitly
    logf = open(work / "G.log", "a")
    tlog = open(work / "timing_G.jsonl", "a")

    def log(msg):
        line = f"{time.strftime('%Y-%m-%dT%H:%M:%S')} {msg}"
        print(line, flush=True)
        logf.write(line + "\n")
        logf.flush()

    def timing(**kw):
        tlog.write(json.dumps({"t": time.strftime('%Y-%m-%dT%H:%M:%S'), **kw}, default=str) + "\n")
        tlog.flush()

    from campaign_front import gpu_used, have_checkpoint, wait_gpu
    from rtt.config import load_machine, load_targets, load_tiers
    from rtt.ledger.store import Ledger
    from rtt.source import lanes as LN
    t, ti, m = load_targets(), load_tiers(), load_machine()
    L = Ledger(Path(a.ledger).expanduser(), run_id=a.run_id)
    L.register_configs(t, ti, m)

    jobs = build_jobs(a.day, a.seed, a.per_job, a.batch, a.g3_jobs, a.g4b_jobs)
    plan = [{"batch_id": j.batch_id, "lane": j.spec.lane, "checkpoint": j.spec.checkpoint,
             "conditioning": dict(j.spec.conditioning), "n": j.n, "seed": j.seed,
             "status": "planned" if have_checkpoint(j.spec.checkpoint) else "SKIPPED_CHECKPOINT_NOT_LOCAL"}
            for j in jobs]
    (work / "lane_plan_G.json").write_text(json.dumps({"plan_summary": LN.plan_summary(jobs), "jobs": plan,
                                                        "per_job": a.per_job, "batch": a.batch}, indent=1))
    log(f"plan: {collections.Counter((p['lane'], p['checkpoint'], p['status']) for p in plan)}")
    if a.phase == "plan":
        return 0

    if a.phase in ("gen", "all"):
        from rtt.source.ingest import ingest_many
        from rtt.source.mattergen import MatterGenRunError, generate
        for j, p in zip(jobs, plan):
            if p["status"] != "planned":
                continue
            waited = wait_gpu(a.gpu_max_mib, log)
            t0 = time.time()
            try:
                raws = generate(j, workdir=work / "S1", machine=m)
            except MatterGenRunError as ex:
                log(f"{j.batch_id} FAILED: {str(ex)[-300:]}")
                timing(kind="mattergen", batch_id=j.batch_id, lane=j.spec.lane, error=str(ex)[-300:])
                continue
            info = raws[0].meta.get("run", {}) if raws else {}
            cids = ingest_many(L, raws, targets=t, tiers=ti)
            log(f"{j.batch_id} {j.spec.checkpoint} {j.spec.conditioning} -> {len(raws)} structures, "
                f"wall {info.get('wall_s')} s{' (reused)' if info.get('reused') else ''}")
            timing(kind="mattergen", batch_id=j.batch_id, lane=j.spec.lane, checkpoint=j.spec.checkpoint,
                   checkpoint_sha256=raws[0].origin.checkpoint_sha256 if raws else None,
                   conditioning=dict(j.spec.conditioning), n=len(raws), wall_s=info.get("wall_s"),
                   reused=bool(info.get("reused")), gpu_wait_s=round(waited, 1), n_cids=len(set(cids)),
                   loop_s=round(time.time() - t0, 1), gpu_used_mib_after=gpu_used())
        from rtt.source.controls import inject_controls
        cc = inject_controls(L, ("Mg2IrH6",), targets=t, tiers=ti)
        s1 = collections.Counter((L.q1("SELECT lane FROM candidates WHERE id=?", (c,))[0],
                                  (L.current_outcome(c, "S1") or {}).get("outcome")) for c in L.candidate_ids())
        log(f"S1: {dict(s1)}; control {cc}")

    if a.phase in ("s2", "all"):
        from rtt.structure.stage import run_s2
        t0 = time.time()
        s2 = run_s2(L, targets=t, tiers=ti, novelty=True, machine=m)
        log(f"S2: {collections.Counter(o.outcome for o in s2.values())} in {time.time() - t0:.0f} s")
        timing(kind="S2", n=len(s2), wall_s=round(time.time() - t0, 1))

    if a.phase in ("s3", "all"):
        from campaign2_lane_p import make_cached_hull
        from rtt.engine.retry import needs_run
        from rtt.mlff.screen import run_s3
        try:   # GPU brief: < 16 GB; 0.6 x 24 GB = 14.4 GB
            import torch
            if torch.cuda.is_available():
                torch.cuda.set_per_process_memory_fraction(float(os.environ.get("RTT_GPU_FRACTION", "0.6")))
        except Exception as ex:  # noqa: BLE001
            log(f"GPU memory cap not set: {ex}")
        cached = make_cached_hull(work, log)
        if a.hull_nonH:
            hull_fn = cached
        else:
            from rtt.source.external import path_of

            def hull_fn(chemsys, *, pressure_GPa=0.0, model="MatterSim-v1.0.0-5M", machine=None, **kw):
                """ML-self hull only for H systems (GPU budget, doc 34 §2). A non-H system uses an
                already-cached hull if one exists; otherwise NO_ML_HULL (run_s3 records the
                exception). This cannot change a decision: for a non-H cell an ML-window
                exceedance is DEFERRED(uncalibrated_window) anyway (doc 31 §7); it changes the
                rank only (neutral hull term). The DFT-reference hull is still computed."""
                cs = "-".join(sorted(set(chemsys.split("-"))))
                if "H" not in cs.split("-"):
                    sub = f"{model}_P{pressure_GPa:g}/{cs}.json"
                    if not any((d / sub).exists() for d in (work / "mlhull", path_of("front_cache", machine) / "mlhull")):
                        raise RuntimeError("ML-self hull not built for a non-H system (GPU budget; doc 34 §2)")
                return cached(chemsys, pressure_GPa=pressure_GPa, model=model, machine=machine, **kw)
        def queued(cid):    # rank-queue entries are re-ranked by rerank_s3, never recomputed
            o = L.current_outcome(cid, "S3") or {}
            return o.get("outcome") == "DEFERRED" and "rank_cut" in (o.get("criterion") or "")
        todo = [cid for cid in L.candidate_ids()
                if (L.current_outcome(cid, "S2") or {}).get("outcome") == "PASS" and needs_run(L, cid, "S3", ti)
                and not queued(cid)]
        lane = {c: L.q1("SELECT lane FROM candidates WHERE id=?", (c,))[0] for c in todo}
        # controls first; then H-containing cells (the target class), then by size
        hasH = {c: "H" in L.candidate_record(c).structure.species for c in todo}
        todo.sort(key=lambda c: (lane[c] != "C", not hasH[c], L.candidate_record(c).structure.n_atoms, c))
        if a.max_n:
            todo = todo[:a.max_n]
        log(f"S3 to compute: {len(todo)} ({sum(hasH[c] for c in todo)} with H)")
        for i in range(0, len(todo), a.chunk):
            chunk = todo[i:i + a.chunk]
            t0 = time.time()
            res = run_s3(L, chunk, targets=t, tiers=ti, machine=m, device="cuda", workdir=work / "S3",
                         phonons=False, alignn=False, close_pool=False, pool_id="S3:campaign2-G",
                         hull_fn=hull_fn)
            dt = time.time() - t0
            c = collections.Counter(o.outcome + (f"({o.criterion})" if o.criterion else "") for o in res.values())
            log(f"S3 chunk {i // a.chunk}: {len(chunk)} in {dt:.0f} s: {dict(c)}")
            timing(kind="S3_chunk", n=len(chunk), wall_s=round(dt, 1),
                   lanes=dict(collections.Counter(lane[x] for x in chunk)))
    if a.phase in ("s3", "all", "rank") and not a.no_close:
        from rtt.mlff.screen import rerank_s3
        fin = rerank_s3(L, targets=t, tiers=ti, rank_cut=a.rank_cut, pool_id="S3:campaign2-G")
        log(f"S3 rank cut {a.rank_cut}: promoted {sum(1 for o in fin.values() if o.outcome == 'PASS')}")

    if a.phase in ("table", "all", "s3", "rank"):
        log(f"table: {table(L, work)}")
    log(f"counts: {json.dumps(L.summary_counts())[:1500]}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
