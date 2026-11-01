"""Campaign front end (S1 lanes G/D/C -> S2 -> S3) on one run ledger (agent E1-campaign).

    export OMP_NUM_THREADS=2 OPENBLAS_NUM_THREADS=2 MKL_NUM_THREADS=2 TORCH_DISABLE_NATIVE_JIT=1
    nice -n 19 ~/.rtt/envs/py/bin/python scripts/campaign_front.py --ledger L.sqlite --work DIR \
        [--gpu-h 2.0] [--batch 16] [--rank-cut 20] [--phase gen|s2|s3|all]

* Lane plan: ``rtt.source.lanes.plan_day`` at reduced GPU hours. G1 draws from the doc 02 S10
  family-1 enumeration (electropositive donor x late 4d/5d metal x H) because
  ``targets.toml [chemsys] list`` points at a file that does not exist (the fallback would also
  include donor-donor-H systems). Jobs whose checkpoint is not in the local HF cache are
  **skipped and recorded** (no downloads).
* Shared GPU: before every GPU job, wait until ``nvidia-smi`` shows < ``--gpu-max-mib`` used.
* S3 day-wide rank cut: chunks run with rank_cut 0 (every non-control survivor enters the rank
  queue as DEFERRED(rank_cut)), then one ``rerank_s3`` with the real cut promotes the top N of
  the whole day. (Calling run_s3 per chunk with the real cut would let the first chunk take all
  slots.)
* Timing per job/lane goes to ``<work>/timing_front.jsonl``.
Resumable: MatterGen outputs are reused, ingest is idempotent, S2/S3 skip finished candidates.
"""
from __future__ import annotations

import argparse
import collections
import json
import os
import subprocess
import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

DONORS = ("Li", "Na", "Mg", "Ca", "Sr", "Ba", "Al", "Sc", "Y", "La")
TMS = ("Rh", "Pd", "Ir", "Pt", "Os", "Au", "Ag")


def gpu_used() -> int | None:
    try:
        cp = subprocess.run(["nvidia-smi", "--query-gpu=memory.used", "--format=csv,noheader,nounits"],
                            capture_output=True, text=True, timeout=30)
        return int(cp.stdout.strip().splitlines()[0])
    except Exception:  # noqa: BLE001
        return None


def wait_gpu(max_mib: int, log) -> float:
    t0 = time.time()
    while True:
        u = gpu_used()
        if u is not None and u < max_mib:
            return time.time() - t0
        log(f"GPU busy ({u} MiB >= {max_mib}); waiting")
        time.sleep(60)


def have_checkpoint(name: str) -> bool:
    hub = Path(os.path.expanduser(os.environ.get("RTT_HF_HUB", "~/.cache/huggingface/hub")))
    for p in hub.glob(f"models--microsoft--mattergen/snapshots/*/checkpoints/{name}/checkpoints/last.ckpt"):
        try:
            if p.resolve().stat().st_size > 1_000_000:
                return True
        except OSError:
            pass
    return False


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--ledger", required=True)
    ap.add_argument("--work", required=True)
    ap.add_argument("--run-id", default="campaign1")
    ap.add_argument("--gpu-h", type=float, default=2.0)
    ap.add_argument("--batch", type=int, default=16)
    ap.add_argument("--day", default="campaign1-day0")
    ap.add_argument("--seed", type=int, default=20261029)
    ap.add_argument("--rank-cut", type=int, default=20)
    ap.add_argument("--chunk", type=int, default=48)
    ap.add_argument("--gpu-max-mib", type=int, default=12000)
    ap.add_argument("--phase", default="all", choices=["plan", "gen", "s2", "s3", "all"])
    ap.add_argument("--d-limit", type=int, default=None)
    a = ap.parse_args()
    work = Path(a.work).expanduser()
    work.mkdir(parents=True, exist_ok=True)
    os.environ["HF_HUB_OFFLINE"] = "1"          # never download a checkpoint implicitly
    logf = open(work / "front.log", "a")
    import faulthandler
    faulthandler.dump_traceback_later(1200, repeat=True, file=open(work / "front_stacks.log", "a"))   # hang diagnosis
    tlog = open(work / "timing_front.jsonl", "a")

    def log(msg):
        line = f"{time.strftime('%Y-%m-%dT%H:%M:%S')} {msg}"
        print(line, flush=True)
        logf.write(line + "\n")
        logf.flush()

    def timing(**kw):
        tlog.write(json.dumps({"t": time.strftime('%Y-%m-%dT%H:%M:%S'), **kw}, default=str) + "\n")
        tlog.flush()

    from rtt.config import load_machine, load_targets, load_tiers
    from rtt.ledger.store import Ledger
    from rtt.source import lanes as LN
    t, ti, m = load_targets(), load_tiers(), load_machine()
    L = Ledger(Path(a.ledger).expanduser(), run_id=a.run_id)
    L.register_configs(t, ti, m)

    fam1 = sorted(dict.fromkeys(LN.norm_chemsys([d, tm, "H"]) for d in DONORS for tm in TMS))
    LN.load_chemsys_list = lambda targets=None: (fam1, ["CAMPAIGN_FAMILY1_DONOR_TM"])   # see docstring
    jobs = LN.plan_day(m, t, tiers=ti, day=a.day, seed=a.seed, gpu_h=a.gpu_h, batch_size=a.batch)
    plan = []
    for j in jobs:
        ok = have_checkpoint(j.spec.checkpoint)
        plan.append({"batch_id": j.batch_id, "lane": j.spec.lane, "checkpoint": j.spec.checkpoint,
                     "conditioning": dict(j.spec.conditioning), "n": j.n, "seed": j.seed,
                     "status": "planned" if ok else "SKIPPED_CHECKPOINT_NOT_LOCAL"})
    (work / "lane_plan.json").write_text(json.dumps({"plan_summary": LN.plan_summary(jobs), "jobs": plan,
                                                      "gpu_h": a.gpu_h, "batch": a.batch}, indent=1))
    log(f"plan: {collections.Counter((p['lane'], p['status']) for p in plan)}")
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
            log(f"{j.batch_id} {j.spec.conditioning} -> {len(raws)} structures, wall {info.get('wall_s')} s"
                f"{' (reused)' if info.get('reused') else ''}")
            timing(kind="mattergen", batch_id=j.batch_id, lane=j.spec.lane, checkpoint=j.spec.checkpoint,
                   conditioning=dict(j.spec.conditioning), n=len(raws), wall_s=info.get("wall_s"),
                   reused=bool(info.get("reused")), gpu_wait_s=round(waited, 1), n_cids=len(set(cids)),
                   loop_s=round(time.time() - t0, 1), gpu_used_mib_at_start=gpu_used())
        # lane D: the G1 systems actually planned + the doc 02 5(c) families
        from rtt.source.runner import run_database_lane
        systems = sorted({p["conditioning"]["chemical_system"] for p in plan
                          if p["lane"] in ("G1", "G4") and "chemical_system" in p["conditioning"]})
        t0 = time.time()
        dc = run_database_lane(L, systems, targets=t, tiers=ti, machine=m, limit=a.d_limit)
        log(f"lane D: {len(dc)} entries from {len(systems)} systems in {time.time() - t0:.0f} s")
        timing(kind="lane_D", systems=systems, n=len(dc), wall_s=round(time.time() - t0, 1))
        from rtt.source.controls import inject_controls
        cc = inject_controls(L, targets=t, tiers=ti)
        log(f"controls: {len(cc)}")

    if a.phase in ("s2", "all"):
        from rtt.structure.stage import run_s2
        t0 = time.time()
        s2 = run_s2(L, targets=t, tiers=ti, novelty=True, machine=m)
        log(f"S2: {collections.Counter(o.outcome for o in s2.values())} in {time.time() - t0:.0f} s")
        timing(kind="S2", n=len(s2), wall_s=round(time.time() - t0, 1))

    if a.phase in ("s3", "all"):
        from rtt.engine.retry import needs_run
        from rtt.mlff.screen import rerank_s3, run_s3
        todo = []
        for cid in L.candidate_ids():
            o2 = L.current_outcome(cid, "S2")
            if o2 and o2["outcome"] == "PASS" and needs_run(L, cid, "S3", ti):
                todo.append(cid)
        # controls first, then by size so a batch waits for similar cells
        lane = {c: L.q1("SELECT lane FROM candidates WHERE id=?", (c,))[0] for c in todo}
        todo.sort(key=lambda c: (lane[c] != "C", L.candidate_record(c).structure.n_atoms))
        log(f"S3 to compute: {len(todo)}")
        for i in range(0, len(todo), a.chunk):
            chunk = todo[i:i + a.chunk]
            waited = wait_gpu(a.gpu_max_mib, log)
            t0 = time.time()
            res = run_s3(L, chunk, targets=t, tiers=ti, machine=m, device="cuda", workdir=work / "S3",
                         alignn=False, rank_cut=0)
            dt = time.time() - t0
            c = collections.Counter(o.outcome + (f"({o.criterion})" if o.criterion else "") for o in res.values())
            log(f"S3 chunk {i // a.chunk}: {len(chunk)} in {dt:.0f} s: {dict(c)}")
            lanes_c = collections.Counter(lane[x] for x in chunk)
            timing(kind="S3_chunk", n=len(chunk), wall_s=round(dt, 1), gpu_wait_s=round(waited, 1),
                   lanes=dict(lanes_c))
        fin = rerank_s3(L, targets=t, tiers=ti, rank_cut=a.rank_cut)
        log(f"S3 day-wide rank cut {a.rank_cut}: promoted {sum(1 for o in fin.values() if o.outcome == 'PASS')}")
    log(f"counts: {json.dumps(L.summary_counts())[:2000]}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
