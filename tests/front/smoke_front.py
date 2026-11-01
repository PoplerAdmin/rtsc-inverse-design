"""Live front-end smoke run (not collected by pytest; needs the GPU, MatterGen venv, LMDB).

One lane job (16 structures) -> S1 ingest -> S2 validate/dedupe/novelty -> S3 MatterSim
relax + hulls + features + phonons -> ledger rows; plus the control lane through S1-S3.

    export OMP_NUM_THREADS=2 OPENBLAS_NUM_THREADS=2 MKL_NUM_THREADS=2 TORCH_DISABLE_NATIVE_JIT=1
    nice -n 19 ~/.rtt/envs/py/bin/python tests/front/smoke_front.py OUTDIR [--chemsys H-Ir-Mg] [--n 16]

Refuses to start if the GPU already has >= 4000 MiB in use (shared-GPU rule).
"""
from __future__ import annotations

import argparse
import collections
import json
import os
import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "src"))


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("outdir")
    ap.add_argument("--chemsys", default="H-Ir-Mg")
    ap.add_argument("--ehull-cond", type=float, default=0.1)
    ap.add_argument("--n", type=int, default=16)
    ap.add_argument("--seed", type=int, default=20261029)
    ap.add_argument("--no-controls", action="store_true")
    ap.add_argument("--mattergen-dir", default=None, help="reuse finished MatterGen runs from here")
    a = ap.parse_args()
    out = Path(a.outdir).expanduser()
    out.mkdir(parents=True, exist_ok=True)
    os.environ.setdefault("RTT_FRONT_CACHE", str(out / "front_cache"))

    from rtt.config import load_machine, load_targets, load_tiers
    from rtt.ledger.store import Ledger
    from rtt.mlff.screen import run_s3
    from rtt.source.controls import inject_controls
    from rtt.source.external import gpu_used_mib
    from rtt.source.ingest import ingest_many
    from rtt.source.lanes import LaneJob, LaneSpec
    from rtt.source.mattergen import generate
    from rtt.structure.stage import run_s2

    used = gpu_used_mib()
    if used is None or used >= 4000:
        print(f"GPU busy or unavailable ({used} MiB used); not starting", flush=True)
        return 2
    t, ti, m = load_targets(), load_tiers(), load_machine()
    L = Ledger(out / "ledger.sqlite", run_id=f"front-smoke-{a.seed}")
    L.register_configs(t, ti, m)
    timing = {}
    t0 = time.time()
    spec = LaneSpec("G1", "chemical_system_energy_above_hull",
                    {"chemical_system": a.chemsys, "energy_above_hull": a.ehull_cond}, batch_size=16,
                    source_id="doc02#5d1")
    job = LaneJob(spec=spec, n=a.n, seed=a.seed, batch_id=f"smoke-G1-{a.chemsys}")
    raws = generate(job, workdir=Path(a.mattergen_dir).expanduser() if a.mattergen_dir else out / "mattergen_runs",
                    machine=m)
    timing["S1_mattergen_s"] = round(time.time() - t0, 1)
    t0 = time.time()
    print(f"generated {len(raws)} in {timing['S1_mattergen_s']} s", flush=True)
    cids = ingest_many(L, raws, targets=t, tiers=ti)
    ctrl = [] if a.no_controls else inject_controls(L, targets=t, tiers=ti)
    timing["S1_ingest_s"] = round(time.time() - t0, 1)  # incl. controls
    t1 = time.time()
    s2 = run_s2(L, targets=t, tiers=ti, novelty=True, machine=m)
    timing["S2_s"] = round(time.time() - t1, 1)
    print("S2", collections.Counter(o.outcome for o in s2.values()), flush=True)
    used = gpu_used_mib()
    if used is None or used >= 4000:
        print(f"GPU became busy ({used} MiB); stopping before S3", flush=True)
        return 3
    t2 = time.time()
    s3 = run_s3(L, targets=t, tiers=ti, machine=m, device="cuda", workdir=out, alignn=True)
    timing["S3_s"] = round(time.time() - t2, 1)
    print("S3", collections.Counter(o.outcome for o in s3.values()), flush=True)

    rows = []
    for cid in cids + ctrl:
        rec = L.candidate_record(cid)
        row = {"cid": cid, "lane": rec.origin.lane, "control": rec.origin.control, "formula": rec.formula,
               "sg": rec.structure.spacegroup, "n": rec.structure.n_atoms, "P": rec.pressure_GPa,
               "cost_class": rec.cost_class}
        for st in ("S1", "S2", "S3"):
            o = L.current_outcome(cid, st)
            row[st] = None if o is None else o["outcome"] + (f"({o['criterion']})" if o.get("criterion") else "")
            if st == "S2" and o is not None and isinstance(o.get("value"), dict):
                row["novel"] = o["value"].get("novel")
            if st == "S3" and o is not None and o.get("tier_run_id"):
                mt = json.loads(L.tier_run(o["tier_run_id"])["metrics_json"] or "{}")
                ph = mt.get("phonons") or {}
                f = mt.get("features") or {}
                row.update(e_hull_self=mt.get("e_hull_ml_self_eV"), e_hull_dftref=mt.get("e_hull_ml_dftref_eV"),
                           min_freq_cm1=ph.get("min_freq_cm1"), omega2_light_K=ph.get("omega2_light_K"),
                           H_f=f.get("H_f"), motif=f.get("motif"), score=mt.get("score_s3"),
                           relax_converged=(mt.get("relax") or {}).get("converged"))
        rows.append(row)
    summ = {"timing": timing, "counts": L.summary_counts(), "rows": rows,
            "lane_spec": spec.to_json(), "seed": a.seed, "n": a.n}
    (out / "smoke_summary.json").write_text(json.dumps(summ, indent=1, default=str))
    L.export_jsonl(out / "export")
    print(json.dumps(summ, indent=1, default=str), flush=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())
