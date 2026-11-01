"""One front-end round (S1 lanes -> S2 -> S3), run by the controller as a child process.

Generalises ``scripts/campaign_front.py`` (campaign 1) and ``scripts/campaign2_lane_p.py``:

* round 1 only: controls (lane C) and the prototype-substitution lane P (doc 32);
* every round with a GPU: the MatterGen lanes of :func:`rtt.source.lanes.plan_day` for
  ``--gpu-h`` GPU hours (a job whose checkpoint is not in the local Hugging Face cache is skipped
  and logged; ``install.sh`` downloads them, nothing is downloaded here);
* lane D (database) for the chemical systems the G1/G4 jobs used (round 1: the first systems
  of the chemsys list), when the Alex-MP reference is installed;
* S2 on everything new, then S3 with one pool-wide rank cut (doc 31 B9).

The controller's ``control`` file is checked between MatterGen jobs, so a drain/stop request
ends the round early; everything already ingested stays in the ledger.
"""
from __future__ import annotations

import argparse
import os
import sys
import time
from pathlib import Path

from .home import machine_path


def log(msg: str) -> None:
    print(f"{time.strftime('%Y-%m-%dT%H:%M:%S')} {msg}", flush=True)


def stop_requested(rd: Path) -> bool:
    try:
        return (rd / "control").read_text().strip().lower() in ("drain", "stop")
    except OSError:
        return False


def have_checkpoint(name: str, machine=None) -> bool:
    from ..source.external import path_of
    hub = Path(path_of("hf_hub", machine))
    for p in hub.glob(f"models--microsoft--mattergen/snapshots/*/checkpoints/{name}/checkpoints/last.ckpt"):
        try:
            if p.resolve().stat().st_size > 1_000_000:
                return True
        except OSError:
            pass
    return False


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(prog="rtt-front")
    ap.add_argument("--run", required=True)
    ap.add_argument("--machine", default=None)
    ap.add_argument("--round", type=int, default=1)
    ap.add_argument("--gpu-h", type=float, default=2.0)
    ap.add_argument("--db-systems", type=int, default=12, help="round 1: database systems from the chemsys list")
    ap.add_argument("--controls-only", nargs="*", default=None,
                    help="smoke test: inject only these controls (e.g. Al), no generation or other lanes")
    a = ap.parse_args(argv)
    rd = Path(a.run).expanduser()
    os.environ["HF_HUB_OFFLINE"] = "1"
    from ..cli import _ctx
    ns = argparse.Namespace(run=str(rd), ledger=str(rd / "ledger.sqlite"), workdir=str(rd / "work"),
                            targets=None, tiers=None, machine=a.machine or str(machine_path()), ranks=None)
    ctx = _ctx(ns)
    L, m = ctx.ledger, ctx.machine
    from ..engine.stages import STAGES
    from ..source import lanes as LN
    from ..source.external import path_of

    if a.controls_only is not None:
        from ..source.controls import inject_controls
        cids = inject_controls(L, tuple(a.controls_only) or None, targets=ctx.targets, tiers=ctx.tiers)
        log(f"smoke: controls {a.controls_only} -> {cids}")
        res = STAGES["S2"].run_batch(ctx, None)
        log(f"S2: {len(res)} candidates decided")
        res = STAGES["S3"].run_batch(ctx, None)
        log(f"S3: {len(res)} candidates decided")
        return 0

    if a.round == 1:
        from ..source.controls import inject_controls
        log(f"controls: {len(inject_controls(L, targets=ctx.targets, tiers=ctx.tiers))} candidates")
        try:
            from ..source.prototype import run_lane_p
            log(f"lane P: {len(run_lane_p(L, targets=ctx.targets, tiers=ctx.tiers))} candidates")
        except Exception as ex:  # noqa: BLE001
            log(f"lane P skipped: {ex!r}")

    systems: set[str] = set()
    if a.gpu_h > 0 and m.get("lanes.G.gpu", False):
        from ..source.runner import run_lane_job
        jobs = LN.plan_day(m, ctx.targets, tiers=ctx.tiers, day=f"{rd.name}-r{a.round}", seed=a.round,
                           gpu_h=a.gpu_h, batch_size=int(m.get("lanes.G.mattergen_batch", 16)))
        log(f"MatterGen plan: {len(jobs)} jobs for {a.gpu_h} GPU-h: {LN.plan_summary(jobs)}")
        for j in jobs:
            if stop_requested(rd):
                log("stop requested: ending generation early")
                break
            if not have_checkpoint(j.spec.checkpoint, m):
                log(f"{j.batch_id} skipped: checkpoint {j.spec.checkpoint} not installed")
                continue
            t0 = time.time()
            try:
                cids = run_lane_job(L, j, targets=ctx.targets, tiers=ctx.tiers, machine=m,
                                    workdir=ctx.workdir / "S1")
            except Exception as ex:  # noqa: BLE001
                log(f"{j.batch_id} FAILED: {str(ex)[-300:]}")
                continue
            cs = dict(j.spec.conditioning).get("chemical_system")
            if cs and j.spec.lane in ("G1", "G4"):
                systems.add(cs)
            log(f"{j.batch_id} {j.spec.lane} {dict(j.spec.conditioning)} -> {len(cids)} candidates "
                f"in {time.time() - t0:.0f} s")
    else:
        log("no usable GPU: MatterGen lanes skipped")

    if a.round == 1 and not systems:
        hyd, _ = LN.load_chemsys_list(ctx.targets)
        systems = set(hyd[: a.db_systems])
    if systems and not stop_requested(rd):
        if Path(path_of("alexmp_lmdb", m)).exists():
            from ..source.runner import run_database_lane
            try:
                dc = run_database_lane(L, sorted(systems), targets=ctx.targets, tiers=ctx.tiers, machine=m)
                log(f"lane D: {len(dc)} entries from {len(systems)} systems")
            except Exception as ex:  # noqa: BLE001
                log(f"lane D failed: {ex!r}")
        else:
            log("lane D skipped: Alex-MP reference not installed (install.sh --with-reference)")

    if stop_requested(rd):
        return 0
    res = STAGES["S2"].run_batch(ctx, None)
    log(f"S2: {len(res)} candidates decided")
    if stop_requested(rd):
        return 0
    res = STAGES["S3"].run_batch(ctx, None)
    log(f"S3: {len(res)} candidates decided")
    return 0


if __name__ == "__main__":
    sys.exit(main())
