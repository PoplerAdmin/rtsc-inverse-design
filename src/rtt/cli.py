"""``rtt`` command line (architecture §5.1): run | step | status | redecide | evidence | export |
controls | rank | decide-check.

Every command works on one run directory (``--run runs/<run_id>``, ledger at
``<run>/ledger.sqlite``). Compute commands respect docs/CONTEXT.md: callers launch them with
``nice -n 19`` and <= 2 threads; ``rtt`` itself sets the thread variables for its children.
"""
from __future__ import annotations

import argparse
import json
import os
import sys
from pathlib import Path

from .config import load_machine, load_targets, load_tiers, thread_env


def _ctx(args):
    from .engine.stage import Context
    from .ledger.store import Ledger
    run = Path(args.run)
    # the live ledger may sit on a Linux file system (SQLite WAL is unreliable on /mnt/<drive>):
    # --ledger / RTT_LEDGER; the run directory keeps the JSONL export (+ a snapshot copy)
    lp = getattr(args, "ledger", None) or os.environ.get("RTT_LEDGER")
    L = Ledger(Path(os.path.expanduser(lp)) if lp else run / "ledger.sqlite", run_id=run.name)
    t, ti, m = load_targets(args.targets), load_tiers(args.tiers), load_machine(args.machine)
    L.register_configs(t, ti, m)
    os.environ.update(thread_env())
    # big files (QE scratch, MatterGen output) belong outside the repo: --workdir / RTT_WORKDIR
    wd = getattr(args, "workdir", None) or os.environ.get("RTT_WORKDIR")
    workdir = Path(os.path.expanduser(wd)) if wd else run / "work"
    launcher = None
    ranks = getattr(args, "ranks", None)
    if ranks:
        from .dft.workflow import default_launcher
        launcher = default_launcher(int(ranks))
    return Context(ledger=L, targets=t, tiers=ti, machine=m, workdir=workdir, launcher=launcher)


def cmd_status(args) -> int:
    ctx = _ctx(args)
    try:
        from .report.summary import daily_summary
        print(daily_summary(ctx.ledger))
    except Exception:  # noqa: BLE001 - summary module optional
        print(json.dumps(ctx.ledger.summary_counts(), indent=1))
    return 0


def cmd_step(args) -> int:
    from .engine.runner import run_job
    from .engine.scheduler import Job
    from .engine.stages import STAGES
    ctx = _ctx(args)
    st = STAGES[args.stage]
    if args.stage in ("S2", "S3"):
        res = st.run_batch(ctx, args.cids or None)
        print(json.dumps({k: v.outcome for k, v in res.items()}, indent=1))
        return 0
    if args.stage == "S1":
        cids = st.run_batch(ctx, database_systems=args.systems or (), controls=args.controls)
        print("\n".join(cids))
        return 0
    rc = 0
    for cid in args.cids:
        if not gate_ok(ctx, cid, args.stage, getattr(args, "override_reason", None)):
            rc = 2
            continue
        o = run_job(Job(cid, args.stage, "", 0.0), ctx)
        print(cid, getattr(o, "outcome", o))
    return rc


def gate_ok(ctx, cid: str, stage: str, override_reason: str | None) -> bool:
    """R30 F5: ``rtt step`` honours the stage's ``requires`` (e.g. S6 needs a passing S5). A
    deliberate bypass needs ``--override-reason``; it is recorded in the ledger (meta key
    ``gate_override:<stage>:<cid>:<n>``, JSON with the unmet requirement, its current outcome and
    the reason) before the stage runs."""
    from .engine.stages import STAGES
    st = STAGES[stage]
    rec = ctx.ledger.candidate_record(cid)
    if st.ready(rec):
        return True
    unmet = {s: (rec.stages[s].outcome if s in rec.stages else None) for s in st.requires if not rec.passed(s)}
    if not override_reason:
        print(f"{cid}: refused: {stage} requires {sorted(unmet)} to pass (current: {unmet}); "
              "use --override-reason '<why>' for a recorded bypass", file=sys.stderr)
        return False
    import datetime as _dt
    n = ctx.ledger.q1("SELECT COUNT(*) FROM meta WHERE key LIKE ?", (f"gate_override:{stage}:{cid}:%",))[0] + 1
    ctx.ledger.set_meta_once(f"gate_override:{stage}:{cid}:{n}", json.dumps(
        {"candidate_id": cid, "stage": stage, "unmet": unmet, "reason": override_reason,
         "time": _dt.datetime.now(_dt.timezone.utc).isoformat(timespec="seconds")}, sort_keys=True))
    print(f"{cid}: GATE OVERRIDE recorded ({stage} without {sorted(unmet)}): {override_reason}", file=sys.stderr)
    return True


def cmd_run(args) -> int:
    from .engine.runner import run_job
    from .engine.scheduler import next_jobs
    ctx = _ctx(args)
    jobs = next_jobs(ctx.ledger, ctx.machine, max_jobs=args.max_jobs, tiers=ctx.tiers)
    for j in jobs:
        print(f"{j.candidate_id} {j.stage} lane={j.lane} est={j.est_core_h:.2f} core-h", flush=True)
        if not args.plan_only:
            o = run_job(j, ctx)
            print("  ->", getattr(o, "outcome", o), flush=True)
    return 0


def cmd_redecide(args) -> int:
    """Re-evaluate decisions from stored records (no compute); appends rows, never overwrites."""
    from .engine.stages import S7Validity, write_outcome
    ctx = _ctx(args)
    n = 0
    if args.stage in ("S6", "S7"):
        for (cid,) in ctx.ledger.qa("SELECT DISTINCT candidate_id FROM decisions WHERE stage IN ('S6','S7')"):
            rec = ctx.ledger.candidate_record(cid)
            o = S7Validity().redecide(rec, ctx)
            print(cid, o.outcome if o else "-")
            n += 1
    elif args.stage == "S5":
        # doc 16: re-apply the current T2-SCREEN rule to stored S5 metrics. By default only
        # current SCREENED_OUT outcomes are re-decided (re-opening removals made by an older
        # rule); --all also re-decides PASS/DEFERRED. A row is appended only if the outcome or
        # criterion changes (append-only ledger: the new row becomes current).
        from .engine.stages import S5Coarse
        st = S5Coarse()
        want = ("SCREENED_OUT",) if not getattr(args, "all", False) else ("SCREENED_OUT", "PASS", "DEFERRED")
        dry = getattr(args, "dry_run", False)
        cids = getattr(args, "cids", None) or [r[0] for r in ctx.ledger.qa(
            "SELECT DISTINCT candidate_id FROM outcomes WHERE stage='S5'")]
        for cid in cids:
            cur = ctx.ledger.current_outcome(cid, "S5")
            if not cur or cur["outcome"] not in want:
                continue
            rec = ctx.ledger.candidate_record(cid)
            o = st.redecide(rec, ctx)
            if o is None:
                print(cid, cur["outcome"], "-> no stored per-mesh metrics, unchanged")
                continue
            same = (o.outcome, o.criterion) == (cur["outcome"], cur.get("criterion"))
            print(cid, cur["outcome"], "->", o.outcome, o.criterion or "",
                  "(unchanged)" if same else ("[dry run]" if dry else ""))
            print("   ", o.reason)
            if not same and not dry:
                write_outcome(ctx, cid, o, rec=rec)
                n += 1
    print(f"redecided {n}")
    return 0


def cmd_evidence(args) -> int:
    from .report.evidence import evidence_statement
    ctx = _ctx(args)
    print(evidence_statement(args.cid, ctx.ledger))
    return 0


def cmd_export(args) -> int:
    ctx = _ctx(args)
    print(ctx.ledger.export_jsonl(args.out or Path(args.run) / "export"))
    return 0


def cmd_controls(args) -> int:
    from .engine.stages import STAGES
    ctx = _ctx(args)
    print("\n".join(STAGES["S1"].run_batch(ctx, controls=True)))
    return 0


def cmd_rank(args) -> int:
    from .score.rank import rows_from_ledger
    ctx = _ctx(args)
    for r in rows_from_ledger(ctx.ledger):
        print(json.dumps(r.to_json()))
    return 0


def cmd_decide_check(args) -> int:
    """Run decide_t on every case of lean/export/decide_cases.json (mirror check)."""
    from .decide.core import decide_t
    from .decide.lean_io import config_from_export, inputs_from_export
    d = json.loads(Path(args.cases).read_text(encoding="utf-8"))
    cfg = config_from_export(d["config"])
    bad = [c["id"] for c in d["cases"] if decide_t(cfg, inputs_from_export(c["inputs"])).outcome
           != c["expected"]["lean_outcome"]]
    print(f"{len(d['cases']) - len(bad)}/{len(d['cases'])} agree" + (f"; disagree: {bad}" if bad else ""))
    return 1 if bad else 0


def cmd_init(args) -> int:
    """Detect the hardware and write $RTT_HOME/machine.toml (portable install, campaign/hw.py)."""
    from .campaign import hw
    from .campaign.home import machine_path, rtt_home
    h = hw.detect(rtt_home())
    mp = machine_path()
    print(json.dumps(h.as_dict(), indent=1))
    for w in hw.minimum_report(h):
        print("warning:", w)
    if mp.exists() and not args.force:
        print(f"{mp} exists (use --force to rewrite)")
        return 0
    mp.parent.mkdir(parents=True, exist_ok=True)
    mp.write_text(hw.machine_toml(h, scratch=str(rtt_home() / "scratch"), max_job_wall_h=args.max_job_wall_h,
                                  gpu_h_per_round=args.gpu_h), encoding="utf-8")
    print(f"wrote {mp}")
    return 0


def cmd_doctor(args) -> int:
    from .campaign.doctor import report
    return report()


def _campaign_rd(args) -> Path:
    from .campaign.home import run_dir
    return run_dir(args.run_id)


def cmd_campaign(args) -> int:
    from .campaign import controller as C
    from .campaign.home import controller_pid, machine_path, read_json
    rd = _campaign_rd(args)
    if args.action == "start":
        if args.foreground:
            return C.main(["--run", str(rd), "--machine", str(machine_path()), "--rounds", str(args.rounds)]
                          + (["--no-front"] if args.no_front else []) + (["--smoke"] if args.smoke else []))
        pid = C.start_detached(rd, rounds=args.rounds, front=not args.no_front, machine=machine_path(),
                               smoke=args.smoke)
        print(f"controller pid {pid} for {rd}" if pid else f"start failed: see {rd / 'controller.log'}")
        return 0 if pid else 1
    if args.action in ("drain", "stop"):
        if not controller_pid(rd):
            print("not running")
            return 0
        C.request(rd, args.action)
        print(f"{args.action} requested for {rd}")
        return 0
    st = read_json(rd / "status.json", {}) or {}
    alive = controller_pid(rd)
    print(f"run {rd}: {'running pid ' + str(alive) if alive else 'not running'}; state={st.get('state')} "
          f"phase={st.get('phase')} round={st.get('round')} ranks={st.get('used_ranks')}/{st.get('max_ranks')}")
    for j in st.get("jobs", []):
        est = j.get("est_h")
        print(f"  {j['key']:28s} ranks={j['ranks']} {j['elapsed_h']:.2f} h / est {f'{est:.2f} h' if est else '-'}")
    for x in st.get("compute_limited", []):
        print(f"  compute-limited {x['key']}: ~{x['est_wall_h']} h")
    return 0


def cmd_ui(args) -> int:
    from .campaign.ui import run
    return run(_campaign_rd(args), rounds=args.rounds)


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(prog="rtt")
    ap.add_argument("--run", default="runs/dev")
    ap.add_argument("--targets")
    ap.add_argument("--tiers")
    ap.add_argument("--machine")
    ap.add_argument("--ledger", help="ledger path (default <run>/ledger.sqlite; env RTT_LEDGER)")
    ap.add_argument("--workdir", help="compute scratch (default <run>/work; env RTT_WORKDIR)")
    ap.add_argument("--ranks", type=int, help="MPI ranks per QE job (default: the dft default, 2)")
    sub = ap.add_subparsers(dest="cmd", required=True)
    sub.add_parser("status").set_defaults(f=cmd_status)
    p = sub.add_parser("step")
    p.add_argument("stage", choices=["S1", "S2", "S3", "S4", "S5", "S6", "S7"])
    p.add_argument("cids", nargs="*")
    p.add_argument("--systems", nargs="*")
    p.add_argument("--controls", action="store_true")
    p.add_argument("--override-reason", default=None,
                   help="run although the stage's requirement is not met (recorded in the ledger)")
    p.set_defaults(f=cmd_step)
    p = sub.add_parser("run")
    p.add_argument("--max-jobs", type=int, default=None)
    p.add_argument("--plan-only", action="store_true")
    p.set_defaults(f=cmd_run)
    p = sub.add_parser("redecide")
    p.add_argument("--stage", default="S6", choices=["S5", "S6", "S7"])
    p.add_argument("cids", nargs="*", help="S5 only: restrict to these candidates")
    p.add_argument("--all", action="store_true", help="S5: also re-decide PASS/DEFERRED")
    p.add_argument("--dry-run", action="store_true", help="S5: print, do not append")
    p.set_defaults(f=cmd_redecide)
    p = sub.add_parser("evidence")
    p.add_argument("cid")
    p.set_defaults(f=cmd_evidence)
    p = sub.add_parser("export")
    p.add_argument("--out")
    p.set_defaults(f=cmd_export)
    sub.add_parser("controls").set_defaults(f=cmd_controls)
    sub.add_parser("rank").set_defaults(f=cmd_rank)
    p = sub.add_parser("decide-check")
    p.add_argument("--cases", default=str(Path(__file__).resolve().parents[2] / "lean" / "export" / "decide_cases.json"))
    p.set_defaults(f=cmd_decide_check)
    # portable runner (src/rtt/campaign): install home $RTT_HOME, runs in $RTT_HOME/runs/<run-id>
    p = sub.add_parser("init", help="detect hardware, write $RTT_HOME/machine.toml")
    p.add_argument("--force", action="store_true")
    p.add_argument("--max-job-wall-h", type=float, default=48.0)
    p.add_argument("--gpu-h", type=float, default=2.0, help="MatterGen GPU hours per round")
    p.set_defaults(f=cmd_init)
    sub.add_parser("doctor", help="check the installation").set_defaults(f=cmd_doctor)
    p = sub.add_parser("campaign", help="start | drain | stop | status the background controller")
    p.add_argument("action", choices=["start", "drain", "stop", "status"])
    p.add_argument("--run-id", default=None, help="run name under $RTT_HOME/runs (default $RTT_RUN or 'default')")
    p.add_argument("--rounds", type=int, default=1, help="generation+verification rounds (0 = until stopped)")
    p.add_argument("--no-front", action="store_true", help="skip S1-S3 (continue verifying the existing ledger)")
    p.add_argument("--foreground", action="store_true")
    p.add_argument("--smoke", action="store_true",
                   help="installation check: the Al control through S1-S5 only (about an hour)")
    p.set_defaults(f=cmd_campaign)
    p = sub.add_parser("ui", help="terminal dashboard")
    p.add_argument("--run-id", default=None)
    p.add_argument("--rounds", type=int, default=1)
    p.set_defaults(f=cmd_ui)
    args = ap.parse_args(argv)
    return int(args.f(args) or 0)


if __name__ == "__main__":
    sys.exit(main())
