"""The campaign controller behind ``rtt campaign start`` (one per run directory).

One *round* = front end (S1 generation and database/prototype/control lanes, S2, S3; one child
process on the GPU) followed by the back end: S4-S6 jobs from :func:`rtt.engine.scheduler.next_jobs`,
each a child ``rtt step <stage> <cid>`` process, run in parallel while their MPI ranks fit in
``[budget] max_ranks``. A round ends when no job is running and none can start. ``--rounds 0``
repeats rounds until stopped.

Every child runs in its own session (process group), so stopping kills QE's whole ``mpirun``
tree. Stop requests come from the ``control`` file (``drain`` / ``stop``), SIGTERM or SIGINT.
``drain`` lets running jobs finish and starts nothing new; ``stop`` terminates every child
(SIGTERM, then SIGKILL after a grace period) and marks their open ledger runs failed, so the
candidates are retried on the next start (a killed run is never recorded as a result).

Jobs whose estimated wall time (est. core-h / ranks) exceeds ``[budget] max_job_wall_h`` are not
started; they are listed as ``compute_limited`` in ``status.json`` and the events log.

Process layout (diagrams and design reasons: docs/RUNNER.md §2-3)::

    rtt ui ──reads──> status.json, events.log, jobs/*.log <──writes── controller
    rtt ui / stop.sh ──writes──> control ──read every 5 s──> controller
    controller ──spawns (own session)──> front child (S1-S3, GPU)
                                      ├─> rtt step S5 <cid> ──> mpirun ──> ph.x x4
                                      └─> rtt step S4 <cid> ──> mpirun ──> pw.x x4
    every process ──appends──> ledger.sqlite (WAL; rows cannot be edited or deleted)

The controller is separate from the dashboard so that a closed terminal never stops days of
compute; each job is a separate process so that one crashing QE run cannot take down the rest.
"""
from __future__ import annotations

import argparse
import collections
import datetime as dt
import os
import signal
import subprocess
import sys
import time
from dataclasses import dataclass, field
from pathlib import Path

from .home import controller_pid, machine_path, read_json, run_dir, write_json_atomic

POLL_S = 5.0
REPLAN_S = 60.0
COUNTS_S = 30.0
KILL_GRACE_S = 30.0


def now_iso() -> str:
    return dt.datetime.now().astimezone().isoformat(timespec="seconds")


@dataclass
class Child:
    key: str                      # "front-<round>" or "<stage>:<cid>"
    kind: str                     # front | S4 | S5 | S6 | S7
    proc: subprocess.Popen
    ranks: int
    started: float
    log: str
    cid: str | None = None
    est_h: float | None = None
    gpu: bool = False


@dataclass
class Controller:
    rd: Path
    machine: Path
    rounds: int = 1
    front: bool = True
    gpu_h: float | None = None
    targets: str | None = None
    tiers: str | None = None
    max_stage: str | None = None          # e.g. "S5": never schedule later stages
    controls_only: list | None = None     # smoke test: only these controls in the front end
    children: dict = field(default_factory=dict)
    state: str = "starting"
    phase: str = "idle"
    round: int = 0
    compute_limited: dict = field(default_factory=dict)
    finished_jobs: collections.deque = field(default_factory=lambda: collections.deque(maxlen=50))
    counts: dict = field(default_factory=dict)
    top: list = field(default_factory=list)
    _stop_signal: str | None = None
    _last_counts: float = 0.0
    _started: str = field(default_factory=now_iso)

    # ------------------------------------------------------------------ setup
    def __post_init__(self):
        self.rd.mkdir(parents=True, exist_ok=True)
        (self.rd / "jobs").mkdir(exist_ok=True)
        self.workdir = self.rd / "work"
        self.ledger_path = self.rd / "ledger.sqlite"
        from ..config import load_machine
        self.m = load_machine(self.machine)
        self.max_ranks = int(self.m.get("budget.max_ranks", 4))
        self.max_wall_h = float(self.m.get("budget.max_job_wall_h", 48.0))
        self.s_ranks = int(self.m.get("lanes.S.mpi_ranks", 4))
        self.f_ranks = int(self.m.get("lanes.F.mpi_ranks", 4))
        self.gpu = bool(self.m.get("lanes.G.gpu", False))
        if self.gpu_h is None:
            self.gpu_h = float(self.m.get("budget.gpu_h_per_round", 2.0))

    def event(self, msg: str) -> None:
        line = f"{now_iso()} {msg}"
        with open(self.rd / "events.log", "a", encoding="utf-8") as fh:
            fh.write(line + "\n")
        print(line, flush=True)

    def _ledger(self):
        from ..config import load_targets, load_tiers
        from ..ledger.store import Ledger
        L = Ledger(self.ledger_path, run_id=self.rd.name)
        t, ti = load_targets(self.targets), load_tiers(self.tiers)
        L.register_configs(t, ti, self.m)
        return L, ti

    def base_argv(self) -> list[str]:
        a = [sys.executable, "-m", "rtt.cli", "--run", str(self.rd), "--ledger", str(self.ledger_path),
             "--workdir", str(self.workdir), "--machine", str(self.machine)]
        if self.targets:
            a += ["--targets", self.targets]
        if self.tiers:
            a += ["--tiers", self.tiers]
        return a

    # ------------------------------------------------------------------ control
    def control_word(self) -> str | None:
        if self._stop_signal:
            return self._stop_signal
        try:
            w = (self.rd / "control").read_text().strip().lower()
        except OSError:
            return None
        return w if w in ("drain", "stop") else None

    def _on_signal(self, signum, _frame):
        self._stop_signal = "stop"

    def used_ranks(self) -> int:
        return sum(c.ranks for c in self.children.values())

    def spawn(self, key: str, kind: str, argv: list[str], ranks: int, *, cid=None, est_h=None, gpu=False) -> None:
        log = self.rd / "jobs" / (key.replace(":", "_") + ".log")
        fh = open(log, "a", encoding="utf-8")
        fh.write(f"\n=== {now_iso()} {' '.join(argv)}\n")
        fh.flush()
        env = dict(os.environ)
        env.setdefault("PYTHONUNBUFFERED", "1")
        p = subprocess.Popen(argv, stdout=fh, stderr=subprocess.STDOUT, stdin=subprocess.DEVNULL,
                             start_new_session=True, env=env, cwd=str(self.rd))
        fh.close()
        self.children[key] = Child(key, kind, p, ranks, time.time(), str(log), cid=cid, est_h=est_h, gpu=gpu)
        self.event(f"start {key} pid={p.pid} ranks={ranks}" + (f" est={est_h:.1f}h" if est_h else ""))

    def reap(self) -> list[Child]:
        done = []
        for key, c in list(self.children.items()):
            rc = c.proc.poll()
            if rc is None:
                continue
            del self.children[key]
            wall = time.time() - c.started
            self.finished_jobs.append({"key": key, "rc": rc, "wall_h": round(wall / 3600, 2), "end": now_iso(),
                                       "log": c.log})
            self.event(f"end {key} rc={rc} wall={wall / 3600:.2f}h")
            done.append(c)
        return done

    def kill_all(self) -> None:
        """End every child and everything it started. Killing the child's process group is not
        enough: the core QE launcher starts ``mpirun`` with ``start_new_session=True``, so the
        MPI ranks live in their own session (found by the stop test on the reference machine).
        So: collect the whole descendant tree first (while the parents still exist, so the
        ancestry is known), signal every process group in it, then sweep for anything still
        running inside this run's work directory."""
        pids = set()
        for c in self.children.values():
            pids |= {c.proc.pid, *descendants(c.proc.pid)}
        _signal_all(pids, signal.SIGTERM)
        t0 = time.time()
        while time.time() - t0 < KILL_GRACE_S and any(_alive(p) for p in pids):
            self.reap()
            time.sleep(1)
        _signal_all({p for p in pids if _alive(p)}, signal.SIGKILL)
        stray = procs_in_dir(self.workdir)
        if stray:
            self.event(f"killing {len(stray)} leftover processes in {self.workdir}")
            _signal_all(stray, signal.SIGKILL)
        time.sleep(1)
        self.reap()
        self.close_open_runs("stopped by user")

    def close_open_runs(self, why: str) -> int:
        """Mark ledger runs left 'started' (by killed children or a crash) as failed, so the
        scheduler retries them; nothing is recorded as a result."""
        try:
            L, _ = self._ledger()
            ids = [r[0] for r in L.qa("SELECT id FROM tier_runs WHERE state='started'")]
            for i in ids:
                L.finish_tier_run(i, "failed", metrics={"interrupted": why, "at": now_iso()})
            if ids:
                self.event(f"closed {len(ids)} open ledger runs as failed ({why})")
            return len(ids)
        except Exception as ex:  # noqa: BLE001
            self.event(f"could not close open runs: {ex!r}")
            return 0

    # ------------------------------------------------------------------ status
    def refresh_counts(self) -> None:
        if time.time() - self._last_counts < COUNTS_S or not self.ledger_path.exists():
            return
        self._last_counts = time.time()
        try:
            L, _ = self._ledger()
            self.counts = L.summary_counts()
            self.top = top_candidates(L)
        except Exception as ex:  # noqa: BLE001
            self.counts = {"error": repr(ex)}

    def write_status(self) -> None:
        self.refresh_counts()
        jobs = [{"key": c.key, "kind": c.kind, "cid": c.cid, "ranks": c.ranks, "pid": c.proc.pid,
                 "elapsed_h": round((time.time() - c.started) / 3600, 2), "est_h": c.est_h, "log": c.log,
                 "gpu": c.gpu} for c in self.children.values()]
        write_json_atomic(self.rd / "status.json", {
            "state": self.state, "phase": self.phase, "round": self.round, "rounds": self.rounds,
            "pid": os.getpid(), "started": self._started, "updated": now_iso(),
            "max_ranks": self.max_ranks, "used_ranks": self.used_ranks(), "gpu": self.gpu,
            "max_job_wall_h": self.max_wall_h, "jobs": jobs, "finished": list(self.finished_jobs)[-15:],
            "compute_limited": list(self.compute_limited.values()), "counts": self.counts, "top": self.top})

    # ------------------------------------------------------------------ phases
    def run_front(self) -> bool:
        self.phase = "front"
        argv = [sys.executable, "-m", "rtt.campaign.front", "--run", str(self.rd), "--machine", str(self.machine),
                "--round", str(self.round), "--gpu-h", str(self.gpu_h if self.gpu else 0.0)]
        if self.controls_only is not None:
            argv += ["--controls-only", *self.controls_only]
        self.spawn(f"front-{self.round}", "front", argv, ranks=1, gpu=self.gpu)
        while self.children:
            w = self.control_word()
            if w == "stop":
                return False
            if w == "drain":
                self.state = "draining"
            self.reap()
            self.write_status()
            time.sleep(POLL_S)
        return self.control_word() is None

    def plan_jobs(self, L, tiers) -> list:
        from ..engine.scheduler import next_jobs
        running = {c.key for c in self.children.values()}
        out = []
        # the per-day lane budgets of next_jobs do not apply here: the rank budget and the wall
        # cap below are the limits of a continuously running controller
        big = {"S": -1e12, "F": -1e12}
        for j in next_jobs(L, self.m, spent_core_h=big, tiers=tiers):
            key = f"{j.stage}:{j.candidate_id}"
            if key in running:
                continue
            if self.max_stage and int(j.stage[1:]) > int(self.max_stage[1:]):
                continue
            ranks =self.f_ranks if j.stage in ("S6", "S7") else self.s_ranks
            est_h = (j.est_core_h or 0.0) / max(1, ranks) / float(self.m.get("derate.contention", 0.85))
            if est_h > self.max_wall_h:
                if key not in self.compute_limited:
                    self.compute_limited[key] = {"key": key, "est_core_h": round(j.est_core_h, 1),
                                                 "est_wall_h": round(est_h, 1), "ranks": ranks}
                    self.event(f"compute-limited {key}: est {est_h:.0f} h on {ranks} ranks > {self.max_wall_h} h cap")
                continue
            out.append((j, key, ranks, est_h))
        return out

    def run_back(self) -> bool:
        self.phase = "back"
        L, tiers = self._ledger()
        last_plan = 0.0
        plan: list = []
        while True:
            w = self.control_word()
            if w == "stop":
                return False
            if w == "drain":
                self.state = "draining"
            finished = self.reap()
            if finished or time.time() - last_plan > REPLAN_S:
                plan = self.plan_jobs(L, tiers) if w is None else []
                last_plan = time.time()
            if w is None:
                for j, key, ranks, est_h in list(plan):
                    if key in self.children:
                        continue
                    if self.used_ranks() + ranks > self.max_ranks:
                        continue
                    argv = self.base_argv() + ["--ranks", str(ranks), "step", j.stage, j.candidate_id]
                    self.spawn(key, j.stage, argv, ranks, cid=j.candidate_id, est_h=est_h)
                    plan.remove((j, key, ranks, est_h))
            self.write_status()
            if not self.children and (w is not None or not plan):
                return w is None
            time.sleep(POLL_S)

    # ------------------------------------------------------------------ main
    def run(self) -> int:
        other = controller_pid(self.rd)
        if other and other != os.getpid():
            print(f"a controller is already running for {self.rd} (pid {other})", file=sys.stderr)
            return 1
        (self.rd / "controller.pid").write_text(str(os.getpid()))
        try:
            (self.rd / "control").unlink()
        except OSError:
            pass
        import threading
        if threading.current_thread() is threading.main_thread():
            signal.signal(signal.SIGTERM, self._on_signal)
            signal.signal(signal.SIGINT, self._on_signal)
        self.state = "running"
        self.event(f"controller start: rounds={self.rounds or 'until stopped'} max_ranks={self.max_ranks} "
                   f"gpu={self.gpu} cap={self.max_wall_h} h")
        self.close_open_runs("left open by an earlier controller")
        rc = 0
        try:
            k = 0
            while self.rounds == 0 or k < self.rounds:
                k += 1
                self.round = k
                self.event(f"round {k}")
                if self.front and not self.run_front():
                    break
                if not self.run_back():
                    break
                if self.control_word():
                    break
            w = self.control_word()
            if w == "stop":
                self.state = "stopping"
                self.write_status()
                self.event("stop requested: terminating jobs")
                self.kill_all()
                self.state = "stopped"
            elif w == "drain":
                self.state = "stopped"
                self.event("drained: all jobs finished")
            else:
                self.state = "finished"
                self.event("all rounds finished")
        except Exception as ex:  # noqa: BLE001
            import traceback
            self.state = "error"
            self.event(f"controller error: {ex!r}\n{traceback.format_exc()}")
            self.kill_all()
            rc = 1
        finally:
            self.phase = "idle"
            self._last_counts = 0.0
            self.write_status()
            try:
                (self.rd / "controller.pid").unlink()
            except OSError:
                pass
        return rc


def _ppid_map() -> dict[int, int]:
    out = {}
    for d in Path("/proc").iterdir():
        if not d.name.isdigit():
            continue
        try:
            stat = (d / "stat").read_text()
            # comm may contain spaces/parentheses: the fields after the last ')' are fixed
            out[int(d.name)] = int(stat.rsplit(")", 1)[1].split()[1])
        except (OSError, ValueError, IndexError):
            pass
    return out


def descendants(pid: int) -> set[int]:
    pm = _ppid_map()
    kids: dict[int, list[int]] = collections.defaultdict(list)
    for p, pp in pm.items():
        kids[pp].append(p)
    out, todo = set(), [pid]
    while todo:
        for k in kids.get(todo.pop(), ()):
            if k not in out:
                out.add(k)
                todo.append(k)
    return out


def procs_in_dir(root: Path) -> set[int]:
    """Processes of this user whose working directory is inside ``root`` (QE runs in work/)."""
    root = Path(root).resolve()
    out = set()
    for d in Path("/proc").iterdir():
        if not d.name.isdigit() or int(d.name) == os.getpid():
            continue
        try:
            cwd = Path(os.readlink(d / "cwd"))
        except OSError:
            continue
        if cwd == root or root in cwd.parents:
            out.add(int(d.name))
    return out


def _alive(pid: int) -> bool:
    try:
        os.kill(pid, 0)
    except ProcessLookupError:
        return False
    except PermissionError:
        return True
    try:       # a zombie (exited, not yet reaped) counts as gone
        return Path(f"/proc/{pid}/stat").read_text().rsplit(")", 1)[1].split()[0] != "Z"
    except (OSError, IndexError):
        return True


def _signal_all(pids, sig) -> None:
    groups = set()
    for p in pids:
        try:
            groups.add(os.getpgid(p))
        except ProcessLookupError:
            pass
    own = os.getpgid(0)
    for g in groups - {own}:
        try:
            os.killpg(g, sig)
        except (ProcessLookupError, PermissionError):
            pass
    for p in pids:
        try:
            os.kill(p, sig)
        except (ProcessLookupError, PermissionError):
            pass


def _formula(species) -> str:
    c = collections.Counter(species)
    return "".join(f"{el}{n if n > 1 else ''}" for el, n in sorted(c.items()))


def top_candidates(L, n: int = 15) -> list[dict]:
    """Candidates that reached S4 or later, furthest stage first, with their latest outcome."""
    order = {"S7": 7, "S6": 6, "S5": 5, "S4": 4}
    rows = []
    for (cid, lane) in L.qa("SELECT DISTINCT o.candidate_id, c.lane FROM outcomes o JOIN candidates c "
                            "ON c.id = o.candidate_id WHERE o.stage IN ('S4','S5','S6','S7')"):
        cur = L.current_outcomes(cid)
        st = max((s for s in cur if s in order), key=lambda s: order[s], default=None)
        if st is None:
            continue
        o = cur[st]
        try:
            formula = _formula(L.candidate_record(cid).structure.species)
        except Exception:  # noqa: BLE001
            formula = "?"
        rows.append({"cid": cid, "formula": formula, "lane": lane, "stage": st, "outcome": o.get("outcome"),
                     "value": o.get("value"), "criterion": o.get("criterion")})
    rows.sort(key=lambda r: (-order[r["stage"]], r["outcome"] != "PASS", r["cid"]))
    return rows[:n]


# ---------------------------------------------------------------------- CLI helpers
def request(rd: Path, word: str) -> None:
    (rd / "control").write_text(word + "\n")


def start_detached(rd: Path, *, rounds: int = 1, front: bool = True, machine: Path | None = None,
                   gpu_h: float | None = None, smoke: bool = False) -> int | None:
    """Start the controller as a detached background process; returns its pid. ``smoke``: the
    installation check (control Al only, stages S1-S5, one round)."""
    if controller_pid(rd):
        return controller_pid(rd)
    rd.mkdir(parents=True, exist_ok=True)
    argv = [sys.executable, "-m", "rtt.campaign.controller", "--run", str(rd), "--rounds", str(rounds),
            "--machine", str(machine or machine_path())]
    if not front:
        argv.append("--no-front")
    if gpu_h is not None:
        argv += ["--gpu-h", str(gpu_h)]
    if smoke:
        argv += ["--smoke"]
    log = open(rd / "controller.log", "a", encoding="utf-8")
    p = subprocess.Popen(argv, stdout=log, stderr=subprocess.STDOUT, stdin=subprocess.DEVNULL,
                         start_new_session=True, cwd=str(rd))
    log.close()
    for _ in range(20):
        if controller_pid(rd):
            break
        if p.poll() is not None:
            return None
        time.sleep(0.25)
    return p.pid


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(prog="rtt-controller")
    ap.add_argument("--run", default=None, help="run directory (default $RTT_HOME/runs/$RTT_RUN or default)")
    ap.add_argument("--machine", default=None)
    ap.add_argument("--rounds", type=int, default=1, help="0 = until stopped")
    ap.add_argument("--no-front", action="store_true")
    ap.add_argument("--gpu-h", type=float, default=None)
    ap.add_argument("--targets")
    ap.add_argument("--tiers")
    ap.add_argument("--max-stage", default=None, help="never schedule stages after this one (e.g. S5)")
    ap.add_argument("--smoke", action="store_true", help="installation check: control Al, S1-S5, one round")
    a = ap.parse_args(argv)
    rd = Path(a.run).expanduser() if a.run else run_dir()
    c = Controller(rd=rd, machine=Path(a.machine).expanduser() if a.machine else machine_path(),
                   rounds=1 if a.smoke else a.rounds, front=not a.no_front, gpu_h=a.gpu_h, targets=a.targets,
                   tiers=a.tiers, max_stage="S5" if a.smoke else a.max_stage,
                   controls_only=["Al"] if a.smoke else None)
    return c.run()


if __name__ == "__main__":
    sys.exit(main())
