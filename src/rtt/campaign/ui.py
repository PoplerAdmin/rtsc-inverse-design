"""Terminal dashboard (``rtt ui``; ``start.sh`` opens it).

It reads only the controller's files in the run directory (``status.json``, ``events.log``,
``jobs/*.log``) and writes only its ``control`` file, so it can be closed and reopened at any time
without touching the running campaign.

Keys: ``s`` start  ``d`` drain (finish running jobs, start nothing new)  ``x`` stop now (ends
every job; asks first)  ``l`` cycle the log pane  ``q`` quit the dashboard (asks whether to leave
the campaign running).
"""
from __future__ import annotations

import os
import shutil
import subprocess
import time
from pathlib import Path

from textual.app import App, ComposeResult
from textual.binding import Binding
from textual.containers import Horizontal, Vertical
from textual.screen import ModalScreen
from textual.widgets import Button, DataTable, Footer, Header, Label, RichLog, Static

from .controller import request, start_detached
from .home import controller_pid, machine_path, read_json

STAGES = ("S1", "S2", "S3", "S4", "S5", "S6", "S7")
PASSISH = {"PASS", "pass", "D_PASS", "undecided", "D_UNDECIDED"}
OUTISH = {"SCREENED_OUT", "rejectedProven", "rejectedAssuming", "REJECTED_PROVEN", "DUPLICATE", "INVALID",
          "EXCLUDED_ASSUMPTION"}


def _cpu_sample():
    try:
        f = Path("/proc/stat").read_text().splitlines()[0].split()[1:]
        v = [int(x) for x in f]
        return sum(v), v[3] + (v[4] if len(v) > 4 else 0)
    except (OSError, ValueError, IndexError):
        return None


def _mem():
    try:
        d = {}
        for ln in Path("/proc/meminfo").read_text().splitlines():
            k, v = ln.split(":", 1)
            d[k] = int(v.split()[0])
        return (d["MemTotal"] - d["MemAvailable"]) / 1024 ** 2, d["MemTotal"] / 1024 ** 2
    except (OSError, KeyError, ValueError):
        return None


def _gpu():
    if not shutil.which("nvidia-smi"):
        return None
    try:
        out = subprocess.run(["nvidia-smi", "--query-gpu=utilization.gpu,memory.used,memory.total",
                              "--format=csv,noheader,nounits"], capture_output=True, text=True, timeout=5).stdout
        u, used, tot = [float(x) for x in out.splitlines()[0].split(",")]
        return u, used / 1024, tot / 1024
    except (OSError, subprocess.SubprocessError, ValueError, IndexError):
        return None


def _bar(frac: float, width: int = 16) -> str:
    frac = max(0.0, min(1.0, frac))
    n = int(round(frac * width))
    return "█" * n + "░" * (width - n)


def _tail(p: Path, n: int = 200) -> list[str]:
    try:
        with open(p, "rb") as fh:
            fh.seek(0, os.SEEK_END)
            size = fh.tell()
            fh.seek(max(0, size - 64 * 1024))
            return fh.read().decode("utf-8", "replace").splitlines()[-n:]
    except OSError:
        return []


class Confirm(ModalScreen[str]):
    """A question with a few answer buttons; dismisses with the chosen button id."""

    def __init__(self, question: str, buttons: list[tuple[str, str, str]]):
        super().__init__()
        self.question, self.buttons = question, buttons

    def compose(self) -> ComposeResult:
        with Vertical(id="dialog"):
            yield Label(self.question, id="question")
            with Horizontal(id="buttons"):
                for bid, text, variant in self.buttons:
                    yield Button(text, id=bid, variant=variant)

    def on_button_pressed(self, ev: Button.Pressed) -> None:
        self.dismiss(ev.button.id)


class Dashboard(App):
    TITLE = "rtsc-targeted"
    CSS = """
    Screen { layout: vertical; }
    #top { height: 12; }
    #machine { width: 44; border: round $accent; padding: 0 1; }
    #funnel { border: round $accent; }
    #jobs { height: 10; border: round $accent; }
    #mid { height: 1fr; min-height: 8; }
    #top_cands { border: round $accent; }
    #limited { width: 44; border: round $accent; padding: 0 1; }
    #log { height: 12; border: round $accent; }
    #state { height: 1; padding: 0 1; }
    Confirm { align: center middle; }
    #dialog { width: 72; height: auto; border: thick $warning; background: $surface; padding: 1 2; }
    #buttons { height: auto; margin-top: 1; }
    #buttons Button { margin-right: 2; }
    """
    BINDINGS = [
        Binding("s", "start", "Start"),
        Binding("d", "drain", "Drain"),
        Binding("x", "stop", "Stop now"),
        Binding("l", "cycle_log", "Log pane"),
        Binding("q", "quit_ui", "Quit"),
    ]

    def __init__(self, rd: Path, *, rounds: int = 1):
        super().__init__()
        self.rd, self.rounds = rd, rounds
        self._cpu_prev = _cpu_sample()
        self._log_choice = 0          # 0 = events, 1.. = running job logs
        self._log_sig = None

    def compose(self) -> ComposeResult:
        yield Header(show_clock=True)
        yield Static("", id="state")
        with Horizontal(id="top"):
            yield Static("", id="machine")
            yield DataTable(id="funnel", zebra_stripes=True)
        yield DataTable(id="jobs", zebra_stripes=True)
        with Horizontal(id="mid"):
            yield DataTable(id="top_cands", zebra_stripes=True)
            yield Static("", id="limited")
        yield RichLog(id="log", max_lines=400, wrap=False)
        yield Footer()

    def on_mount(self) -> None:
        self.sub_title = f"run {self.rd.name}  ({self.rd})"
        f = self.query_one("#funnel", DataTable)
        f.add_columns("stage", "pass", "out", "deferred", "error", "other")
        j = self.query_one("#jobs", DataTable)
        j.add_columns("job", "ranks", "elapsed h", "est h", "progress", "pid")
        t = self.query_one("#top_cands", DataTable)
        t.add_columns("candidate", "formula", "lane", "stage", "outcome", "value")
        self.query_one("#funnel").border_title = "funnel (current outcome per candidate)"
        self.query_one("#jobs").border_title = "running jobs"
        self.query_one("#top_cands").border_title = "furthest candidates"
        self.query_one("#limited").border_title = "compute-limited"
        self.query_one("#machine").border_title = "machine"
        self.refresh_all()
        self.set_interval(2.0, self.refresh_all)

    # ------------------------------------------------------------------ refresh
    def refresh_all(self) -> None:
        st = read_json(self.rd / "status.json", {}) or {}
        alive = controller_pid(self.rd) is not None
        state = st.get("state", "not started") if alive else (
            f"not running (last: {st.get('state')})" if st else "not started")
        jobs = st.get("jobs", [])
        self.query_one("#state", Static).update(
            f"[b]state:[/b] {state}   [b]phase:[/b] {st.get('phase', '-') if alive else '-'}   "
            f"[b]round:[/b] {st.get('round', 0)}/{st.get('rounds') or '∞'}   "
            f"[b]ranks:[/b] {st.get('used_ranks', 0)}/{st.get('max_ranks', '?')}   "
            f"[b]updated:[/b] {st.get('updated', '-')}")
        self._machine(st)
        self._funnel(st.get("counts", {}).get("outcomes", {}))
        self._jobs(jobs)
        self._top(st.get("top", []))
        lim = st.get("compute_limited", [])
        self.query_one("#limited", Static).update(
            "\n".join(f"{x['key']}\n  ~{x['est_wall_h']} h on {x['ranks']} ranks" for x in lim[:8])
            or f"none (cap {st.get('max_job_wall_h', 48)} h per job)")
        self._log(jobs)

    def _machine(self, st: dict) -> None:
        lines = []
        cur = _cpu_sample()
        if cur and self._cpu_prev:
            dt, di = cur[0] - self._cpu_prev[0], cur[1] - self._cpu_prev[1]
            frac = 1 - di / dt if dt > 0 else 0.0
            lines.append(f"CPU  {_bar(frac)} {frac * 100:5.1f}%")
        self._cpu_prev = cur
        m = _mem()
        if m:
            lines.append(f"RAM  {_bar(m[0] / m[1])} {m[0]:.1f}/{m[1]:.0f} GB")
        g = _gpu()
        if g:
            lines.append(f"GPU  {_bar(g[0] / 100)} {g[0]:5.1f}%")
            lines.append(f"VRAM {_bar(g[1] / g[2])} {g[1]:.1f}/{g[2]:.0f} GB")
        else:
            lines.append("GPU  none")
        mr = st.get("max_ranks")
        if mr:
            lines.append(f"QE   {_bar(st.get('used_ranks', 0) / mr)} {st.get('used_ranks', 0)}/{mr} ranks")
        lines.append(f"cores {os.cpu_count()} logical")
        self.query_one("#machine", Static).update("\n".join(lines))

    def _funnel(self, counts: dict) -> None:
        t = self.query_one("#funnel", DataTable)
        t.clear()
        for s in STAGES:
            c = counts.get(s, {})
            p = sum(v for k, v in c.items() if k in PASSISH)
            o = sum(v for k, v in c.items() if k in OUTISH)
            d = sum(v for k, v in c.items() if k == "DEFERRED")
            e = sum(v for k, v in c.items() if k == "ERROR")
            rest = sum(c.values()) - p - o - d - e
            t.add_row(s, str(p), str(o), str(d), str(e), str(rest))

    def _jobs(self, jobs: list) -> None:
        t = self.query_one("#jobs", DataTable)
        t.clear()
        for j in jobs:
            est = j.get("est_h")
            prog = (_bar(j["elapsed_h"] / est) + f" {100 * j['elapsed_h'] / est:3.0f}%") if est else "-"
            t.add_row(j["key"], str(j["ranks"]), f"{j['elapsed_h']:.2f}", f"{est:.1f}" if est else "-", prog,
                      str(j["pid"]))

    def _top(self, rows: list) -> None:
        t = self.query_one("#top_cands", DataTable)
        t.clear()
        for r in rows:
            t.add_row(r["cid"], r["formula"], r["lane"], r["stage"], str(r["outcome"]), str(r.get("value") or "")[:24])

    def _log(self, jobs: list) -> None:
        logs = [("events", self.rd / "events.log")] + [(j["key"], Path(j["log"])) for j in jobs]
        name, path = logs[self._log_choice % len(logs)]
        lines = _tail(path, 200)
        sig = (name, len(lines), lines[-1] if lines else "")
        w = self.query_one("#log", RichLog)
        w.border_title = f"log: {name}  (l to switch)"
        if sig != self._log_sig:
            self._log_sig = sig
            w.clear()
            for ln in lines:
                w.write(ln)

    # ------------------------------------------------------------------ actions
    def action_cycle_log(self) -> None:
        self._log_choice += 1
        self._log_sig = None
        self.refresh_all()

    def action_start(self) -> None:
        if controller_pid(self.rd):
            self.notify("already running")
            return
        pid = start_detached(self.rd, rounds=self.rounds, machine=machine_path())
        self.notify(f"campaign started (pid {pid})" if pid else "start failed: see controller.log",
                    severity="information" if pid else "error")

    def action_drain(self) -> None:
        if not controller_pid(self.rd):
            self.notify("not running")
            return
        request(self.rd, "drain")
        self.notify("draining: running jobs finish, nothing new starts")

    def action_stop(self) -> None:
        if not controller_pid(self.rd):
            self.notify("not running")
            return

        def done(choice: str | None) -> None:
            if choice == "stop":
                request(self.rd, "stop")
                self.notify("stopping: every job is being terminated")

        self.push_screen(Confirm("Stop now? Every running job is terminated; its partial results are "
                                 "discarded and it is retried on the next start.",
                                 [("stop", "Stop now", "error"), ("cancel", "Cancel", "default")]), done)

    def action_quit_ui(self) -> None:
        if not controller_pid(self.rd):
            self.exit()
            return

        def done(choice: str | None) -> None:
            if choice == "leave":
                self.exit()
            elif choice == "drain":
                request(self.rd, "drain")
                self.exit()
            elif choice == "stop":
                request(self.rd, "stop")
                self.exit()

        self.push_screen(Confirm("The campaign is running. Close the dashboard and…",
                                 [("leave", "Keep running", "primary"), ("drain", "Drain", "warning"),
                                  ("stop", "Stop now", "error"), ("cancel", "Cancel", "default")]), done)


def run(rd: Path, rounds: int = 1) -> int:
    rd.mkdir(parents=True, exist_ok=True)
    Dashboard(rd, rounds=rounds).run()
    return 0
