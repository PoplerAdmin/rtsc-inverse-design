"""Portable runner (src/rtt/campaign): sizing, machine.toml, controller start/drain/stop.

The controller tests replace the scheduler with a fixed job list and the child command with a
short ``sleep``, so no QE, GPU or MatterGen is needed."""
from __future__ import annotations

import os
import sys
import threading
import time
from pathlib import Path

import pytest

from rtt.campaign import hw
from rtt.campaign import controller as C
from rtt.campaign.home import controller_pid, read_json
from rtt.config import load_machine
from rtt.engine.scheduler import Job

pytestmark = pytest.mark.skipif(os.name != "posix", reason="process groups need POSIX")


def _hw(cores=8, ram=15.0, gpu=24.0):
    return hw.Hardware(logical_cpus=2 * cores, physical_cores=cores, ram_GB=ram,
                       gpus=[{"name": "GPU", "mem_GB": gpu}] if gpu else [], disk_free_GB=500.0)


def test_sizing_reference_machine():
    s = hw.sizing(_hw())
    assert s == {"max_ranks": 7, "s_ranks": 4, "s_jobs": 1, "f_ranks": 7, "gpu": True, "mattergen_batch": 16}


def test_sizing_scales_and_respects_ram():
    big = hw.sizing(_hw(cores=32, ram=128, gpu=80))
    assert big["max_ranks"] == 31 and big["s_jobs"] == 7 and big["f_ranks"] == 8 and big["mattergen_batch"] == 32
    lowram = hw.sizing(_hw(cores=32, ram=6, gpu=0))
    assert lowram["max_ranks"] == 4 and lowram["gpu"] is False


def test_minimum_report_flags_small_machines():
    assert hw.minimum_report(_hw()) == []
    w = " ".join(hw.minimum_report(_hw(cores=4, ram=8, gpu=0)))
    assert "physical cores" in w and "RAM" in w and "no NVIDIA GPU" in w


def test_machine_toml_loads(tmp_path):
    p = tmp_path / "machine.toml"
    p.write_text(hw.machine_toml(_hw(), scratch=str(tmp_path / "s")))
    m = load_machine(p)
    assert m.get("budget.max_ranks") == 7 and m.get("lanes.S.mpi_ranks") == 4 and m.get("lanes.G.gpu") is True


def _controller(tmp_path, jobs, sleep_s=1.0, max_ranks=8):
    mp = tmp_path / "machine.toml"
    txt = hw.machine_toml(_hw(cores=max_ranks + 1), scratch=str(tmp_path / "s"))
    mp.write_text(txt)
    c = C.Controller(rd=tmp_path / "run", machine=mp, rounds=1, front=False)
    c.jobs_queue = list(jobs)

    def plan_jobs(L, tiers):
        out = []
        for j in c.jobs_queue:
            key = f"{j.stage}:{j.candidate_id}"
            if key in c.children:
                continue
            out.append((j, key, c.s_ranks, 0.01))
        return out

    def spawn(key, kind, argv, ranks, **kw):
        c.jobs_queue = [j for j in c.jobs_queue if f"{j.stage}:{j.candidate_id}" != key]
        return C.Controller.spawn(c, key, kind, [sys.executable, "-c", f"import time; time.sleep({sleep_s})"],
                                  ranks, **kw)

    c.plan_jobs = plan_jobs
    c.spawn = spawn
    return c


def _jobs(n):
    return [Job(f"c-{i:03d}", "S5", "S", 1.0, 4.0) for i in range(n)]


def test_controller_runs_all_jobs_within_rank_budget(tmp_path):
    c = _controller(tmp_path, _jobs(3), sleep_s=0.5)
    peak = []
    orig = C.Controller.write_status

    def ws(self=c):
        peak.append(self.used_ranks())
        orig(self)
    c.write_status = ws
    assert c.run() == 0
    st = read_json(c.rd / "status.json")
    assert st["state"] == "finished" and len(st["finished"]) == 3
    assert max(peak) <= c.max_ranks
    assert controller_pid(c.rd) is None


def test_controller_stop_kills_children(tmp_path):
    c = _controller(tmp_path, _jobs(2), sleep_s=60)
    t = threading.Thread(target=c.run)
    t.start()
    for _ in range(100):
        if c.children:
            break
        time.sleep(0.1)
    pids = [ch.proc.pid for ch in c.children.values()]
    assert pids
    C.request(c.rd, "stop")
    t.join(timeout=60)
    assert not t.is_alive()
    st = read_json(c.rd / "status.json")
    assert st["state"] == "stopped" and st["jobs"] == []
    for p in pids:
        with pytest.raises(ProcessLookupError):
            os.kill(p, 0)


def test_controller_drain_finishes_running_and_starts_nothing(tmp_path):
    c = _controller(tmp_path, _jobs(4), sleep_s=1.5)
    t = threading.Thread(target=c.run)
    t.start()
    for _ in range(100):
        if c.children:
            break
        time.sleep(0.1)
    C.request(c.rd, "drain")
    t.join(timeout=60)
    st = read_json(c.rd / "status.json")
    assert st["state"] == "stopped"
    assert 0 < len(st["finished"]) < 4 and all(f["rc"] == 0 for f in st["finished"])


def test_compute_limited_jobs_are_reported_not_started(tmp_path, monkeypatch):
    mp = tmp_path / "machine.toml"
    mp.write_text(hw.machine_toml(_hw(), scratch=str(tmp_path / "s"), max_job_wall_h=1.0))
    c = C.Controller(rd=tmp_path / "run", machine=mp, rounds=1, front=False)
    import rtt.engine.scheduler as S
    monkeypatch.setattr(S, "next_jobs", lambda *a, **k: [Job("c-big", "S6", "F", 1.0, 1e4)])
    out = c.plan_jobs(None, None)
    assert out == [] and "S6:c-big" in c.compute_limited


def test_cli_init_and_campaign_status(tmp_path, monkeypatch, capsys):
    from rtt import cli
    monkeypatch.setenv("RTT_HOME", str(tmp_path))
    assert cli.main(["init"]) == 0
    assert (tmp_path / "machine.toml").exists()
    load_machine(tmp_path / "machine.toml")
    assert cli.main(["campaign", "status", "--run-id", "x"]) == 0
    assert "not running" in capsys.readouterr().out
    assert cli.main(["campaign", "stop", "--run-id", "x"]) == 0


def test_dashboard_renders_status_and_quits(tmp_path, monkeypatch):
    textual = pytest.importorskip("textual")
    import asyncio
    from rtt.campaign.home import write_json_atomic
    from rtt.campaign.ui import Dashboard
    monkeypatch.setenv("RTT_HOME", str(tmp_path))
    rd = tmp_path / "runs" / "t"
    (rd / "jobs").mkdir(parents=True)
    (rd / "events.log").write_text("2026-11-01T10:00:00 round 1\n")
    write_json_atomic(rd / "status.json", {
        "state": "running", "phase": "back", "round": 1, "rounds": 1, "max_ranks": 7, "used_ranks": 4,
        "jobs": [{"key": "S5:c-1", "kind": "S5", "cid": "c-1", "ranks": 4, "pid": 1, "elapsed_h": 1.0,
                  "est_h": 4.0, "log": str(rd / "jobs" / "S5_c-1.log"), "gpu": False}],
        "counts": {"outcomes": {"S1": {"PASS": 10, "SCREENED_OUT": 2}, "S5": {"SCREENED_OUT": 1}}},
        "top": [{"cid": "c-1", "formula": "H6IrMg2", "lane": "C", "stage": "S5", "outcome": "PASS",
                 "value": "1639"}],
        "compute_limited": [{"key": "S6:c-2", "est_wall_h": 160.0, "ranks": 7}], "max_job_wall_h": 48})

    async def go():
        app = Dashboard(rd)
        async with app.run_test(size=(140, 50)) as pilot:
            await pilot.pause()
            assert app.query_one("#funnel").row_count == 7
            assert app.query_one("#jobs").row_count == 1
            assert app.query_one("#top_cands").row_count == 1
            await pilot.press("l")
            await pilot.press("q")       # controller not running: quits without asking
            await pilot.pause()
        assert app.return_code in (0, None)
    asyncio.run(go())


def test_stop_kills_grandchildren_in_their_own_session(tmp_path):
    """The core QE launcher starts mpirun with start_new_session=True: a process-group kill
    of the job alone leaves the MPI ranks running. Stop must end the whole tree."""
    gpid = tmp_path / "grandchild.pid"
    code = ("import subprocess, sys, time; "
            "p = subprocess.Popen([sys.executable, '-c', 'import time; time.sleep(120)'], start_new_session=True); "
            f"open({str(gpid)!r}, 'w').write(str(p.pid)); time.sleep(120)")
    c = _controller(tmp_path, _jobs(1))
    c.spawn = lambda key, kind, argv, ranks, **kw: (
        setattr(c, "jobs_queue", []),
        C.Controller.spawn(c, key, kind, [sys.executable, "-c", code], ranks, **kw))
    t = threading.Thread(target=c.run)
    t.start()
    for _ in range(200):
        if gpid.exists() and gpid.read_text():
            break
        time.sleep(0.1)
    g = int(gpid.read_text())
    assert C._alive(g)
    C.request(c.rd, "stop")
    t.join(timeout=90)
    assert not t.is_alive()
    time.sleep(0.5)
    assert not C._alive(g)
