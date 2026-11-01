"""CLI run-context options (campaign 1): --ledger / --workdir / --ranks keep big files and the
live SQLite ledger outside the repo run directory."""
from __future__ import annotations

import argparse
from pathlib import Path


def _args(tmp_path: Path, **kw) -> argparse.Namespace:
    d = dict(run=str(tmp_path / "runs" / "c1"), targets=None, tiers=None, machine=None,
             ledger=None, workdir=None, ranks=None)
    d.update(kw)
    return argparse.Namespace(**d)


def test_defaults_inside_run_dir(tmp_path, monkeypatch):
    from rtt.cli import _ctx
    monkeypatch.delenv("RTT_LEDGER", raising=False)
    monkeypatch.delenv("RTT_WORKDIR", raising=False)
    (tmp_path / "runs" / "c1").mkdir(parents=True)
    ctx = _ctx(_args(tmp_path))
    assert Path(ctx.ledger.path) == tmp_path / "runs" / "c1" / "ledger.sqlite"
    assert ctx.workdir == tmp_path / "runs" / "c1" / "work"
    assert ctx.launcher is None


def test_ledger_workdir_ranks_options(tmp_path, monkeypatch):
    from rtt.cli import _ctx
    monkeypatch.delenv("RTT_LEDGER", raising=False)
    (tmp_path / "runs" / "c1").mkdir(parents=True)
    ctx = _ctx(_args(tmp_path, ledger=str(tmp_path / "live" / "L.sqlite"), workdir=str(tmp_path / "scratch"),
                     ranks=3))
    assert Path(ctx.ledger.path) == tmp_path / "live" / "L.sqlite"
    assert ctx.ledger.meta("run_id") == "c1"
    assert ctx.workdir == tmp_path / "scratch"
    assert ctx.launcher.nprocs == 3


def test_env_fallbacks(tmp_path, monkeypatch):
    from rtt.cli import _ctx
    (tmp_path / "runs" / "c1").mkdir(parents=True)
    monkeypatch.setenv("RTT_LEDGER", str(tmp_path / "envL.sqlite"))
    monkeypatch.setenv("RTT_WORKDIR", str(tmp_path / "envW"))
    ctx = _ctx(_args(tmp_path))
    assert Path(ctx.ledger.path) == tmp_path / "envL.sqlite"
    assert ctx.workdir == tmp_path / "envW"
