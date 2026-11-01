"""Install home and run directory layout.

``RTT_HOME`` (default ``~/.rtt``) holds everything ``install.sh`` sets up::

    $RTT_HOME/env.sh            environment for every rtt process (written by install.sh)
    $RTT_HOME/machine.toml      resources of this machine (``rtt init``)
    $RTT_HOME/runs/<run_id>/    one campaign: ledger.sqlite, work/, controller files

Controller files in a run directory::

    controller.pid      pid of the running controller (its process group = its children)
    control             one word written by the UI/CLI: ``drain`` (finish running jobs, start no
                        new ones) or ``stop`` (end every job now)
    status.json         snapshot written by the controller every few seconds
    events.log          one line per event (job start/finish, phase changes, errors)
    jobs/<id>.log       stdout/stderr of each child job
"""
from __future__ import annotations

import json
import os
from pathlib import Path


def rtt_home() -> Path:
    return Path(os.path.expanduser(os.environ.get("RTT_HOME", "~/.rtt")))


def runs_dir() -> Path:
    return rtt_home() / "runs"


def run_dir(run_id: str | None = None) -> Path:
    rid = run_id or os.environ.get("RTT_RUN") or "default"
    return runs_dir() / rid


def machine_path() -> Path:
    return rtt_home() / "machine.toml"


def read_json(p: Path, default=None):
    try:
        return json.loads(p.read_text(encoding="utf-8"))
    except (OSError, ValueError):
        return default


def write_json_atomic(p: Path, obj) -> None:
    tmp = p.with_suffix(p.suffix + ".tmp")
    tmp.write_text(json.dumps(obj, indent=1, sort_keys=True, default=str), encoding="utf-8")
    os.replace(tmp, p)


def pid_alive(pid: int | None) -> bool:
    if not pid:
        return False
    try:
        os.kill(int(pid), 0)
    except ProcessLookupError:
        return False
    except PermissionError:
        return True
    return True


def controller_pid(rd: Path) -> int | None:
    try:
        pid = int((rd / "controller.pid").read_text().strip())
    except (OSError, ValueError):
        return None
    return pid if pid_alive(pid) else None
