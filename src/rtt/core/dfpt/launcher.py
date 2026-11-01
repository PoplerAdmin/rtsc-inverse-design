"""Launcher abstraction and a single-step runner with timeout and error capture.

Command forms produced (input passed with QE's ``-in`` flag, handled by ``input_from_file`` in
pw.x, ph.x, q2r.x and matdyn.x; stdout to ``<step>.out``, stderr to ``<step>.err``)::

    plain : [nice -n N] pw.x [-nk K] -in step.in
    mpirun: [nice -n N] mpirun -np N [extra] pw.x [-nk K] -in step.in
    srun  : srun -n N [extra] pw.x [-nk K] -in step.in

``-nk`` (pools) is only passed to pw.x, ph.x and epw.x (EPW requires one process per pool:
set npool = nprocs for epw.x, otherwise elphel2_shuffle stops with "only one proc per pool"). q2r.x / matdyn.x run serially by default
(``serial_codes``), they are cheap.
"""
from __future__ import annotations

import os
import re
import shlex
import signal
import subprocess
import time
from dataclasses import asdict, dataclass, field
from pathlib import Path
from typing import Any

POOL_CODES = frozenset({"pw.x", "ph.x", "epw.x"})  # epw.x needs npool == nprocs


@dataclass
class Launcher:
    kind: str = "plain"                     # "plain" | "mpirun" | "srun"
    nprocs: int = 1
    npool: int | None = None
    mpirun: str = "mpirun"
    srun: str = "srun"
    mpi_args: list[str] = field(default_factory=list)   # e.g. ["--allow-run-as-root", "--bind-to", "none"]
    bin_dir: str | None = None
    nice: int | None = None                 # prepend `nice -n <nice>` (not with srun)
    env: dict[str, str] = field(default_factory=dict)
    serial_codes: tuple[str, ...] = ("q2r.x", "matdyn.x")
    omp_threads: int | None = 1

    def __post_init__(self) -> None:
        if self.kind not in ("plain", "mpirun", "srun"):
            raise ValueError(f"unknown launcher kind {self.kind!r}")
        if self.nprocs < 1:
            raise ValueError("nprocs must be >= 1")

    @classmethod
    def from_config(cls, cfg: dict[str, Any] | None) -> "Launcher":
        cfg = dict(cfg or {})
        if isinstance(cfg.get("mpi_args"), str):
            cfg["mpi_args"] = shlex.split(cfg["mpi_args"])
        if "serial_codes" in cfg:
            cfg["serial_codes"] = tuple(cfg["serial_codes"])
        return cls(**cfg)

    def exe(self, code: str) -> str:
        return str(Path(self.bin_dir) / code) if self.bin_dir else code

    def command(self, code: str, input_file: str) -> list[str]:
        serial = code in self.serial_codes or self.kind == "plain"
        prog = [self.exe(code)]
        if code in POOL_CODES and self.npool and self.npool > 1 and not serial:
            prog += ["-nk", str(int(self.npool))]
        prog += ["-in", input_file]
        if serial:
            cmd = prog
        elif self.kind == "mpirun":
            cmd = [self.mpirun, "-np", str(int(self.nprocs)), *self.mpi_args, *prog]
        else:
            cmd = [self.srun, "-n", str(int(self.nprocs)), *self.mpi_args, *prog]
        if self.nice is not None and self.kind != "srun":
            cmd = ["nice", "-n", str(int(self.nice)), *cmd]
        return cmd

    def environment(self) -> dict[str, str]:
        e = dict(os.environ)
        if self.omp_threads is not None:
            e["OMP_NUM_THREADS"] = str(int(self.omp_threads))
        e.update({k: str(v) for k, v in self.env.items()})
        return e

    def to_dict(self) -> dict[str, Any]:
        return asdict(self)


@dataclass
class RunRecord:
    step: str
    command: list[str]
    cwd: str
    returncode: int | None
    elapsed_s: float
    timed_out: bool
    job_done: bool
    stdout_file: str
    stderr_file: str
    error_message: str = ""          # QE '%%%%' error block, CRASH file, or stderr tail
    started_at: float = 0.0

    @property
    def ok(self) -> bool:
        return (not self.timed_out) and self.returncode == 0 and self.job_done

    def to_dict(self) -> dict[str, Any]:
        d = asdict(self)
        d["ok"] = self.ok
        return d


_QE_ERR = re.compile(r"%{20,}\s*\n(.*?)\n\s*%{20,}", re.S)


def extract_qe_error(text: str) -> str:
    """Return the last QE ``%%%%`` error block (``Error in routine ...``), if any."""
    blocks = _QE_ERR.findall(text)
    return " | ".join(" ".join(b.split()) for b in blocks[-1:]) if blocks else ""


def tail(path: Path, n: int = 30) -> str:
    try:
        with open(path, "r", errors="replace") as fh:
            lines = fh.readlines()
        return "".join(lines[-n:])
    except OSError:
        return ""


def run_command(step: str, cmd: list[str], cwd: Path, *, env: dict[str, str] | None,
                timeout_s: float | None, out_name: str | None = None) -> RunRecord:
    """Run one QE step; never raises on QE failure (the record carries the diagnosis)."""
    cwd = Path(cwd)
    out = cwd / (out_name or f"{step}.out")
    err = cwd / f"{step}.err"
    crash = cwd / "CRASH"
    if crash.exists():
        crash.unlink()
    t0 = time.time()
    timed_out = False
    with open(out, "w") as fo, open(err, "w") as fe:
        try:
            proc = subprocess.Popen(cmd, cwd=str(cwd), stdout=fo, stderr=fe,
                                    stdin=subprocess.DEVNULL, env=env, start_new_session=True)
        except OSError as exc:  # binary not found etc.
            fe.write(f"failed to start {cmd!r}: {exc}\n")
            return RunRecord(step, cmd, str(cwd), None, 0.0, False, False, str(out), str(err),
                             f"failed to start: {exc}", t0)
        try:
            rc = proc.wait(timeout=timeout_s)
        except subprocess.TimeoutExpired:
            timed_out = True
            try:
                os.killpg(proc.pid, signal.SIGTERM)
                proc.wait(timeout=20)
            except Exception:
                try:
                    os.killpg(proc.pid, signal.SIGKILL)
                except Exception:
                    pass
                proc.wait()
            rc = proc.returncode
    elapsed = time.time() - t0
    text = out.read_text(errors="replace") if out.exists() else ""
    job_done = "JOB DONE" in text
    msg = ""
    if timed_out:
        msg = f"timeout after {timeout_s} s"
    if not job_done or rc != 0:
        parts = [msg] if msg else []
        qe = extract_qe_error(text)
        if qe:
            parts.append("QE: " + qe)
        if crash.exists():
            parts.append("CRASH: " + " ".join(crash.read_text(errors="replace").split())[:2000])
        et = tail(err, 15).strip()
        if et:
            parts.append("stderr: " + et[-2000:])
        if not parts:
            parts.append(f"exit code {rc}; no 'JOB DONE' in output; last lines: "
                         + tail(out, 8).strip()[-1500:])
        msg = " || ".join(parts)
    return RunRecord(step, cmd, str(cwd), rc, elapsed, timed_out, job_done, str(out), str(err),
                     msg, t0)
