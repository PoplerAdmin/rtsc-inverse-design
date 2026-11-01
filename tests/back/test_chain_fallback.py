"""Campaign 1: one sticky robust-solver retry for known pw.x failures (La nscf Davidson, YH3 BFGS)."""
from __future__ import annotations

import json
from pathlib import Path

from rtt.dft import workflow as W


IN = """ &control
    calculation = 'vc-relax'
 /
 &electrons
    conv_thr = 1e-10
 /
 &ions
    ion_dynamics = 'bfgs'
 /
 &cell
    cell_dynamics = 'bfgs'
 /
"""


def test_fallback_text():
    t = W.pw_fallback_text(IN)
    assert "diagonalization = 'cg'" in t and "'damp'" in t and "'damp-w'" in t and "bfgs" not in t
    assert "conv_thr = 1e-10" in t


class _R:
    def __init__(self, ok, err=""):
        self.ok, self.error_message, self.returncode = ok, err, 0 if ok else 1

    def to_dict(self):
        return {"ok": self.ok}


class _L:
    nprocs, nice = 1, 19

    def command(self, code, inp):
        return ["true"]

    def exe(self, code):
        return "true"

    def environment(self):
        return {}


def test_chain_retries_once_with_marker(tmp_path, monkeypatch):
    import rtt.core.dfpt.launcher as LA
    calls = []

    def fake_run(name, cmd, w, env=None, timeout_s=None):
        txt = (Path(w) / f"{name}.in").read_text()
        calls.append(txt)
        if "'cg'" in txt:
            (Path(w) / f"{name}.out").write_text("JOB DONE")
            return _R(True)
        (Path(w) / f"{name}.out").write_text("Error in routine bfgs (1):\n dE0s is positive which should never happen")
        return _R(False, "QE: Error in routine bfgs (1): dE0s is positive")

    monkeypatch.setattr(LA, "run_command", fake_run)
    chk = lambda w: ((w / "relax.out").exists() and "JOB DONE" in (w / "relax.out").read_text(), "x")  # noqa: E731
    st = W.ChainStep("relax", "pw.x", lambda: IN, chk)
    ch = W.Chain(tmp_path, [st], _L(), material_id="m", profile="s4", policy=W.RunPolicy(max_machine_ranks=None))
    ok, msg = ch.run()
    assert ok, msg
    assert len(calls) == 2 and "bfgs" in calls[0] and "'damp'" in calls[1]
    assert json.loads((tmp_path / "relax.fallback").read_text())["error"] == "dE0s is positive"
    # second run adopts (sticky marker -> same token), no QE call
    ch2 = W.Chain(tmp_path, [st], _L(), material_id="m", profile="s4", policy=W.RunPolicy(max_machine_ranks=None))
    assert ch2.run()[0] and len(calls) == 2


def test_unknown_error_not_retried(tmp_path, monkeypatch):
    import rtt.core.dfpt.launcher as LA
    n = []

    def fake_run(name, cmd, w, env=None, timeout_s=None):
        n.append(1)
        return _R(False, "segfault")

    monkeypatch.setattr(LA, "run_command", fake_run)
    st = W.ChainStep("scf", "pw.x", lambda: IN, lambda w: (False, "no"))
    ch = W.Chain(tmp_path, [st], _L(), material_id="m", profile="s4", policy=W.RunPolicy(max_machine_ranks=None))
    assert not ch.run()[0] and len(n) == 1 and not (tmp_path / "scf.fallback").exists()


def test_stale_dyn_files_moved_aside(tmp_path):
    """Campaign 1: stale dyn files of a killed ph.x run made the per-q pruner delete live wfc."""
    for n in ("dyn0", "dyn1", "dyn2"):
        (tmp_path / n).write_text("old")
    (tmp_path / "dyna2F").write_text("keep")
    moved = W._clear_stale_ph_outputs(tmp_path)
    assert moved == ["dyn0", "dyn1", "dyn2"]
    assert not (tmp_path / "dyn0").exists() and (tmp_path / "dyna2F").exists()
    assert len(list(tmp_path.glob("stale_ph_*/dyn*"))) == 3
