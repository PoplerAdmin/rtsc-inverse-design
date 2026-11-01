"""Doc 31 §4-§5: S5 single-mesh admission failure -> DEFERRED (helper for S5Coarse), and one
recorded vc-relax continuation when pw.x hits nstep (campaign 1: NaIrH3 q 3^3, BaIrH6 S4)."""
from __future__ import annotations

import json
from pathlib import Path
from types import SimpleNamespace

from rtt.dft import workflow as W
from rtt.dft.admission import admission_deferral, mesh_status, s5_admission_outcome
from rtt.dft.coarse import CoarseMetrics, SigmaRow
from rtt.dft.collect import Collected

ADM = ("smearing admission: no smearing has non-negative per-mode alpha^2F (max negative fraction per "
       "smearing {'1': 0.0027, '2': 3e-06})")


def _run(ok=True, q=(3, 3, 3)):
    return SimpleNamespace(ok=ok, q=list(q), error=None if ok else "ph.x crashed")


def _col(errors=(), spectra=("sp",), report=None):
    return Collected(spectra=list(spectra) if not errors else [], artifacts=[], eph=None,
                     admission=report or {"rule": "mode_positive", "tol_rel": 1e-6,
                                          "max_mode_negative_fraction": {"1": 0.0027, "2": 3e-06}},
                     errors=list(errors), flags=["SMEARING_ADMISSION_FAIL"] if errors else [])


def _met(q, ok=True):
    rows = [SigmaRow(index=i, sigma_Ry=0.01 * i, admitted=ok, A_K=100.0 + i if ok else None, on_plateau=ok and i >= 3)
            for i in range(1, 6)]
    return CoarseMetrics(rows=rows, plateau_sigma_min_Ry=0.03, min_freq_matdyn_cm1=10.0, q=list(q),
                         flags=[] if ok else ["SMEARING_ADMISSION_FAIL"], admission_ok=ok)


def test_single_mesh_admission_failure_is_deferred():
    runs = [_run(q=(3, 3, 3)), _run(q=(4, 4, 4))]
    o = s5_admission_outcome([[3, 3, 3], [4, 4, 4]], runs, [_col([ADM]), _col()], tier_run_id=7)
    assert o.outcome == "DEFERRED" and o.criterion == "admission_failed" and o.tier_run_id == 7
    assert "ADMISSION_FAILED:q3x3x3" in o.tags and "SINGLE_MESH_ADMITTED" in o.tags
    assert o.value["meshes"][0]["max_mode_negative_fraction"]["1"] == 0.0027
    assert "both meshes" in o.reason


def test_real_failures_stay_error():
    runs = [_run(q=(3, 3, 3)), _run(ok=False, q=(4, 4, 4))]
    assert s5_admission_outcome(None, runs, [_col([ADM]), None]) is None              # crash on q4
    other = _col(["a2F.dos1 changed while it was parsed"])
    assert mesh_status(_run(), other) == "error"
    assert s5_admission_outcome(None, [_run(), _run()], [_col([ADM, "a2F.dos1 changed"]), _col()]) is None
    assert s5_admission_outcome(None, [_run(), _run()], [_col(), _col()]) is None       # nothing failed


def test_admission_deferral_from_stage_result_payload():
    """The S5Coarse payload carries CoarseMetrics per mesh and the joined error text."""
    runs = [_run(q=(3, 3, 3)), _run(q=(4, 4, 4))]
    res = SimpleNamespace(ok=False, tier_run_id=11, error=f"q[3, 3, 3]: {ADM}",
                          payload={"runs": runs, "metrics_by_mesh": [_met((3, 3, 3), ok=False), _met((4, 4, 4))]})
    o = admission_deferral(res)
    assert o is not None and o.outcome == "DEFERRED" and o.tier_run_id == 11
    assert o.value["meshes"][1]["A_mid_K"] is not None                                # admitted mesh kept
    res.error = f"q[3, 3, 3]: {ADM}; q[3, 3, 3]: run not done: pw.x"                    # mixed -> ERROR path
    assert admission_deferral(res) is None
    assert admission_deferral(SimpleNamespace(ok=True, payload={})) is None


# ------------------------------------------------------------------ nstep continuation
IN = """ &control
    calculation = 'vc-relax'
    nstep = 100
 /
 &system
    ibrav = 0
    nat = 2
 /
 &ions
    ion_dynamics = 'bfgs'
 /
CELL_PARAMETERS angstrom
      3.00000000000000     0.00000000000000     0.00000000000000
      0.00000000000000     3.00000000000000     0.00000000000000
      0.00000000000000     0.00000000000000     3.00000000000000
ATOMIC_POSITIONS crystal
 Mg      0.00000000000000     0.00000000000000     0.00000000000000
 H       0.50000000000000     0.50000000000000     0.50000000000000
K_POINTS automatic
 4 4 4 0 0 0
"""

OUT_NSTEP = """
CELL_PARAMETERS (alat=  5.66918640)
   1.010000000   0.000000000   0.000000000
   0.000000000   1.010000000   0.000000000
   0.000000000   0.000000000   1.020000000

ATOMIC_POSITIONS (crystal)
Mg            0.0000000000        0.0000000000        0.0000000000
H             0.5000000000        0.5000000000        0.4900000000

     The maximum number of steps has been reached.

     End of BFGS Geometry Optimization
JOB DONE.
"""


def test_pw_continue_text():
    t, info = W.pw_continue_text(IN, OUT_NSTEP)
    assert "nstep = 200" in t and info["nstep_old"] == 100 and info["nstep_new"] == 200
    assert "0.4900000000" in t and "ATOMIC_POSITIONS (crystal)" in t
    assert info["geometry"] == "last printed cell and positions"
    a = 5.66918640 * W.BOHR_A * 1.01
    assert f"{a:20.14f}".strip() in t
    assert "K_POINTS automatic" in t and "ion_dynamics = 'bfgs'" in t                   # physics unchanged


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


def test_chain_nstep_retry_once_recorded(tmp_path, monkeypatch):
    import rtt.core.dfpt.launcher as LA
    calls = []

    def fake_run(name, cmd, w, env=None, timeout_s=None):
        txt = (Path(w) / f"{name}.in").read_text()
        calls.append(txt)
        if "nstep = 200" in txt:
            (Path(w) / f"{name}.out").write_text("Begin final coordinates\nEnd final coordinates\nJOB DONE")
        else:
            (Path(w) / f"{name}.out").write_text(OUT_NSTEP)
        return _R(True)

    monkeypatch.setattr(LA, "run_command", fake_run)
    chk = lambda w: ((w / "relax.out").exists() and "End final coordinates" in (w / "relax.out").read_text(),  # noqa: E731
                     "no final coords")
    st = W.ChainStep("relax", "pw.x", lambda: IN, chk)
    ch = W.Chain(tmp_path, [st], _L(), material_id="m", profile="s4", policy=W.RunPolicy(max_machine_ranks=None))
    ok, msg = ch.run()
    assert ok, msg
    assert len(calls) == 2 and "nstep = 100" in calls[0] and "0.4900000000" in calls[1]
    rec = json.loads((tmp_path / "relax.nstep_retry").read_text())
    assert rec["nstep_new"] == 200 and rec["first_output"] == "relax.nstep1.out"
    assert (tmp_path / "relax.nstep1.out").read_text() == OUT_NSTEP and ch.retries[0]["kind"] == "nstep_retry"
    # sticky: a second chain adopts the finished continuation without running QE
    ch2 = W.Chain(tmp_path, [st], _L(), material_id="m", profile="s4", policy=W.RunPolicy(max_machine_ranks=None))
    assert ch2.run()[0] and len(calls) == 2


def test_chain_nstep_retry_only_once(tmp_path, monkeypatch):
    import rtt.core.dfpt.launcher as LA
    calls = []

    def fake_run(name, cmd, w, env=None, timeout_s=None):
        calls.append(1)
        (Path(w) / f"{name}.out").write_text(OUT_NSTEP)
        return _R(True)

    monkeypatch.setattr(LA, "run_command", fake_run)
    chk = lambda w: ((w / "relax.out").exists() and "End final coordinates" in (w / "relax.out").read_text(),  # noqa: E731
                     "no final coords")
    st = W.ChainStep("relax", "pw.x", lambda: IN, chk)
    ch = W.Chain(tmp_path, [st], _L(), material_id="m", profile="s4", policy=W.RunPolicy(max_machine_ranks=None))
    assert not ch.run()[0] and len(calls) == 2


def test_s5coarse_decide_uses_admission_helper():
    """Wired in S5Coarse.decide (doc 31 §4): NaIrH3-like q3 admission failure is DEFERRED, not ERROR."""
    from rtt.config import load_targets, load_tiers
    from rtt.engine.stages import S5Coarse
    from rtt.records import StageResult
    runs = [_run(q=(3, 3, 3)), _run(q=(4, 4, 4))]
    res = StageResult("S5", "c", False, {}, 5, error=f"q[3, 3, 3]: {ADM}",
                      payload={"runs": runs, "collected_by_mesh": [_col([ADM]), _col()],
                               "metrics_by_mesh": [_met((3, 3, 3), ok=False), _met((4, 4, 4))]})
    o = S5Coarse().decide(res, load_tiers(), load_targets())
    assert o.outcome == "DEFERRED" and o.criterion == "admission_failed"
    res.payload["runs"][1] = _run(ok=False, q=(4, 4, 4))
    assert S5Coarse().decide(res, load_tiers(), load_targets()).outcome == "ERROR"
