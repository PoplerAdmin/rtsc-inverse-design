"""Doc 31 §6-§7: S5 cost refit, `rtt step` stage gate (R30 F5), stage_outcomes migration."""
from __future__ import annotations

import argparse
import json
import sqlite3

import pytest

from rtt.config import load_machine, load_targets, load_tiers
from rtt.engine.scheduler import S5_COST_FIT, estimate_s5_core_h, s5_ecut_Ry
from rtt.ledger.store import Ledger
from rtt.records import Origin, Outcome, StructureData

#: doc 30 §5.2/§5.3 measured S5 core-h (q 3^3 + q 4^3, 60 Ry)
MEASURED = {"ScH2": (3, 1.67), "LaH2": (3, 1.35), "NaIrH3": (5, 3.30), "BaPtH3": (5, 6.47), "Mg2IrH6": (9, 29.7)}


def test_s5_refit_bounds_the_measurements():
    for name, (n, ch) in MEASURED.items():
        est = estimate_s5_core_h(n, {"H"}, ecut_Ry=60.0)
        assert ch <= est <= 2.5 * ch, name                  # upper envelope, within x2.5
    assert estimate_s5_core_h(1, {"Al"}, ecut_Ry=60.0) == S5_COST_FIT["floor_core_h"]
    # the old C1 value 0.525 was ~10x low for 5-atom 5d hydrides
    assert estimate_s5_core_h(5, {"H", "Ir", "Na"}, ecut_Ry=60.0) > 10 * 0.525


def test_s5_cutoff_scaling():
    t = load_tiers()
    assert s5_ecut_Ry({"H", "Ir"}, t) == 80.0 and s5_ecut_Ry({"Al"}, t) == 60.0
    r = estimate_s5_core_h(5, {"H"}, tiers=t) / estimate_s5_core_h(5, {"H"}, ecut_Ry=60.0)
    assert r == pytest.approx((80 / 60) ** 1.5)


def _ctx(tmp_path):
    from rtt.engine.stage import Context
    L = Ledger(tmp_path / "L.sqlite")
    t, ti, m = load_targets(), load_tiers(), load_machine()
    L.register_configs(t, ti, m)
    return Context(ledger=L, targets=t, tiers=ti, machine=m, workdir=tmp_path / "w")


def test_step_gate_refuses_and_records_override(tmp_path, monkeypatch):
    from rtt import cli
    from rtt.engine.stages import write_outcome
    ctx = _ctx(tmp_path)
    st = StructureData(lattice=[[3, 0, 0], [0, 3, 0], [0, 0, 3]], species=["Na", "Ir", "H"],
                       frac_coords=[[0, 0, 0], [.5, .5, .5], [.5, .5, 0]])
    ctx.ledger.insert_candidate("x", "a" * 64, st, Origin(lane="G1"), cost_class="C1")
    for s in ("S1", "S2", "S3", "S4"):
        write_outcome(ctx, "x", Outcome(s, "PASS"))
    write_outcome(ctx, "x", Outcome("S5", "ERROR", reason="admission"))
    monkeypatch.setattr(cli, "_ctx", lambda args: ctx)
    ran = []
    import rtt.engine.runner as R
    monkeypatch.setattr(R, "run_job", lambda job, c: ran.append(job.stage) or Outcome("S6", "ERROR", reason="fake"))
    args = argparse.Namespace(stage="S6", cids=["x"], override_reason=None, systems=None, controls=False)
    assert cli.cmd_step(args) == 2 and ran == []
    args.override_reason = "first real S6 chain test (doc 30 §5.6)"
    assert cli.cmd_step(args) == 0 and ran == ["S6"]
    rows = ctx.ledger.qa("SELECT value FROM meta WHERE key LIKE 'gate_override:S6:x:%'")
    rec = json.loads(rows[0][0])
    assert rec["unmet"] == {"S5": "ERROR"} and "doc 30" in rec["reason"]
    ctx.ledger.close()


def test_old_ledger_migrates_s1_deferred(tmp_path):
    """A ledger created before doc 31 (36 fixed pairs, guard at 36) gains (S1, DEFERRED) once."""
    p = tmp_path / "old.sqlite"
    L = Ledger(p)
    L.close()
    c = sqlite3.connect(p)
    c.executescript("DROP TRIGGER stage_outcomes_no_insert; DROP TRIGGER stage_outcomes_no_delete;"
                    "DELETE FROM stage_outcomes WHERE stage='S1' AND outcome='DEFERRED';"
                    "CREATE TRIGGER stage_outcomes_no_insert BEFORE INSERT ON stage_outcomes "
                    "WHEN (SELECT COUNT(*) FROM stage_outcomes) >= 36 BEGIN SELECT RAISE(ABORT, 'fixed'); END;"
                    "DELETE FROM meta WHERE key LIKE 'migration:%';")
    c.commit()
    c.close()
    L = Ledger(p)
    assert L.q1("SELECT COUNT(*) FROM stage_outcomes")[0] == 37
    assert L.q1("SELECT COUNT(*) FROM meta WHERE key LIKE 'migration:stage_outcomes:%'")[0] == 1
    with pytest.raises(sqlite3.DatabaseError):
        L.x("INSERT INTO stage_outcomes VALUES ('S9','PASS')")              # guard back in place
    L.close()
    L = Ledger(p)                                                           # idempotent
    assert L.q1("SELECT COUNT(*) FROM stage_outcomes")[0] == 37
    L.close()


def test_rank_pools_append_only(tmp_path):
    L = Ledger(tmp_path / "r.sqlite")
    L.insert_rank_pool({"pool_id": "p", "cut": 1, "free_slots": 1, "ranking": [], "ranking_sha256": "0" * 64})
    with pytest.raises(sqlite3.DatabaseError):
        L.x("UPDATE rank_pools SET cut = 2")
    assert L.rank_pools()[0]["pool_id"] == "p"
    L.close()
