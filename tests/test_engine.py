"""Engine glue, scheduler, scoring/ranking, config loader, import-graph contract (R1b/R1c)."""
from __future__ import annotations

import ast
import json
from pathlib import Path

import pytest

from rtt.config import ConfigError, load_machine, load_targets, load_tiers, validate_sources
from rtt.engine.scheduler import next_jobs, next_stage, rebuild_queue  # noqa: F401
from rtt.engine.stage import Context, Stage
from rtt.engine.stages import STAGES, combine, record_dfpt_decision, store_verdict, write_outcome
from rtt.ledger.artifacts import register
from rtt.ledger.store import Ledger
from rtt.records import Origin, Outcome, Prov, StructureData
from rtt.score.rank import RankRow, final_rank
from rtt.score.s3 import score_s3

ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "src" / "rtt"
DATA = Path(__file__).parent / "data"


@pytest.fixture
def ctx(tmp_path):
    L = Ledger(tmp_path / "ledger.sqlite")
    t, ti, m = load_targets(), load_tiers(), load_machine()
    L.register_configs(t, ti, m)
    yield Context(ledger=L, targets=t, tiers=ti, machine=m, workdir=tmp_path / "work")
    L.close()


def add_al(ctx, cid="al", lane="C", T_s=None):
    st = StructureData(lattice=[[0, 2.0, 2.0], [2.0, 0, 2.0], [2.0, 2.0, 0]], species=["Al"], frac_coords=[[0, 0, 0]])
    ctx.ledger.insert_candidate(cid, cid * 10 + "0" * (64 - 10 * len(cid)), st,
                                Origin(lane=lane, control="Al" if lane == "C" else None, T_star_surrogate_K=T_s),
                                cost_class="C1")
    return ctx.ledger.candidate_record(cid)


# ------------------------------------------------------------------------------------ combination rule
def test_combine():
    assert combine("S6", "rejectedProven", "N1-ISO", "REJECTED_PROVEN") == ("rejectedProven", None)
    assert combine("S6", "rejectedProven", "N1-ISO", "UNDECIDED") == ("undecided", "CROSSCHECK_DISAGREE")
    assert combine("S6", "rejectedProven", "N1-ISO", None) == ("undecided", "CROSSCHECK_DISAGREE")
    assert combine("S6", "rejectedProven", "N1-KAPPA", "UNDECIDED") == ("rejectedProven", None)
    assert combine("S6", "screenedOut", None, None) == ("undecided", "WOULD_FIRE_LOWER_TIER")
    assert combine("S6", "pass", None, "UNDECIDED") == ("pass", None)


# ------------------------------------------------------------------------------------ S6 chain on real spectra
def _s6_chain(ctx, core_override=None, tamper=False):
    from rtt.core.dfpt import collect_qe_eph
    from rtt.certify import certify, certify_config
    from rtt.decide.inputs import envelope_from_ladder, tags_from_spectra
    from rtt.records import DecideInputs, Validity
    from rtt.units import EV_TO_K
    rec = add_al(ctx)
    L = ctx.ledger
    tr = L.start_tier_run("al", "S6", "dfpt", profile="full", params={"fixture": "qe76_al_bench"})
    d = DATA / "qe76_al_bench"
    refs = {p.name: register(L, tr, p) for p in sorted(d.glob("a2F.dos*"))}
    res = collect_qe_eph(d, "Al", el_ph_nsigma=10, smearing_admission="mode_positive")
    sp = []
    for s in res.spectra:
        r = refs[Path(s.provenance["source_file"]).name]
        s.provenance["source_artifact_id"] = r.id
        s.provenance["source_sha256"] = ("0" * 64) if tamper else r.sha256
        sp.append(s.with_envelope(0.10))
    L.finish_tier_run(tr, "done", metrics={"fixture": True})
    cc = certify_config(300.0, 0.15, mode="decide")
    v = certify("al", sp, cc)
    vid = store_verdict(ctx, "al", tr, "full", v, sp, cc, ctx.workdir)
    fine = {"A_K": [s.A() for s in sp], "lambda": [s.lambda_() for s in sp],
            "m_K2": [2 * sum(o * w for o, w in zip(s.omega_K, s.weight_K)) for s in sp]}
    p = Prov("fullDFPT", tr)
    iso, eps = envelope_from_ladder(fine, fine, p, omega_log_by_q=[1.0, 1.0])
    val = Validity(EF=11.7 * EV_TO_K, omegaMax=max(s.omega_max() for s in sp), imagAbsMaxK=0.0, containsHLiBe=False,
                   anharmonicCheck=False, anharmShiftAbs=0.0, alpha2FRecomputedAnharmonic=False,
                   muStarNegativeSuspected=False, admittedSmearings=len(sp), qMeshes=2, nSheets=1)
    di = DecideInputs(tier="fullDFPT", val=val, iso=iso, eps=eps, tags=tags_from_spectra(sp),
                      provenance={"tier_run_ids": [tr]})
    status = core_override or v.status.value
    return record_dfpt_decision(ctx, rec, "S6", tr, di, vid, status), vid


def test_s6_chain_rejects_al_at_300K(ctx):
    (o, out), vid = _s6_chain(ctx)
    assert o.outcome == "rejectedProven" and o.rule == "N1-ISO"
    assert out.outcome == "rejectedProven"
    cur = ctx.ledger.current_outcome("al", "S6")
    assert cur["outcome"] == "rejectedProven" and cur["decision_id"]
    v = ctx.ledger.verdict(vid)
    assert v["status"] == "REJECTED_PROVEN" and v["recheck_ok"] == 1 and v["spectrum_tier"] == "full"
    rec = ctx.ledger.candidate_record("al")
    assert rec.stages["S6"].decide["rule"] == "N1-ISO" and "M-ISO" in rec.stages["S6"].decide["assumptions"]
    s = ctx.ledger.summary_counts()
    assert s["core_rejected_proven_by_tier"] == {"full": 1}


def test_s6_chain_crosscheck_disagree(ctx):
    (o, out), _ = _s6_chain(ctx, core_override="UNDECIDED")
    assert o.outcome == "rejectedProven" and out.outcome == "undecided" and out.reason == "CROSSCHECK_DISAGREE"


def test_s6_chain_tampered_hash_refused(ctx):
    import sqlite3
    with pytest.raises(sqlite3.IntegrityError, match="hashed"):
        _s6_chain(ctx, tamper=True)


# ------------------------------------------------------------------------------------ scheduler
def test_scheduler_queue_and_budget(ctx):
    add_al(ctx, "al", lane="C", T_s=0.5875)
    st = StructureData(lattice=[[3, 0, 0], [0, 3, 0], [0, 0, 3]], species=["Mg", "H", "H"],
                       frac_coords=[[0, 0, 0], [0.5, 0.5, 0], [0.5, 0, 0.5]])
    ctx.ledger.insert_candidate("g1", "g" * 64, st, Origin(lane="G1"), cost_class="C3")
    for cid in ("al", "g1"):
        rec = ctx.ledger.candidate_record(cid)
        for stg in ("S1", "S2", "S3"):
            write_outcome(ctx, cid, Outcome(stg, "PASS"), rec=rec)
    ctx.ledger.insert_label("g1", "score_s3", 0.7, "S3")
    assert next_stage(ctx.ledger, "al") == "S4" and rebuild_queue(ctx.ledger) == 2
    jobs = next_jobs(ctx.ledger, ctx.machine)
    assert [j.candidate_id for j in jobs] == ["al", "g1"] and jobs[0].control
    write_outcome(ctx, "g1", Outcome("S4", "EXCLUDED_ASSUMPTION", criterion="H1", value=0.3, threshold=0.0,
                                     source_id="doc02#5a:H1"))
    assert next_stage(ctx.ledger, "g1") is None
    # control lane outcomes use the surrogate target (registered in inputs)
    assert ctx.ledger.current_outcome("al", "S3")["targets_sha256"] != ctx.targets.sha256


def test_stage_protocol():
    for sid, st in STAGES.items():
        assert isinstance(st, Stage) and st.id == sid


# ------------------------------------------------------------------------------------ scoring / ranking
def test_score_s3_forms():
    m = {"features": {"H_f": 0.6, "light_fraction": 0.7, "M_eff_amu": 2.0}, "phonons": {"omega2_light_K": 1500.0},
         "e_hull_ml_self_eV": 0.05}
    a = score_s3(m, load_targets())
    b = score_s3(m["features"], m["phonons"], None)
    assert 0 < b <= 1 and 0 < a <= 1
    worse = score_s3({**m, "e_hull_ml_self_eV": 0.4})
    assert worse < score_s3(m)
    assert score_s3({**m, "features": {**m["features"], "motif": "H2_molecular"}}) < score_s3(m)


def test_final_rank_order():
    rows = final_rank([
        {"candidate_id": "a", "decide_outcome": "rejectedProven", "T_low_K": 400.0, "A_mid_K": 3000},
        {"candidate_id": "b", "decide_outcome": "undecided", "T_low_K": 100.0, "A_mid_K": 1000},
        {"candidate_id": "c", "decide_outcome": "undecided", "T_low_K": 150.0, "A_mid_K": 900, "T_up_K": 1e4},
        {"candidate_id": "d", "decide_outcome": "pass", "T_low_K": None, "A_mid_K": 500},
        {"candidate_id": "e", "decide_outcome": "undecided", "T_low_K": 150.0, "A_mid_K": 1200},
    ])
    assert [r.candidate_id for r in rows] == ["d", "e", "c", "b", "a"] and rows[0].rank == 1


# ------------------------------------------------------------------------------------ config (R7)
def test_unsourced_threshold_refused(tmp_path):
    with pytest.raises(ConfigError, match="R7"):
        validate_sources({"S3": {"ehull_self_max_eV": 0.42}})
    with pytest.raises(ConfigError, match="R7"):
        validate_sources({"S5": {"margin_M_c": {"value": 2.0}}})
    validate_sources({"S5": {"margin_M_c": {"value": 2.0, "source": "doc11#4"}}})
    p = tmp_path / "t.toml"
    p.write_text((SRC / "configs" / "tiers.toml").read_text(encoding="utf-8").replace(
        'source = "doc11#3.1"\n', ""), encoding="utf-8")
    with pytest.raises(ConfigError):
        load_tiers(p)


def test_threshold_lookup():
    ti = load_tiers()
    th = ti.threshold("S5.margin_M_c")
    assert th.value == 2.0 and th.provisional and "doc11#4" in th.source_id
    assert ti.threshold("S5.q").source_id == "doc11#3.1"
    t = load_targets()
    assert float(t.floor) == pytest.approx(1299.444435)
    assert float(t.with_T_star(80.0).floor) == pytest.approx(4.33148145 * 80)


# ------------------------------------------------------------------------------------ import graph (R1b, R1c)
def _imports(path: Path) -> set[str]:
    tree = ast.parse(path.read_text(encoding="utf-8"))
    names = set()
    for n in ast.walk(tree):
        if isinstance(n, ast.ImportFrom):
            names.update(f"{n.module}.{a.name}" for a in n.names)
        elif isinstance(n, ast.Import):
            names.update(a.name for a in n.names)
    return names


def _src_files():
    return [p for p in SRC.rglob("*.py") if p.relative_to(SRC).parts[0] != "core"]


def test_only_certify_run_calls_core_engine():
    bad = []
    for p in _src_files():
        rel = p.relative_to(SRC).as_posix()
        imps = _imports(p)
        hits = {i for i in imps if i.endswith((".certify_spectra", ".classify", ".apply_input_policy"))
                and ("rtt.core." in i or i.startswith(("core.", ".core.", "..core.")))}
        if hits and rel != "certify/run.py":
            bad.append((rel, hits))
    assert not bad, bad


def test_no_verdict_or_status_construction_outside_core():
    """R1b: Verdict objects are built only by the core engine."""
    bad = []
    for p in _src_files():
        tree = ast.parse(p.read_text(encoding="utf-8"))
        for n in ast.walk(tree):
            if isinstance(n, ast.Call):
                f = n.func
                name = f.id if isinstance(f, ast.Name) else (f.attr if isinstance(f, ast.Attribute) else "")
                if name in ("Verdict", "from_dict") and isinstance(f, ast.Name | ast.Attribute):
                    if name == "Verdict" or (isinstance(f, ast.Attribute) and getattr(f.value, "id", "") == "Verdict"):
                        bad.append(p.relative_to(SRC).as_posix())
    assert not bad, bad


def test_decide_outcomes_only_from_decide_core():
    """R1c: decisions rows are written only with rtt.decide.core.decide_t output (origin string
    used in exactly one place outside records/ledger)."""
    users = [p.relative_to(SRC).as_posix() for p in _src_files()
             if "insert_decision(" in p.read_text(encoding="utf-8")]
    assert set(users) <= {"engine/stages.py", "ledger/store.py"}, users


def test_cli_decide_check():
    from rtt.cli import main
    if not (ROOT / "lean" / "export" / "decide_cases.json").exists():
        pytest.skip("no Lean export")
    assert main(["decide-check"]) == 0


# ------------------------------------------------------------------------------------ RC1 F12/F13/N6
def _pass_to_s3(ctx, cid):
    rec = ctx.ledger.candidate_record(cid)
    for stg in ("S1", "S2", "S3"):
        write_outcome(ctx, cid, Outcome(stg, "PASS"), rec=rec)
    return rec


def test_retry_error_and_backoff(ctx):
    import datetime as dt
    from rtt.engine.retry import backoff_ready, exhausted, failed_attempts, needs_run
    add_al(ctx, "g", lane="D")
    _pass_to_s3(ctx, "g")
    for i in range(3):
        assert next_stage(ctx.ledger, "g", ctx.tiers) == "S4"
        write_outcome(ctx, "g", Outcome("S4", "ERROR", reason="spin scf not converged", retry=i + 1))
        assert failed_attempts(ctx.ledger, "g", "S4") == i + 1
    assert not needs_run(ctx.ledger, "g", "S4", ctx.tiers) and next_stage(ctx.ledger, "g", ctx.tiers) is None
    assert exhausted(ctx.ledger, ctx.tiers) == [("g", "S4", "ERROR", 3)]
    assert not backoff_ready(ctx.ledger, "g", "S4", ctx.tiers)
    later = dt.datetime.now(dt.timezone.utc) + dt.timedelta(hours=25)
    assert backoff_ready(ctx.ledger, "g", "S4", ctx.tiers, now=later)


def test_rank_cut_deferral_is_not_a_failure(ctx):
    from rtt.engine.retry import failed_attempts, needs_run
    add_al(ctx, "g", lane="D")
    rec = ctx.ledger.candidate_record("g")
    for stg in ("S1", "S2"):
        write_outcome(ctx, "g", Outcome(stg, "PASS"), rec=rec)
    for _ in range(5):
        write_outcome(ctx, "g", Outcome("S3", "DEFERRED", reason="rank_cut(30, 25)"))
    assert failed_attempts(ctx.ledger, "g", "S3") == 0 and needs_run(ctx.ledger, "g", "S3", ctx.tiers)


def test_terminal_removal_stops(ctx):
    add_al(ctx, "g", lane="D")
    _pass_to_s3(ctx, "g")
    write_outcome(ctx, "g", Outcome("S4", "EXCLUDED_ASSUMPTION", criterion="H1", value=1.0, threshold=0.0,
                                    source_id="doc02#5a:H1"))
    assert next_stage(ctx.ledger, "g", ctx.tiers) is None


def test_s4_decide_uses_sourced_thresholds(ctx, tmp_path):
    """RC1 F1: the engine hands decide_s4 Threshold objects; an insulating magnet is excluded."""
    from rtt.dft import electronic as E
    from rtt.dft.workflow import DFPTRun
    from rtt.records import StageResult
    run = DFPTRun(material_id="x", profile={"name": "s4", "kind": "electronic"}, workdir=str(tmp_path), status="done")
    r = E.ElectronicResult(run=run, nat=4)
    r.gap_eV, r.moment_muB, r.abs_moment_muB = 2.0, 3.0, 3.0
    o = STAGES["S4"].decide(StageResult("S4", "x", True, {}, None, payload=r), ctx.tiers, ctx.targets)
    assert o.outcome == "EXCLUDED_ASSUMPTION" and o.criterion.startswith("H1")
    th = ctx.targets.threshold("hard.H1_metallic")
    assert th.value == 0.0 and th.source_id == "doc02#5a:H1"


def test_control_surrogate_target_reaches_decisions(ctx):
    """N6: a control's stages are decided at its surrogate T*; a generated candidate of the same
    structure is evaluated at the real T* (not deduplicated against the control)."""
    rec = add_al(ctx, "al", lane="C", T_s=0.5875)
    t = ctx.targets_for(rec)
    assert float(t.T_star) == 0.5875 and float(ctx.targets.T_star) == 300.0
    from rtt.records import Origin as O
    from rtt.structure.stage import run_s2
    st = rec.structure
    ctx.ledger.insert_candidate("gen", "h" * 64, st, O(lane="G1"))
    for cid in ("al", "gen"):
        write_outcome(ctx, cid, Outcome("S1", "PASS"), rec=ctx.ledger.candidate_record(cid))
    out = run_s2(ctx.ledger, targets=ctx.targets, tiers=ctx.tiers, novelty=False)
    assert out["al"].outcome == "PASS" and out["gen"].outcome == "PASS"
    assert ctx.targets_for(ctx.ledger.candidate_record("gen")).sha256 == ctx.targets.sha256


# ------------------------------------------------------------------------------------ B3 wiring (2026-10-29)
def test_s5_decide_gets_s4_sheets(ctx):
    """S5 decide passes S4's n_sheets and sheet_w_min to decide_coarse (multiband A/w_min rule)."""
    from rtt.dft.coarse import CoarseMetrics, SigmaRow
    from rtt.records import StageResult
    rows = [SigmaRow(index=i + 3, sigma_Ry=0.01 * (i + 3), admitted=True, A_K=300.0, on_plateau=True) for i in range(3)]
    met = CoarseMetrics(rows=rows, plateau_sigma_min_Ry=0.03, min_freq_matdyn_cm1=50.0, q=[4, 4, 4],
                        elements=["Al"], ecutwfc_Ry=60.0)     # doc 16: unknown elements block removal
    s5 = STAGES["S5"]
    o1 = s5.decide(StageResult("S5", "x", True, {"n_sheets": 1}, None, payload={"metrics": met}), ctx.tiers, ctx.targets)
    o0 = s5.decide(StageResult("S5", "x", True, {}, None, payload={"metrics": met}), ctx.tiers, ctx.targets)
    o2 = s5.decide(StageResult("S5", "x", True, {"n_sheets": 2, "sheet_w_min": 0.1}, None, payload={"metrics": met}),
                   ctx.tiers, ctx.targets)
    assert o1.outcome == "SCREENED_OUT"          # single sheet, A far below floor / M_c
    assert o0.outcome == "DEFERRED"              # sheets unknown -> no removal
    assert o2.outcome != "SCREENED_OUT"          # A / w_min = 3000 K >= floor / M_c


def test_s6_min_freq_reaches_decide_inputs():
    from rtt.decide.inputs import build_inputs
    from rtt.records import CandidateRecord, Origin, StageEntry, StructureData
    st = StructureData(lattice=[[3, 0, 0], [0, 3, 0], [0, 0, 3]], species=["Al"], frac_coords=[[0, 0, 0]])
    rec = CandidateRecord("c", "s", "Al", st, Origin(lane="D"))
    rec.stages["S4"] = StageEntry("PASS", metrics={"E_F_eV": 11.7, "n_sheets": 1})
    for mf, want in ((-50.0, 50.0), (120.0, 0.0)):
        rec.stages["S6"] = StageEntry("pending", tier_run_id=7, metrics={"min_freq_K": mf, "max_imag_mode_K": max(0, -mf)})
        assert build_inputs(rec, "fullDFPT").val.imagAbsMaxK == want
    rec.stages["S6"] = StageEntry("pending", tier_run_id=7, metrics={"min_freq_K": -80.0})
    assert build_inputs(rec, "fullDFPT").val.imagAbsMaxK == 80.0


def test_s5_cost_counts_both_meshes():
    from rtt.engine.scheduler import COST_CORE_H
    assert COST_CORE_H["S5"]["C1"] == pytest.approx(0.062 * 125 * (80 / 60) ** 1.5)   # doc 31 §6 refit
    assert load_tiers().get("S5.q_meshes") == [[3, 3, 3], [4, 4, 4]]


def test_s5_decide_two_mesh_t2_screen(ctx):
    from rtt.dft.coarse import CoarseMetrics, SigmaRow
    from rtt.records import ALLOWED_OUTCOMES, StageResult

    def m(q, A):
        rows = [SigmaRow(index=i + 3, sigma_Ry=0.01 * (i + 3), admitted=True, A_K=A, lambda_=0.5, omega_log_K=300.0,
                         on_plateau=True) for i in range(3)]
        return CoarseMetrics(rows=rows, plateau_sigma_min_Ry=0.03, min_freq_matdyn_cm1=50.0, q=q)
    m3, m4 = m([3, 3, 3], 290.0), m([4, 4, 4], 300.0)
    res = StageResult("S5", "x", True, {"n_sheets": 1}, None, payload={"metrics": m4, "metrics_by_mesh": [m3, m4]})
    o = STAGES["S5"].decide(res, ctx.tiers, ctx.targets)
    assert o.outcome in ALLOWED_OUTCOMES["S5"]
    assert o.outcome != "SCREENED_OUT" or "envelope" in (o.criterion or "")
