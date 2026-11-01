"""Ledger triggers (architecture §6.2, contract R1b/R1c/R2-R5): forbidden rows are rejected by
SQLite itself, independent of the Python layer."""
from __future__ import annotations

import json
import sqlite3

import pytest

from rtt.config import load_machine, load_targets, load_tiers
from rtt.ledger.store import Ledger
from rtt.records import (ALLOWED_OUTCOMES, ORIGIN_CERTIFY, ORIGIN_DECIDE, DecideInputs, EpsData, IsoData, Origin,
                         Outcome, Prov, StructureData, Validity)

T, TI, M = load_targets(), load_tiers(), load_machine()
SHA = "a" * 64


@pytest.fixture
def L(tmp_path):
    led = Ledger(tmp_path / "ledger.sqlite")
    led.register_configs(T, TI, M)
    st = StructureData(lattice=[[3, 0, 0], [0, 3, 0], [0, 0, 3]], species=["Al"], frac_coords=[[0, 0, 0]])
    led.insert_candidate("c1", "s" * 64, st, Origin(lane="C", control="Al"))
    led.insert_candidate("c2", "t" * 64, st, Origin(lane="D"))
    yield led
    led.close()


def ins_outcome(L, cid, o: Outcome, decision_id=None, origin="rtt.engine.stages"):
    return L.insert_outcome(cid, o, targets_sha256=T.sha256, tiers_sha256=TI.sha256, origin=origin,
                            decision_id=decision_id)


def raw_outcome(L, cid, stage, outcome, **kw):
    cols = dict(candidate_id=cid, stage=stage, outcome=outcome, targets_sha256=T.sha256, tiers_sha256=TI.sha256,
                origin="test", created="now", provisional=0)
    cols.update(kw)
    q = f"INSERT INTO outcomes({','.join(cols)}) VALUES ({','.join('?' * len(cols))})"
    return L.x(q, list(cols.values()))


def dfpt_run(L, cid="c1", profile="full", kind="dfpt", done=True, stage="S6"):
    tr = L.start_tier_run(cid, stage, kind, profile=profile, params={"q": [4, 4, 4]})
    a_spec = L.insert_artifact(tr, "/x/a2F.dos3", SHA, 100)
    a_cert = L.insert_artifact(tr, "/x/cert.json", "b" * 64, 10)
    if done:
        L.finish_tier_run(tr, "done", wall_s=1.0)
    return tr, a_spec, a_cert


def spectra_json(a_spec, sha=SHA):
    return [{"alt_index": 0, "spectrum_digest": "d" * 64, "source_artifact_id": a_spec, "sha256": sha}]


def verdict(L, tr, a_cert, a_spec, cid="c1", status="REJECTED_PROVEN", spectrum_tier="full", recheck_ok=True,
            origin=ORIGIN_CERTIFY, spectra=None):
    return L.insert_verdict(cid, tr, spectrum_tier, status=status, T_low_K=1.0, T_up_K=9.0, T_star_K=300.0,
                            mu_hi=0.13, rel_envelope=0.1, certificate_artifact_id=a_cert, recheck_ok=recheck_ok,
                            spectra=spectra if spectra is not None else spectra_json(a_spec), origin=origin)


def decide_input(L, tier="fullDFPT", tr_ids=(), cid="c1"):
    p = Prov("fullDFPT" if tier in ("fullDFPT", "converged") else tier, 1)
    di = DecideInputs(tier=tier, val=Validity(), iso=IsoData(provA=p, provLam=p, provM=p), eps=EpsData(prov=p),
                      provenance={"tier_run_ids": list(tr_ids)})
    return L.insert_decide_input(cid, di)


ASSUME = ["ME-MODEL", "CERT-H", "ENV-SOUND", "M-ISO"]


# ------------------------------------------------------------------------------------ append-only
@pytest.mark.parametrize("table", ["outcomes", "decisions", "verdicts", "verdict_spectra", "artifacts",
                                   "candidates", "inputs", "structures", "decide_inputs", "meta", "stage_outcomes"])
def test_append_only(L, table):
    tr, a_spec, a_cert = dfpt_run(L)
    vid = verdict(L, tr, a_cert, a_spec)
    di = decide_input(L, tr_ids=[tr])
    L.insert_decision("c1", "S6", "fullDFPT", di, "pass")
    ins_outcome(L, "c1", Outcome("S1", "PASS"))
    with pytest.raises(sqlite3.IntegrityError):
        L.x(f"DELETE FROM {table}")
    with pytest.raises(sqlite3.IntegrityError):
        L.x(f"UPDATE {table} SET rowid = rowid")


def test_queue_is_mutable(L):
    L.queue_put("c1", "S3", 1.0)
    L.queue_put("c1", "S3", 2.0)
    assert L.queue_items("S3") == [("c1", "S3", 2.0)]
    L.queue_clear()


def test_stage_outcomes_fixed(L):
    with pytest.raises(sqlite3.IntegrityError):
        L.x("INSERT INTO stage_outcomes VALUES ('S3','EXCLUDED_ASSUMPTION')")
    rows = set(L.qa("SELECT stage, outcome FROM stage_outcomes"))
    assert rows == {(s, o) for s, os in ALLOWED_OUTCOMES.items() for o in os}


# ------------------------------------------------------------------------------------ tier_runs
def test_tier_run_lifecycle(L):
    with pytest.raises(sqlite3.IntegrityError):
        L.x("INSERT INTO tier_runs(candidate_id, stage, kind, profile, attempt, state, params_json, versions_json)"
            " VALUES ('c1','S6','dfpt','full',1,'done','{}','{}')")
    tr = L.start_tier_run("c1", "S6", "dfpt", profile="full")
    with pytest.raises(sqlite3.IntegrityError):
        L.x("UPDATE tier_runs SET profile='coarse' WHERE id=?", (tr,))
    L.finish_tier_run(tr, "done")
    with pytest.raises(sqlite3.IntegrityError):
        L.finish_tier_run(tr, "failed")
    tr2 = L.start_tier_run("c1", "S6", "dfpt", profile="full")
    L.finish_tier_run(tr2, "failed")
    with pytest.raises(sqlite3.IntegrityError):
        L.insert_artifact(tr2, "/x", SHA, 1)


# ------------------------------------------------------------------------------------ outcomes_allowed
ALL_OUTCOMES = ["PASS", "SCREENED_OUT", "EXCLUDED_ASSUMPTION", "DEFERRED", "INVALID", "DUPLICATE", "ERROR",
                "rejectedProven", "rejectedAssuming", "screenedOut", "undecided", "pass"]
FORBIDDEN = [(s, o) for s in ["S1", "S2", "S3", "S4", "S5", "S6", "S7", "S8", "S9"] for o in ALL_OUTCOMES
             if o not in ALLOWED_OUTCOMES.get(s, frozenset())]


@pytest.mark.parametrize("stage,outcome", FORBIDDEN)
def test_forbidden_stage_outcome_rejected(L, stage, outcome):
    with pytest.raises(sqlite3.IntegrityError):
        raw_outcome(L, "c1", stage, outcome, criterion="coarse-A", value="1", threshold="2", source_id="x",
                    reason="r")


def test_examples_from_spec(L):
    for stage, outcome in [("S3", "EXCLUDED_ASSUMPTION"), ("S5", "rejectedProven"), ("S5", "rejectedAssuming"),
                           ("S1", "ERROR"), ("S6", "SCREENED_OUT"), ("S6", "screenedOut")]:
        with pytest.raises(sqlite3.IntegrityError):
            raw_outcome(L, "c1", stage, outcome, criterion="H1", value="1", threshold="0", source_id="x")


def test_screened_out_needs_criterion_value_threshold_source(L):
    for missing in ("criterion", "value", "threshold", "source_id"):
        kw = dict(criterion="H7-ml", value="0.5", threshold="0.42", source_id="doc02#5a:H7")
        kw.pop(missing)
        with pytest.raises(sqlite3.IntegrityError, match="criterion"):
            raw_outcome(L, "c1", "S3", "SCREENED_OUT", **kw)
    raw_outcome(L, "c1", "S3", "SCREENED_OUT", criterion="H7-ml", value="0.5", threshold="0.42", source_id="doc02")
    with pytest.raises(sqlite3.IntegrityError):
        raw_outcome(L, "c1", "S3", "DEFERRED")


def test_excluded_assumption_criteria(L):
    with pytest.raises(sqlite3.IntegrityError, match="R3"):
        raw_outcome(L, "c1", "S4", "EXCLUDED_ASSUMPTION", criterion="H7", value="1", threshold="0", source_id="x")
    raw_outcome(L, "c1", "S4", "EXCLUDED_ASSUMPTION", criterion="H1", value="0.3", threshold="0", source_id="x")
    raw_outcome(L, "c1", "S7", "EXCLUDED_ASSUMPTION", criterion="H7-zpe", value="0.5", threshold="0.35", source_id="x")


def test_s5_screen_only_coarse_A(L):
    with pytest.raises(sqlite3.IntegrityError, match="coarse-A"):
        raw_outcome(L, "c1", "S5", "SCREENED_OUT", criterion="imag-modes", value="1", threshold="0", source_id="x")
    raw_outcome(L, "c1", "S5", "SCREENED_OUT", criterion="coarse-A", value="300", threshold="649", source_id="doc11#4")


def test_unregistered_config_hash(L):
    with pytest.raises(sqlite3.IntegrityError, match="R7"):
        raw_outcome(L, "c1", "S1", "PASS", targets_sha256="0" * 64)


# ------------------------------------------------------------------------------------ verdicts_bind
def test_valid_verdict_accepted(L):
    tr, a_spec, a_cert = dfpt_run(L)
    vid = verdict(L, tr, a_cert, a_spec)
    assert L.q1("SELECT COUNT(*) FROM verdict_spectra WHERE verdict_id=?", (vid,))[0] == 1


def test_verdict_without_dfpt_tier_run(L):
    tr = L.start_tier_run("c1", "S3", "ml")
    a = L.insert_artifact(tr, "/x", SHA, 1)
    L.finish_tier_run(tr, "done")
    with pytest.raises(sqlite3.IntegrityError, match="verdicts_bind"):
        verdict(L, tr, a, a)


def test_verdict_run_not_done(L):
    tr, a_spec, a_cert = dfpt_run(L, done=False)
    with pytest.raises(sqlite3.IntegrityError, match="verdicts_bind"):
        verdict(L, tr, a_cert, a_spec)


def test_verdict_profile_tier_mismatch(L):
    tr, a_spec, a_cert = dfpt_run(L, profile="coarse", stage="S5")
    with pytest.raises(sqlite3.IntegrityError, match="verdicts_bind"):
        verdict(L, tr, a_cert, a_spec, spectrum_tier="full")
    verdict(L, tr, a_cert, a_spec, spectrum_tier="coarse")       # coarse evidence is fine as coarse


def test_verdict_other_candidate_run(L):
    tr, a_spec, a_cert = dfpt_run(L, cid="c2")
    with pytest.raises(sqlite3.IntegrityError, match="verdicts_bind"):
        verdict(L, tr, a_cert, a_spec, cid="c1")


def test_verdict_certificate_of_other_run(L):
    tr, a_spec, _ = dfpt_run(L)
    tr2, _, a_cert2 = dfpt_run(L)
    with pytest.raises(sqlite3.IntegrityError, match="certificate"):
        verdict(L, tr, a_cert2, a_spec)


@pytest.mark.parametrize("status", ["REJECTED_PROVEN", "CERTIFIED_CANDIDATE"])
def test_proven_verdict_needs_recheck(L, status):
    tr, a_spec, a_cert = dfpt_run(L)
    with pytest.raises(sqlite3.IntegrityError, match="recheck"):
        verdict(L, tr, a_cert, a_spec, status=status, recheck_ok=False)


def test_undecided_verdict_without_recheck_ok_allowed(L):
    tr, a_spec, a_cert = dfpt_run(L)
    verdict(L, tr, a_cert, a_spec, status="UNDECIDED", recheck_ok=False)


def test_verdict_origin(L):
    tr, a_spec, a_cert = dfpt_run(L)
    with pytest.raises(sqlite3.IntegrityError, match="origin"):
        verdict(L, tr, a_cert, a_spec, origin="rtt.core.certify.classify")


def test_verdict_artifact_hash_must_match(L):
    tr, a_spec, a_cert = dfpt_run(L)
    with pytest.raises(sqlite3.IntegrityError, match="hashed"):
        verdict(L, tr, a_cert, a_spec, spectra=spectra_json(a_spec, sha="e" * 64))
    with pytest.raises(sqlite3.IntegrityError, match="hashed"):
        verdict(L, tr, a_cert, a_spec, spectra=[])
    with pytest.raises(sqlite3.IntegrityError, match="hashed"):
        verdict(L, tr, a_cert, a_spec, spectra=[{"alt_index": 0, "source_artifact_id": a_spec, "sha256": SHA}])


def test_verdict_spectra_rows_must_match(L):
    tr, a_spec, a_cert = dfpt_run(L)
    vid = verdict(L, tr, a_cert, a_spec)
    with pytest.raises(sqlite3.IntegrityError, match="verdict_spectra_bind"):
        L.x("INSERT INTO verdict_spectra VALUES (?,?,?,?)", (vid, 1, "f" * 64, a_spec))


# ------------------------------------------------------------------------------------ decisions_guard
def full_setup(L, status="REJECTED_PROVEN"):
    tr, a_spec, a_cert = dfpt_run(L)
    vid = verdict(L, tr, a_cert, a_spec, status=status, recheck_ok=status != "UNDECIDED")
    di = decide_input(L, "fullDFPT", [tr])
    return tr, vid, di


def test_valid_rejection_chain(L):
    tr, vid, di = full_setup(L)
    did = L.insert_decision("c1", "S6", "fullDFPT", di, "rejectedProven", firing_rule="N1-ISO",
                            assumptions=ASSUME, crosscheck_verdict_id=vid)
    ins_outcome(L, "c1", Outcome("S6", "rejectedProven", tier_run_id=tr), decision_id=did)
    assert L.current_outcome("c1", "S6")["outcome"] == "rejectedProven"


def test_rejection_requires_proof_grade_tier(L):
    tr, a_spec, a_cert = dfpt_run(L, profile="coarse", stage="S5")
    verdict(L, tr, a_cert, a_spec, spectrum_tier="coarse")
    di = decide_input(L, "coarseDFPT", [tr])
    with pytest.raises(sqlite3.IntegrityError, match="decisions_guard"):
        L.insert_decision("c1", "S5", "coarseDFPT", di, "rejectedProven", firing_rule="NM-G", assumptions=ASSUME)
    L.insert_decision("c1", "S5", "coarseDFPT", di, "screenedOut")


def test_ml_tier_never_rejects(L):
    tr = L.start_tier_run("c1", "S3", "ml")
    L.finish_tier_run(tr, "done")
    di = decide_input(L, "ml", [tr])
    for o in ("rejectedProven", "rejectedAssuming"):
        with pytest.raises(sqlite3.IntegrityError):
            L.insert_decision("c1", "S6", "ml", di, o, firing_rule="N1-KAPPA", assumptions=ASSUME)


def test_decision_origin(L):
    tr, vid, di = full_setup(L)
    with pytest.raises(sqlite3.IntegrityError, match="origin"):
        L.insert_decision("c1", "S6", "fullDFPT", di, "pass", origin="somebody")


def test_decision_tier_must_match_input_and_runs(L):
    tr, vid, di = full_setup(L)
    with pytest.raises(sqlite3.IntegrityError, match="decide"):
        L.insert_decision("c1", "S6", "converged", di, "pass")         # input says fullDFPT
    di_bad = decide_input(L, "fullDFPT", [])                             # references no full run
    with pytest.raises(sqlite3.IntegrityError, match="tier_runs"):
        L.insert_decision("c1", "S6", "fullDFPT", di_bad, "pass")
    tr2, _, _ = dfpt_run(L, cid="c2")
    di_other = decide_input(L, "fullDFPT", [tr, tr2])                    # run of another candidate
    with pytest.raises(sqlite3.IntegrityError, match="tier_runs"):
        L.insert_decision("c1", "S6", "fullDFPT", di_other, "pass")


def test_rejection_without_rechecked_verdict(L):
    tr, a_spec, a_cert = dfpt_run(L)
    di = decide_input(L, "fullDFPT", [tr])
    with pytest.raises(sqlite3.IntegrityError, match="rechecked verdict"):
        L.insert_decision("c1", "S6", "fullDFPT", di, "rejectedProven", firing_rule="N1-KAPPA", assumptions=ASSUME)
    verdict(L, tr, a_cert, a_spec, status="UNDECIDED", recheck_ok=False)
    with pytest.raises(sqlite3.IntegrityError, match="rechecked verdict"):
        L.insert_decision("c1", "S6", "fullDFPT", di, "rejectedProven", firing_rule="N1-KAPPA", assumptions=ASSUME)


def test_n1_iso_needs_crosscheck(L):
    tr, vid, di = full_setup(L)
    with pytest.raises(sqlite3.IntegrityError, match="crosscheck"):
        L.insert_decision("c1", "S6", "fullDFPT", di, "rejectedProven", firing_rule="N1-ISO", assumptions=ASSUME)
    # crosscheck verdict must be REJECTED_PROVEN
    tr2, vid_u, di2 = full_setup(L, status="UNDECIDED")
    _, vid_ok, _ = full_setup(L)
    with pytest.raises(sqlite3.IntegrityError):
        L.insert_decision("c1", "S6", "fullDFPT", di2, "rejectedProven", firing_rule="N1-ISO",
                          assumptions=ASSUME, crosscheck_verdict_id=vid_u)
    # multiband rules rest on decide_t alone (need a rechecked verdict, no crosscheck)
    L.insert_decision("c1", "S6", "fullDFPT", di, "rejectedProven", firing_rule="N1-KAPPA", assumptions=ASSUME)


def test_rejection_needs_rule_and_assumptions(L):
    tr, vid, di = full_setup(L)
    for rule, assume in [(None, ASSUME), ("RANK-AIM", ASSUME), ("N1-KAPPA", None), ("N1-KAPPA", [])]:
        with pytest.raises(sqlite3.IntegrityError):
            L.insert_decision("c1", "S6", "fullDFPT", di, "rejectedProven", firing_rule=rule,
                              assumptions=assume if assume is not None else [])
    with pytest.raises(sqlite3.IntegrityError, match="rejectedAssuming"):
        L.insert_decision("c1", "S6", "fullDFPT", di, "rejectedAssuming", firing_rule="NM-G", assumptions=ASSUME)


def test_screened_out_decision_at_full_tier_allowed_but_not_as_s6_removal(L):
    tr, vid, di = full_setup(L)
    did = L.insert_decision("c1", "S6", "fullDFPT", di, "screenedOut")
    with pytest.raises(sqlite3.IntegrityError, match="outcomes_allowed"):
        ins_outcome(L, "c1", Outcome.__new__(Outcome), decision_id=did) if False else \
            raw_outcome(L, "c1", "S6", "screenedOut", decision_id=did)
    raw_outcome(L, "c1", "S6", "undecided", decision_id=did, reason="WOULD_FIRE_LOWER_TIER")


# ------------------------------------------------------------------------------------ outcomes bound to decisions
def test_decide_outcome_without_decision(L):
    tr, vid, di = full_setup(L)
    with pytest.raises(sqlite3.IntegrityError, match="decisions row"):
        raw_outcome(L, "c1", "S6", "rejectedProven", tier_run_id=tr)
    did = L.insert_decision("c1", "S6", "fullDFPT", di, "pass")
    with pytest.raises(sqlite3.IntegrityError, match="decisions row"):
        raw_outcome(L, "c1", "S6", "rejectedProven", tier_run_id=tr, decision_id=did)   # outcome != decision
    with pytest.raises(sqlite3.IntegrityError, match="decisions row"):
        raw_outcome(L, "c2", "S6", "pass", decision_id=did)                             # other candidate


def test_rejection_outcome_needs_full_run(L):
    tr, vid, di = full_setup(L)
    did = L.insert_decision("c1", "S6", "fullDFPT", di, "rejectedProven", firing_rule="N1-KAPPA", assumptions=ASSUME)
    trc, _, _ = dfpt_run(L, profile="coarse", stage="S5")
    with pytest.raises(sqlite3.IntegrityError, match="R4"):
        raw_outcome(L, "c1", "S6", "rejectedProven", tier_run_id=trc, decision_id=did)
    with pytest.raises(sqlite3.IntegrityError, match="R4"):
        raw_outcome(L, "c1", "S6", "rejectedProven", decision_id=did)
    raw_outcome(L, "c1", "S6", "rejectedProven", tier_run_id=tr, decision_id=did)


def test_combination_downgrade_allowed(L):
    """CROSSCHECK_DISAGREE: decision rejectedProven, stage outcome undecided (§3.2)."""
    tr, vid, di = full_setup(L)
    did = L.insert_decision("c1", "S6", "fullDFPT", di, "rejectedProven", firing_rule="N1-ISO",
                            assumptions=ASSUME, crosscheck_verdict_id=vid)
    raw_outcome(L, "c1", "S6", "undecided", decision_id=did, reason="CROSSCHECK_DISAGREE")


# ------------------------------------------------------------------------------------ store helpers
def test_current_outcome_latest_and_record(L, tmp_path):
    ins_outcome(L, "c1", Outcome("S3", "DEFERRED", reason="rank_cut(30, 25)"))
    ins_outcome(L, "c1", Outcome("S3", "PASS"))
    assert L.current_outcome("c1", "S3")["outcome"] == "PASS"
    rec = L.candidate_record("c1")
    assert rec.stages["S3"].outcome == "PASS" and rec.formula == "Al"
    man = L.export_jsonl(tmp_path / "export", verify_artifacts=False)
    assert man.exists() and (tmp_path / "export" / "outcomes.jsonl").read_text().count("\n") == 2
    s = L.summary_counts()
    assert s["outcomes"]["S3"] == {"PASS": 1} and s["candidates"] == 2


def test_python_outcome_validation():
    with pytest.raises(ValueError):
        Outcome("S3", "EXCLUDED_ASSUMPTION", criterion="H1", source_id="x")
    with pytest.raises(ValueError):
        Outcome("S5", "SCREENED_OUT")
    with pytest.raises(ValueError):
        Outcome("S3", "DEFERRED")
