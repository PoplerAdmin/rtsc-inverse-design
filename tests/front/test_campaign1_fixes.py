"""Doc 31 (campaign-1 fixes, R30 items): pool-wide S3 rank cut (B9), lane minimum shares,
metallicity prior + valence parity, chemsys list (B11), S1 generator-domain deferral, S3 window
domain, S2 representative order. CPU only, fake MatterSim hooks."""
from __future__ import annotations

import json
import logging

import pytest

from rtt.mlff.screen import close_rank_pool, decide_s3, rank_pool, run_s3
from rtt.records import DEFERRED, PASS, Outcome
from rtt.score.metallicity import metallicity_prior, parse_formula, valence_electrons
from rtt.score.s3 import score_s3
from rtt.structure.stage import run_s2

from test_s3 import _fake_hull, _fake_relax  # noqa: E402 (rootdir-prepended tests/front)


def _ingest(ledger, targets, tiers, species, lane="G1", a=3.0):
    from rtt.records import Origin, RawStructure, StructureData
    from rtt.source.ingest import ingest_raw
    n = len(species)
    fr = [[0.5 * (i % 2), 0.5 * ((i // 2) % 2), (0.37 * i) % 1] for i in range(n)]
    return ingest_raw(ledger, RawStructure(structure=StructureData(
        lattice=[[a, 0, 0], [0, a, 0], [0, 0, a * (1 + 0.15 * n)]], species=list(species), frac_coords=fr),
        origin=Origin(lane=lane)), targets=targets, tiers=tiers)


# ------------------------------------------------------------------ B9: one cut per pool
def _items(spec):
    return [(cid, lane, Outcome(stage="S3", outcome=PASS), s) for cid, lane, s in spec]


def test_rank_pool_deterministic_and_order_free():
    spec = [("a", "D", 0.9), ("b", "D", 0.8), ("c", "G1", 0.7), ("d", "G1", 0.7), ("e", "G4", 0.1), ("k", "C", 0.0)]
    o1, r1 = rank_pool(_items(spec), 3, "src")
    o2, r2 = rank_pool(_items(list(reversed(spec))), 3, "src")
    assert {k: v.outcome for k, v in o1.items()} == {k: v.outcome for k, v in o2.items()}
    assert r1 == r2                                                        # input order irrelevant
    assert [c for c, o in o1.items() if o.outcome == PASS] == ["a", "b", "c", "k"]   # tie c/d -> id
    assert o1["d"].reason == "rank_cut(4, 3)" and "CONTROL_EXEMPT_RANK_CUT" in o1["k"].tags


def test_lane_min_share_reserves_slots():
    spec = [(f"d{i}", "D", 0.9 - 0.01 * i) for i in range(10)] + [(f"g{i}", "G1", 0.5 - 0.01 * i) for i in range(5)]
    out, rows = rank_pool(_items(spec), 5, "cut", min_share={"G": 0.4}, share_src="share", pool_id="p1")
    passed = sorted(c for c, o in out.items() if o.outcome == PASS)
    assert passed == ["d0", "d1", "d2", "g0", "g1"]                       # floor(0.4*5) = 2 G slots
    assert "LANE_MIN_SHARE:G" in out["g0"].tags and out["g0"].source_id == "share"
    assert "LANE_MIN_SHARE:G" not in out["d0"].tags and "RANK_POOL:p1" in out["d0"].tags
    assert {r["cid"]: r["via"] for r in rows}["g1"] == "lane_min_share:G"
    # a share cannot exceed the cut, and an empty lane leaves its slots to the global order
    out, _ = rank_pool(_items(spec[:10]), 5, "cut", min_share={"G": 0.4})
    assert sum(o.outcome == PASS for o in out.values()) == 5


def test_tiers_carry_share_and_weight(tiers):
    from rtt.mlff.screen import lane_min_share
    shares, src = lane_min_share(tiers)
    assert shares == {"G": 0.40} and "CAMPAIGN_GENERATED_QUOTA" in src
    assert tiers.threshold("S3.metallicity_weight").value == 0.20


def test_chunked_run_cuts_once_over_the_pool(ledger, targets, tiers, tmp_path):
    """B9: with chunk_size 2, the first chunk must not take every slot."""
    cids = [_ingest(ledger, targets, tiers, ["Mg"] + ["H"] * (k + 1), a=3.0 + 0.3 * k) for k in range(5)]
    run_s2(ledger, targets=targets, tiers=tiers, novelty=False)
    out = run_s3(ledger, targets=targets, tiers=tiers, workdir=tmp_path, alignn=False, dftref=False, phonons=False,
                 relax_fn=_fake_relax, hull_fn=_fake_hull, rank_cut=2, chunk_size=2, pool_id="day1")
    cur = {c: ledger.current_outcome(c, "S3") for c in cids}
    passed = [c for c in cids if cur[c]["outcome"] == PASS]
    assert len(passed) == 2
    # the two passes are the top-2 by stored score over all 5 (not the first chunk's 2)
    scores = {c: json.loads(ledger.tier_run(cur[c]["tier_run_id"])["metrics_json"])["score_s3"] for c in cids}
    top2 = sorted(cids, key=lambda c: (-scores[c], c))[:2]
    assert sorted(passed) == sorted(top2)
    pools = ledger.rank_pools()
    assert len(pools) == 1 and pools[0]["pool_id"] == "day1" and pools[0]["free_slots"] == 2
    assert {r["cid"] for r in pools[0]["ranking"]} == set(cids)
    assert all("RANK_POOL:day1" in (cur[c]["tags"] or []) for c in passed)
    assert set(out) == set(cids)


def test_queue_only_then_single_close(ledger, targets, tiers, tmp_path):
    cids = [_ingest(ledger, targets, tiers, ["Mg"] + ["H"] * (k + 1), a=3.0 + 0.3 * k) for k in range(4)]
    run_s2(ledger, targets=targets, tiers=tiers, novelty=False)
    for half in (cids[:2], cids[2:]):
        run_s3(ledger, half, targets=targets, tiers=tiers, workdir=tmp_path, alignn=False, dftref=False,
               phonons=False, relax_fn=_fake_relax, hull_fn=_fake_hull, close_pool=False, pool_id="d")
    assert all(ledger.current_outcome(c, "S3")["outcome"] == DEFERRED for c in cids)
    assert all("rank_pending" in ledger.current_outcome(c, "S3")["reason"] for c in cids)
    prom = close_rank_pool(ledger, targets=targets, tiers=tiers, rank_cut=3, pool_id="d")
    assert len(prom) == 3 and all("PROMOTED_FROM_RANK_CUT" in o.tags for o in prom.values())
    assert len(ledger.rank_pools()) == 1


# ------------------------------------------------------------------ metallicity prior (+ parity)
#: campaign 1 S4 outcomes (doc 30 §5.1), primitive-cell formulas: True = S4 metal
CAMPAIGN1_S4 = {"Mg2IrH6": True, "NaIrH3": True, "BaPtH3": True, "LaRhH2": True, "BaRhH6": True,
                "La2Rh2H7": True, "LaRh2H7": True, "Sr2Pt2H8": True, "SrPtH6": False, "BaPtH2": False,
                "Ba2Pt2H8": False, "NaAuH2": False, "LaH2": True, "ScH2": True, "LaH3": True, "YH3": True,
                "Y2H6": True, "Ba2IrH6": True, "BaRh2H7": True, "La4H10": True, "Al2H6": False,
                "BaPtH6": False}


def test_prior_on_campaign1_s4_outcomes():
    """Ranking prior only. The check is in-sample (doc 31 §2): every S4 insulator gets the
    'insulator-likely' level, and 17 of 22 overall are on the correct side of 0.5."""
    right = 0
    for f, metal in CAMPAIGN1_S4.items():
        p = metallicity_prior(parse_formula(f))["p_metal"]
        if not metal:
            assert p < 0.5, f
        right += (p > 0.5) == metal
    assert right >= 17


def test_parity_and_rules():
    assert valence_electrons({"Na": 1, "Ir": 1, "H": 3}) == 13
    assert metallicity_prior({"Na": 1, "Ir": 1, "H": 3})["rule"] == "odd_valence_electron_count"
    assert metallicity_prior({"Mg": 2, "Fe": 1, "H": 6})["rule"] == "closed_shell_18_electron"
    assert metallicity_prior({"Al": 1, "H": 3})["rule"] == "charge_balanced_saline_hydride"
    assert metallicity_prior({"Mg": 1, "B": 2})["p_metal"] == 0.5               # no rule for non-hydrides
    assert metallicity_prior({"Pb": 1})["rule"] == "all_metal_composition"
    # parity needs the cell: from a reduced formula it is not used
    assert metallicity_prior({"La": 1, "H": 2}, cell=False)["rule"] == "excess_electrons"


def test_score_orders_metal_before_isomorphic_insulator(tiers):
    base = {"e_hull_ml_self_eV": 0.05, "phonons": {}}
    ins = {**base, "features": {"H_f": 0.75, "light_fraction": 0.75, "composition": {"Al": 2, "H": 6}}}
    met = {**base, "features": {"H_f": 0.75, "light_fraction": 0.75, "composition": {"Sc": 2, "H": 4}}}
    assert score_s3(met, tiers=tiers) > score_s3(ins, tiers=tiers)
    # ranking only: the S3 decision does not read it
    assert decide_s3({**ins, "features": {**ins["features"], "n_H": 6}}, tiers).outcome == "PASS"


# ------------------------------------------------------------------ B11: chemsys list
def test_shipped_chemsys_list_from_profile(targets):
    from rtt.source.lanes import default_chemsys_list, load_chemsys_list
    hyd, flags = load_chemsys_list(targets)
    assert hyd == default_chemsys_list() and len(hyd) == 70 and flags[0].startswith("CHEMSYS_FILE")
    assert all(any(t in c.split("-") for t in ("Rh", "Pd", "Ir", "Pt", "Os", "Au", "Ag")) for c in hyd)
    nh, _ = load_chemsys_list(targets, section="nonhydride")
    assert "B-C-Li" in nh and all("H" not in c.split("-") for c in nh)


def test_missing_chemsys_list_is_loud(tmp_path, caplog):
    from rtt.config import load_targets
    from rtt.source.lanes import load_chemsys_list, plan_day
    from rtt.config import CONFIG_DIR
    txt = (CONFIG_DIR / "targets.toml").read_text().replace('list = "targets/chemsys_list.txt"',
                                                            'list = "targets/does_not_exist.txt"')
    p = tmp_path / "targets.toml"
    p.write_text(txt)
    t = load_targets(p)
    with caplog.at_level(logging.WARNING, logger="rtt.source.lanes"):
        hyd, flags = load_chemsys_list(t)
    assert flags[0].startswith("DEFAULT_CHEMSYS_FALLBACK") and "does_not_exist" in flags[0]
    assert any("fallback" in r.message for r in caplog.records)
    jobs = [j for j in plan_day(None, t, gpu_h=2.0, day="x") if j.spec.lane == "G1"]
    assert jobs and all("DEFAULT_CHEMSYS_FALLBACK" in j.spec.source_id for j in jobs)


# ------------------------------------------------------------------ R30 F1: domains
def test_s1_size_limit_defers_database_entries(tiers):
    from rtt.records import Origin, RawStructure, StructureData
    from rtt.source.prefilter import prefilter
    n = 21
    sd = StructureData(lattice=[[10, 0, 0], [0, 10, 0], [0, 0, 10]], species=["H"] * (n - 1) + ["Mg"],
                       frac_coords=[[0.04 * i, 0.03 * i, 0.05 * i] for i in range(n)])
    info = {"n_atoms_primitive": n, "elements": ["H", "Mg"], "cost_class": "C3"}
    o = prefilter(RawStructure(structure=sd, origin=Origin(lane="D")), tiers, info=info)
    assert o.outcome == DEFERRED and o.criterion == "cost_class" and o.reason
    o = prefilter(RawStructure(structure=sd, origin=Origin(lane="G1")), tiers, info=info)
    assert o.outcome == "SCREENED_OUT" and o.criterion == "n_atoms_primitive"


def test_s1_deferred_allowed_in_ledger(ledger, targets, tiers):
    cid = _ingest(ledger, targets, tiers, ["Mg", "H", "H"], lane="D")
    from rtt.source.ingest import write_outcome
    write_outcome(ledger, cid, Outcome(stage="S1", outcome=DEFERRED, criterion="cost_class", reason="x"),
                  origin_obj=ledger.candidate_record(cid).origin, targets=targets, tiers=tiers, origin="t")
    assert ledger.current_outcome(cid, "S1")["outcome"] == DEFERRED


def test_s3_window_does_not_remove_non_hydrides(tiers):
    thr = float(tiers.threshold("S3.ehull_self_max_eV").value)
    o = decide_s3({"e_hull_ml_self_eV": thr + 0.03, "features": {"n_H": 0}}, tiers)
    assert o.outcome == DEFERRED and o.criterion == "uncalibrated_window" and o.provisional
    from rtt.engine.retry import is_retryable
    assert not is_retryable({"outcome": DEFERRED, "criterion": o.criterion, "reason": o.reason})
    assert decide_s3({"e_hull_ml_self_eV": thr + 0.03, "features": {"n_H": 3}}, tiers).outcome == "SCREENED_OUT"


def test_s2_database_copy_is_representative(ledger, targets, tiers):
    """R30 F3: a generated copy of a database entry is the DUPLICATE, tagged as a rediscovery."""
    g = _ingest(ledger, targets, tiers, ["Mg", "H", "H"], lane="G1")
    d = _ingest(ledger, targets, tiers, ["Mg", "H", "H"], lane="D")
    out = run_s2(ledger, targets=targets, tiers=tiers, novelty=False)
    assert out[d].outcome == PASS
    assert out[g].outcome == "DUPLICATE" and out[g].value == d and "GENERATED_REDISCOVERY" in out[g].tags
