"""rtt.report.evidence (§3.6 template: every item checked or NOT CHECKED) and summary."""
from __future__ import annotations

import pytest

from rtt.report import evidence as EV
from rtt.report.summary import control_recall, daily_summary


def test_empty_candidate_everything_not_checked():
    txt = EV.render(EV.EvidenceInputs(candidate_id="c-x"))
    for item in EV.TEMPLATE_ITEMS:
        assert f"**{item}**" in txt, item
    lines = [l for l in txt.splitlines() if l.startswith("- **")]
    assert len(lines) == len(EV.TEMPLATE_ITEMS)
    nc = [l for l in lines if "NOT CHECKED" in l]
    # only "provisional margins" can read checked (none used) with no data
    assert len(nc) == len(lines) - 1


def test_checked_items_carry_values_and_ids():
    ev = EV.EvidenceInputs(
        candidate_id="c-1", formula="Mg2IrH6", structure={"spacegroup": 225, "n_atoms": 9},
        structure_sha256="ab" * 32, origin={"lane": "C", "control": "Mg2IrH6"}, pressure_GPa=0.0,
        stages={"S4": {"outcome": "PASS", "tier_run_id": 140, "metrics": {"gap_eV": 0.0, "moment_muB": 0.01, "n_sheets": 2}},
                "S5": {"outcome": "PASS", "tier_run_id": 188, "provisional": True,
                       "tags": ["PROVISIONAL_MARGIN"], "metrics": {}},
                "S6": {"outcome": "undecided", "tier_run_id": 230, "metrics": {"min_freq_cm1": 85.0, "A_mid_K": 900.0}}},
        decisions=[{"id": 9, "stage": "S6", "tier": "fullDFPT", "outcome": "undecided",
                    "firing_rule": "G-VALID", "decide_input_id": 90}],
        decide_inputs={"validity": {"EF": 0.35, "omegaMax": 2500.0, "maxImagModeK": 0.0,
                                    "containsHLiBe": True, "anharmonicCheck": False,
                                    "muStarNegativeSuspected": False, "admittedSmearings": 7,
                                    "qMeshes": 2, "nSheets": 2, "bandResolvedAvailable": False,
                                    "dosFractionsKnown": True},
                       "envelope": {"Anom": 900.0, "lamNom": 1.5, "mNom": 1e6, "dA": 0.1,
                                    "dLam": 0.2, "dM": 0.1, "epsSys": 0.1, "epsOmega": 0.02,
                                    "tier": "fullDFPT"}},
        verdicts=[{"id": 61, "spectrum_tier": "full", "status": "UNDECIDED", "T_low_K": 40.0,
                   "T_up_K": 400.0, "mu_hi": 0.13, "rule_ids": "[]", "assumptions": "[]",
                   "flags": '["ANHARMONIC_H"]', "recheck_ok": 1}],
        control_recall={"S4": {"passed": 1, "total": 1, "removed": 0}})
    d = dict(EV.items(ev))
    assert d["H1 metallic (S4)"].startswith("checked (gap_eV = 0")
    assert "tier_run #140" in d["H1 metallic (S4)"]
    assert d["G-VALID 2 stable/anharmonic"].startswith("checked")
    assert "anharmonicCheck=False" in d["G-VALID 2 stable/anharmonic"]
    assert d["G-VALID 4 convergence data"].startswith("checked")
    assert d["core verdict"].startswith("checked (UNDECIDED [40, 400] K")
    assert "verdict #61" in d["core verdict"]
    assert "NOT CHECKED at [0.1, 0.15]" in d["mu* sensitivity"]
    assert d["anharmonic check"].startswith("NOT CHECKED")
    assert d["ZPE-corrected E_hull (H7)"].startswith("NOT CHECKED")
    assert "S5: PASS" in d["provisional margins used upstream"]
    assert d["decide outcome and firing rule"].startswith("checked (undecided at tier fullDFPT")


def test_missing_validity_field_is_not_checked():
    ev = EV.EvidenceInputs(candidate_id="c", decisions=[{"id": 1, "decide_input_id": 2}],
                           decide_inputs={"validity": {"EF": None, "omegaMax": 10.0}, "envelope": None})
    d = dict(EV.items(ev))
    assert d["G-VALID 1 non-adiabatic"].startswith("NOT CHECKED (missing: EF")


def test_gather_from_ledger_and_summary():
    from rtt.config import load_machine, load_targets, load_tiers
    from rtt.ledger.store import Ledger
    from rtt.records import Origin, Outcome
    t, ti, m = load_targets(), load_tiers(), load_machine()
    L = Ledger(":memory:")
    L.register_configs(t, ti, m)
    st = {"lattice": [[1, 0, 0], [0, 1, 0], [0, 0, 1]], "species": ["Al"], "frac_coords": [[0, 0, 0]]}
    L.insert_candidate("c-al", "a" * 64, st, Origin(lane="C", control="Al", T_star_surrogate_K=0.6))
    tr = L.start_tier_run("c-al", "S4", "dft", profile="s4")
    L.finish_tier_run(tr, "done", wall_s=10.0, core_h=0.01,
                      metrics={"gap_eV": 0.0, "moment_muB": 0.0, "n_sheets": 1})
    L.insert_outcome("c-al", Outcome(stage="S4", outcome="PASS", tier_run_id=tr),
                     targets_sha256=t.sha256, tiers_sha256=ti.sha256, origin="test")
    tr5 = L.start_tier_run("c-al", "S5", "dfpt", profile="coarse")
    L.finish_tier_run(tr5, "done", wall_s=60.0, core_h=0.03, metrics={"A_K_by_sigma": [60.0]})
    L.insert_outcome("c-al", Outcome(stage="S5", outcome="SCREENED_OUT", criterion="coarse-A",
                                     value=66.0, threshold=649.7, source_id="doc11#4",
                                     provisional=True, tier_run_id=tr5, tags=["PROVISIONAL_MARGIN"]),
                     targets_sha256=t.sha256, tiers_sha256=ti.sha256, origin="test")
    txt = EV.evidence_statement("c-al", L)
    assert "**H1 metallic (S4)**: checked" in txt
    assert "S5: SCREENED_OUT (coarse-A)" in txt
    assert "**core verdict**: NOT CHECKED" in txt
    cr = control_recall(L)
    assert cr["S5"] == {"passed": 0, "total": 1, "removed": 1}
    s = daily_summary(L)
    assert "| S5 | SCREENED_OUT | 1 |" in s and "recall < 1" in s
    assert "coarse-A (doc11#4) PROVISIONAL_MARGIN" in s


def test_novelty_from_s2_outcome_value_and_s5_line():
    """Campaign 1: novelty lives in the S2 outcome value; the S5 coarse A is shown in the header."""
    ev = EV.EvidenceInputs(candidate_id="c-x", formula="H2LaRh",
                           stages={"S2": {"outcome": "PASS", "metrics": {"novel": True, "alexmp_match": None}},
                                   "S5": {"outcome": "DEFERRED", "criterion": "multiband", "tier_run_id": 7,
                                          "metrics": {"A_mid_K": 210.0, "A_lo_K": 190.0}}})
    txt = EV.render(ev)
    assert "novel=True" in txt
    assert "A_coarse,mid = 210 K" in txt and "DEFERRED (multiband)" in txt
