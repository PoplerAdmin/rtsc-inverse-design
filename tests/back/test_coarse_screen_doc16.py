"""Doc 16 (coarse-screen fix): all-sigma statistic, sigma -> 0 extrapolation, H cutoff floor,
re-decision of stored S5 outcomes. Ladder numbers are the campaign-1 Mg2IrH6 coarse runs
(60 Ry; ~/rtsc_targeted_work/E1-campaign/work/c-9a8fa100a45b/S5, doc 16 §2)."""
from __future__ import annotations

import json
from fractions import Fraction

import pytest

from rtt.config import Threshold, load_targets, load_tiers, load_machine, validate_sources
from rtt.dft import coarse as C
from rtt.dft.profiles import COARSE, apply_ecut_floor, ecut_floor_Ry

FLOOR = Fraction(433148145, 100000000) * 300            # 1299.44 K

# campaign-1 Mg2IrH6, sigma = 0.01..0.10 Ry: A (K), omega_log (K), admitted
Q4_A = [904.7, 640.7, 518.3, 439.3, 385.8, 351.4, 330.6, 318.2, 310.0, 303.6]
Q4_W = [575.3, 574.9, 574.2, 573.8, 575.4, 579.8, 585.6, 591.4, 596.3, 599.9]
Q3_A = [888.3, 670.6, 551.5, 468.6, 411.6, 374.8, 352.7, 339.7, 331.3, 325.1]
Q3_W = [602.4, 608.9, 610.1, 610.2, 612.3, 617.9, 625.7, 633.5, 640.1, 644.9]
Q3_ADM = [False] * 5 + [True] * 5


def _m(As, ws, adm, q, elements=("H", "Ir", "Mg"), ecut=80.0):
    rows = [C.SigmaRow(index=i + 1, sigma_Ry=0.01 * (i + 1), admitted=a, A_K=A, omega_log_K=w,
                       on_plateau=a and 0.01 * (i + 1) >= 0.03 - 1e-12)
            for i, (A, w, a) in enumerate(zip(As, ws, adm))]
    return C.CoarseMetrics(rows=rows, plateau_sigma_min_Ry=0.03, min_freq_matdyn_cm1=0.0, q=[q] * 3,
                           elements=list(elements), ecutwfc_Ry=ecut)


def _th():
    return dict(eps_sys=Threshold(0.10, "doc01#0.2b"), e_c=Threshold(0.5, "doc16#5", provisional=True),
                margin_M_c=Threshold(2.0, "doc13#5", provisional=True))


FLOORS = {"H": 80.0, "source": "doc15#1", "provisional": True}


def test_a_top_extrapolates_and_keeps_excluded_rows():
    q3 = _m(Q3_A, Q3_W, Q3_ADM, 3)
    # excluded sigma 1-5 are used: 2*888.3 - 670.6 = 1106.0
    assert C.a_top(q3) == pytest.approx(2 * 888.3 - 670.6)
    flat = _m([66.0, 66.0, 66.0], [300.0] * 3, [True] * 3, 4)
    assert C.a_top(flat) == 66.0                       # flat ladder: extrapolation = value
    rising = _m([36.0, 37.0, 38.0], [300.0] * 3, [True] * 3, 4)
    assert C.a_top(rising) == 38.0                     # max dominates when A rises with sigma


def test_mg2irh6_doc13_statistic_undercovers_and_doc16_does_not():
    """doc 16 §2: at A_target = the reference's A_hi lower bound (1822 K, P1 80 Ry), the doc 13
    statistic would remove Mg2IrH6, the all-sigma statistic keeps it with margin."""
    fine, coarse = _m(Q4_A, Q4_W, [True] * 10, 4), _m(Q3_A, Q3_W, Q3_ADM, 3)
    old = C.coarse_envelope(fine, coarse, 0.10)
    assert 1.5 * old["A_hi_K"] < 1822.0                              # 1639 K: unsound
    new = C.coarse_envelope_all(fine, coarse, 0.10)
    assert new["A_nom_K"] == pytest.approx(2 * 904.7 - 640.7)
    assert 1.5 * new["A_hi_K"] > 1.9 * 1822.0                         # 3558 K
    o = C.decide_coarse(fine, metrics_coarser=coarse, floor_K=Fraction(1822), n_sheets=1,
                        ecut_floors=FLOORS, **_th())
    assert o.outcome == "PASS"


def test_counterfactual_admission_cut_on_fine_mesh():
    """If admission also cut sigma 1-5 on q 4^3 (as it did on q 3^3), the doc 13 statistic falls
    to ~0.41 x the reference A_hi; the doc 16 statistic is unchanged (all rows are used)."""
    fine = _m(Q4_A, Q4_W, Q3_ADM, 4)
    coarse = _m(Q3_A, Q3_W, Q3_ADM, 3)
    old = C.coarse_envelope(fine, coarse, 0.10)
    assert 1.5 * old["A_hi_K"] / 1822.0 < 0.45
    full = C.coarse_envelope_all(_m(Q4_A, Q4_W, [True] * 10, 4), coarse, 0.10)
    assert C.coarse_envelope_all(fine, coarse, 0.10)["A_hi_K"] == pytest.approx(full["A_hi_K"])


def test_ecut_floor_blocks_hydride_removal():
    low = [40.0, 38.0, 36.0, 35.0, 34.0]
    w = [300.0] * 5
    for ecut, want in ((60.0, "DEFERRED"), (80.0, "SCREENED_OUT"), (None, "DEFERRED")):
        f, c = _m(low, w, [True] * 5, 4, ecut=ecut), _m(low, w, [True] * 5, 3, ecut=ecut)
        o = C.decide_coarse(f, metrics_coarser=c, floor_K=FLOOR, n_sheets=1, ecut_floors=FLOORS, **_th())
        assert o.outcome == want, (ecut, o)
        if want == "DEFERRED":
            assert "BLOCKED:ecut_below_floor" in o.tags
    # no H: no floor applies at 60 Ry
    f, c = (_m(low, w, [True] * 5, q, elements=("Nb",), ecut=60.0) for q in (4, 3))
    assert C.decide_coarse(f, metrics_coarser=c, floor_K=FLOOR, n_sheets=1, ecut_floors=FLOORS,
                           **_th()).outcome == "SCREENED_OUT"
    # elements unknown with a floor table -> blocked
    f, c = (_m(low, w, [True] * 5, q, elements=(), ecut=80.0) for q in (4, 3))
    o = C.decide_coarse(f, metrics_coarser=c, floor_K=FLOOR, n_sheets=1, ecut_floors=FLOORS, **_th())
    assert o.outcome == "DEFERRED" and "BLOCKED:ecut_unknown" in o.tags


def test_ecut_floor_profile_helpers_and_tiers():
    assert ecut_floor_Ry(["H", "Mg"], FLOORS) == 80.0
    assert ecut_floor_Ry(["Nb"], FLOORS) is None and ecut_floor_Ry(["H"], None) is None
    assert apply_ecut_floor(COARSE, ["Mg", "Ir", "H"], FLOORS).ecutwfc_Ry == 80.0
    assert ("ecutwfc_Ry", 80.0) in apply_ecut_floor(COARSE, ["H"], FLOORS).overrides
    assert apply_ecut_floor(COARSE, ["Al"], FLOORS) is COARSE
    hi = COARSE.with_overrides(ecutwfc_Ry=100.0)
    assert apply_ecut_floor(hi, ["H"], FLOORS).ecutwfc_Ry == 100.0         # never lowers
    t = load_tiers()
    validate_sources(t.raw)
    tab = t.get("S5.ecut_floor_Ry")
    assert tab["H"] == 80.0 and tab.get("source") and tab.get("provisional") is True


def test_from_metrics_roundtrip():
    m = _m(Q4_A, Q4_W, [True] * 10, 4)
    d = json.loads(json.dumps(m.to_metrics()))
    r = C.CoarseMetrics.from_metrics(d)
    assert [x.A_K for x in r.rows] == Q4_A and r.ecutwfc_Ry == 80.0 and r.elements == ["H", "Ir", "Mg"]
    # pre-doc-16 record: no ecut / elements stored -> filled from arguments
    for k in ("ecutwfc_Ry", "elements", "admission_ok"):
        d.pop(k)
    r = C.CoarseMetrics.from_metrics(d, elements=["H", "La"], ecutwfc_Ry=60.0)
    assert r.ecutwfc_Ry == 60.0 and r.elements == ["H", "La"] and r.admission_ok


def test_redecide_reopens_old_screened_out(tmp_path):
    """A campaign-1 style SCREENED_OUT (60 Ry, H) re-decided under doc 16 becomes
    DEFERRED(ecut_below_floor): appended, the old row stays (append-only)."""
    from rtt.engine.stage import Context
    from rtt.engine.stages import S5Coarse, write_outcome
    from rtt.ledger.store import Ledger
    from rtt.records import Origin, Outcome, StructureData
    L = Ledger(tmp_path / "l.sqlite")
    t, ti, mc = load_targets(), load_tiers(), load_machine()
    L.register_configs(t, ti, mc)
    ctx = Context(ledger=L, targets=t, tiers=ti, machine=mc, workdir=tmp_path / "work")
    st = StructureData(lattice=[[0, 2.8, 2.8], [2.8, 0, 2.8], [2.8, 2.8, 0]], species=["La", "H", "H"],
                       frac_coords=[[0, 0, 0], [0.25, 0.25, 0.25], [0.75, 0.75, 0.75]])
    L.insert_candidate("lah2", "a" * 64, st, Origin(lane="D"), cost_class="C1")
    low = [39.2, 37.0, 36.0, 35.3, 34.6, 33.6, 32.4, 31.6, 32.4, 36.9]
    adm = [False] * 5 + [True] * 5
    meshes = [_m(low, [224.0] * 10, adm, q, elements=(), ecut=None) for q in (3, 4)]
    by_mesh = []
    for q, m in zip((3, 4), meshes):
        d = m.to_metrics()
        for k in ("ecutwfc_Ry", "elements", "admission_ok"):   # as stored before doc 16
            d.pop(k)
        by_mesh.append({"q": [q] * 3, **d})
    trs = [L.start_tier_run("lah2", "S5", "dfpt", profile="coarse", params={"q": [q] * 3}) for q in (3, 4)]
    L.finish_tier_run(trs[0], "done", metrics={"q": [3, 3, 3]})
    L.finish_tier_run(trs[1], "done", metrics={**by_mesh[-1], "by_mesh": by_mesh, "n_sheets": 1})
    # the run's own record gives the cutoff (60 Ry)
    for q in (3, 4):
        wd = tmp_path / "work" / "lah2" / "S5" / f"q{q}x{q}x{q}"
        wd.mkdir(parents=True)
        (wd / "dfpt_run.json").write_text(json.dumps({"settings": {"ecutwfc_Ry": 60.0}}))
    write_outcome(ctx, "lah2", Outcome("S5", "SCREENED_OUT", criterion="coarse-A envelope", value=65.2,
                                       threshold=1299.44, source_id="doc13#5", tier_run_id=trs[1]))
    rec = L.candidate_record("lah2")
    o = S5Coarse().redecide(rec, ctx)
    assert o.outcome == "DEFERRED" and o.criterion == "ecut_below_floor"
    assert "REDECIDED_DOC16" in o.tags and "supersedes outcome" in o.reason
    write_outcome(ctx, "lah2", o, rec=rec)
    assert L.current_outcome("lah2", "S5")["outcome"] == "DEFERRED"
    assert len(L.qa("SELECT id FROM outcomes WHERE candidate_id='lah2' AND stage='S5'")) == 2
    from rtt.engine.retry import is_retryable
    assert is_retryable(L.current_outcome("lah2", "S5"))          # S5 is re-run (at 80 Ry)
    L.close()
