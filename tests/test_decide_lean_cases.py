"""decide_t against every case exported by Lean (lean/export/decide_cases.json, decideT at K = Q):
outcome, rule, assumptions, inputs-used provenance, gate results, rule floors/provenance and the
derived envelope numbers must all agree exactly."""
from __future__ import annotations

import json
from fractions import Fraction
from pathlib import Path

import pytest

from rtt.decide import core as C
from rtt.decide.lean_io import EXPORT_INPUT_NAMES, config_from_export, inputs_from_export

ROOT = Path(__file__).resolve().parents[1]
CASES = ROOT / "lean" / "export" / "decide_cases.json"
DATA = json.loads(CASES.read_text(encoding="utf-8")) if CASES.exists() else {"cases": [], "config": None}
CFG = config_from_export(DATA["config"]) if DATA["config"] else None
pytestmark = pytest.mark.skipif(not CASES.exists(), reason="lean/export/decide_cases.json missing")


def test_case_count():
    assert len(DATA["cases"]) == DATA["n_cases"] >= 561


@pytest.mark.parametrize("case", DATA["cases"], ids=[c["id"] for c in DATA["cases"]])
def test_lean_case(case):
    x = inputs_from_export(case["inputs"])
    exp = case["expected"]
    o = C.decide_t(CFG, x)
    assert o.outcome == exp["lean_outcome"]
    assert (o.outcome == "screenedOut") == bool(exp["would_fire_ranking_only"])
    if exp["lean_outcome"] == "rejectedProven":
        assert o.rule == exp["rule"]
        assert o.assumptions == exp["assumptions"]
        used = [(u["name"], u["prov"]) for u in o.inputs_used]
        want = [(EXPORT_INPUT_NAMES[u["name"]], {"tier": {"ML": "ml"}.get(u["prov"]["tier"], u["prov"]["tier"]),
                                                 "run": u["prov"]["run"]}) for u in exp["inputs_used"]]
        assert used == want
    else:
        assert o.assumptions == exp.get("assumptions", []) == []
    # gates and per-rule booleans
    X = C.exact_view(x)
    g = case["gates"]
    assert C.wf(X, CFG) == g["well_formed"]
    assert C.blocked(x.tags) == g["blocked"]
    if g["well_formed"]:
        assert C.gvalid(X) == g["g_valid"]
        for r, d in case["rules"].items():
            assert C.floor_r(CFG, X, r) == d["floor"], r
            assert C.prov_ok(X, r) == d["prov_ok"], r
        der = case["derived"]
        q = lambda s: None if s is None else Fraction(s)   # noqa: E731
        assert C.Ahi(X) == q(der["A_hi"]) and C.lamHi(X) == q(der["lambda_hi"]) and C.mHi(X) == q(der["m_hi"])
        assert CFG.floor == q(der["floor_A"])
        assert 4 * CFG.floor ** 2 / C.lamHi(X) == q(der["nm_cert_floor"])
        assert CFG.gUp(C.lamHi(X)) == q(der["g_up"])
        if der.get("kappa_hi") is not None:
            assert C.kappaHi(X.mat, X) == q(der["kappa_hi"])
            assert C.rowMaxHi(X.mat, X) == q(der["rowmax_hi"])
        if der.get("wmin_lo") is not None:
            assert C.wminLo(X.dos) == q(der["wmin_lo"])


def test_export_config_equals_rtt_config():
    """The Lean export's a0 and certified NM table equal what rtt loads from targets.toml and
    docs/theory/certificates/nm_*.json."""
    from rtt.targets.spec import load_spec
    mine = C.DecideConfig.from_spec(load_spec())
    assert (mine.a0, mine.Tstar, mine.tbl, mine.gInf) == (CFG.a0, CFG.Tstar, CFG.tbl, CFG.gInf)
