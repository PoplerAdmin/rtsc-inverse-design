"""Lane P (prototype substitution, doc 32): electron count, prototypes, enumeration, S1 ingest."""
from __future__ import annotations

from fractions import Fraction

import pytest

from rtt.records import Origin
from rtt.source import prototype as P


def test_electron_count_rule():
    # closed shell 18: Mg2FeH6, K2PtH6, SrPtH6-like [PtH6]2- with a 2+ cation
    assert P.complex_count(("Mg", "Mg"), ("Fe",)) == 18
    assert P.complex_count(("K", "K"), ("Pt",)) == 18
    # Mg2IrH6 control: 19 = one electron in the eg sigma* band
    assert P.complex_count(("Mg", "Mg"), ("Ir",)) == 19
    f = P.filling(Fraction(19))
    assert f["class"] == "eg_sigma_star_electrons" and f["eg_filling_fraction"] == 0.25
    # isoelectronic alloys (rock-salt A order) and half-electron alloys
    assert P.complex_count(("Na", "Mg"), ("Pt",)) == 19
    assert P.complex_count(("Mg", "Sc"), ("Ru",)) == 19
    assert P.complex_count(("Mg",) * 3 + ("Li",), ("Ir", "Ir")) == Fraction(37, 2)
    assert P.complex_count(("Mg",) * 4, ("Ir", "Pt")) == Fraction(39, 2)
    assert P.filling(Fraction(17))["class"] == "t2g_holes"


def test_parity_matches_campaign1_s4():
    # R30 F1(c): odd cell counts were metals, even-count 18e complexes were insulators
    assert P.cell_valence(("Mg", "Mg"), ("Ir",), 6) == 19          # metal
    assert P.cell_valence(("Ba", "Ba"), ("Ir",), 6) == 19          # metal
    assert P.cell_valence(("Mg", "Mg"), ("Fe",), 6) % 2 == 0       # Mg2FeH6, 18 e insulator


@pytest.mark.parametrize("builder,args,n,sg", [
    (P.build_a2mh6, ("Mg", "Ir"), 9, 225),
    (P.build_aamh6, ("Na", "Mg", "Pt"), 9, 216),
    (P.build_a3a_m2h12, ("Mg", "Li", "Ir"), 18, None),
    (P.build_a2mm_h12, ("Mg", "Ir", "Pt"), 18, None),
])
def test_prototype_cells(builder, args, n, sg):
    from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
    st = builder(*args)
    assert len(st) == n
    comp = st.composition
    assert comp["H"] == 6 * (n // 9)
    if sg is not None:
        assert SpacegroupAnalyzer(st, symprec=0.01).get_space_group_number() == sg
    # no unphysical contacts: M-H about 1.6-1.7 A, H-H > 2 A
    import numpy as np
    d = st.distance_matrix
    H = [i for i, s in enumerate(st) if s.specie.symbol == "H"]
    assert min(d[i, j] for i in H for j in H if i != j) > 2.0
    assert 1.5 < float(np.min(d[np.nonzero(d)])) < 1.8


def test_mg2irh6_geometry_guess():
    a, x = P.lattice_guess(("Mg", "Mg"), ("Ir",))
    assert abs(a - 6.65) < 0.03 and abs(x - 0.2556) < 0.003


def test_enumeration_targets_and_determinism():
    c1 = P.enumerate_lane_p()
    c2 = P.enumerate_lane_p()
    assert [c.priority for c in c1] == [c.priority for c in c2]
    assert all(c.N in P.DEFAULT_N_TARGETS for c in c1)
    assert all(c.N != 18 for c in c1)
    assert all("Tc" not in c.m_sites for c in c1)
    first = c1[0]
    assert first.formula == "Mg2IrH6" and first.prototype == "A2MH6"   # the control's filling first
    alloys = [c for c in c1 if c.prototype != "A2MH6"]
    assert all(Fraction(37, 2) <= c.N <= Fraction(39, 2) for c in alloys)
    with pytest.raises(ValueError):
        P.enumerate_lane_p(M=("Tc",))


def test_config_plan_and_provenance():
    cfg = P.load_lane_p_config()
    plan = P.plan_lane_p(cfg)
    caps = {k: v for k, v in cfg.get("caps").items() if k != "source"}
    for proto, cap in caps.items():
        assert sum(c.prototype == proto for c in plan) <= cap
    raw = P.to_raw(plan[0])
    assert raw.origin.lane == "P" and isinstance(raw.origin, Origin)
    cond = raw.origin.conditioning
    assert cond["prototype"] == "A2MH6" and cond["substitution"] == {"A": "Mg", "M": "Ir"}
    assert cond["N_per_complex"] == 19.0 and cond["valence_parity"] == "odd"
    assert "PROTOTYPE_SUBSTITUTION" in raw.origin.flags


def test_ingest_lane_p(tmp_path):
    from rtt.config import load_targets, load_tiers
    from rtt.ledger.store import Ledger
    t, ti = load_targets(), load_tiers()
    L = Ledger(tmp_path / "l.sqlite", run_id="t")
    cands = [c for c in P.enumerate_lane_p() if set(c.a_sites) | set(c.m_sites) <= {"Mg", "Ir", "Pd", "Li", "Na", "Pt"}][:3]
    cids = P.run_lane_p(L, targets=t, tiers=ti, cands=cands)
    assert len(cids) == 3
    for cid in cids:
        rec = L.candidate_record(cid)
        assert rec.origin.lane == "P" and rec.origin.conditioning["prototype"]
        assert L.current_outcome(cid, "S1") is not None
