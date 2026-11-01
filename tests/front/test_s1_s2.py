"""S1 prefilter / controls / ingestion and S2 validity / dedupe (CPU)."""
from __future__ import annotations

import numpy as np
import pytest

from rtt.records import ALLOWED_OUTCOMES, Origin, RawStructure, StructureData
from rtt.source.controls import CONTROLS, control, control_raws, inject_controls
from rtt.source.ingest import ingest_raw
from rtt.source.prefilter import prefilter, primitive_info
from rtt.structure.canonical import canonicalize, structure_sha256, to_pymatgen
from rtt.structure.dedupe import DedupeIndex, is_novel
from rtt.structure.stage import run_s2
from rtt.structure.validate import validate


def _raw(sd: dict, lane="G5", p=0.0) -> RawStructure:
    return RawStructure(structure=StructureData(**sd), origin=Origin(lane=lane, flags=["HYPOTHETICAL_STRUCTURE"]),
                        pressure_GPa=p)


def _cubic(a, species, frac):
    return {"lattice": [[a, 0, 0], [0, a, 0], [0, 0, a]], "species": species, "frac_coords": frac}


# ------------------------------------------------------------------ controls
def test_controls_structures():
    names = [c.name for c in CONTROLS]
    assert names == ["Mg2IrH6", "H3S", "MgB2", "Nb", "Pb", "Al"]
    st = to_pymatgen(control("Mg2IrH6").raw().structure.to_json())
    assert dict(st.composition.reduced_composition.as_dict()) == {"Mg": 2, "Ir": 1, "H": 6}
    info = primitive_info(control("Mg2IrH6").raw())
    assert info["spacegroup"] == 225 and info["n_atoms_primitive"] == 9 and info["cost_class"] == "C2"
    ir = [i for i, s in enumerate(st) if s.specie.symbol == "Ir"][0]
    dIrH = min(st.get_distance(ir, j) for j, s in enumerate(st) if s.specie.symbol == "H")
    assert abs(dIrH - 1.7288) < 1e-3
    h3s = control("H3S").raw()
    assert h3s.pressure_GPa == 200.0 and h3s.origin.pressure_override_GPa == 200.0
    assert "CALIBRATION_ONLY" in h3s.origin.flags
    assert "BLOCKED_PH_CRASH" in control("MgB2").raw().origin.flags
    for r in control_raws():
        assert r.origin.lane == "C" and r.origin.T_star_surrogate_K > 0


def test_inject_controls_idempotent_and_pass_s1(ledger, targets, tiers):
    ids = inject_controls(ledger, targets=targets, tiers=tiers)
    assert len(ids) == 6
    assert inject_controls(ledger, targets=targets, tiers=tiers) == ids
    assert ledger.q1("SELECT COUNT(*) FROM candidates")[0] == 6
    for cid in ids:
        o = ledger.current_outcome(cid, "S1")
        assert o["outcome"] == "PASS"
        # recorded against the surrogate target, registered in inputs (R7 trigger accepted it)
        assert o["targets_sha256"] != targets.sha256
    rec = ledger.candidate_record(ids[1])
    assert rec.pressure_GPa == 200.0 and rec.structure.spacegroup == 229


# ------------------------------------------------------------------ prefilter
def test_prefilter_criteria(tiers):
    ok = _raw(_cubic(3.0, ["Mg", "H"], [[0, 0, 0], [0.5, 0.5, 0.5]]))
    assert prefilter(ok, tiers).outcome == "PASS"
    nop = _raw(_cubic(4.0, ["Cs", "H"], [[0, 0, 0], [0.5, 0.5, 0.5]]))
    o = prefilter(nop, tiers)
    assert o.outcome == "SCREENED_OUT" and o.criterion == "no_pseudo" and o.value == ["Cs"] and o.source_id
    xe = _raw(_cubic(4.0, ["Xe", "H"], [[0, 0, 0], [0.5, 0.5, 0.5]]))
    o = prefilter(xe, tiers)
    assert o.criterion == "outside_mattergen_domain"
    rng = np.random.default_rng(0)
    big = _raw(_cubic(9.0, ["Mg"] * 7 + ["H"] * 14 + ["Li"], rng.random((22, 3)).tolist()))
    o = prefilter(big, tiers)
    assert o.outcome == "SCREENED_OUT" and o.criterion == "n_atoms_primitive" and o.threshold == 20
    for oo in (o,):
        assert oo.outcome in ALLOWED_OUTCOMES["S1"]


def test_supercell_is_primitive_for_size_rule(tiers):
    from pymatgen.core import Lattice, Structure
    st = Structure.from_spacegroup("Fm-3m", Lattice.cubic(4.05), ["Al"], [[0, 0, 0]]) * (3, 3, 3)   # 108 atoms
    sd = {"lattice": st.lattice.matrix.tolist(), "species": [s.specie.symbol for s in st],
          "frac_coords": st.frac_coords.tolist()}
    assert prefilter(_raw(sd), tiers).outcome == "PASS"


# ------------------------------------------------------------------ canonical / validate / dedupe
def test_hash_invariant_under_supercell_and_permutation():
    from pymatgen.core import Lattice, Structure
    a = Structure.from_spacegroup("Pm-3m", Lattice.cubic(3.2), ["Mg", "H"], [[0, 0, 0], [0.5, 0.5, 0.5]])
    b = a.copy()
    b.make_supercell([2, 1, 1])
    c = Structure(b.lattice, list(reversed([s.specie for s in b])), list(reversed(b.frac_coords)))
    c.translate_sites(list(range(len(c))), [0.1, 0.2, 0.3])
    h = {structure_sha256(canonicalize(x)) for x in (a, b, c)}
    assert len(h) == 1


def test_validate():
    assert validate(_cubic(3.0, ["Mg", "H"], [[0, 0, 0], [0.5, 0.5, 0.5]])).outcome == "PASS"
    o = validate(_cubic(3.0, ["Mg", "H"], [[0, 0, 0], [0.05, 0, 0]]))
    assert o.outcome == "INVALID" and o.criterion == "min_distance_A"
    o = validate(_cubic(30.0, ["Mg"], [[0, 0, 0]]))
    assert o.outcome == "INVALID" and o.criterion == "volume_per_atom_A3"
    o = validate({"lattice": [[1, 0, 0], [0, 1, 0], [0, 0, 0]], "species": ["H"], "frac_coords": [[0, 0, 0]]})
    assert o.outcome == "INVALID" and o.criterion == "finite"
    # H2 molecule distance is valid (only penalised in ranking)
    assert validate(_cubic(5.0, ["H", "H"], [[0, 0, 0], [0.148, 0, 0]])).outcome == "PASS"


def test_dedupe_index_matches_strained_copy():
    base = _cubic(3.2, ["Mg", "H"], [[0, 0, 0], [0.5, 0.5, 0.5]])
    idx = DedupeIndex()
    idx.add("c-1", base, 0.0)
    strained = _cubic(3.25, ["H", "Mg"], [[0.51, 0.5, 0.5], [0, 0, 0]])
    assert idx.find(strained, 0.0) == "c-1"
    assert idx.find(strained, 10.0) is None                       # other pressure
    other = _cubic(3.2, ["Mg", "H", "H"], [[0, 0, 0], [0.5, 0.5, 0.5], [0.5, 0, 0]])
    assert idx.find(other, 0.0) is None
    nov, rid = is_novel(strained, [("ref-1", base)])
    assert not nov and rid == "ref-1"


def test_run_s2_duplicates_and_controls(ledger, targets, tiers):
    cids = inject_controls(ledger, targets=targets, tiers=tiers)
    # a generated copy of fcc Al (conventional cell, slightly strained) and an invalid structure
    al = _cubic(4.07, ["Al"] * 4, [[0, 0, 0], [0.5, 0.5, 0], [0.5, 0, 0.5], [0, 0.5, 0.5]])
    c_al = ingest_raw(ledger, _raw(al), targets=targets, tiers=tiers)
    bad = _cubic(3.0, ["Mg", "H"], [[0, 0, 0], [0.02, 0, 0]])
    c_bad = ingest_raw(ledger, _raw(bad), targets=targets, tiers=tiers)
    # a second generated copy (other lane, slightly different strain) is a real duplicate of c_al
    al2 = _cubic(4.08, ["Al"] * 4, [[0, 0, 0], [0.5, 0.5, 0], [0.5, 0, 0.5], [0, 0.5, 0.5]])
    c_al2 = ingest_raw(ledger, _raw(al2, lane="G2"), targets=targets, tiers=tiers)
    out = run_s2(ledger, targets=targets, tiers=tiers, novelty=False)
    assert all(out[c].outcome == "PASS" for c in cids)
    # RC1 F16: a generated copy of a control is linked, not dropped (evaluated at the real T*)
    assert out[c_al].outcome == "PASS" and "SAME_AS_CONTROL" in out[c_al].tags
    assert out[c_al].value["same_structure_as_control"] == cids[-1]
    assert out[c_al2].outcome == "DUPLICATE" and out[c_al2].value == c_al
    assert out[c_bad].outcome == "INVALID"
    assert ledger.q1("SELECT duplicate_of FROM duplicates WHERE candidate_id=?", (c_al2,))[0] == c_al
    assert run_s2(ledger, targets=targets, tiers=tiers, novelty=False) == {}       # nothing left to do
    for o in out.values():
        assert o.outcome in ALLOWED_OUTCOMES["S2"]


def test_s1_screened_candidate_skips_s2(ledger, targets, tiers):
    c = ingest_raw(ledger, _raw(_cubic(4.0, ["Cs", "H"], [[0, 0, 0], [0.5, 0.5, 0.5]])), targets=targets, tiers=tiers)
    assert ledger.current_outcome(c, "S1")["outcome"] == "SCREENED_OUT"
    assert c not in run_s2(ledger, targets=targets, tiers=tiers, novelty=False)


def test_ledger_trigger_rejects_disallowed_s1(ledger, targets, tiers):
    import sqlite3
    from rtt.records import Outcome
    c = ingest_raw(ledger, _raw(_cubic(3.0, ["Mg", "H"], [[0, 0, 0], [0.5, 0.5, 0.5]])), targets=targets, tiers=tiers)
    with pytest.raises((ValueError, sqlite3.DatabaseError)):
        ledger.insert_outcome(c, Outcome(stage="S1", outcome="INVALID"), targets_sha256=targets.sha256,
                              tiers_sha256=tiers.sha256, origin="test")
