"""Campaign 1 fix: S4 computes E_form vs elemental references and passes it to decide_s4, so the
H7-heur screen is applied instead of being NOT CHECKED on every candidate."""
from __future__ import annotations

from types import SimpleNamespace

import pytest

from rtt.config import Threshold
from rtt.dft import element_refs as ER
from rtt.dft.refhull import RY_EV


def test_formation_energy_and_missing_ref():
    mu = {"Mg": -10.0, "H": -15.0}
    e_tot_eV = 2 * -10.0 + 4 * -15.0 - 6 * 0.2                 # 0.2 eV/atom above elements... below
    assert ER.formation_energy(e_tot_eV / RY_EV, {"Mg": 2, "H": 4}, mu) == pytest.approx(-0.2)
    assert ER.formation_energy(e_tot_eV / RY_EV, {"Mg": 2, "Ir": 1}, mu) is None
    assert ER.formation_energy(None, {"Mg": 2}, mu) is None


def test_elemental_structure_picks_lowest_and_h2():
    rows = [{"chemsys": "Mg", "n_atoms": 2, "energy_per_atom_eV": -1.5, "entry_id": "a",
             "structure": {"lattice": [[3, 0, 0], [0, 3, 0], [0, 0, 5]], "species": ["Mg", "Mg"],
                           "frac_coords": [[0, 0, 0], [0.5, 0.5, 1.25]]}},
            {"chemsys": "Mg", "n_atoms": 1, "energy_per_atom_eV": -1.4, "entry_id": "b",
             "structure": {"lattice": [[3, 0, 0], [0, 3, 0], [0, 0, 3]], "species": ["Mg"], "frac_coords": [[0, 0, 0]]}},
            {"chemsys": "H-Mg", "n_atoms": 3, "energy_per_atom_eV": -9.0, "entry_id": "c", "structure": {}}]
    st, ref = ER.elemental_structure("Mg", query=lambda el, machine=None: {"rows": rows})
    assert ref == "alexmp:a" and st["frac_coords"][1][2] == pytest.approx(0.25)
    st, ref = ER.elemental_structure("H")
    assert st["molecule"] and ref == "H2_box10A"


def test_s4_stage_wires_eform_into_decision(tmp_path):
    from rtt.dft import electronic as E
    from rtt.engine.stages import S4Electronic
    from rtt.records import StageResult
    st = S4Electronic()
    res = SimpleNamespace(composition={"Mg": 2, "H": 4}, e_total_Ry=(2 * -10.0 + 4 * -15.0 + 6 * 0.9) / RY_EV)
    ctx = SimpleNamespace(workdir=tmp_path, launcher=None, machine=None,
                          extra={"element_refs_fn": lambda els, prof, root, **kw: ({"Mg": -10.0, "H": -15.0}, {})})
    out = st.formation_energy(res, None, ctx)
    assert out["e_form_eV"] == pytest.approx(0.9)
    # pure element: no E_form
    assert st.formation_energy(SimpleNamespace(composition={"Al": 1}, e_total_Ry=-1.0), None, ctx)["e_form_eV"] is None
    # failing references: None, never raised
    bad = SimpleNamespace(workdir=tmp_path, launcher=None, machine=None,
                          extra={"element_refs_fn": lambda *a, **k: (_ for _ in ()).throw(RuntimeError("x"))})
    assert st.formation_energy(res, None, bad)["e_form_eV"] is None

    # decide() receives it: 0.9 eV/atom > 0.56 window -> SCREENED_OUT(H7-heur)
    from tests.back.test_electronic_belli import _res  # noqa: F401 - reuse the fixture builder
    er = _res()
    sr = StageResult("S4", "c-x", True, {"e_form_eV": 0.9}, 1, payload=er)
    th = SimpleNamespace(threshold=lambda k: Threshold(0.56, "doc22#RC1-c", provisional=True))
    tgt = SimpleNamespace(threshold=lambda k: {"hard.H1_metallic": Threshold(0.0, "doc02#5a:H1"),
                                                "hard.H2_nonmagnetic": Threshold(0.05, "doc02#5a:H2", provisional=True)}[k])
    o = st.decide(sr, th, tgt)
    assert (o.outcome, o.criterion) == ("SCREENED_OUT", "H7-heur")
    sr2 = StageResult("S4", "c-x", True, {"e_form_eV": None}, 1, payload=er)
    assert "H7_HEUR_NOT_CHECKED" in st.decide(sr2, th, tgt).tags
