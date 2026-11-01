"""rtt.dft.refhull: hull arithmetic (pymatgen) and the per-profile cache key."""
from __future__ import annotations

import json

import pytest

from rtt.dft import refhull as R
from rtt.dft.profiles import PROFILES


def _p(name, comp, epa, zpe=None):
    n = sum(comp.values())
    return R.PhaseEnergy(name=name, composition=comp, nat=n, e_total_Ry=epa * n / R.RY_EV,
                         e_per_atom_eV=epa, profile="ref_hull", cache_key=name, zpe_meV_atom=zpe)


def test_e_above_hull_and_zpe():
    pytest.importorskip("pymatgen")
    h = R.RefHull(chemsys=["Mg", "H"], profile="ref_hull",
                  phases=[_p("Mg", {"Mg": 1}, 0.0, 0.0), _p("H2", {"H": 2}, 0.0, 130.0),
                          _p("MgH2", {"Mg": 1, "H": 2}, -0.3, 80.0)])
    assert h.element_refs() == {"Mg": 0.0, "H": 0.0}
    # on the hull
    assert h.e_above_hull({"Mg": 1, "H": 2}, -0.3) == pytest.approx(0.0, abs=1e-9)
    # 0.1 eV/atom above MgH2
    assert h.e_above_hull({"Mg": 1, "H": 2}, -0.2) == pytest.approx(0.1, abs=1e-9)
    # ZPE-corrected: candidate ZPE 100 meV vs MgH2 80 -> +0.02
    assert h.e_above_hull({"Mg": 1, "H": 2}, -0.3, zpe_meV_atom=100.0, with_zpe=True) == \
        pytest.approx(0.02, abs=1e-9)
    with pytest.raises(ValueError):
        h.e_above_hull({"Mg": 1, "Li": 1}, 0.0)
    h.phases[0].zpe_meV_atom = None
    with pytest.raises(ValueError):
        h.e_above_hull({"Mg": 1, "H": 2}, -0.3, zpe_meV_atom=1.0, with_zpe=True)


def test_h2_molecule():
    m = R.h2_molecule()
    assert m["molecule"] and m["species"] == ["H", "H"]
    dx = (m["frac_coords"][1][0] - m["frac_coords"][0][0]) * 10.0
    assert dx == pytest.approx(0.7414)


def test_phase_energy_cache(tmp_path, monkeypatch):
    pdir = tmp_path / "pseudo"
    pdir.mkdir()
    (pdir / "Al_ONCV_PBE-1.2.upf").write_text('<PP_HEADER pseudo_type="NC" relativistic="scalar"/>')
    monkeypatch.setenv("RTT_PSEUDO_DIR", str(pdir))
    st = {"lattice": [[0, 2.02, 2.02], [2.02, 0, 2.02], [2.02, 2.02, 0]], "species": ["Al"],
          "frac_coords": [[0, 0, 0]], "masses_amu": {"Al": 26.98}}
    prof = PROFILES["ref_hull"]
    cache = tmp_path / "cache"
    # first call would need QE; emulate a cached result under the key phase_energy computes
    import hashlib
    from rtt.core.dfpt.structure import to_celldata
    from rtt.dft.profiles import pseudo_md5s, sg15_pseudos
    cell = to_celldata(st)
    key = hashlib.sha256("|".join([R._canon(cell), json.dumps(
        pseudo_md5s(sg15_pseudos(["Al"]), ["Al"]), sort_keys=True), prof.sha256(),
        "False"]).encode()).hexdigest()
    cache.mkdir()
    pe = _p("Al", {"Al": 1}, -1883.2)
    pe.cache_key = key
    (cache / f"{key}.json").write_text(json.dumps(pe.to_dict()))
    got = R.phase_energy("Al", st, prof, tmp_path / "w", cache_dir=cache)
    assert got.e_per_atom_eV == -1883.2 and not (tmp_path / "w").exists()
    # a different profile is a different key -> would run (we only check the key differs)
    assert prof.with_overrides(ecutwfc_Ry=60.0).sha256() != prof.sha256()
