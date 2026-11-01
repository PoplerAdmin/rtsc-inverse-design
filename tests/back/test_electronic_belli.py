"""rtt.dft.electronic band analysis + S4 decision; rtt.dft.belli percolation."""
from __future__ import annotations

import math

import numpy as np
import pytest

from rtt.dft import belli as B
from rtt.dft import electronic as E
from rtt.dft.workflow import DFPTRun

HA = E.HA_EV


def _xml(path, eigs_ha, weights, ef_ha, nat=2, lsda=False, nbnd_up=None):
    ks = "".join(f'<ks_energies><k_point weight="{w}">0 0 0</k_point><npw>10</npw>'
                 f'<eigenvalues size="{len(e)}">{" ".join(map(str, e))}</eigenvalues></ks_energies>'
                 for e, w in zip(eigs_ha, weights))
    extra = f"<nbnd_up>{nbnd_up}</nbnd_up>" if lsda else ""
    path.write_text(f'<?xml version="1.0"?><qes:espresso xmlns:qes="http://www.quantum-espresso.org/ns/qes/qes-1.0">'
                    f'<output><atomic_structure nat="{nat}"/><band_structure><lsda>{str(lsda).lower()}</lsda>'
                    f'{extra}<fermi_energy>{ef_ha}</fermi_energy>{ks}</band_structure></output></qes:espresso>')


def test_metal_bands(tmp_path):
    # band 0 fully below, band 1 crosses (0.20..0.30 Ha) with E_Fermi 0.25, band 2 above
    _xml(tmp_path / "d.xml", [[0.0, 0.20, 0.5], [0.05, 0.30, 0.6]], [1.0, 1.0], 0.25)
    b = E.parse_pw_xml(tmp_path / "d.xml")
    assert b.eig_eV.shape == (1, 2, 3) and abs(b.weights.sum() - 1) < 1e-12
    a = E.analyse_bands(b, dos_sigma_eV=0.1)
    assert a.n_sheets == 1 and a.gap_eV == 0.0
    # E_F (G-VALID): min(E_F - bottom, width/2) = min(0.05, 0.05) Ha
    assert a.E_F_eV == pytest.approx(0.05 * HA)
    assert a.E_F_hole_eV == pytest.approx(0.05 * HA)
    assert a.N_EF_per_eV_atom == pytest.approx(a.N_EF_per_eV_cell / 2)


def test_insulator_gap(tmp_path):
    _xml(tmp_path / "d.xml", [[0.0, 0.1, 0.4], [0.02, 0.12, 0.45]], [1.0, 3.0], 0.25)
    a = E.analyse_bands(E.parse_pw_xml(tmp_path / "d.xml"))
    assert a.n_sheets == 0 and a.gap_eV == pytest.approx((0.4 - 0.12) * HA)
    assert a.E_F_eV is None and a.N_EF_per_eV_cell < 1e-6


def test_lsda_sheets_count_per_spin(tmp_path):
    # 2 bands per spin; spin up band 1 crosses, spin down none
    _xml(tmp_path / "d.xml", [[0.0, 0.2, 0.0, 0.4], [0.05, 0.3, 0.05, 0.5]], [1, 1], 0.25,
         lsda=True, nbnd_up=2)
    a = E.analyse_bands(E.parse_pw_xml(tmp_path / "d.xml"))
    assert a.n_sheets_by_spin == [1, 0] and a.n_sheets == 1


def test_parse_magnetization():
    t = ("total magnetization       =     2.10 Bohr mag/cell\nabsolute magnetization    =     2.30 Bohr mag/cell\n"
         "total magnetization       =    -0.01 Bohr mag/cell\nabsolute magnetization    =     0.80 Bohr mag/cell\n")
    assert E.parse_magnetization(t) == (0.01, 0.80)
    assert E.parse_total_energy_Ry("!    total energy              =     -39.5 Ry\n") == -39.5


def _res(**kw):
    run = DFPTRun(material_id="x", profile={"name": "s4"}, workdir="/tmp", status="done")
    r = E.ElectronicResult(run=run, gap_eV=0.0, moment_muB=0.0, abs_moment_muB=0.0, nat=2,
                           H_f=0.0, n_sheets=1, N_EF_per_eV_atom=0.3)
    for k, v in kw.items():
        setattr(r, k, v)
    return r


def test_decide_s4():
    from rtt.config import Threshold
    gap = Threshold(0.0, "doc02#5a:H1")
    mom = Threshold(0.05, "doc02#5a:H2", provisional=True)
    hull = Threshold(0.42, "doc02#5a:H7", provisional=True)
    o = E.decide_s4(_res(), gap_max=gap, moment_tol=mom, ehull_heur_max=hull, e_form_eV=0.1)
    assert o.outcome == "PASS" and o.stage == "S4" and o.provisional
    assert "PROVISIONAL_THRESHOLD(H2)" in o.tags
    ins = dict(n_sheets=0, N_EF_per_eV_atom=0.0)
    o = E.decide_s4(_res(gap_eV=1.0, **ins), gap_max=gap, moment_tol=mom)
    assert (o.outcome, o.criterion, o.source_id) == ("EXCLUDED_ASSUMPTION", "H1", "doc02#5a:H1")
    o = E.decide_s4(_res(gap_eV=0.2, **ins), gap_max=gap, moment_tol=mom)
    assert (o.outcome, o.criterion) == ("DEFERRED", "H1-marginal")
    # magnetic ground state (E_spin lower by 6.8 eV/2 atoms) -> excluded; per-atom moment 0.25
    o = E.decide_s4(_res(moment_muB=0.5, abs_moment_muB=0.5, e_total_Ry=-10.0,
                         e_total_spin_Ry=-10.5), gap_max=gap, moment_tol=mom)
    assert (o.outcome, o.criterion, o.provisional) == ("EXCLUDED_ASSUMPTION", "H2", True)
    # metastable magnetic solution (higher than NM) -> PASS with a flag
    o = E.decide_s4(_res(moment_muB=0.5, abs_moment_muB=0.5, e_total_Ry=-10.5,
                         e_total_spin_Ry=-10.0), gap_max=gap, moment_tol=mom)
    assert o.outcome == "PASS" and any(t.startswith("H2_METASTABLE") for t in o.tags)
    # energies missing -> DEFERRED; marginal moment (<= 2 tol per atom) -> DEFERRED
    o = E.decide_s4(_res(moment_muB=0.5, abs_moment_muB=0.5), gap_max=gap, moment_tol=mom)
    assert (o.outcome, o.criterion) == ("DEFERRED", "H2-energy-unknown")
    o = E.decide_s4(_res(moment_muB=0.15, abs_moment_muB=0.15, e_total_Ry=-10.0,
                         e_total_spin_Ry=-10.5), gap_max=gap, moment_tol=mom)
    assert (o.outcome, o.criterion) == ("DEFERRED", "H2-marginal")
    o = E.decide_s4(_res(moment_muB=0.5), gap_max=gap, moment_tol=mom, enable_h2=False)
    assert o.outcome == "PASS" and "H2_NOT_CHECKED" in o.tags
    o = E.decide_s4(_res(), gap_max=gap, moment_tol=mom, ehull_heur_max=hull, e_form_eV=0.5)
    assert (o.outcome, o.criterion) == ("SCREENED_OUT", "H7-heur")
    assert "PROVISIONAL_THRESHOLD(H7)" in o.tags
    o = E.decide_s4(_res(), gap_max=gap, moment_tol=mom)
    assert "H7_HEUR_NOT_CHECKED" in o.tags
    # missing thresholds are visible, never a silent "checked"
    o = E.decide_s4(_res(), gap_max=gap, moment_tol=None)
    assert "H2_NOT_CHECKED" in o.tags
    bad = _res(error="boom")
    assert E.decide_s4(bad, gap_max=gap, moment_tol=mom).outcome == "ERROR"


def test_toml_inline_threshold_tables():
    """targets.toml [hard] uses named keys, not 'value' (RC1 F1)."""
    from rtt.config import load_targets
    hard = load_targets().get("hard")
    v, s, p = E.th_parts(hard["H1_metallic"])
    assert v == 0.0 and s
    v, s, p = E.th_parts(hard["H2_nonmagnetic"])
    assert v is not None and s and p
    with pytest.raises(ValueError):
        E.th_parts({"foo": 1, "source": "x"})
    o = E.decide_s4(_res(gap_eV=2.0, n_sheets=0, N_EF_per_eV_atom=0.0),
                    gap_max=hard["H1_metallic"], moment_tol=hard["H2_nonmagnetic"])
    assert (o.outcome, o.criterion) == ("EXCLUDED_ASSUMPTION", "H1")


def test_parse_dosx(tmp_path):
    p = tmp_path / "pwscf.dos"
    rows = "\n".join(f"{e:.3f} {1.0 + 0.1 * e:.4f} {2.0 * e:.4f}" for e in np.arange(0, 10.01, 0.5))
    p.write_text("#  E (eV)   dos(E)     Int dos(E) EFermi =    5.250 eV\n" + rows + "\n")
    d = E.parse_dosx(p)
    assert d["EFermi_eV"] == 5.25
    assert d["N_EF_per_eV_cell"] == pytest.approx(1.525)
    assert d["int_dos_at_EF"] == pytest.approx(10.5)


def test_sheet_dos_fractions(tmp_path):
    # two crossing bands, one with many states at E_F, one barely crossing
    _xml(tmp_path / "d.xml", [[0.20, 0.249], [0.30, 0.26], [0.25, 0.4]], [1, 1, 1], 0.25)
    a = E.analyse_bands(E.parse_pw_xml(tmp_path / "d.xml"), dos_sigma_eV=0.1)
    assert a.n_sheets == 2
    fr = [c["dos_fraction"] for c in a.crossing_bands]
    assert sum(fr) == pytest.approx(1.0, abs=1e-6) or sum(fr) <= 1.0
    assert a.sheet_w_min == pytest.approx(min(min(c["dos_fraction"], c["dos_fraction_2sigma"])
                                              for c in a.crossing_bands))


def test_formation_energy_lower_bound():
    # 1 Ry total for A1B1 with refs 0 -> 6.8 eV/atom
    e = E.formation_energy_per_atom_eV(1.0, {"A": 1, "B": 1}, {"A": 0.0, "B": 0.0})
    assert e == pytest.approx(E.RY_EV / 2)


def test_networking_value_channels():
    n = 12
    elf = np.zeros((n, n, n))
    elf[:, 5, 5] = 0.8          # line along x
    elf[5, :, 5] = 0.8          # line along y
    elf[5, 5, :] = 0.6          # line along z: weakest link
    assert B.networking_value(elf, tol=1e-3) == pytest.approx(0.6, abs=2e-3)
    elf2 = np.zeros((n, n, n))
    elf2[:, :, 5] = 0.9         # a slab percolates in 2 directions only
    assert B.networking_value(elf2) == 0.0


def test_h_fraction_and_tc_phi():
    cell = np.eye(3) * 4.0
    pos = np.array([[0.0, 0, 0], [2.0, 2.0, 2.0]])
    dens = np.zeros((8, 8, 8))
    dens[4, 4, 4] = 1.0         # at (2,2,2): inside the H cell
    dens[0, 0, 0] = 3.0
    assert B.h_fraction_voronoi(dens, cell, pos, ["Mg", "H"]) == pytest.approx(0.25)
    assert B.tc_phi(0.729, 0.75, 0.664) == pytest.approx(750 * 0.729 * 0.75 * 0.664 ** (1 / 3) - 85)
    assert B.tc_phi(0.5, 0.0, 0.0) is None


def test_score_s4():
    from rtt.score.s4 import components, score_s4
    noh = {"N_EF_per_eV_atom": 0.3, "H_f": 0.0}
    assert score_s4(noh) == pytest.approx(0.5)
    h = {"N_EF_per_eV_atom": 0.3, "H_f": 0.75, "Tc_phi_K": 300.0, "H_DOS": 0.5}
    c = components(h)
    assert c == {"s_phi": 1.0, "s_dos": 0.5, "s_hdos": 1.0}
    assert score_s4(h) == pytest.approx(0.6 + 0.125 + 0.15)
    assert score_s4({}) == 0.0
    assert score_s4({"H_f": 0.5, "Tc_phi_K": -50.0}) == 0.0


def test_h1_gaussian_tail_does_not_override_tetrahedron():
    """Campaign 1 (NaAuH2): PBE gap 2.04 eV, tetrahedron N(E_F) = 0, 2-sigma Gaussian 0.017 > nef_min.
    The tetrahedron value decides: EXCLUDED(H1), not PASS."""
    from rtt.config import Threshold
    gap = Threshold(0.0, "doc02#5a:H1")
    mom = Threshold(0.05, "doc02#5a:H2", provisional=True)
    r = _res(gap_eV=2.04, n_sheets=0, N_EF_per_eV_atom=0.0, N_EF_gauss_per_eV_atom=0.0002,
             N_EF_2sigma_per_eV_atom=0.017)
    o = E.decide_s4(r, gap_max=gap, moment_tol=mom)
    assert (o.outcome, o.criterion) == ("EXCLUDED_ASSUMPTION", "H1")
    # without a tetrahedron value the Gaussian fallback still applies
    r = _res(gap_eV=0.0, n_sheets=0, N_EF_per_eV_atom=None, N_EF_gauss_per_eV_atom=0.2)
    assert E.decide_s4(r, gap_max=gap, moment_tol=mom).outcome == "PASS"
