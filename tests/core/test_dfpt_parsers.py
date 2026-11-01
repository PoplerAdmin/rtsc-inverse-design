"""Tests for rtt.core.dfpt parsers, unit conversion and the atomic-measure discretisation.

Two kinds of test data are used, and they are labelled:

* SYNTHETIC files written here with the exact Fortran WRITE formats copied from the QE / EPW
  source (PHonon/PH/matdyn.f90 SUBROUTINE a2Fdos; EPW/src/supercond.f90 evaluate_a2f_lambda).
  Their alpha^2F is an analytic toy function; it is not data for any material.
* REAL output of a tiny, deliberately unconverged fcc-Al QE 6.7 run in tests/data_dfpt/
  qe67_al_smoke (see PROVENANCE.txt there). Used only to check parsing of genuine files.
"""
from __future__ import annotations

import math
import sys
from pathlib import Path

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from rtt.core.dfpt import (collect_qe_eph, epw_a2f_to_spectra, grid_to_atoms, parse_epw_a2f,  # noqa: E402
                       parse_qe_a2f_dos, parse_qe_lambda_file, qe_a2f_dos_to_spectrum,
                       read_two_column, spectrum_lambda, spectrum_omega_log, check_epw_spectrum)
from rtt.core.dfpt import units  # noqa: E402
from rtt.core.dfpt.parsers import ffloat  # noqa: E402
from rtt.core.spectrum import Spectrum  # noqa: E402

REAL = Path(__file__).parent / "data_dfpt" / "qe67_al_smoke"


# ----------------------------------------------------------------------------- Fortran emulation
def fortran_E(x: float, w: int = 16, d: int = 6) -> str:
    """gfortran Ew.d output: 0.dddddd E+xx (mantissa in [0.1, 1))."""
    if x == 0:
        s = "0." + "0" * d + "E+00"
    else:
        sign = "-" if x < 0 else ""
        a = abs(x)
        e = math.floor(math.log10(a)) + 1
        m = round(a / 10 ** e, d)
        if m >= 1.0:
            m /= 10
            e += 1
        s = f"{sign}0.{int(round(m * 10 ** d)):0{d}d}E{e:+03d}"
    return s.rjust(w)


def write_qe_a2f_dos(path: Path, Emin_Ry: float, delta_Ry: float, a2f_modes: np.ndarray,
                     synthetic_note: bool = True) -> float:
    """SYNTHETIC a2F.dos in the format of matdyn.f90 a2Fdos. Returns QE's footer lambda."""
    ndos, nmodes = a2f_modes.shape
    lines = ["", " # Eliashberg function a2F (per both spin)", " #  frequencies in Rydberg  ",
             " # DOS normalized to E in Rydberg: a2F_total, a2F(mode) ", ""]
    if synthetic_note:
        lines.insert(1, " # SYNTHETIC test data written by tests/test_dfpt_parsers.py")
    lam = 0.0
    for n in range(1, ndos + 1):
        E = Emin_Ry + (n - 1) * delta_Ry + 0.5 * delta_Ry
        tot = float(np.sum(a2f_modes[n - 1]))
        lam += 2.0 * tot / E * delta_Ry
        lines.append("   " + "".join(fortran_E(v) for v in [E, tot, *a2f_modes[n - 1]]))
    lines.append(f"  lambda =  {lam!r}         Delta =    {delta_Ry:.16E}")
    path.write_text("\n".join(lines) + "\n")
    return lam


def toy_a2f(E: np.ndarray, peaks=((0.0015, 0.0002, 0.6), (0.003, 0.0003, 0.3))) -> np.ndarray:
    """Analytic toy a2F(E) (E in Ry): sum of Gaussians, zero for E <= 0. SYNTHETIC."""
    out = np.zeros_like(E)
    for c, s, h in peaks:
        out += h * np.exp(-0.5 * ((E - c) / s) ** 2)
    return np.where(E > 0, out, 0.0)


# ----------------------------------------------------------------------------- units
def test_unit_factors_match_qe_constants():
    # values derived from QE constants.f90 (exact SI + CODATA-2018 Hartree)
    assert units.RY_TO_K == pytest.approx(157887.512402, rel=1e-11)
    assert units.RY_TO_CMM1 == pytest.approx(109737.315682, rel=1e-11)
    assert units.RY_TO_THZ == pytest.approx(3289.84196025, rel=1e-11)
    assert units.MEV_TO_K == pytest.approx(11.6045181216, rel=1e-11)
    assert units.THZ_TO_K == pytest.approx(47.9924307337, rel=1e-11)
    assert units.CMM1_TO_K == pytest.approx(1.4387768775, rel=1e-10)
    # internal consistency: 1 Ry expressed in cm^-1 and THz maps to the same kelvin
    assert units.RY_TO_CMM1 * units.CMM1_TO_K == pytest.approx(units.RY_TO_K, rel=1e-13)
    assert units.RY_TO_THZ * units.THZ_TO_K == pytest.approx(units.RY_TO_K, rel=1e-13)
    assert units.to_kelvin_factor("cm^-1") == units.CMM1_TO_K
    with pytest.raises(ValueError):
        units.to_kelvin_factor("furlong")


def test_ffloat_fortran_variants():
    assert ffloat("0.123456E-02") == pytest.approx(1.23456e-3)
    assert ffloat("1.5D+01") == 15.0
    assert ffloat("0.100000-100") == pytest.approx(1e-101)


# ----------------------------------------------------------------------------- QE a2F.dos
def test_qe_a2f_dos_synthetic_roundtrip(tmp_path):
    delta = 4e-5
    ndos = 120
    E = 0.0 + (np.arange(1, ndos + 1) - 0.5) * delta
    a = toy_a2f(E)
    modes = np.stack([0.5 * a, 0.3 * a, 0.2 * a], axis=1)
    lam_qe = write_qe_a2f_dos(tmp_path / "a2F.dos3", 0.0, delta, modes)
    tab = parse_qe_a2f_dos(tmp_path / "a2F.dos3")
    assert tab.isig == 3 and tab.nmodes == 3 and len(tab.E_Ry) == ndos
    assert tab.delta_Ry == pytest.approx(delta, rel=1e-15)
    assert tab.lambda_reported == pytest.approx(lam_qe, rel=1e-15)
    sp = qe_a2f_dos_to_spectrum(tab, label="SYNTHETIC toy")
    assert isinstance(sp, Spectrum)
    # A in kelvin = sum a2F * Delta * Ry->K ; lambda dimensionless, unchanged by units
    A_expected = float(np.sum(a)) * delta * units.RY_TO_K
    assert sp.A() == pytest.approx(A_expected, rel=2e-6)       # E16.6 rounding of a2F
    assert spectrum_lambda(sp) == pytest.approx(lam_qe, rel=2e-6)
    assert sp.lambda_() == pytest.approx(spectrum_lambda(sp), rel=1e-12)
    assert sp.omega_log() == pytest.approx(spectrum_omega_log(sp), rel=1e-12)
    assert max(sp.omega_K) == pytest.approx(E[a > 0].max() * units.RY_TO_K, rel=1e-6)
    assert sp.provenance["dropped_nonpositive_omega"]["count"] == 0
    assert sp.provenance["warnings"] == []


def test_qe_a2f_dos_negative_bins_dropped_and_recorded(tmp_path):
    # Emin < 0 as matdyn sets when imaginary modes exist; put SYNTHETIC weight below zero
    delta = 5e-5
    Emin = -10 * delta
    ndos = 100
    E = Emin + (np.arange(1, ndos + 1) - 0.5) * delta
    a = toy_a2f(E) + np.where(E < 0, 0.05, 0.0)
    lam_qe = write_qe_a2f_dos(tmp_path / "a2F.dos1", Emin, delta, a[:, None])
    tab = parse_qe_a2f_dos(tmp_path / "a2F.dos1")
    sp = qe_a2f_dos_to_spectrum(tab, label="SYNTHETIC negative-bin test")
    d = sp.provenance["dropped_nonpositive_omega"]
    assert d["count"] == 10
    assert d["abs_mass_K"] == pytest.approx(10 * 0.05 * delta * units.RY_TO_K, rel=1e-5)
    assert any("omega <= 0" in w and w.startswith("WARNING") for w in sp.provenance["warnings"])
    assert min(sp.omega_K) > 0
    # QE footer includes the negative bins (with negative 2 a2F/E contributions)
    from rtt.core.dfpt import check_qe_spectrum
    cs = check_qe_spectrum(sp, tab, None)
    assert cs[0].name == "lambda_vs_a2F_footer" and cs[0].ok
    assert tab.lambda_all_bins() == pytest.approx(lam_qe, rel=2e-6)
    assert lam_qe < spectrum_lambda(sp)


def test_qe_a2f_dos_rejects_non_qe_file(tmp_path):
    p = tmp_path / "a2F.dos1"
    p.write_text("# meV a2F\n1 0.1\n2 0.2\n")
    with pytest.raises(ValueError, match="Rydberg"):
        parse_qe_a2f_dos(p)


def test_negative_a2f_policy(tmp_path):
    delta = 4e-5
    E = (np.arange(1, 81) - 0.5) * delta
    a = toy_a2f(E)
    a[40] = -0.2      # SYNTHETIC large negative value at omega > 0
    write_qe_a2f_dos(tmp_path / "a2F.dos1", 0.0, delta, a[:, None])
    tab = parse_qe_a2f_dos(tmp_path / "a2F.dos1")
    with pytest.raises(ValueError, match="negative"):
        qe_a2f_dos_to_spectrum(tab)
    sp = qe_a2f_dos_to_spectrum(tab, negative_policy="clip")
    assert sp.provenance["clipped_negative_mass_K"] == pytest.approx(0.2 * delta * units.RY_TO_K,
                                                                     rel=1e-5)
    assert any("clipped" in w for w in sp.provenance["warnings"])


# ----------------------------------------------------------------------------- lambda file
def test_qe_lambda_file_format(tmp_path):
    # SYNTHETIC, format '(" Broadening ",F8.4," lambda ",F12.4," dos(Ef)",F8.4," omega_ln [K]",F12.4)'
    rows = [(0.005, 1.2345, 12.3456, 412.5), (0.010, 0.9876, 1234.5678, 400.0)]
    txt = "\n  Electron-phonon coupling constant, lambda \n\n"
    for b, lam, dos, w in rows:
        dos_s = f"{dos:8.4f}" if dos < 1000 else "********"
        txt += f" Broadening {b:8.4f} lambda {lam:12.4f} dos(Ef){dos_s} omega_ln [K]{w:12.4f}\n"
    p = tmp_path / "lambda"
    p.write_text(txt)
    ent = parse_qe_lambda_file(p)
    assert len(ent) == 2
    assert ent[0].broadening_Ry == 0.005 and ent[0].lambda_ == 1.2345
    assert ent[0].dos_ef == 12.3456 and ent[0].omega_ln_K == 412.5
    assert ent[1].dos_ef is None   # Fortran overflow '********'


# ----------------------------------------------------------------------------- real QE output
@pytest.mark.skipif(not REAL.exists(), reason="real QE fixture missing")
def test_real_qe67_al_smoke_files():
    """REAL QE 6.7 output (smoke test, unconverged): parse + internal consistency."""
    res = collect_qe_eph(REAL, "Al-smoke-REAL-unconverged",
                         base_provenance={"synthetic": False, "note": "smoke test only"})
    assert len(res.spectra) == 10
    assert len(res.lambda_entries) == 10
    assert [e.broadening_Ry for e in res.lambda_entries] == pytest.approx(
        [0.01 * i for i in range(1, 11)])
    for sp in res.spectra:
        tab = parse_qe_a2f_dos(REAL / f"a2F.dos{sp.provenance['smearing_index']}")
        # histogram lambda must agree with QE's own footer to print precision
        assert spectrum_lambda(sp) == pytest.approx(tab.lambda_reported, rel=2e-6)
        # first bin centre is Delta/2 (Emin = 0, no imaginary modes in this run)
        assert tab.E_Ry[0] == pytest.approx(0.5 * tab.delta_Ry, rel=1e-5)
        assert sp.provenance["synthetic"] is False
    assert res.checks_ok, [c for c in res.checks if not c.ok]
    names = {c.name for c in res.checks}
    assert names == {"lambda_vs_a2F_footer", "lambda_vs_lambda_file", "omega_log_vs_lambda_file"}
    assert res.spread["n_alternatives"] == 10 and res.spread["lambda"]["rel_spread"] > 0
    # check fails loudly when tolerance is tight: lambda file is a direct q-sum, not a histogram
    from rtt.core.dfpt import Tolerances
    tight = collect_qe_eph(REAL, "x", tol=Tolerances(lambda_file_rel=1e-6, lambda_file_abs=0))
    assert not tight.checks_ok
    assert any(f.startswith("DFPT_CHECK_FAIL:lambda_vs_lambda_file") for f in tight.flags)


# ----------------------------------------------------------------------------- EPW
def write_epw_a2f(path: Path, w_meV: np.ndarray, a2f: np.ndarray, smear_meV, dw: float) -> list[float]:
    """SYNTHETIC <prefix>.a2f in the format of EPW/src/supercond.f90 (qe-6.7.0)."""
    nsm = a2f.shape[1]
    cum = np.cumsum(2.0 * a2f / w_meV[:, None] * dw, axis=0)
    lines = [f" w[meV] a2f and integrated 2*a2f/w for {nsm:4d} smearing values"]
    for i in range(len(w_meV)):
        lines.append(f"{w_meV[i]:12.7f}" + "".join(f"{v:12.7f}" for v in [*a2f[i], *cum[i]]))
    lines.append(" Integrated el-ph coupling")
    lines.append("  #         " + "".join(f"{v:12.7f}" for v in cum[-1]))
    lines.append(" Phonon smearing (meV)")
    lines.append("  #         " + "".join(f"{v:12.7f}" for v in smear_meV))
    lines.append(f"Electron smearing (eV){0.1:12.7f}")
    lines.append(f"Fermi window (eV){0.4:12.7f}")
    lines.append(f"Summed el-ph coupling {float(cum[-1, 0]):12.7f}")
    path.write_text("\n".join(lines) + "\n")
    return list(cum[-1])


def test_epw_a2f_synthetic(tmp_path):
    nq, wmax = 200, 100.0
    dw = wmax / nq
    w = dw * np.arange(1, nq + 1)                                # wsph(i) = i*dw
    a = np.stack([0.5 * np.exp(-0.5 * ((w - 40) / s) ** 2) for s in (4.0, 6.0)], axis=1)
    lam = write_epw_a2f(tmp_path / "x.a2f", w, a, [0.05, 0.10], dw)
    tab = parse_epw_a2f(tmp_path / "x.a2f")
    assert tab.a2f.shape == (nq, 2) and tab.phonon_smearing_meV == [0.05, 0.10]
    assert tab.electron_smearing_eV == 0.1 and tab.fermi_window_eV == 0.4
    sps = epw_a2f_to_spectra(tab, label_prefix="SYNTHETIC:")
    assert len(sps) == 2
    for j, sp in enumerate(sps):
        assert spectrum_lambda(sp) == pytest.approx(lam[j], rel=1e-5)
        assert sp.A() == pytest.approx(float(np.sum(a[:, j])) * dw * units.MEV_TO_K, rel=1e-5)
        assert all(c.ok for c in check_epw_spectrum(sp))


# ----------------------------------------------------------------------------- generic reader
@pytest.mark.parametrize("unit,factor", [("meV", units.MEV_TO_K), ("THz", units.THZ_TO_K),
                                         ("cm-1", units.CMM1_TO_K), ("Ry", units.RY_TO_K),
                                         ("K", 1.0)])
def test_two_column_units(tmp_path, unit, factor):
    # SYNTHETIC Einstein-like peak at 300 K, expressed in `unit`
    wK = np.linspace(5.0, 600.0, 400)
    aK = 1.2 * np.exp(-0.5 * ((wK - 300.0) / 20.0) ** 2)
    x = wK / factor
    p = tmp_path / f"a2f_{unit}.dat"
    np.savetxt(p, np.column_stack([x, aK]), header=f"SYNTHETIC omega[{unit}] a2F")
    sp = read_two_column(p, unit, label="SYNTHETIC")
    A_exact = 1.2 * 20.0 * math.sqrt(2 * math.pi)                # Gaussian integral, kelvin
    assert sp.A() == pytest.approx(A_exact, rel=1e-6)
    from scipy.integrate import quad
    f = lambda w: 1.2 * math.exp(-0.5 * ((w - 300.0) / 20.0) ** 2)  # noqa: E731
    l0 = quad(lambda w: f(w) / w, 5.0, 600.0, points=[300.0])[0]
    l1 = quad(lambda w: f(w) * math.log(w) / w, 5.0, 600.0, points=[300.0])[0]
    assert spectrum_omega_log(sp) == pytest.approx(math.exp(l1 / l0), rel=1e-5)
    assert spectrum_lambda(sp) == pytest.approx(2 * l0, rel=1e-5)
    assert sp.provenance["input_unit"] == units.canonical_unit(unit)


def test_rules_integrate_polynomial():
    x = np.linspace(0.0, 2.0, 201)
    y = x ** 2
    for rule in ("bins", "trapezoid", "midpoint"):
        pos, m = grid_to_atoms(x, y, rule)
        # all three are consistent quadratures of int_0^2 x^2 = 8/3 (bins extends half a cell)
        tol = 2e-2 if rule == "bins" else 1e-4
        assert float(np.sum(m)) == pytest.approx(8.0 / 3.0, rel=tol), rule
    pos, m = grid_to_atoms(x, y, "midpoint")
    assert len(pos) == 200 and pos[0] == pytest.approx(0.005)
    with pytest.raises(ValueError):
        grid_to_atoms(x, y, "simpson")


def test_two_column_drops_nonpositive_with_record(tmp_path):
    x = np.linspace(-2.0, 10.0, 13)         # meV, SYNTHETIC with weight at omega <= 0
    y = np.full_like(x, 0.1)
    p = tmp_path / "a.dat"
    np.savetxt(p, np.column_stack([x, y]))
    sp = read_two_column(p, "meV", rule="trapezoid")
    d = sp.provenance["dropped_nonpositive_omega"]
    assert d["count"] == 3                   # -2, -1, 0
    assert d["abs_mass_K"] == pytest.approx(0.1 * 2.5 * units.MEV_TO_K)
    assert min(sp.omega_K) > 0


# ----------------------------------------------------------------------------- real failure modes
EXP50 = Path(__file__).parent / "data_dfpt" / "qe67_al_expanded_a5.0"
EXP56 = Path(__file__).parent / "data_dfpt" / "qe67_al_expanded_a5.6"


@pytest.mark.skipif(not EXP50.exists(), reason="fixture missing")
def test_real_qe_nan_a2f_is_an_error_not_dropped():
    """REAL QE 6.7 output for a dynamically unstable (over-expanded) Al lattice: NaN a2F."""
    tab = parse_qe_a2f_dos(EXP50 / "a2F.dos1")
    assert np.isnan(tab.a2f_total).all() and tab.E_Ry[0] < 0
    with pytest.raises(ValueError, match="non-finite"):
        qe_a2f_dos_to_spectrum(tab)
    res = collect_qe_eph(EXP50, "Al-a5.0-REAL-unstable", el_ph_nsigma=1)
    assert not res.ok and res.spectra == []
    assert "DFPT_SPECTRUM_ERROR:sigma1" in res.flags
    assert any(f.startswith("IMAGINARY_MODES") for f in res.flags)
    assert any("imaginary modes present" in e for e in res.errors)


@pytest.mark.skipif(not EXP56.exists(), reason="fixture missing")
def test_real_qe_negative_grid_zero_a2f():
    """REAL QE 6.7 output: grid starting at E < 0 (Emin = min freq) and a2F identically 0."""
    tab = parse_qe_a2f_dos(EXP56 / "a2F.dos1")
    assert tab.E_Ry[0] < 0 < tab.E_Ry[-1]
    assert np.all(tab.a2f_total == 0) and tab.lambda_reported == 0.0
    with pytest.raises(ValueError, match="identically zero"):
        qe_a2f_dos_to_spectrum(tab)
    res = collect_qe_eph(EXP56, "Al-a5.6-REAL-unstable", el_ph_nsigma=1)
    assert not res.ok and res.lambda_entries[0].lambda_ == 0.0
    assert res.lambda_entries[0].omega_ln_K is None or math.isnan(res.lambda_entries[0].omega_ln_K)
