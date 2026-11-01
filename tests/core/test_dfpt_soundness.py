"""Soundness fixes requested by the audit (truncated inputs, clipping, EPW bound conditions,
BandCoupling semantics). All files written here are SYNTHETIC unless a test says REAL.
"""
from __future__ import annotations

import struct
import sys
from pathlib import Path

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from rtt.core.dfpt import (BandCoupling, IncompleteFileError, SpectrumCheckError,  # noqa: E402
                       band_coupling_from_epw_lambda, load_attached_spectra, parse_epw_a2f,
                       parse_epw_lambda_fs, parse_qe_a2f_dos, qe_a2f_dos_to_spectrum,
                       read_two_column, epw_omega_max_meV)
from rtt.core.dfpt.band_coupling import check_epw_frequencies, read_epw_freq  # noqa: E402

REAL_EPW = Path(__file__).parent / "data_dfpt" / "epw53_pb_smoke"


def qe_file(path: Path, truncate_at_Ry: float | None = None, footer_lambda_scale: float = 1.0,
            extra_after: str | None = None) -> float:
    """SYNTHETIC QE a2F.dos (same construction as the audit test)."""
    d = 1.0e-4
    E = d * (np.arange(200) + 0.5)
    a2f = 7.6 * np.exp(-0.5 * ((E - 0.012) / 0.001) ** 2)
    lam = float(np.sum(2.0 * a2f / E * d))
    lines = ["", " # Eliashberg function a2F (per both spin)", " #  frequencies in Rydberg  ",
             " # DOS normalized to E in Rydberg: a2F_total, a2F(mode) ", " # SYNTHETIC", ""]
    for e, a in zip(E, a2f):
        if truncate_at_Ry is not None and e > truncate_at_Ry:
            break
        lines.append(f"   {e:16.6E}{a:16.6E}{a / 2:16.6E}{a / 2:16.6E}")
    else:
        lines.append(f"  lambda =   {lam * footer_lambda_scale:.8f}   Delta =   {d:.8E}")
    if extra_after:
        lines.append(extra_after)
    path.write_text("\n".join(lines) + "\n")
    return lam


# ----------------------------------------------------------------------------- (1) completeness
def test_truncated_qe_file_is_an_error(tmp_path):
    full, trunc = tmp_path / "a2F.dos1", tmp_path / "t" / "a2F.dos1"
    trunc.parent.mkdir()
    qe_file(full)
    qe_file(trunc, truncate_at_Ry=0.0075)
    assert parse_qe_a2f_dos(full).complete
    with pytest.raises(IncompleteFileError, match="footer"):
        parse_qe_a2f_dos(trunc)
    with pytest.raises(IncompleteFileError):
        load_attached_spectra(trunc, "qe_a2f_dos")
    # diagnostics path: accepted but flagged at assumption level
    tab = parse_qe_a2f_dos(trunc, allow_incomplete=True)
    assert not tab.complete
    sp = qe_a2f_dos_to_spectrum(tab)
    assert "INCOMPLETE_INPUT_FILE" in sp.provenance["flags"]
    assert any(a.startswith("A-FILE-INCOMPLETE") for a in sp.provenance["assumptions"])
    sp2 = load_attached_spectra(trunc, "qe_a2f_dos", allow_incomplete=True)[0]
    assert "INCOMPLETE_INPUT_FILE" in sp2.provenance["flags"]
    assert sp2.A() < load_attached_spectra(full, "qe_a2f_dos")[0].A()


def test_content_after_footer_is_an_error(tmp_path):
    p = tmp_path / "a2F.dos1"
    qe_file(p, extra_after="   1.0E-02   1.0E+00   0.5E+00   0.5E+00")
    with pytest.raises(ValueError, match="after the"):
        parse_qe_a2f_dos(p)


def write_epw(path: Path, nrows: int | None = None, block: bool = True, bad_block: bool = False):
    """SYNTHETIC EPW .a2f (EPW/src/supercond.f90 format), optionally truncated."""
    nq, dw = 100, 0.5
    w = dw * np.arange(1, nq + 1)
    a = np.stack([0.4 * np.exp(-0.5 * ((w - 20) / s) ** 2) for s in (2.0, 3.0)], axis=1)
    cum = np.cumsum(2 * a / w[:, None] * dw, axis=0)
    lines = [" w[meV] a2f and integrated 2*a2f/w for    2 smearing values"]
    for i in range(nrows or nq):
        lines.append(f"{w[i]:12.7f}" + "".join(f"{v:12.7f}" for v in [*a[i], *cum[i]]))
    if block:
        last = cum[-1] * (1.01 if bad_block else 1.0)
        lines += [" Integrated el-ph coupling", "  #         " + "".join(f"{v:12.7f}" for v in last),
                  " Phonon smearing (meV)", "  #         " + "".join(f"{v:12.7f}" for v in (0.5, 0.55))]
    path.write_text("\n".join(lines) + "\n")


def test_epw_completeness(tmp_path):
    ok, trunc, bad = tmp_path / "ok.a2f", tmp_path / "t.a2f", tmp_path / "b.a2f"
    write_epw(ok)
    write_epw(trunc, nrows=60, block=False)
    write_epw(bad, bad_block=True)
    assert parse_epw_a2f(ok).complete
    sps = load_attached_spectra(ok, "epw")
    assert len(sps) == 2 and all(s.provenance["checks"][0]["ok"] for s in sps)
    with pytest.raises(IncompleteFileError, match="Integrated"):
        parse_epw_a2f(trunc)
    with pytest.raises(IncompleteFileError):
        load_attached_spectra(trunc, "epw")
    with pytest.raises(ValueError, match="cumulative"):
        parse_epw_a2f(bad)
    diag = load_attached_spectra(trunc, "epw", allow_incomplete=True)
    assert all("INCOMPLETE_INPUT_FILE" in s.provenance["flags"] for s in diag)


# ----------------------------------------------------------------------------- (2) entry point
def test_load_attached_qe_checks(tmp_path):
    good, tampered = tmp_path / "a2F.dos1", tmp_path / "x" / "a2F.dos1"
    tampered.parent.mkdir()
    qe_file(good)
    qe_file(tampered, footer_lambda_scale=1.5)       # rows do not match the footer
    sp = load_attached_spectra(good, "qe_a2f_dos")[0]
    assert sp.provenance["assumptions"] == []
    assert sp.provenance["checks"][0]["name"] == "lambda_vs_a2F_footer"
    with pytest.raises(SpectrumCheckError, match="lambda_vs_a2F_footer"):
        load_attached_spectra(tampered, "qe_a2f_dos")
    with pytest.raises(ValueError, match="Rydberg"):
        load_attached_spectra(good, "qe_a2f_dos", unit="meV")
    with pytest.raises(ValueError, match="unknown kind"):
        load_attached_spectra(good, "vasp")


def test_load_attached_qe_with_lambda_file(tmp_path):
    p = tmp_path / "a2F.dos1"
    lam = qe_file(p)
    lf = tmp_path / "lambda"
    lf.write_text(f" Broadening {0.005:8.4f} lambda {lam:12.4f} dos(Ef){3.0:8.4f} omega_ln [K]"
                  f"{1000.0:12.4f}\n")     # SYNTHETIC; omega_ln deliberately inconsistent
    with pytest.raises(SpectrumCheckError, match="omega_log"):
        load_attached_spectra(p, "qe_a2f_dos", lambda_file=lf)


def test_load_attached_two_column(tmp_path):
    p = tmp_path / "a.dat"
    x = np.linspace(1, 60, 60)
    np.savetxt(p, np.column_stack([x, 0.3 * np.exp(-0.5 * ((x - 30) / 5) ** 2)]),
               header="SYNTHETIC omega[meV] a2F")
    with pytest.raises(ValueError, match="unit"):
        load_attached_spectra(p, "two_column")
    sp = load_attached_spectra(p, "two_column", "meV")[0]
    assert any(a.startswith("A-FILE-COMPLETE(unverified)") for a in sp.provenance["assumptions"])
    lam = sp.lambda_()
    sp2 = load_attached_spectra(p, "two_column", "meV", expected_lambda=lam)[0]
    assert sp2.provenance["assumptions"] == []
    with pytest.raises(SpectrumCheckError):
        load_attached_spectra(p, "two_column", "meV", expected_lambda=1.2 * lam)
    with pytest.raises(SpectrumCheckError):
        load_attached_spectra(p, "two_column", "meV", expected_A_K=2 * sp.A())


# ----------------------------------------------------------------------------- (3) clip / drop
def test_clip_and_drop_are_assumptions(tmp_path):
    p = tmp_path / "a.dat"
    np.savetxt(p, np.array([[-1.0, 0.2], [0.0, 0.2], [1.0, 0.5], [2.0, -0.3], [3.0, 0.5]]))
    with pytest.raises(ValueError, match="negative"):
        read_two_column(p, "meV")
    sp = read_two_column(p, "meV", negative_policy="clip")
    a = sp.provenance["assumptions"]
    assert any(s.startswith("A-NEG-CLIP(") and "CERTIFIED" in s for s in a)
    assert any(s.startswith("A-DROP-NONPOS(") and "REJECTED" in s for s in a)


# ----------------------------------------------------------------------------- (4) EPW bound
def write_freq(path: Path, xq, wf_Ry, nqf=(2, 2, 2)) -> None:
    """SYNTHETIC EPW ephmat/freq in gfortran sequential format (4-byte markers)."""
    def rec(b: bytes) -> bytes:
        return struct.pack("<i", len(b)) + b + struct.pack("<i", len(b))
    nmodes = len(wf_Ry[0])
    out = rec(struct.pack("<5i", len(xq), *nqf, nmodes))
    for q, w in zip(xq, wf_Ry):
        out += rec(struct.pack("<3d", *q))
        for v in w:
            out += rec(struct.pack("<d", v))
    path.write_bytes(out)


def _states(tmp_path):
    p = tmp_path / "x.lambda_FS"
    p.write_text("# SYNTHETIC\n" + "".join(
        f"{0.1 * i:12.6f}{0.0:12.6f}{0.0:12.6f}{1:8d}{0.01 * i:12.6f}{0.5 + 0.1 * i:24.15f}\n"
        for i in range(5)))
    return parse_epw_lambda_fs(p)


def test_epw_bound_mode_checks(tmp_path):
    st = _states(tmp_path)
    xq = [(0, 0, 0), (0.5, 0, 0)]
    good, imag, excl = tmp_path / "good", tmp_path / "imag", tmp_path / "excl"
    write_freq(good, xq, [[1e-9, 1e-9, 1e-9], [5e-4, 6e-4, 9e-4]])
    write_freq(imag, xq, [[1e-9, 1e-9, 1e-9], [-2e-4, 6e-4, 9e-4]])
    write_freq(excl, xq, [[1e-9, 1e-9, 1e-9], [5e-6, 6e-4, 9e-4]])     # 5e-6 Ry ~ 1.1 cm^-1
    wmax_meV = 9e-4 * 13605.693122994 * 1.1
    bc = band_coupling_from_epw_lambda(st, omega_max=wmax_meV, freq_file=good, eps_acustic_cm1=2.0)
    a = bc.provenance["assumptions"]
    assert len(a) == 3 and all("(verified:" in s for s in a)
    assert bc.provenance["mode_check"]["n_excluded_gamma"] == 3
    with pytest.raises(ValueError, match="imaginary"):
        band_coupling_from_epw_lambda(st, omega_max=wmax_meV, freq_file=imag, eps_acustic_cm1=2.0)
    with pytest.raises(ValueError, match="eps_acustic"):
        band_coupling_from_epw_lambda(st, omega_max=wmax_meV, freq_file=excl, eps_acustic_cm1=2.0)
    with pytest.raises(ValueError, match="omega_max"):
        band_coupling_from_epw_lambda(st, omega_max=1.0, freq_file=good, eps_acustic_cm1=2.0)
    loose = band_coupling_from_epw_lambda(st, omega_max=wmax_meV, freq_file=imag,
                                          eps_acustic_cm1=2.0, strict_modes=False)
    assert loose.provenance["mode_problems"]
    assert not any("NO-IMAGINARY" in s and "verified" in s for s in loose.provenance["assumptions"])
    # without a freq file: conditions are unverified assumptions
    plain = band_coupling_from_epw_lambda(st, omega_max=wmax_meV)
    assert {s.split("(")[0].split(":")[0] for s in plain.provenance["assumptions"]} == \
        {"EPW-OMEGA-MAX", "EPW-NO-IMAGINARY", "EPW-NO-EXCLUDED"}
    assert not any("verified" in s for s in plain.provenance["assumptions"])
    # damaged binary
    (tmp_path / "broken").write_bytes(good.read_bytes()[:-3])
    with pytest.raises(ValueError):
        read_epw_freq(tmp_path / "broken")


@pytest.mark.skipif(not (REAL_EPW / "pb.ephmat.freq").exists(), reason="REAL fixture missing")
def test_real_epw_freq_verifies_conditions():
    """REAL EPW 5.3 freq file: grid end of .a2f is 1.1 x max fine-grid frequency."""
    f = read_epw_freq(REAL_EPW / "pb.ephmat.freq")
    assert f["nmodes"] == 3 and f["nqf"] == (6, 6, 6)
    wmax = epw_omega_max_meV(REAL_EPW / "pb.a2f")
    rep = check_epw_frequencies(f, omega_max_K=wmax * 11.604518121550082, eps_acustic_cm1=2.0)
    assert rep["omega_max_K"] == pytest.approx(1.1 * rep["max_freq_K"], rel=1e-6)
    assert rep["n_imaginary"] == 0 and rep["n_excluded_nongamma"] == 0
    bc = band_coupling_from_epw_lambda(parse_epw_lambda_fs(REAL_EPW / "pb.lambda_FS"),
                                       omega_max=wmax, freq_file=REAL_EPW / "pb.ephmat.freq",
                                       eps_acustic_cm1=2.0)
    assert all("(verified:" in s for s in bc.provenance["assumptions"])


# ----------------------------------------------------------------------------- (5) BandCoupling
def test_band_coupling_row_sum_semantics():
    with pytest.raises(ValueError, match="must not carry A_K"):
        BandCoupling(("a",), "row_sum_upper_bound", [[1.0]], row_sum_upper_K=(2.0,))
    # exact comparison: 0.1 + 0.2 > 0.3 in binary (by 2^-54 relative), must be refused
    with pytest.raises(ValueError, match="exceeds"):
        BandCoupling(("a", "b"), "computed", [[0.1, 0.2], [0.0, 0.0]], row_sum_upper_K=(0.3, 0.0))
    ok = BandCoupling(("a", "b"), "computed", [[0.1, 0.2], [0.0, 0.0]],
                      row_sum_upper_K=(0.30000000000000004, 0.0))
    assert ok.max_row_sum_K() == pytest.approx(0.3)
