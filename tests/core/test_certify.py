"""Tests for rtt.core.certify (team CERT).  All spectra are SYNTHETIC.

Run: nice python3 -m pytest -q tests/core/test_certify.py
"""
from __future__ import annotations

import math
import sys
from fractions import Fraction
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
sys.path.insert(0, str(Path(__file__).resolve().parent))

from synthetic import einstein, float_tc_reference, multipeak  # noqa: E402
from rtt.core.constants import A0_FRACTION  # noqa: E402
from rtt.core.spectrum import Spectrum  # noqa: E402
from rtt.core.status import Status, Verdict  # noqa: E402
from rtt.core.certify import (AltEvidence, ArbSpectrum, CertifyParams, check_cell,  # noqa: E402
                          check_rayleigh, check_univ, classify, classify_evidence, classifyAlt,
                          recheck)


# ------------------------------------------------------------------ Lean mirror / contract
def test_classify_empty_is_undecided():               # ERRATA F1
    v = classify([], CertifyParams())
    assert v.status == Status.UNDECIDED
    assert classify_evidence([]) == Status.UNDECIDED


@pytest.mark.parametrize("u,c,r,exp", [
    (True, False, False, Status.REJECTED_PROVEN), (True, True, True, Status.REJECTED_PROVEN),
    (False, True, False, Status.REJECTED_PROVEN), (False, True, True, Status.REJECTED_PROVEN),
    (False, False, True, Status.CERTIFIED_CANDIDATE), (False, False, False, Status.UNDECIDED)])
def test_classifyAlt_truth_table(u, c, r, exp):
    assert classifyAlt(AltEvidence(u, c, r)) == exp


def test_classify_conjunction():
    R = AltEvidence(True, False, False)
    Cc = AltEvidence(False, False, True)
    U = AltEvidence(False, False, False)
    assert classify_evidence([R, R]) == Status.REJECTED_PROVEN
    assert classify_evidence([Cc, Cc]) == Status.CERTIFIED_CANDIDATE
    assert classify_evidence([R, Cc]) == Status.UNDECIDED
    assert classify_evidence([R, U]) == Status.UNDECIDED


# ------------------------------------------------------------------ R-UNIV exactness
def test_univ_strict_and_exact():
    a0 = A0_FRACTION
    T = 300.0
    A_eq = a0 * Fraction(T)                     # exactly a0 * T*  (not a double in general)
    A_lo = math.nextafter(float(A_eq), 0.0)
    while Fraction(A_lo) >= A_eq:
        A_lo = math.nextafter(A_lo, 0.0)
    A_hi = math.nextafter(float(A_eq), math.inf)
    while Fraction(A_hi) <= A_eq:
        A_hi = math.nextafter(A_hi, math.inf)
    assert check_univ(Spectrum((100.0,), (A_lo,)), T) is True
    assert check_univ(Spectrum((100.0,), (A_hi,)), T) is False
    # envelope enters as (1+e) A
    assert check_univ(Spectrum((100.0,), (A_lo,), 0.01), T) is False


# ------------------------------------------------------------------ end to end
@pytest.fixture(scope="module")
def einstein_case():
    s = einstein(500.0, 200.0)                  # lambda = 0.8, Tc ~ 34 K  (synthetic)
    v = classify([s], CertifyParams(T_star=30.0, mu_hi=0.0))
    return s, v


def test_einstein_bracket_contains_float_tc(einstein_case):
    s, v = einstein_case
    ref = float_tc_reference(s, N=400)
    assert v.T_low_K <= ref <= v.T_up_K
    assert (v.T_up_K - v.T_low_K) / ref < 1e-3
    assert v.status == Status.CERTIFIED_CANDIDATE       # T_low >= T* = 30 K
    assert any("CUTOFF-BLOCK" in a for a in v.assumptions)   # ERRATA F3


def test_recheck_roundtrip(einstein_case):
    s, v = einstein_case
    w = Verdict.from_json(v.to_json())
    r = recheck(w, [s])
    assert r["ok"], r["errors"]


def _tampered(v, f):
    w = Verdict.from_json(v.to_json())
    f(w)
    return w


def test_recheck_rejects_tampering(einstein_case):
    s, v = einstein_case
    alt = lambda w: w.certificate["alternatives"][0]

    def bump_u(w):
        c = alt(w)["cells"]["cells"][0]
        c["u"][len(c["u"]) // 2] *= 2.0

    def bump_eta(w):
        alt(w)["cells"]["cells"][0]["eta"] *= 1e4

    def stretch(w):
        c = alt(w)["cells"]["cells"][0]
        T1 = float.fromhex(c["T1"]) * 0.995
        c["T1"] = T1.hex()
        alt(w)["cells"]["T_up"] = T1.hex()
        w.T_up_K = T1

    def ray_up(w):
        r = alt(w)["rayleigh"]
        T = float.fromhex(r["T"]) * 1.01
        r["T"] = T.hex()
        w.T_low_K = T

    def flip(w):
        w.status = Status.REJECTED_PROVEN

    for f in (bump_u, bump_eta, stretch, ray_up, flip):
        r = recheck(_tampered(v, f), [s])
        assert not r["ok"], f.__name__
    # other spectrum, same certificate
    s2 = Spectrum(s.omega_K, (s.weight_K[0] * 1.001,), s.rel_envelope)
    assert not recheck(Verdict.from_json(v.to_json()), [s2])["ok"]


def test_check_cell_structural_rejections(einstein_case):
    s, v = einstein_case
    asp = ArbSpectrum(s.omega_K, s.weight_K, s.rel_envelope)
    c0 = v.certificate["alternatives"][0]["cells"]["cells"][0]
    assert check_cell(asp, c0)["ok"]
    bad = []
    for k, val in [("eta", 0.0), ("C", -1.0), ("N_e", c0["N"]), ("M_v", c0["N_e"]),
                   ("variant", "xx"), ("p", 4)]:
        c = dict(c0)
        c[k] = val
        bad.append(c)
    c = dict(c0)
    c["u"] = list(c0["u"])
    c["u"][0] = -c["u"][0]
    bad.append(c)
    c = dict(c0)
    c["T2"] = (float.fromhex(c0["T1"]) * 0.5).hex()
    bad.append(c)
    for c in bad:
        assert not check_cell(asp, c)["ok"]


def test_rayleigh_above_tc_fails(einstein_case):
    s, v = einstein_case
    asp = ArbSpectrum(s.omega_K, s.weight_K, s.rel_envelope)
    rc = dict(v.certificate["alternatives"][0]["rayleigh"])
    assert check_rayleigh(asp, rc)["ok"]
    rc["T"] = (float.fromhex(rc["T"]) * 1.02).hex()
    from rtt.core.certify.rayleigh import certified_n_c
    rc["N_c"] = certified_n_c(float.fromhex(rc["T"]), rc["omega_c"])
    assert not check_rayleigh(asp, rc)["ok"]


def test_mu_and_envelope_bracket():
    s = einstein(300.0, 150.0, e=0.05)          # synthetic
    v = classify([s], CertifyParams(T_star=300.0, mu_hi=0.1))
    wc = 10 * 300.0
    for scale, mu in [(1.0, 0.0), (1.0, 0.1), (0.95, 0.1), (1.05, 0.0)]:
        ref = float_tc_reference(s, mu=mu, omega_c=wc, N=300, scale=scale)
        assert v.T_low_K <= ref <= v.T_up_K, (scale, mu, ref, v.T_low_K, v.T_up_K)
    assert v.status == Status.REJECTED_PROVEN   # A_hi < a0 * 300
    assert "R-UNIV" in v.rule_ids


def test_two_alternatives_s5():
    s1 = einstein(500.0, 200.0)                  # Tc ~ 34 K
    s2 = einstein(500.0, 80.0)                   # Tc ~ 10 K
    v = classify([s1, s2], CertifyParams(T_star=30.0, mu_hi=0.0))
    assert v.status == Status.UNDECIDED          # certified for s1, rejected for s2
    t1 = float_tc_reference(s1, N=300)
    t2 = float_tc_reference(s2, N=300)
    assert v.T_low_K <= min(t1, t2) and v.T_up_K >= max(t1, t2)     # ERRATA F7
    assert recheck(Verdict.from_json(v.to_json()), [s1, s2])["ok"]


def test_multipeak_small():
    s = multipeak(5, 1.2, n_pts=300, om_lo=30.0, om_hi=900.0)
    v = classify([s], CertifyParams(T_star=300.0, mu_hi=0.0))
    ref = float_tc_reference(s, N=300)
    assert v.T_low_K <= ref <= v.T_up_K
    assert v.status == Status.REJECTED_PROVEN


def test_G8_cform_cells_rejected_unless_declared(einstein_case):
    s, v = einstein_case
    w = Verdict.from_json(v.to_json())
    for c in w.certificate["alternatives"][0]["cells"]["cells"]:
        c["variant"] = "cform"                     # row_check of the run is "mhat"
    r = recheck(w, [s])
    assert not r["ok"]
    assert any("variant" in ch for a in r["alternatives"] for ch in a["checks"])


def test_expected_params_binding(einstein_case):
    s, v = einstein_case
    w = Verdict.from_json(v.to_json())
    assert recheck(w, [s], expected_params=CertifyParams(T_star=30.0, mu_hi=0.0))["ok"]
    for bad in (CertifyParams(T_star=31.0, mu_hi=0.0), CertifyParams(T_star=30.0, mu_hi=0.1),
                CertifyParams(T_star=30.0, mu_hi=0.0, omega_c=1234.0),
                CertifyParams(T_star=30.0, mu_hi=0.0, row_check="cform"),
                CertifyParams(T_star=30.0, mu_hi=0.0, rejection_model="multiband")):
        assert not recheck(w, [s], expected_params=bad)["ok"]
    assert any(a.startswith("TARGET(T*=30.0") for a in v.assumptions)
