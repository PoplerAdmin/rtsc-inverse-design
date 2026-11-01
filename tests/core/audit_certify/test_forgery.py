"""Referee: certificate-forgery attacks on certify.recheck.

Each test builds a genuine Verdict, tampers with it, and asserts that recheck REJECTS the
tampered verdict.  A failing test is a demonstrated flaw (see the docstring of each test).
Run: nice python -m pytest -q tests/audit_certify/test_forgery.py
"""
from __future__ import annotations

import copy
import json
import sys

import pytest


from rtt.core.spectrum import Spectrum  # noqa: E402
from rtt.core.status import Status, Verdict  # noqa: E402
from rtt.core.certify import classify, CertifyParams, recheck  # noqa: E402
from rtt.core.certify.cells import check_cell  # noqa: E402
from rtt.core.certify.arbcore import ArbSpectrum  # noqa: E402

E1000 = Spectrum((1000.0,), (1250.0,), 0.0, "E1000")


@pytest.fixture(scope="module")
def rej():
    v = classify([E1000], CertifyParams(mode="bracket", time_budget_s=30, mu_hi=0.0))
    assert v.status == Status.REJECTED_PROVEN
    assert len(v.certificate["alternatives"][0]["cells"]["cells"]) >= 2
    return v


@pytest.fixture(scope="module")
def cert200():
    # same spectrum, T* = 200 K: Tc ~ 247.6 K -> CERTIFIED with mu_hi = 0
    v = classify([E1000], CertifyParams(T_star=200.0, mode="bracket", time_budget_s=30, mu_hi=0.0))
    assert v.status == Status.CERTIFIED_CANDIDATE
    return v


def clone(v):
    return Verdict.from_json(v.to_json())


def test_roundtrip_ok(rej):
    assert recheck(clone(rej), [E1000])["ok"]


def test_gap_between_cells(rej):
    v = clone(rej)
    cells = v.certificate["alternatives"][0]["cells"]["cells"]
    t = float.fromhex(cells[1]["T1"])
    import math
    cells[1]["T1"] = math.nextafter(t, 1e9).hex()          # one-ulp gap
    assert not recheck(v, [E1000])["ok"]


def test_missing_top_cell(rej):
    v = clone(rej)
    v.certificate["alternatives"][0]["cells"]["cells"].pop()
    assert not recheck(v, [E1000])["ok"]


def test_negative_head_entry(rej):
    v = clone(rej)
    v.certificate["alternatives"][0]["cells"]["cells"][0]["u"][0] = -1.0
    assert not recheck(v, [E1000])["ok"]


def test_inconsistent_tail_constant(rej):
    v = clone(rej)
    c = v.certificate["alternatives"][0]["cells"]["cells"][0]
    c["C"] = c["C"] * 1e-6
    assert not recheck(v, [E1000])["ok"]


def test_stored_check_booleans_ignored(rej):
    v = clone(rej)
    c = v.certificate["alternatives"][0]["cells"]["cells"][0]
    c["u"] = [1.0] * len(c["u"])                         # garbage vector
    c["check"] = {"ok": True}
    v.certificate["alternatives"][0]["evidence"] = {"univFires": True}
    assert not recheck(v, [E1000])["ok"]


def test_wrong_spectrum(rej):
    other = Spectrum((1000.0,), (1250.0000001,), 0.0)
    assert not recheck(clone(rej), [other])["ok"]


def test_univ_from_other_spectrum():
    small = Spectrum((1000.0,), (100.0,), 0.0)
    v = classify([small], CertifyParams(mode="decide"))
    assert v.status == Status.REJECTED_PROVEN
    big = Spectrum((1000.0,), (5000.0,), 0.0)
    vv = clone(v)
    vv.certificate["alternatives"][0]["spectrum_sha256"] = __import__(
        "rtt.core.certify.classify", fromlist=["x"]).spectrum_digest(big)
    assert not recheck(vv, [big])["ok"]


def test_lower_T_up(rej):
    v = clone(rej)
    blk = v.certificate["alternatives"][0]["cells"]
    blk["T_up"] = (float.fromhex(blk["T_up"]) * 0.9).hex()
    v.T_up_K = float.fromhex(blk["T_up"])
    assert not recheck(v, [E1000])["ok"]


# ------------------------------------------------------------------ binding of run parameters
def test_T_star_rebinding(cert200):
    """FLAW (binding): the run parameters live only inside the certificate; recheck takes no
    expected T*/mu_hi/omega_c.  A CERTIFIED verdict obtained at T* = 200 K is turned into a
    REJECTED_PROVEN verdict by editing certificate.params.T_star (no proof is forged, but the
    re-checked verdict no longer says what the run asked).  A consumer that trusts
    recheck()['ok'] and the run's configured T* gets a 'rejected' material whose Tc is ~247 K
    >= the configured T* = 200 K."""
    v = clone(cert200)
    v.certificate["params"]["T_star"] = 300.0
    v.status = Status.REJECTED_PROVEN
    # the stored cells reach below 300 K, so the rechecked status is REJECTED as well
    r = recheck(v, [E1000])
    # recheck has no way of knowing the run used T* = 200: report what it says
    assert not r["ok"], f"recheck accepted a verdict re-targeted to T*=300: {r}"


def test_mu_hi_rebinding():
    """FLAW (binding): mu_hi / omega_c for CERTIFIED are taken from certificate.params, not from
    the caller; the assumptions text (U-MU, CUTOFF-BLOCK) is not compared at all."""
    sp = Spectrum((1000.0,), (1700.0,), 0.0)
    v = classify([sp], CertifyParams(mode="bracket", time_budget_s=30, mu_hi=0.0))
    assert v.status == Status.CERTIFIED_CANDIDATE
    vv = clone(v)
    vv.assumptions = [a.replace("U-MU[0,0.0]", "U-MU[0,0.15]") for a in vv.assumptions]
    r = recheck(vv, [sp])
    assert not r["ok"], "recheck accepted a CERTIFIED verdict whose assumptions claim mu_hi=0.15"


def test_cform_assumption_stripped():
    """FLAW (provenance): recheck accepts cform cells (hand proof, NOT Lean) even when the verdict
    does not carry the 'R-CELLS[cform]' assumption and params.row_check says 'mhat'.
    NOTE: in this instance the verdict is carried by R-UNIV (the cform search found no cells),
    so the test shows only that assumption stripping is undetected; the cform-acceptance part
    is by code inspection (classify.evidence_for_alternative accepts any variant in VARIANTS)."""
    v = classify([E1000], CertifyParams(mode="decide", row_check="cform", mu_hi=0.0))
    assert v.status == Status.REJECTED_PROVEN
    vv = clone(v)
    vv.assumptions = [a for a in vv.assumptions if "cform" not in a]
    vv.certificate["params"]["row_check"] = "mhat"
    r = recheck(vv, [E1000])
    assert not r["ok"], "recheck accepted cform cells presented as Lean-backed mhat cells"


def test_singleband_flag_stripped(rej):
    """FLAW (provenance): MULTIBAND_RISK / M-ISO-SINGLEBAND can be stripped; recheck still ok."""
    vv = clone(rej)
    vv.flags = []
    vv.assumptions = [a for a in vv.assumptions if "SINGLEBAND" not in a]
    assert not recheck(vv, [E1000])["ok"]
