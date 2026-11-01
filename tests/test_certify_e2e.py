"""End-to-end certify on the QE 7.6 benchmark spectra (fcc Al and Pb, 0 GPa; tests/data/
qe76_{al,pb}_bench, from an earlier benchmark run of this pipeline).

* T* = 300 K: core certify_spectra gives REJECTED_PROVEN (rechecked), and decide_t on inputs
  built from the same spectra gives rejectedProven by N1-ISO; the §3.2 combination rule agrees.
* Low T*: the core bracket must equal the pre-registered validation results
  (Al CERTIFIED at T* = 0.5875 K, Pb CERTIFIED at T* = 3.595 K, μ*_hi = 0.13, bracket mode;
  earlier validation run, 2026-10-28).

Fixture assumptions for the decide_t part (NOT DFT results, stated so the test is honest): the
benchmark has one q mesh, so the q ladder is declared as two meshes with δ^q = 0 and q-spread 0;
n_sheets = 1 is declared (Al and Pb really have several Fermi sheets; the M-ISO assumption is
what the outcome then records); E_F = 11.7 eV (Al) / 9.5 eV (Pb), free-electron values.
"""
from __future__ import annotations

from pathlib import Path

import pytest

from rtt.core.dfpt import collect_qe_eph
from rtt.core.status import Status
from rtt.certify import certify, certify_config, recheck_or_raise
from rtt.certify.alternatives import envelope_e, freq_scaled
from rtt.decide.core import DecideConfig, decide_t
from rtt.decide.inputs import envelope_from_ladder, tags_from_spectra
from rtt.records import DecideInputs, Prov, Validity
from rtt.targets.spec import load_spec
from rtt.units import EV_TO_K

DATA = Path(__file__).parent / "data"
OLD = {  # material: (T*, mu_hi, status, T_low, T_up) from the pre-registered validation runs
    "al": (0.5875, 0.13, "CERTIFIED_CANDIDATE", 0.8028157119046102, 10.302066039359831),
    "pb": (3.595, 0.13, "CERTIFIED_CANDIDATE", 4.515087885251041, 9.423668041719516),
}
EF_EV = {"al": 11.7, "pb": 9.5}


def spectra(m: str):
    res = collect_qe_eph(DATA / f"qe76_{m}_bench", f"{m.capitalize()}-fcc@0GPa", el_ph_nsigma=10,
                         smearing_admission="mode_positive")
    assert res.ok, res.errors
    return [s if s.rel_envelope >= 0.10 else s.with_envelope(0.10) for s in res.spectra]


def decide_inputs(m: str, sp) -> DecideInputs:
    fine = {"A_K": [s.A() for s in sp], "lambda": [s.lambda_() for s in sp],
            "m_K2": [2.0 * sum(o * w for o, w in zip(s.omega_K, s.weight_K)) for s in sp]}
    p = Prov("fullDFPT", 1)
    iso, eps = envelope_from_ladder(fine, fine, p, omega_log_by_q=[1.0, 1.0])
    val = Validity(EF=EF_EV[m] * EV_TO_K, omegaMax=max(s.omega_max() for s in sp), imagAbsMaxK=0.0,
                   containsHLiBe=False, anharmonicCheck=False, anharmShiftAbs=0.0,
                   alpha2FRecomputedAnharmonic=False, muStarNegativeSuspected=False,
                   admittedSmearings=len(sp), qMeshes=2, nSheets=1)
    return DecideInputs(tier="fullDFPT", val=val, iso=iso, eps=eps, tags=tags_from_spectra(sp))


@pytest.mark.parametrize("m", ["al", "pb"])
def test_reject_at_300K(m):
    sp = spectra(m)
    cc = certify_config(300.0, 0.15, mode="decide")
    v = certify(f"{m}-bench", sp, cc)
    assert v.status is Status.REJECTED_PROVEN
    assert "R-UNIV" in v.rule_ids
    recheck_or_raise(v, sp, cc)
    o = decide_t(DecideConfig.from_spec(load_spec()), decide_inputs(m, sp))
    assert o.outcome == "rejectedProven" and o.rule == "N1-ISO"
    assert "M-ISO" in o.assumptions
    # §3.2 combination rule: an N1-ISO rejection needs the core REJECTED_PROVEN -> agrees here


@pytest.mark.parametrize("m", ["al", "pb"])
def test_reject_at_300K_with_frequency_alternatives(m):
    sp = spectra(m)
    e = envelope_e(0.0, 0.10)
    alts = freq_scaled([s.with_envelope(e) for s in sp], 0.02)
    assert len(alts) == 3 * len(sp)
    v = certify(f"{m}-bench-fs", alts, certify_config(300.0, 0.15, mode="decide"))
    assert v.status is Status.REJECTED_PROVEN


@pytest.mark.parametrize("m", ["al", "pb"])
def test_low_T_matches_old_validation(m):
    T, mu, status, T_low, T_up = OLD[m]
    sp = spectra(m)
    cc = certify_config(T, mu, mode="bracket")
    v = certify(f"{m}-bench", sp, cc)
    assert v.status.value == status
    assert v.T_low_K == pytest.approx(T_low, rel=1e-6)
    assert v.T_up_K == pytest.approx(T_up, rel=1e-6)
    recheck_or_raise(v, sp, cc)


def test_overrides_refused():
    sp = spectra("al")
    with pytest.raises(ValueError, match="not allowed"):
        certify("x", sp, {**certify_config(300.0, 0.15, mode="decide"), "allow_unverified_rejection": True})


def test_freq_scaled_hopfield():
    sp = spectra("al")[:1]
    lo, hi = freq_scaled(sp, 0.1)[1:]
    s = sp[0]
    assert lo.A() == pytest.approx(s.A() / 0.9) and hi.A() == pytest.approx(s.A() / 1.1)
    assert lo.lambda_() == pytest.approx(s.lambda_() / 0.81)
    m = lambda x: sum(o * w for o, w in zip(x.omega_K, x.weight_K))   # noqa: E731
    assert m(lo) == pytest.approx(m(s))
    assert envelope_e(0.0, 0.10) == pytest.approx(0.10) and envelope_e(0.5, 0.5) is None
