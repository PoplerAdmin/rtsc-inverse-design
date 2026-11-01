"""R-MULTI tests (ARCH 4c; pending audit, behind rejection_model="multiband").  SYNTHETIC data."""
from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
sys.path.insert(0, str(Path(__file__).resolve().parent))

from synthetic import einstein  # noqa: E402
from rtt.core.dfpt.band_coupling import BandCoupling  # noqa: E402
from rtt.core.status import Status, Verdict  # noqa: E402
from rtt.core.certify import CertifyParams, classify, recheck  # noqa: E402
from rtt.core.certify.multiband import (bc_payload, check_route, propose_rmulti,  # noqa: E402
                                    rmulti_evidence, rmulti_for)


def two_band(scale=1.0, e=0.0, kind="computed", asym=False):
    A = [[600.0 * scale, 150.0 * scale], [75.0 * scale, 300.0 * scale]]   # N1 A12 = N2 A21
    if asym:
        A[1][0] *= 1.01
    return BandCoupling(labels=("s", "p"), kind=kind, A_K=A, N=(1.0, 2.0), rel_envelope=e,
                        label="synthetic two-band")


def test_eig_route_and_rowsum():
    bc = two_band()
    blk = rmulti_for(bc, 300.0)
    assert blk["route"]["route"] == "eig"
    # kappa: K = [[600, 150/sqrt2],[150/sqrt2, 300]] -> 633.7 K; row sum bound = 750
    assert 633.0 < blk["kappa_hi"] < 634.5          # (900 + sqrt(135000))/2 = 633.71
    fires, info = rmulti_evidence(blk, 300.0)
    assert fires and info["T_bound"] < 300.0
    fires2, _ = rmulti_evidence(blk, 0.99 * blk["kappa_hi"] / 4.33148145)
    assert not fires2


def test_rmulti_strict_and_envelope():
    bc = two_band(e=0.1)
    blk = rmulti_for(bc, 300.0)
    assert blk["kappa_hi"] >= 1.1 * 633.7
    p = bc_payload(bc)
    bad = dict(blk["route"])
    bad["kappa_hi"] = blk["kappa_hi"] / 1.1          # forgets the envelope
    assert not check_route(p, bad)["ok"]
    bad = dict(blk["route"])
    bad["u"] = [1.0, 1e-3]                            # wrong vector
    assert not check_route(p, bad)["ok"]


def test_symmetrized_assumption_and_rowsum_kind():
    blk = rmulti_for(two_band(asym=True), 300.0)
    _, info = rmulti_evidence(blk, 300.0)
    assert any("symmetrized" in a for a in info["assumptions"])
    bc = BandCoupling(labels=("a", "b"), kind="row_sum_upper_bound", row_sum_upper_K=(900.0, 1000.0))
    blk = rmulti_for(bc, 300.0)
    assert blk["route"]["route"] == "rowsum" and blk["kappa_hi"] >= 1000.0
    f, info = rmulti_evidence(blk, 300.0)
    assert f


def test_classify_multiband_mode():
    s = einstein(500.0, 200.0)                       # single band: Tc ~ 34 K  (synthetic)
    bc = two_band()
    # isotropic mode: single-band rejection carries M-ISO-SINGLEBAND + MULTIBAND_RISK
    v = classify([s], CertifyParams(T_star=300.0, mu_hi=0.0, mode="decide"))
    assert v.status == Status.REJECTED_PROVEN
    assert any("M-ISO-SINGLEBAND" in a for a in v.assumptions) and "MULTIBAND_RISK" in v.flags
    # multiband mode without band data -> UNDECIDED
    v = classify([s], CertifyParams(T_star=300.0, mu_hi=0.0, mode="decide",
                                    rejection_model="multiband"))
    assert v.status == Status.UNDECIDED
    # with band data -> R-MULTI rejection, recheckable
    v = classify([s], CertifyParams(T_star=300.0, mu_hi=0.0, mode="decide",
                                    rejection_model="multiband"), band_couplings=[bc])
    assert v.status == Status.REJECTED_PROVEN and "R-MULTI" in v.rule_ids
    assert any("PSD" in a for a in v.assumptions)
    w = Verdict.from_json(v.to_json())
    assert recheck(w, [s], band_couplings=[bc])["ok"]
    assert not recheck(w, [s], band_couplings=[two_band(scale=1.5)])["ok"]
    assert not recheck(w, [s])["ok"]
    # a band coupling with kappa above a0 T* -> UNDECIDED
    v = classify([s], CertifyParams(T_star=100.0, mu_hi=0.0, mode="decide",
                                    rejection_model="multiband"), band_couplings=[bc])
    assert v.status == Status.UNDECIDED


def test_contribution_normalised_data_refused():
    bc = BandCoupling(labels=("s", "p"), kind="computed", A_K=[[600.0, 150.0], [75.0, 300.0]],
                      N=(1.0, 2.0), provenance={"normalisation": "contribution (sums to isotropic a2F)"})
    blk = rmulti_for(bc, 300.0)
    assert blk["route"] is None and "refused" in blk
    f, info = rmulti_evidence(blk, 300.0)
    assert not f
    # tampering with the stored refusal breaks the digest
    blk2 = dict(blk)
    blk2["payload"] = dict(blk["payload"], refusal=None)
    f2, info2 = rmulti_evidence(blk2, 300.0)
    assert not f2 and "digest" in info2["reason"]


def test_cw_route_without_dos():
    bc = BandCoupling(labels=("s", "p"), kind="entrywise_upper_bound",
                      A_K=[[600.0, 150.0], [75.0, 300.0]])          # no N given
    blk = rmulti_for(bc, 300.0)
    assert blk["route"]["route"] == "cw"
    assert 633.0 < blk["kappa_hi"] < 634.5                          # rho(A) = 633.71
    f, info = rmulti_evidence(blk, 300.0)
    assert f and any("upper bound" in a for a in info["assumptions"])
    bad = dict(blk["route"])
    bad["v"] = [1.0, 1.0]
    from rtt.core.certify.multiband import check_route
    assert not check_route(bc_payload(bc), bad)["ok"]


def test_G12_empty_band_matrix_never_fires():
    for payload in ({"kind": "computed", "A_K": [], "N": [], "row_sum_upper_K": None},
                    {"kind": "row_sum_upper_bound", "A_K": None, "N": None, "row_sum_upper_K": []}):
        payload.update(rel_envelope=0.0, label="", labels=[], data_assumptions=[], refusal=None)
        from rtt.core.certify.multiband import bc_digest
        for rt in ({"route": "rowsum", "source": "A_K", "kappa_hi": 1.0},
                   {"route": "rowsum", "source": "row_sum_upper_K", "kappa_hi": 1.0},
                   {"route": "eig", "kappa_hi": 1.0, "u": []}, {"route": "cw", "kappa_hi": 1.0, "v": []}):
            assert not check_route(payload, rt)["ok"]
            blk = {"payload": payload, "digest": bc_digest(payload), "route": rt}
            assert rmulti_evidence(blk, 300.0)[0] is False


def test_G7_undecided_multiband_carries_singleband_tag():
    s = einstein(500.0, 200.0)
    v = classify([s], CertifyParams(T_star=300.0, mu_hi=0.0, mode="decide",
                                    rejection_model="multiband"))
    assert v.status == Status.UNDECIDED and v.T_up_K is not None
    assert sum(a.startswith("M-ISO-SINGLEBAND") for a in v.assumptions) == 1
