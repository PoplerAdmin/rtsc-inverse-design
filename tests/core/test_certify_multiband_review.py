"""Independent review tests of rtt/core/certify/multiband.py (team CERT) against the Lean obligations
of the Lean library MultiBand (lean/MultiBand).  ALL DATA SYNTHETIC.

Obligations checked: symmetric W (W-bar + recorded assumption), Collatz-Wielandt vector checked in
Arb (sub_of_K), strict kappa_hi < a0 T*, envelope (sub_env), empty-list rule (ERRATA F1, via
combine_multiband / classify), propagation of M-BAND-ISO; plus a float soundness spot check with
the wMB float oracle mb.py (set RTT_WMB_DIR; skipped if absent).

Tests marked xfail(strict=True) document defects found in review (see the reason strings); they
turn into failures once the defect is fixed, so the marker can then be removed.
"""
from __future__ import annotations

import copy
import importlib.util
import math
import sys
from fractions import Fraction
from pathlib import Path

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
sys.path.insert(0, str(Path(__file__).resolve().parent))

from synthetic import einstein  # noqa: E402
from rtt.core.constants import A0, A0_FRACTION  # noqa: E402
from rtt.core.dfpt.band_coupling import (BandCoupling, EPWLambdaStates,  # noqa: E402
                                     band_coupling_from_epw_lambda)
from rtt.core.status import Status, Verdict  # noqa: E402
from rtt.core.certify import CertifyParams, classify, recheck  # noqa: E402
from rtt.core.certify.classify import combine_multiband  # noqa: E402
from rtt.core.certify.multiband import (bc_digest, bc_payload, check_route,  # noqa: E402
                                    rmulti_evidence, rmulti_for)

A0_EXACT = Fraction(A0_FRACTION)
MB_PATH = next((Path(c) / "mb.py" for c in [__import__("os").environ.get("RTT_WMB_DIR"),
                str(Path(__file__).resolve().parents[2] / "certificates" / "wMB")]
                if c and (Path(c) / "mb.py").exists()), Path("certificates/wMB/mb.py"))


# ----------------------------------------------------------------------------- helpers
def bcm(A, N=None, e=0.0, kind="computed", R=None, label="synthetic", prov=None):
    B = len(A) if A is not None else len(R)
    return BandCoupling(tuple(f"b{i}" for i in range(B)), kind, A, N, R, {}, e, label,
                        prov if prov is not None else {"synthetic": True})


def kappa_ref(A, N):
    """Reference float lambda_max(N^-1/2 Wbar N^-1/2)."""
    A, N = np.asarray(A, float), np.asarray(N, float)
    W = N[:, None] * A
    W = 0.5 * (W + W.T)
    s = 1 / np.sqrt(N)
    return float(np.linalg.eigvalsh(s[:, None] * W * s[None, :])[-1])


def run(bc, T):
    blk = rmulti_for(bc, T)
    f, info = rmulti_evidence(blk, T)
    return blk, f, info


def conj(bcs, T):
    """S5 conjunction exactly as classify combines it (multiband model, single-band UNDECIDED)."""
    fires = [run(bc, T)[1] for bc in bcs]
    return combine_multiband(Status.UNDECIDED, fires, "multiband") == Status.REJECTED_PROVEN


def load_mb():
    if not MB_PATH.exists():
        pytest.skip("wMB oracle mb.py not present")
    spec = importlib.util.spec_from_file_location("wmb_mb_oracle", MB_PATH)
    mb = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mb)
    return mb


# ----------------------------------------------------------------------------- basic cases
def test_decoupled_bands_kappa_is_max_band_A():
    bc = bcm([[500.0, 0.0, 0.0], [0.0, 800.0, 0.0], [0.0, 0.0, 120.0]], N=(1.0, 0.3, 2.0))
    blk, f, info = run(bc, 300.0)
    assert 800.0 <= blk["kappa_hi"] <= 800.0 * (1 + 1e-8)
    assert f and info["T_bound"] >= 800.0 / A0
    assert not run(bc, 800.0 / A0 * (1 - 1e-9))[1]


def test_two_band_detailed_balance_eig_route():
    A, N = [[600.0, 150.0], [75.0, 300.0]], (1.0, 2.0)        # N1 A12 = N2 A21 exactly
    kap = (900.0 + math.sqrt(135000.0)) / 2
    blk, f, info = run(bcm(A, N), 300.0)
    rt = blk["route"]
    assert rt["route"] == "eig" and not rt["check"]["symmetrized"]
    assert kap <= blk["kappa_hi"] <= kap * (1 + 1e-8)
    assert blk["kappa_hi"] < 750.0                             # tighter than max row sum
    assert all(x > 0 for x in rt["u"])                         # Collatz-Wielandt vector > 0
    assert not any("symmetrized" in a for a in info["assumptions"])


def test_equal_row_sums():
    bc = bcm([[300.0, 200.0], [100.0, 400.0]], N=(1.0, 2.0))  # W symmetric, rows sum to 500
    blk, f, _ = run(bc, 200.0)
    assert 500.0 <= blk["kappa_hi"] <= 500.0 * (1 + 1e-8)
    assert f == (500.0 < A0 * 200.0)


def test_reducible_K():
    N = (1.0, 1.0, 2.0)
    for A in ([[400.0, 0.0, 0.0], [0.0, 300.0, 50.0], [0.0, 25.0, 200.0]],     # isolated band wins
              [[100.0, 0.0, 0.0], [0.0, 300.0, 200.0], [0.0, 100.0, 150.0]]):  # coupled block wins
        k = kappa_ref(A, N)
        blk, _, _ = run(bcm(A, N), 300.0)
        assert k <= blk["kappa_hi"] <= k * (1 + 1e-8)
    p = bc_payload(bcm(A, N))
    rt = rmulti_for(bcm(A, N), 300.0)["route"]
    assert rt["route"] == "eig" and check_route(p, rt)["ok"]


def test_zero_entries_and_zero_matrix():
    blk, f, _ = run(bcm([[0.0, 0.0], [0.0, 700.0]], N=(3.0, 1.0)), 300.0)
    assert 700.0 <= blk["kappa_hi"] <= 700.0 * (1 + 1e-8) and f
    blk, f, info = run(bcm([[0.0, 0.0], [0.0, 0.0]], N=(1.0, 1.0)), 1e-3)
    assert f and 0.0 <= blk["kappa_hi"] < 1e-200


def test_no_N_uses_row_sums_and_row_sum_kind():
    blk, f, _ = run(bcm([[600.0, 150.0], [75.0, 300.0]]), 300.0)
    # (CERT, after coordinator note (a)): without N the direct Collatz-Wielandt route on A is
    # used; it certifies rho(A) = 633.71 <= max row sum 750.
    rho = (900.0 + math.sqrt(135000.0)) / 2
    assert blk["route"]["route"] in ("cw", "rowsum") and rho <= blk["kappa_hi"] <= 750.0 * (1 + 1e-10)
    blk, f, info = run(bcm(None, kind="row_sum_upper_bound", R=(900.0, 1000.0), e=0.05), 300.0)
    assert blk["route"]["source"] == "row_sum_upper_K"
    assert 1050.0 <= blk["kappa_hi"] <= 1050.0 * (1 + 1e-10) and f


# ----------------------------------------------------------------------------- Lean obligations
def test_symmetric_W_obligation():
    """Asymmetric W -> kappa of Wbar (>= rho(A), PROOFS 5b) and the assumption is recorded."""
    A, N = [[600.0, 150.0], [75.75, 300.0]], (1.0, 2.0)
    blk, f, info = run(bcm(A, N), 300.0)
    k = kappa_ref(A, N)
    assert blk["route"]["route"] == "eig" and blk["route"]["check"]["symmetrized"]
    assert k <= blk["kappa_hi"] <= k * (1 + 1e-8)
    assert blk["kappa_hi"] >= max(abs(np.linalg.eigvals(np.asarray(A))))
    assert any("symmetrized" in a for a in info["assumptions"])


def test_collatz_wielandt_vector_is_checked_in_arb():
    bc = bcm([[600.0, 150.0], [75.0, 300.0]], N=(1.0, 2.0))
    p = bc_payload(bc)
    rt = rmulti_for(bc, 300.0)["route"]
    kap = (900.0 + math.sqrt(135000.0)) / 2
    for bad in ({"u": [1.0, 1e-3]}, {"u": [1.0, -1.0]}, {"u": [1.0]},
                {"kappa_hi": kap * (1 - 1e-12)},              # below the Perron root: impossible
                {"u": [float("nan"), 1.0]}, {"route": "eigx"}):
        r2 = dict(rt, **bad)
        assert not check_route(p, r2)["ok"], bad
    # Arb, not floats: a kappa_hi a few ulps above the certified value still passes, the exact
    # Perron root (as nearest double, below the true value) fails for its own u
    assert check_route(p, dict(rt, kappa_hi=math.nextafter(rt["kappa_hi"], math.inf)))["ok"]


def test_strict_inequality_at_the_boundary():
    """kappa_hi < a0 T* is strict and exact: with kappa_hi the double just above / below the
    exact rational a0 T*, only the latter fires."""
    T = 300.0
    target = A0_EXACT * Fraction(T)                            # 1299.444435 exactly
    lo = float(target)
    if Fraction(lo) > target:
        lo = math.nextafter(lo, 0.0)
    hi = math.nextafter(lo, math.inf)
    assert Fraction(lo) < target < Fraction(hi)
    for kh, expect in ((lo, True), (hi, False)):
        bc = bcm([[kh]], N=(1.0,))
        blk = {"payload": bc_payload(bc), "route": {"route": "rowsum", "source": "A_K",
                                                     "kappa_hi": kh}}
        blk["digest"] = bc_digest(blk["payload"])
        f, info = rmulti_evidence(blk, T)
        assert info["route_check"]["ok"] and f is expect
        assert Fraction(info["T_bound"]) >= Fraction(kh) / A0_EXACT


def test_envelope_obligation():
    A, N = [[600.0, 150.0], [75.0, 300.0]], (1.0, 2.0)
    kap = (900.0 + math.sqrt(135000.0)) / 2
    for e in (0.0, 0.1, 0.5):
        blk, _, info = run(bcm(A, N, e=e), 300.0)
        assert (1 + e) * kap <= blk["kappa_hi"] <= (1 + e) * kap * (1 + 1e-8)
        assert any(a.startswith(f"E-ENV-MB({e})") for a in info["assumptions"])
        # dropping the envelope factor from kappa_hi must fail the check (sub_env)
        if e > 0:
            p = bc_payload(bcm(A, N, e=e))
            assert not check_route(p, dict(blk["route"], kappa_hi=blk["kappa_hi"] / (1 + e)))["ok"]
    blk, _, _ = run(bcm(None, kind="row_sum_upper_bound", R=(1000.0,), e=0.2), 300.0)
    assert blk["kappa_hi"] >= 1200.0


def test_empty_list_rule():
    assert combine_multiband(Status.UNDECIDED, [], "multiband") == Status.UNDECIDED
    assert combine_multiband(Status.REJECTED_PROVEN, [], "multiband") == Status.UNDECIDED
    s = einstein(500.0, 200.0)
    v = classify([s], CertifyParams(T_star=300.0, mu_hi=0.0, mode="decide",
                                    rejection_model="multiband"), band_couplings=[])
    assert v.status == Status.UNDECIDED


def test_M_BAND_ISO_propagated():
    st = EPWLambdaStates("synthetic", np.zeros((4, 3)), np.array([1, 1, 2, 2]),
                         np.array([0.0, 0.01, -0.01, 0.0]), np.array([0.5, 0.7, 0.9, 0.4]))
    bc = band_coupling_from_epw_lambda(st, omega_max=100.0, omega_unit="meV",
                                       reduction="dos_average", degaussw_eV=0.05,
                                       groups=lambda k, b, e: f"sheet{b}")
    _, f, info = run(bc, 300.0)
    assert f and any(a.startswith("M-BAND-ISO") for a in info["assumptions"])
    v = classify([einstein(500.0, 200.0)], CertifyParams(T_star=300.0, mu_hi=0.0, mode="decide",
                                                         rejection_model="multiband"),
                 band_couplings=[bc])
    assert v.status == Status.REJECTED_PROVEN
    assert any(a.startswith("M-BAND-ISO") for a in v.assumptions)


# ----------------------------------------------------------------------------- mixed kinds, S5
def test_mixture_of_kinds_conjunction():
    comp = bcm([[600.0, 150.0], [75.0, 300.0]], N=(1.0, 2.0), label="computed")
    ub = BandCoupling.from_dict({"labels": ["s", "p"], "lambda": [[1.0, 0.25], [0.125, 0.5]],
                                 "omega_max": 1000.0, "N": [1.0, 2.0], "unit": "K"})
    assert ub.kind == "entrywise_upper_bound"
    rs = bcm(None, kind="row_sum_upper_bound", R=(700.0, 650.0))
    k_ub = kappa_ref(ub.A_K, ub.N)
    blk, _, info = run(ub, 300.0)
    assert k_ub <= blk["kappa_hi"] <= k_ub * (1 + 1e-8)
    assert any("entrywise_upper_bound" in a for a in info["assumptions"])
    assert conj([comp, ub, rs], 300.0)
    hot = bcm(None, kind="row_sum_upper_bound", R=(2000.0,))
    assert not conj([comp, ub, hot], 300.0)                    # one alternative blocks (S5)
    s = einstein(500.0, 200.0)
    prm = CertifyParams(T_star=300.0, mu_hi=0.0, mode="decide", rejection_model="multiband")
    assert classify([s], prm, band_couplings=[comp, ub, rs]).status == Status.REJECTED_PROVEN
    assert classify([s], prm, band_couplings=[comp, ub, hot]).status == Status.UNDECIDED


# ----------------------------------------------------------------------------- tampering
def _verdict(bcs, T=300.0):
    s = einstein(500.0, 200.0)
    v = classify([s], CertifyParams(T_star=T, mu_hi=0.0, mode="decide",
                                    rejection_model="multiband"), band_couplings=bcs)
    assert v.status == Status.REJECTED_PROVEN
    return s, Verdict.from_json(v.to_json())


def test_tampered_certificates_fail_recheck():
    bc = bcm([[600.0, 150.0], [75.0, 300.0]], N=(1.0, 2.0), e=0.05)
    rs = bcm(None, kind="row_sum_upper_bound", R=(900.0, 1000.0))
    s, w = _verdict([bc, rs])
    assert recheck(w, [s], band_couplings=[bc, rs])["ok"]

    def lower_kappa(m):
        m[0]["route"]["kappa_hi"] = m[0]["route"]["kappa_hi"] / 1.05

    def wrong_u(m):
        m[0]["route"]["u"] = [1.0, 1e-3]

    def low_rowsum(m):
        m[1]["route"]["kappa_hi"] = 950.0

    def drop_route(m):
        m[0]["route"] = None

    def eig_route_on_rowsum_kind(m):
        m[1]["route"] = copy.deepcopy(m[0]["route"])

    def payload_scaled(m):                                     # digest recomputed consistently
        p = m[0]["payload"]
        p["A_K"] = [[x / 2 for x in r] for r in p["A_K"]]
        m[0]["digest"] = bc_digest(p)

    for fn in (lower_kappa, wrong_u, low_rowsum, drop_route, eig_route_on_rowsum_kind,
               payload_scaled):
        w2 = copy.deepcopy(w)
        fn(w2.certificate["multiband"])
        assert not recheck(w2, [s], band_couplings=[bc, rs])["ok"], fn.__name__
    w2 = copy.deepcopy(w)
    w2.certificate["multiband"].pop()
    assert not recheck(w2, [s], band_couplings=[bc])["ok"]
    assert not recheck(w, [s], band_couplings=[rs, bc])["ok"]                     # order
    assert not recheck(w, [s], band_couplings=[bcm([[600.0, 150.0], [75.0, 300.0]],
                                                   N=(1.0, 2.0)), rs])["ok"]      # envelope
    assert not recheck(w, [s])["ok"]


# ----------------------------------------------------------------------------- soundness (oracle)
def _tc_float(mb, atoms, N, T_hi, nmax=160):
    f = lambda T: mb.topeig([(Om / T, np.asarray(As) / T) for Om, As in atoms], N, nmax)
    grid = T_hi * np.geomspace(1.0, 1e-3, 60)
    top = next((T for T in grid if f(T) >= 0), None)
    if top is None:
        return 0.0, f
    lo, hi = top, top * grid[0] / grid[1]
    while hi / lo - 1 > 1e-11:
        mid = math.sqrt(lo * hi)
        lo, hi = (mid, hi) if f(mid) >= 0 else (lo, mid)
    return lo, f


def _random_model(rng, B, n_atoms):
    N = rng.uniform(0.2, 2.0, B)
    atoms = []
    for _ in range(n_atoms):
        W = rng.uniform(0, 1, (B, B)) * (rng.uniform(0, 1, (B, B)) < 0.7)
        W = W + W.T                                            # detailed balance
        atoms.append((float(rng.uniform(100, 1500)), W / N[:, None] * rng.uniform(50, 400)))
    A = sum(a for _, a in atoms)
    return atoms, tuple(float(x) for x in N), [[float(x) for x in r] for r in A]


@pytest.mark.parametrize("seed", range(8))
def test_soundness_vs_wMB_oracle(seed):
    """Whenever R-MULTI fires at T*, the wMB float multiband Tc (mu*=0) is < T* and T* is
    float-stable.  Random SYNTHETIC detailed-balanced models, 2-4 bands, 1-3 atoms."""
    mb = load_mb()
    rng = np.random.default_rng(4242 + seed)
    atoms, N, A = _random_model(rng, 2 + seed % 3, 1 + seed % 3)
    bc = bcm(A, N, e=0.02 * (seed % 2))
    _, _, info = run(bc, 1.0)
    tb = info["T_bound"]
    tc, f = _tc_float(mb, atoms, N, tb * 1.5)
    assert 0 < tc < tb
    for fac in (1 + 1e-9, 1.05, 2.0):
        T_star = tb * fac
        fires = run(bc, T_star)[1]
        assert fires
        assert tc < T_star and f(T_star) < 0
    # no fire below the certified bound (the rule is a pure kappa_hi test)
    assert not run(bc, tb * (1 - 1e-9))[1]


def test_soundness_near_sharp_einstein():
    """Near-extremal case (wMB 4.6): one Einstein atom, equal row sums kappa = 1000 K,
    Omega/Tc = 7.5748661: the float Tc sits within ~1e-8 of kappa/a0 and must stay below it."""
    mb = load_mb()
    N = (1.0, 3.0)
    W = np.array([[400.0, 600.0], [600.0, 2400.0]])
    A = (W / np.asarray(N)[:, None]).tolist()
    atoms = [(7.5748661 * 1000.0 / 4.33148149, A)]
    _, _, info = run(bcm(A, N), 1.0)
    tb = info["T_bound"]
    assert 1000.0 / A0 <= tb <= 1000.0 / A0 * (1 + 1e-8)
    tc, f = _tc_float(mb, atoms, N, tb * 1.01, nmax=200)
    assert tb * (1 - 1e-6) < tc < tb
    assert not run(bcm(A, N), tc)[1]
    assert run(bcm(A, N), tb * (1 + 1e-10))[1] and f(tb * (1 + 1e-10)) < 0


# ----------------------------------------------------------------------------- review findings
def _epw_max_bc():
    st = EPWLambdaStates("synthetic", np.zeros((4, 3)), np.array([1, 1, 2, 2]),
                         np.array([0.0, 0.01, -0.01, 0.0]), np.array([0.5, 0.7, 0.9, 0.4]))
    return band_coupling_from_epw_lambda(st, omega_max=100.0, omega_unit="meV", reduction="max",
                                         eps_acustic_cm1=5.0)


def test_R1_data_caveats_propagated():
    bc = _epw_max_bc()
    assert bc.provenance["caveats"]
    _, f, info = run(bc, 300.0)
    assert f
    text = " ".join(info["assumptions"])
    assert "eps_acustic" in text and "omega_max" in text


def test_R2_stripped_assumption_detected_by_recheck():
    st = EPWLambdaStates("synthetic", np.zeros((4, 3)), np.array([1, 1, 2, 2]),
                         np.array([0.0, 0.01, -0.01, 0.0]), np.array([0.5, 0.7, 0.9, 0.4]))
    bc = band_coupling_from_epw_lambda(st, omega_max=100.0, omega_unit="meV",
                                       reduction="dos_average", degaussw_eV=0.05)
    s, w = _verdict([bc])
    w.certificate["multiband"][0]["payload"]["data_assumptions"] = []
    w.assumptions = [a for a in w.assumptions if not a.startswith("M-BAND-ISO")]
    assert not recheck(w, [s], band_couplings=[bc])["ok"]


def test_R3_row_sum_kind_uses_only_declared_bounds():
    # (coordinator) BandCoupling now forbids A_K for kind "row_sum_upper_bound" (ambiguous
    # meaning, root cause of R3); the mixed object can no longer be constructed.
    with pytest.raises(ValueError):
        BandCoupling(("a", "b"), "row_sum_upper_bound", [[600.0, 150.0], [75.0, 300.0]],
                     (1.0, 2.0), (2000.0, 2000.0))
    bc = BandCoupling(("a", "b"), "row_sum_upper_bound", None, (1.0, 2.0), (2000.0, 2000.0))
    blk, f, _ = run(bc, 300.0)
    assert blk["route"]["source"] == "row_sum_upper_K" and blk["kappa_hi"] >= 2000.0
    assert not f


@pytest.mark.parametrize("T", [math.inf, math.nan])
def test_R4_nonfinite_T_star(T):
    f, _ = rmulti_evidence(rmulti_for(bcm([[600.0]], N=(1.0,)), 300.0), T)
    assert f is False
