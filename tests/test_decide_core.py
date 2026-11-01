"""decide_t (mirror of Lean TargetLogic.Decide v2): rule-by-rule truth tests (doc 01 §0.4/§0.6,
review RL1 counterexamples, R01 two-band counterexample, provenance, tags)."""
from __future__ import annotations

import dataclasses
import random
from fractions import Fraction
from pathlib import Path

import pytest

from rtt.decide import core as C
from rtt.decide.core import DecideConfig, decide_t, exact_view
from rtt.decide.inputs import mat_from_matrix, tag_from_string
from rtt.records import (D_PASS, D_REJECTED_PROVEN, D_SCREENED_OUT, D_UNDECIDED, DecideInputs, DosData, EpsData,
                         InputTag, IsoData, MatData, Prov, Validity)
from rtt.targets.spec import load_spec

SPEC = load_spec()
CFG = DecideConfig.from_spec(SPEC)
FULL = Prov("fullDFPT", 1)
FLOOR = float(CFG.floor)

VAL_KEYS = {f.name for f in dataclasses.fields(Validity)}
ISO_KEYS = {"A", "lam", "m", "dA", "dLam", "dM"}
EPS_KEYS = {"epsSys", "epsOmega", "qSpread"}


def base(prov: Prov = FULL, tier: str = "fullDFPT", mat=None, dos=None, tags=(), iso_prov=None, eps_prov=None,
         **kw) -> DecideInputs:
    """Doc-conformant single-sheet proof-grade input that passes (no rule fires)."""
    v = dict(EF=1.0e5, omegaMax=2000.0, imagAbsMaxK=0.0, containsHLiBe=False, anharmonicCheck=False,
             anharmShiftAbs=0.0, alpha2FRecomputedAnharmonic=False, muStarNegativeSuspected=False,
             admittedSmearings=5, qMeshes=2, nSheets=1)
    i = dict(A=2000.0, lam=2.5, m=8.0e6, dA=0.05, dLam=0.05, dM=0.05)
    e = dict(epsSys=0.10, epsOmega=0.02, qSpread=0.01)
    for k, val in kw.items():
        (v if k in VAL_KEYS else i if k in ISO_KEYS else e)[k] = val
    ip = iso_prov or prov
    return DecideInputs(tier=tier, val=Validity(**v), iso=IsoData(**i, provA=ip, provLam=ip, provM=ip),
                        eps=EpsData(**e, prov=eps_prov or prov), mat=mat, dos=dos, tags=list(tags))


def out(x):
    return decide_t(CFG, x)


def two_sheet_dos(w=(0.5, 0.5), dW=0.0, prov=FULL):
    return DosData(w=list(w), dW=dW, prov=prov)


def diag_mat(a11, a22=0.0, n=(0.5, 0.5), dCalA=0.0, prov=FULL):
    return mat_from_matrix([[a11, 0.0], [0.0, a22]], list(n), dCalA, prov)


# ------------------------------------------------------------------------------------ basics
def test_floor_and_table():
    assert CFG.floor == Fraction(433148145, 10**8) * 300
    assert FLOOR == pytest.approx(1299.444435)
    gt = SPEC.g_table                     # certified m0 = g_lo/T^2 from docs/theory/certificates/SUMMARY.json
    m0 = dict(zip(gt.lams, gt.m0))
    assert CFG.tbl[-1][0] == 100 and CFG.gInf == gt.m0_inf * 90000
    assert CFG.gUp(Fraction(3)) == m0[Fraction(3)] * 90000
    assert CFG.gUp(Fraction("3.1")) == m0[Fraction("3.25")] * 90000       # lambda_up = 3.25
    assert 29.9 < float(gt.m0_inf) < 29.96                                 # g(inf) ~ 2.6955e6 K^2
    assert CFG.gUp(Fraction(101)) == CFG.gInf


def test_base_passes():
    o = out(base())
    assert o.outcome == D_PASS and o.rule is None and o.would_fire == []


# ------------------------------------------------------------------------------------ N1-ISO
def test_n1_iso_rejects_with_margin():
    o = out(base(A=500.0))
    assert o.outcome == D_REJECTED_PROVEN and o.rule == "N1-ISO"
    assert o.assumptions == ["ME-MODEL", "CERT-H", "ENV-SOUND", "M-ISO"]
    assert [u["name"] for u in o.inputs_used] == ["isoA", "eps"]


def test_n1_iso_margin_clause_and_k_resolved():
    x = base(A=1000.0)
    a = float(C.Ahi(exact_view(x)))
    assert 0.75 * FLOOR < a < FLOOR
    assert out(x).outcome == D_PASS                      # margin < 25 %, no k-resolved data
    o = out(base(A=1000.0, mat=diag_mat(5000.0, 5000.0)))   # band/k-resolved data present
    assert o.outcome == D_REJECTED_PROVEN and o.rule == "N1-ISO"


def test_n1_iso_needs_single_sheet():
    o = out(base(A=500.0, nSheets=2))
    assert o.outcome == D_UNDECIDED and o.rule == "G-VALID" and "GVALID:5:multiband_without_data" in o.reasons


# ------------------------------------------------------------------------------------ multiband rules
def test_n1_dos_fires_and_uses_lower_w():
    assert float(C.Ahi(exact_view(base(A=500.0)))) / 0.5 < FLOOR
    o = out(base(A=500.0, nSheets=2, dos=two_sheet_dos((0.5, 0.5), 0.0)))
    assert o.outcome == D_REJECTED_PROVEN and o.rule == "N1-DOS"
    assert o.assumptions[-2:] == ["M-SHEET-ISO", "SHEETS-CORRECT"]
    o2 = out(base(A=500.0, nSheets=2, dos=two_sheet_dos((0.5, 0.5), 0.3)))   # w_lo = 0.35
    assert o2.outcome == D_PASS


@pytest.mark.parametrize("nS,d", [(3, DosData([0.5, 0.5], 0.0, FULL)),            # sheet count mismatch
                                  (3, DosData([0.3, 0.3, 0.3], 0.0, FULL)),       # sum 0.9
                                  (3, DosData([0.4, 0.4, 0.2], 0.6, FULL)),       # dW > 1/2
                                  (2, DosData([0.0, 1.0], 0.0, FULL))])           # w not > 0
def test_dos_wf(nS, d):
    o = out(base(nSheets=nS, dos=d))
    assert o.outcome == D_UNDECIDED and o.rule == "WELL-FORMED"


def test_n1_kappa_and_rowmax():
    o = out(base(nSheets=2, mat=diag_mat(1000.0, 900.0)))
    assert o.outcome == D_REJECTED_PROVEN and o.rule == "N1-KAPPA" and "M-DISCRETE" in o.assumptions
    # a poor CW test vector makes the kappa bound loose; the row-sum bound still fires
    m = MatData(calA=[[500.0, 500.0], [500.0, 500.0]], n=[0.5, 0.5], v=[1.0, 1e-6], dCalA=0.0, prov=FULL)
    o = out(base(nSheets=2, mat=m))
    assert o.outcome == D_REJECTED_PROVEN and o.rule == "N1-ROWMAX"


@pytest.mark.parametrize("m", [
    MatData([[100.0, 10.0], [0.0, 100.0]], [0.5, 0.5], [1.0, 1.0], 0.0, FULL),      # N_i A_ij not symmetric
    MatData([[-1.0, 0.0], [0.0, 1.0]], [0.5, 0.5], [1.0, 1.0], 0.0, FULL),        # negative entry
    MatData([[1.0, 0.0], [0.0, 1.0]], [0.5, 0.5], [1.0, 0.0], 0.0, FULL),         # v not > 0
    MatData([[1.0, 0.0], [0.0, 1.0]], [0.5, 0.5], [1.0, 1.0], 0.7, FULL),         # dCalA > 1/2
    MatData([[1.0, 0.0]], [0.5, 0.5], [1.0, 1.0], 0.0, FULL)])                    # shape
def test_mat_wf(m):
    o = out(base(nSheets=2, mat=m))
    assert o.outcome == D_UNDECIDED and o.rule == "WELL-FORMED"


# ------------------------------------------------------------------------------------ NM rules
def test_nm_g_rejects_with_num_assumptions():
    x = base(A=1150.0, lam=3.0, m=2.0e6)
    X = exact_view(x)
    assert C.Ahi(X) > CFG.floor and X.A ** 2 <= X.lam * X.m / 4          # consistent data; N1 silent
    o = out(x)
    assert o.outcome == D_REJECTED_PROVEN and o.rule == "NM-G"
    assert o.assumptions == ["ME-MODEL", "CERT-H", "ENV-SOUND", "M-ISO", "CERT-NM"]
    assert [u["name"] for u in o.inputs_used] == ["isoM", "isoLam", "eps"]


def test_nm_cert_when_A_envelope_loose():
    # consistent nominal moments (A^2 <= lam m / 4), but a wide A ladder (dA = 0.5): N1-ISO silent,
    # lam_hi m_hi < 4 floor^2 -> the certified moment floor fires
    x = base(A=1000.0, dA=0.5, lam=1.0, dLam=0.0, m=4.0e6, dM=0.0)
    X = exact_view(x)
    assert C.Ahi(X) >= CFG.floor and C.lamHi(X) * C.mHi(X) < 4 * CFG.floor ** 2
    o = out(x)
    assert o.outcome == D_REJECTED_PROVEN and o.rule == "NM-CERT"


# ------------------------------------------------------------------------------------ R01 counterexample
def r01(**kw):
    """Decoupled two-band case (doc 01 §5.1): w = (1/2, 1/2); band 1 Einstein λ = 3, Ω = 1075.3 K,
    band 2 uncoupled. Band average: A = 806 K, λ = 1.5, m = 1.73e6 K², κ = calA_11 = 1613 K."""
    return base(**{"A": 806.0, "lam": 1.5, "m": 1.73e6, "dA": 0.0, "dLam": 0.0, "dM": 0.0, "nSheets": 2, **kw})


def test_r01_isotropic_data_only_undecided():
    o = out(r01())
    assert o.outcome == D_UNDECIDED


def test_r01_with_dos_fractions_not_rejected():
    o = out(r01(dos=two_sheet_dos((0.5, 0.5))))
    assert o.outcome == D_PASS            # 806·1.1/0.98/0.5 = 1809 >= 1299.4


def test_r01_with_band_matrix_not_rejected():
    m = mat_from_matrix([[1613.0, 0.0], [0.0, 0.0]], [0.5, 0.5], 0.0, FULL)
    X = exact_view(r01(mat=m))
    assert float(C.kappaHi(X.mat, X)) >= 1613.0 * 1.1 / 0.98 * (1 - 1e-12)
    o = out(r01(mat=m))
    assert o.outcome == D_PASS


def test_r01_declared_single_band_carries_m_iso():
    o = out(r01(nSheets=1))
    assert o.outcome == D_REJECTED_PROVEN and o.rule == "N1-ISO" and "M-ISO" in o.assumptions


# ------------------------------------------------------------------------------------ G-VALID
@pytest.mark.parametrize("kw,reason", [
    (dict(EF=None), "GVALID:1:E_F_missing"),
    (dict(EF=100.0), "GVALID:1:nonadiabatic"),
    (dict(imagAbsMaxK=31.0), "GVALID:2:imaginary_modes"),
    (dict(containsHLiBe=True), "GVALID:2:anharmonic_check_missing"),
    (dict(containsHLiBe=True, anharmonicCheck=True, anharmShiftAbs=0.2, epsOmega=0.2),
     "GVALID:2:anharmonic_shift_not_recomputed"),
    (dict(muStarNegativeSuspected=True), "GVALID:3"),
    (dict(admittedSmearings=2), "GVALID:4:admitted_smearings"),
    (dict(qMeshes=1), "GVALID:4:q_meshes"),
    (dict(epsOmega=0.015), "GVALID:env:eps_omega"),
    (dict(qSpread=0.05), "GVALID:env:eps_omega"),
    (dict(anharmShiftAbs=0.05), "GVALID:env:eps_omega"),
])
def test_gvalid_failures(kw, reason):
    o = out(base(**{"A": 500.0, **kw}))
    assert o.outcome == D_UNDECIDED and o.rule == "G-VALID"
    assert any(r.startswith(reason) for r in o.reasons), o.reasons


@pytest.mark.parametrize("field", ["omegaMax", "imagAbsMaxK", "containsHLiBe", "muStarNegativeSuspected",
                                   "admittedSmearings", "nSheets", "A", "dA", "epsSys", "qSpread",
                                   "anharmShiftAbs"])
def test_missing_field_is_undecided(field):
    o = out(base(**{"A": 500.0, field: None}))
    assert o.outcome == D_UNDECIDED and o.rule == "MISSING"


def test_hydride_with_anharmonic_check():
    o = out(base(A=500.0, containsHLiBe=True, anharmonicCheck=True, anharmShiftAbs=0.2,
                 alpha2FRecomputedAnharmonic=True, epsOmega=0.2))
    assert o.outcome == D_REJECTED_PROVEN
    assert out(base(A=500.0, containsHLiBe=True, anharmonicCheck=False)).outcome == D_UNDECIDED


# ------------------------------------------------------------------------------------ RL1 counterexamples
@pytest.mark.parametrize("kw", [
    dict(A=-1.0),                                                       # ce_negA
    dict(imagAbsMaxK=-500.0),                                           # ce_sign: signed imaginary mode
    dict(containsHLiBe=True, anharmonicCheck=True, anharmShiftAbs=-0.30),   # ce_sign: softening
    dict(A=1250.0, dA=0.03, epsSys=0.0, epsOmega=0.0),                  # ce_eps0
    dict(lam=-1.0), dict(omegaMax=-1.0), dict(m=-1.0), dict(epsSys=0.05), dict(dA=0.6), dict(dA=-0.1),
    dict(qSpread=-0.01), dict(epsOmega=1.0), dict(nSheets=0), dict(EF=-5.0), dict(A=float("nan")),
    dict(admittedSmearings=3.5), dict(containsHLiBe=1),
])
def test_rl1_malformed_inputs_never_reject(kw):
    o = out(base(**{"A": 500.0, **kw}))
    assert o.outcome == D_UNDECIDED


def test_rl1_ce_eps0_with_doc_values():
    o = out(base(A=1250.0, dA=0.03, mat=diag_mat(5000.0, 5000.0)))
    assert o.outcome == D_PASS


@pytest.mark.parametrize("tier", ["ml", "coarseDFPT", "external"])
def test_rl1_ce_mix(tier):
    """Full-DFPT isotropic data, n_sheets = 2, band matrix of lower tier: the κ floor holds but the
    rule may not fire -> screenedOut (doc UNDECIDED, 'would fire, ranking only')."""
    o = out(base(nSheets=2, mat=diag_mat(800.0, 800.0, prov=Prov(tier, 1))))
    assert o.outcome == D_SCREENED_OUT and "N1-KAPPA" in o.would_fire and o.rule is None


def test_mixed_runs_do_not_combine():
    o = out(base(A=500.0, nSheets=2, dos=two_sheet_dos(prov=Prov("fullDFPT", 2))))
    assert o.outcome == D_SCREENED_OUT and "N1-DOS" in o.would_fire
    o = out(base(A=500.0, eps_prov=Prov("fullDFPT", 2)))
    assert o.outcome == D_SCREENED_OUT


@pytest.mark.parametrize("tier", ["ml", "coarseDFPT", "external"])
def test_low_tier_inputs_screen_only(tier):
    p = Prov(tier, 1)
    assert out(base(prov=p, A=500.0)).outcome == D_SCREENED_OUT
    assert out(base(prov=p)).outcome == D_PASS


# ------------------------------------------------------------------------------------ random properties
def random_inputs(rng: random.Random, low_tier_only: bool = False, sane: bool = False) -> DecideInputs:
    def maybe(v, p=0.05):
        return None if rng.random() < p else v
    if sane:
        return _sane_inputs(rng)

    def P():
        tiers = ["coarseDFPT", "ml", "external"] if low_tier_only else ["fullDFPT", "fullDFPT", "coarseDFPT", "ml"]
        return Prov(rng.choice(tiers), rng.choice([1, 1, 2]))
    ns = rng.choice([1, 1, 1, 2, 3])
    mat = dos = None
    if ns > 1 and rng.random() < 0.5:
        a = [[rng.uniform(0, 1500) for _ in range(ns)] for _ in range(ns)]
        a = [[(a[i][j] + a[j][i]) / 2 for j in range(ns)] for i in range(ns)]
        mat = mat_from_matrix(a, [1.0 / ns] * ns, rng.uniform(0, 0.3), P())
    if ns > 1 and rng.random() < 0.5:
        dos = DosData(w=[1.0 / ns] * ns, dW=rng.uniform(0, 0.3), prov=P())
    kw = dict(EF=maybe(rng.choice([50.0, 1e4, 1e5])), omegaMax=maybe(rng.uniform(-100, 3000)),
              imagAbsMaxK=maybe(rng.choice([0.0, 10.0, 40.0, -50.0])), containsHLiBe=maybe(rng.random() < 0.4),
              anharmonicCheck=rng.random() < 0.5, anharmShiftAbs=rng.choice([0.0, 0.0, 0.05, 0.2, -0.2]),
              alpha2FRecomputedAnharmonic=rng.random() < 0.5, muStarNegativeSuspected=rng.random() < 0.1,
              admittedSmearings=rng.choice([2, 3, 7]), qMeshes=rng.choice([1, 2, 2]), nSheets=maybe(ns),
              A=maybe(rng.uniform(-50, 3000)), lam=rng.uniform(-0.5, 5), m=rng.uniform(-1e5, 1e7),
              dA=rng.uniform(-0.1, 0.6), dLam=rng.uniform(0, 0.6), dM=rng.uniform(0, 0.6),
              epsSys=rng.choice([0.0, 0.1, 0.2]), epsOmega=rng.choice([0.0, 0.02, 0.1, 0.25]),
              qSpread=rng.uniform(0, 0.05))
    tags = [InputTag("A-NEG-CLIP(3K)", "negClip", massK=3.0)] if rng.random() < 0.2 else []
    return base(prov=P(), iso_prov=P() if rng.random() < 0.3 else None, mat=mat, dos=dos, tags=tags, **kw)


def test_low_tier_never_reaches_rejection():
    """Lean decide_not_proven_of_lowTier: no fullDFPT input at all -> never REJECTED_PROVEN."""
    rng = random.Random(1234)
    seen = set()
    for _ in range(4000):
        o = out(random_inputs(rng, low_tier_only=True))
        seen.add(o.outcome)
        assert o.outcome != D_REJECTED_PROVEN
    assert {D_SCREENED_OUT, D_UNDECIDED, D_PASS} <= seen


def test_rejection_implies_lean_preconditions():
    """Lean decide_rejectedProven: gates passed, not blocked, the rule's floor holds, all used
    inputs fullDFPT from one run, assumptions recorded."""
    rng = random.Random(7)
    n = 0
    for k in range(6000):
        x = random_inputs(rng, sane=k % 2 == 0)
        o = out(x)
        if o.outcome == D_REJECTED_PROVEN:
            n += 1
            X = exact_view(x)
            assert C.wf(X) and C.gvalid(X) and not o.blocked
            assert C.floor_r(CFG, X, o.rule)
            provs = [C.prov_of(X, k) for k in C.RULE_USES[o.rule]]
            assert all(p.tier == "fullDFPT" for p in provs) and len({p.run for p in provs}) == 1
            assert o.assumptions == list(C.COMMON_ASSUMPTIONS) + list(C.RULE_ASSUMPTIONS[o.rule])
            assert X.A > 0 and X.lam > 0 and X.m > 0
    assert n > 20


def test_decide_t_agrees_with_core_b():
    rng = random.Random(99)
    for _ in range(4000):
        x = random_inputs(rng)
        o = out(x)
        try:
            X = exact_view(x)
        except Exception:
            assert o.outcome == D_UNDECIDED
            continue
        g = C.wf(X) and C.gvalid(X)
        first = g and C.first_rule(CFG, X) is not None
        would = g and C.would_fire(CFG, X)
        assert o.outcome == C.core_b(g, first, would, C.blocked(x.tags))


# ------------------------------------------------------------------------------------ tags
BLOCKING = [InputTag("A-FILE-INCOMPLETE(x)", "fileIncomplete"),
            InputTag("A-FILE-COMPLETE(unverified)", "fileCompleteUnverified"),
            InputTag("A-DROP-NONPOS(2.5K)", "dropNonpos", massK=2.5),
            InputTag("A-SIGMA-EXCLUDED(sigma1:lam=1:A=2K)", "sigmaExcluded", blocksRejection=True),
            InputTag("A-SIGMA-EXCLUDED(sigma1:lam=1:A=3K)", "sigmaExcluded", blocksRejection=None),
            InputTag("A-TYPO-INCOMPLET", "unknown"), InputTag("x", "notAKind")]
NONBLOCKING = [InputTag("A-NEG-CLIP(4K)", "negClip", massK=4.0), InputTag("A-DROP-NONPOS(0K)", "dropNonpos", massK=0.0),
               InputTag("A-SIGMA-EXCLUDED(sigma2:lam=1:A=2K)", "sigmaExcluded", blocksRejection=False),
               InputTag("M-ISO-SINGLEBAND", "informational")]


@pytest.mark.parametrize("tag", BLOCKING, ids=lambda t: t.name)
def test_blocking_tag_downgrades(tag):
    o = out(base(A=500.0, tags=[tag]))
    assert o.outcome == D_UNDECIDED and o.core_outcome == D_REJECTED_PROVEN and o.blocked
    assert out(base(prov=Prov("coarseDFPT", 1), A=500.0, tags=[tag])).outcome == D_UNDECIDED
    assert out(base(tags=[tag])).outcome == D_PASS


@pytest.mark.parametrize("tag", NONBLOCKING, ids=lambda t: t.name)
def test_nonblocking_tag(tag):
    assert out(base(A=500.0, tags=[tag])).outcome == D_REJECTED_PROVEN


def test_tag_addition_is_downgrade_only():
    """Lean decide_downgrade / decide_downgrade_rank."""
    rng = random.Random(5)
    pool = BLOCKING + NONBLOCKING
    for k in range(3000):
        x = random_inputs(rng, sane=k % 2 == 0)
        o1 = out(x)
        y = dataclasses.replace(x, tags=list(x.tags) + rng.sample(pool, rng.randint(1, 3)))
        o2 = out(y)
        assert o2.outcome == o1.outcome or o2.outcome == D_UNDECIDED
        assert C.RANK[o2.outcome] <= C.RANK[o1.outcome]


def test_tag_mapper():
    assert tag_from_string("A-FILE-INCOMPLET").kind == "unknown"
    assert tag_from_string("A-NEG-CLIP(1.5K)").kind == "negClip"
    t = tag_from_string("A-SIGMA-EXCLUDED(sigma1:lam=0.4:A=10K)", ["A-SIGMA-EXCLUDED(sigma1:lam=0.4:A=10K) [x]"])
    assert t.kind == "sigmaExcluded" and t.blocksRejection is True
    assert tag_from_string("A-SIGMA-EXCLUDED(sigma1:lam=0.4:A=10K)", []).blocksRejection is False
    assert tag_from_string("CUTOFF(block...)").kind == "informational"


# ------------------------------------------------------------------------------------ monotonicity
def test_envelope_growth_never_creates_rejection():
    """Larger δ, ε_sys or ε_ω never turn a non-rejection into REJECTED_PROVEN (RL1 #7c)."""
    rng = random.Random(11)
    for _ in range(2000):
        x = random_inputs(rng, sane=rng.random() < 0.5)
        try:
            X = exact_view(x)
        except Exception:
            continue
        if out(x).outcome == D_REJECTED_PROVEN or not (C.wf(X) and C.gvalid(X)):
            continue   # the gates are validity conditions, not envelope monotonicity
        y = dataclasses.replace(
            x, iso=dataclasses.replace(x.iso, dA=x.iso.dA + rng.uniform(0, 0.1), dLam=x.iso.dLam + rng.uniform(0, 0.1),
                                       dM=x.iso.dM + rng.uniform(0, 0.1)),
            eps=dataclasses.replace(x.eps, epsSys=x.eps.epsSys + rng.uniform(0, 0.1),
                                    epsOmega=min(0.9, x.eps.epsOmega + rng.uniform(0, 0.1))))
        assert out(y).outcome != D_REJECTED_PROVEN


def _sane_inputs(rng: random.Random) -> DecideInputs:
    """Mostly valid inputs near the decision boundaries (to exercise rejections)."""
    def P():
        return Prov(rng.choice(["fullDFPT"] * 4 + ["coarseDFPT"]), rng.choice([1, 1, 1, 2]))
    ns = rng.choice([1, 1, 2])
    mat = dos = None
    if ns == 2 and rng.random() < 0.5:
        mat = diag_mat(rng.uniform(300, 2000), rng.uniform(0, 2000), prov=P())
    if ns == 2 and rng.random() < 0.7:
        w = rng.uniform(0.2, 0.8)
        dos = two_sheet_dos((w, 1 - w), rng.uniform(0, 0.2), prov=P())
    tags = [rng.choice(BLOCKING + NONBLOCKING)] if rng.random() < 0.2 else []
    p = P()
    return base(prov=p, iso_prov=P() if rng.random() < 0.2 else p, mat=mat, dos=dos, tags=tags, nSheets=ns,
                A=rng.uniform(100, 2500), lam=rng.uniform(0.5, 4), m=rng.uniform(5e5, 1e7),
                dA=rng.uniform(0, 0.3), dLam=rng.uniform(0, 0.3), dM=rng.uniform(0, 0.3),
                epsSys=rng.choice([0.1, 0.15]), epsOmega=rng.choice([0.02, 0.05]), qSpread=rng.uniform(0, 0.02))


# ------------------------------------------------------------------------------------ pending switches
def test_moment_consistency_gate():
    x = base(A=2000.0, lam=1.0, m=1.0e6)                 # A^2 > lam m / 4 (1 + 1/100): impossible moments
    o = out(x)
    assert o.outcome == D_UNDECIDED and o.rule == "WELL-FORMED"
    assert decide_t(dataclasses.replace(CFG, moment_consistency=False), x).rule == "NM-CERT"   # switch off


def test_w_lo_sum_gate():
    x = base(A=500.0, nSheets=2, dos=two_sheet_dos((0.504, 0.504), 0.0))    # sum w (1 - dW) = 1.008 > 1
    assert out(x).outcome == D_UNDECIDED and out(x).rule == "WELL-FORMED"
    assert decide_t(dataclasses.replace(CFG, w_lo_sum_check=False), x).outcome == D_REJECTED_PROVEN
    assert out(base(A=500.0, nSheets=2, dos=two_sheet_dos((0.504, 0.504), 0.01))).outcome == D_REJECTED_PROVEN


def test_switch_nmg_uncertified_assumptions():
    o = decide_t(dataclasses.replace(CFG, nmg_certified=False), base(A=1150.0, lam=3.0, m=2.0e6))
    assert o.rule == "NM-G" and o.assumptions[-3:] == ["M-ISO", "NUM-EINSTEIN-MIN", "NUM-G-MONOTONE"]


def test_switches_in_targets_toml():
    assert SPEC.decide_flags == {"w_lo_sum_check": True, "moment_consistency": True, "nmg_certified": True}


def test_gtable_from_certificates():
    import json
    d = Path(__file__).resolve().parents[1] / "docs" / "theory" / "certificates"
    gt = SPEC.g_table
    for f in d.glob("nm_*.json"):
        c = json.loads(f.read_text(encoding="utf-8"))
        m0 = Fraction(c["m0_over_T2"])
        if c["lambda"] == "inf":
            assert gt.m0_inf == m0
        else:
            assert gt.m0[gt.lams.index(Fraction(c["lambda"]))] == m0
