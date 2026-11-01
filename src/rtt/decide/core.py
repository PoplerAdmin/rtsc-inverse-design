"""Python mirror of Lean ``TargetLogic.Decide.decideT`` (lean/TargetLogic/Decide.lean, v2 after
review RL1; doc 01 §0.4 YAML and §0.6).

    decideCore = if wf x && gvalid x then
                   match firstRule (floorR && provOK, order KAPPA, ROWMAX, DOS, ISO, NM-CERT, NM-G)
                   | some r => rejectedProven r
                   | none   => if wouldFire (some floorR on any tier) then screenedOut else pass
                 else undecided
    decideT    = applyPolicy(blocked x, decideCore)      # removal (rank >= 1) + blocking tag -> undecided

Every Lean definition used is transcribed one-to-one below (same names in snake_case), over exact
rationals: inputs are converted with ``fractions.Fraction`` (exact IEEE value), constants are
exact (``1/10``, ``1/50``, ``3/4``, ...).

Python-only, on top of Lean (documented in docs/design/20_build_core.md):
* A ``None`` in a field Lean requires (anything but ``val.EF``) cannot form a Lean ``Inputs``;
  the result is ``undecided`` with rule ``MISSING`` (doc 01 §6.1: missing data -> undecided).
  Shape errors of ``mat``/``dos`` (non-square matrix, length mismatch) count as WELL-FORMED
  failures.
* ``core_b(wfgv, first_fires, would_fire, blocked)`` is the Bool skeleton (Lean ``decideB``).
"""
from __future__ import annotations

import math
from dataclasses import dataclass, field
from fractions import Fraction
from typing import Any, Optional, Sequence

from ..records import (D_PASS, D_REJECTED_PROVEN, D_SCREENED_OUT, D_UNDECIDED, DecideInputs, DecideOutcome,
                       InputTag, Prov)

ORIGIN = "rtt.decide.core.decide_t"

HALF = Fraction(1, 2)
TENTH = Fraction(1, 10)
FIFTIETH = Fraction(1, 50)
HUNDREDTH = Fraction(1, 100)
THREE_QUARTERS = Fraction(3, 4)

RULES = ("N1-KAPPA", "N1-ROWMAX", "N1-DOS", "N1-ISO", "NM-CERT", "NM-G")      # Lean allRules order
LEAN_RULE_NAMES = {"N1-KAPPA": "kappa", "N1-ROWMAX": "rowmax", "N1-DOS": "dos", "N1-ISO": "iso",
                   "NM-CERT": "nmCert", "NM-G": "nmG"}
COMMON_ASSUMPTIONS = ("ME-MODEL", "CERT-H", "ENV-SOUND")
RULE_ASSUMPTIONS = {
    "N1-KAPPA": ("M-DISCRETE",),
    "N1-ROWMAX": ("M-DISCRETE",),
    "N1-DOS": ("M-SHEET-ISO", "SHEETS-CORRECT"),
    "N1-ISO": ("M-ISO",),
    "NM-CERT": ("M-ISO",),
    "NM-G": ("M-ISO", "CERT-NM"),
}
NMG_UNCERTIFIED_ASSUMPTIONS = ("M-ISO", "NUM-EINSTEIN-MIN", "NUM-G-MONOTONE")
#: Lean Rule.uses
RULE_USES = {
    "N1-KAPPA": ("mat", "eps"),
    "N1-ROWMAX": ("mat", "eps"),
    "N1-DOS": ("isoA", "eps", "dos"),
    "N1-ISO": ("isoA", "eps"),
    "NM-CERT": ("isoM", "isoLam", "eps"),
    "NM-G": ("isoM", "isoLam", "eps"),
}
CERTIFIED = {r: True for r in RULES}                 # Lean Rule.certified: all [THM]/[THM-C] now


# ------------------------------------------------------------------------------------ config
@dataclass(frozen=True)
class DecideConfig:
    """Lean ``Config``: a0, T*, the §0.3 table ``tbl = [(λ_t, g_t at 300 K)]`` sorted by λ (last
    λ = 100) and ``gInf``."""
    a0: Fraction
    Tstar: Fraction
    tbl: tuple[tuple[Fraction, Fraction], ...]
    gInf: Fraction
    gUp_source: str = ""
    # pending doc 01 §0.6a / Lean additions, switchable (targets.toml [decide]):
    w_lo_sum_check: bool = True       # dosWF also requires sum_i w_i (1 - dW) <= 1
    moment_consistency: bool = True   # wf also requires A^2 <= lam m / 4 (1 + 1/100) (nominal moments)
    nmg_certified: bool = True        # NM-G records [M-ISO, CERT-NM]; False: the older NUM-* steps

    @property
    def floor(self) -> Fraction:
        return self.a0 * self.Tstar

    def gUp(self, lh: Fraction) -> Fraction:
        """Lean ``Config.gUp``: (T*/300)^2 * tableLookup tbl gInf lh."""
        return (self.Tstar / 300) ** 2 * table_lookup(self.tbl, self.gInf, lh)

    @classmethod
    def from_spec(cls, spec) -> "DecideConfig":
        gt = spec.g_table
        tbl = tuple((lk, m0 * 90000) for lk, m0 in zip(gt.lams, gt.m0))
        flags = dict(getattr(spec, "decide_flags", {}) or {})
        return cls(a0=spec.a0, Tstar=spec.T_star, tbl=tbl, gInf=gt.m0_inf * 90000, gUp_source=gt.source, **flags)

    def record(self) -> dict:
        return {"a0": f"{self.a0.numerator}/{self.a0.denominator}", "Tstar": str(self.Tstar),
                "gUp_source": self.gUp_source}


def table_lookup(tbl: Sequence[tuple[Fraction, Fraction]], gInf: Fraction, lh: Fraction) -> Fraction:
    """Lean ``tableLookup``: first (λ_t, g_t) with lh <= λ_t, else gInf."""
    for lt, gt in tbl:
        if lh <= lt:
            return gt
    return gInf


def F(v) -> Fraction:
    """Exact rational of an int/float/str/Fraction (bool and non-finite rejected)."""
    if isinstance(v, bool) or v is None:
        raise TypeError(f"not a number: {v!r}")
    if isinstance(v, Fraction):
        return v
    if isinstance(v, float) and not math.isfinite(v):
        raise ValueError("non-finite")
    return Fraction(v)


# ------------------------------------------------------------------------------------ exact view
class _Missing(Exception):
    pass


class _Malformed(Exception):
    pass


@dataclass
class _X:
    """Exact (Fraction) view of DecideInputs, the Lean ``Inputs ℚ``."""
    EF: Optional[Fraction]
    omegaMax: Fraction
    imagAbsMaxK: Fraction
    containsHLiBe: bool
    anharmonicCheck: bool
    anharmShiftAbs: Fraction
    alpha2FRecomputedAnharmonic: bool
    muStarNegativeSuspected: bool
    admittedSmearings: int
    qMeshes: int
    nSheets: int
    A: Fraction
    lam: Fraction
    m: Fraction
    dA: Fraction
    dLam: Fraction
    dM: Fraction
    provA: Prov
    provLam: Prov
    provM: Prov
    epsSys: Fraction
    epsOmega: Fraction
    qSpread: Fraction
    epsProv: Prov
    mat: Optional[dict]
    dos: Optional[dict]
    tags: list = field(default_factory=list)


def _req(obj, name, kind="num"):
    v = getattr(obj, name)
    if v is None:
        raise _Missing(name)
    if kind == "bool":
        if not isinstance(v, bool):
            raise _Malformed(f"{name} not a bool")
        return v
    if kind == "nat":
        if isinstance(v, bool) or not isinstance(v, int) or v < 0:
            raise _Malformed(f"{name}={v!r} not a natural number")
        return v
    if kind == "prov":
        if not isinstance(v, Prov):
            raise _Malformed(f"{name} not a Prov")
        return v
    try:
        return F(v)
    except (TypeError, ValueError) as e:
        raise _Malformed(f"{name}: {e}")


def exact_view(x: DecideInputs) -> _X:
    v, i, e = x.val, x.iso, x.eps
    if v is None or i is None or e is None:
        raise _Missing("val/iso/eps")
    EF = None
    if v.EF is not None:
        try:
            EF = F(v.EF)
        except (TypeError, ValueError) as err:
            raise _Malformed(f"EF: {err}")
    mat = None
    if x.mat is not None:
        md = x.mat
        B = len(md.calA)
        if B != len(md.n) or B != len(md.v) or any(len(r) != B for r in md.calA):
            raise _Malformed("mat: shape mismatch")
        try:
            mat = {"B": B, "calA": [[F(a) for a in r] for r in md.calA], "n": [F(a) for a in md.n],
                   "v": [F(a) for a in md.v], "dCalA": F(md.dCalA), "prov": md.prov}
        except (TypeError, ValueError) as err:
            raise _Malformed(f"mat: {err}")
        if not isinstance(md.prov, Prov):
            raise _Malformed("mat.prov")
    dos = None
    if x.dos is not None:
        dd = x.dos
        try:
            dos = {"n": len(dd.w), "w": [F(a) for a in dd.w], "dW": F(dd.dW), "prov": dd.prov}
        except (TypeError, ValueError) as err:
            raise _Malformed(f"dos: {err}")
        if not isinstance(dd.prov, Prov):
            raise _Malformed("dos.prov")
    return _X(EF=EF, omegaMax=_req(v, "omegaMax"), imagAbsMaxK=_req(v, "imagAbsMaxK"),
              containsHLiBe=_req(v, "containsHLiBe", "bool"), anharmonicCheck=_req(v, "anharmonicCheck", "bool"),
              anharmShiftAbs=_req(v, "anharmShiftAbs"),
              alpha2FRecomputedAnharmonic=_req(v, "alpha2FRecomputedAnharmonic", "bool"),
              muStarNegativeSuspected=_req(v, "muStarNegativeSuspected", "bool"),
              admittedSmearings=_req(v, "admittedSmearings", "nat"), qMeshes=_req(v, "qMeshes", "nat"),
              nSheets=_req(v, "nSheets", "nat"),
              A=_req(i, "A"), lam=_req(i, "lam"), m=_req(i, "m"), dA=_req(i, "dA"), dLam=_req(i, "dLam"),
              dM=_req(i, "dM"), provA=_req(i, "provA", "prov"), provLam=_req(i, "provLam", "prov"),
              provM=_req(i, "provM", "prov"),
              epsSys=_req(e, "epsSys"), epsOmega=_req(e, "epsOmega"), qSpread=_req(e, "qSpread"),
              epsProv=_req(e, "prov", "prov"), mat=mat, dos=dos, tags=list(x.tags))


# ------------------------------------------------------------------------------------ derived quantities
def fin_max(vals) -> Fraction:
    """Lean finMax: max(0, max_i f i)."""
    return max([Fraction(0)] + list(vals))


def fin_min1(vals) -> Fraction:
    """Lean finMin1: min(1, min_i f i)."""
    return min([Fraction(1)] + list(vals))


def Ahi(x: _X) -> Fraction:
    return x.A * (1 + x.dA) * (1 + x.epsSys) / (1 - x.epsOmega)


def lamHi(x: _X) -> Fraction:
    return x.lam * (1 + x.dLam) * (1 + x.epsSys) / (1 - x.epsOmega) ** 2


def mHi(x: _X) -> Fraction:
    return x.m * (1 + x.dM) * (1 + x.epsSys)


def calAhi(md: dict, x: _X, i: int, j: int) -> Fraction:
    return md["calA"][i][j] * (1 + md["dCalA"]) * (1 + x.epsSys) / (1 - x.epsOmega)


def kappaHi(md: dict, x: _X) -> Fraction:
    B = md["B"]
    return fin_max(sum(calAhi(md, x, i, j) * md["v"][j] for j in range(B)) / md["v"][i] for i in range(B))


def rowAhi(md: dict, x: _X, i: int) -> Fraction:
    return sum((calAhi(md, x, i, j) for j in range(md["B"])), Fraction(0))


def rowMaxHi(md: dict, x: _X) -> Fraction:
    return fin_max(rowAhi(md, x, i) for i in range(md["B"]))


def wminLo(dd: dict) -> Fraction:
    return fin_min1(w * (1 - dd["dW"]) for w in dd["w"])


def mb_present(x: _X) -> bool:
    return x.mat is not None or x.dos is not None


def k_resolved(x: _X) -> bool:
    return x.mat is not None


# ------------------------------------------------------------------------------------ gates
def in0half(d: Fraction) -> bool:
    return 0 <= d and d <= HALF


def mat_wf(md: Optional[dict]) -> bool:
    if md is None:
        return True
    B, A, n, v = md["B"], md["calA"], md["n"], md["v"]
    rng = range(B)
    return (1 <= B and all(0 <= A[i][j] for i in rng for j in rng) and all(0 < n[i] for i in rng)
            and all(0 < v[i] for i in rng) and in0half(md["dCalA"])
            and all(abs(n[i] * A[i][j] - n[j] * A[j][i]) <= HUNDREDTH * max(n[i] * A[i][j], n[j] * A[j][i])
                    for i in rng for j in rng))


def dos_wf(nS: int, dd: Optional[dict], lo_sum_check: bool = True) -> bool:
    if dd is None:
        return True
    w = dd["w"]
    s = sum(w, Fraction(0))
    ok = (dd["n"] == nS and all(0 < wi and wi <= 1 for wi in w)
          and (Fraction(99, 100) <= s and s <= Fraction(101, 100)) and in0half(dd["dW"]))
    if ok and lo_sum_check:
        ok = sum((wi * (1 - dd["dW"]) for wi in w), Fraction(0)) <= 1
    return ok


def wf(x: _X, cfg: Optional["DecideConfig"] = None) -> bool:
    """Lean ``wf`` (WELL-FORMED, §0.6a), plus the switchable pending checks of ``cfg``."""
    lo_sum = True if cfg is None else cfg.w_lo_sum_check
    if (cfg is None or cfg.moment_consistency) and not (x.A ** 2 <= x.lam * x.m / 4 * (1 + HUNDREDTH)):
        return False
    return (0 < x.A and 0 < x.lam and 0 < x.m and 0 < x.omegaMax
            and (x.EF is None or 0 < x.EF)
            and 0 <= x.imagAbsMaxK and 0 <= x.anharmShiftAbs and 1 <= x.nSheets
            and in0half(x.dA) and in0half(x.dLam) and in0half(x.dM) and in0half(x.qSpread)
            and TENTH <= x.epsSys and 0 <= x.epsOmega and x.epsOmega < 1
            and mat_wf(x.mat) and dos_wf(x.nSheets, x.dos, lo_sum))


def gvalid(x: _X) -> bool:
    """Lean ``gvalid`` (G-VALID, revised)."""
    return ((x.EF is not None and 0 < x.EF and lamHi(x) * x.omegaMax <= x.EF / 10)
            and x.imagAbsMaxK <= 30
            and (not x.containsHLiBe or x.anharmonicCheck)
            and (not (TENTH < x.anharmShiftAbs) or x.alpha2FRecomputedAnharmonic)
            and not x.muStarNegativeSuspected
            and 3 <= x.admittedSmearings and 2 <= x.qMeshes
            and TENTH <= x.epsSys
            and max(FIFTIETH, max(x.anharmShiftAbs, x.qSpread)) <= x.epsOmega
            and (not (1 < x.nSheets) or mb_present(x)))


def gvalid_reasons(x: _X) -> list[str]:
    """Which G-VALID conjuncts fail (diagnostics only; ``gvalid`` decides)."""
    r = []
    if x.EF is None:
        r.append("GVALID:1:E_F_missing")
    elif not (0 < x.EF and lamHi(x) * x.omegaMax <= x.EF / 10):
        r.append("GVALID:1:nonadiabatic")
    if not x.imagAbsMaxK <= 30:
        r.append("GVALID:2:imaginary_modes")
    if x.containsHLiBe and not x.anharmonicCheck:
        r.append("GVALID:2:anharmonic_check_missing")
    if TENTH < x.anharmShiftAbs and not x.alpha2FRecomputedAnharmonic:
        r.append("GVALID:2:anharmonic_shift_not_recomputed")
    if x.muStarNegativeSuspected:
        r.append("GVALID:3:mu_star_sign")
    if x.admittedSmearings < 3:
        r.append("GVALID:4:admitted_smearings<3")
    if x.qMeshes < 2:
        r.append("GVALID:4:q_meshes<2")
    if not TENTH <= x.epsSys:
        r.append("GVALID:env:eps_sys<0.10")
    if not max(FIFTIETH, max(x.anharmShiftAbs, x.qSpread)) <= x.epsOmega:
        r.append("GVALID:env:eps_omega_below_minimum")
    if 1 < x.nSheets and not mb_present(x):
        r.append("GVALID:5:multiband_without_data")
    return r


# ------------------------------------------------------------------------------------ rules
def floor_r(cfg: DecideConfig, x: _X, r: str) -> bool:
    """Lean ``floorR``: the floor test of rule r with its own pre-conditions."""
    fl = cfg.floor
    if r == "N1-KAPPA":
        return x.mat is not None and kappaHi(x.mat, x) < fl
    if r == "N1-ROWMAX":
        return x.mat is not None and rowMaxHi(x.mat, x) < fl
    if r == "N1-DOS":
        if x.dos is None:
            return False
        w = wminLo(x.dos)
        return 0 < w and Ahi(x) / w < fl
    if r == "N1-ISO":
        a = Ahi(x)
        return x.nSheets == 1 and a < fl and (a < THREE_QUARTERS * fl or k_resolved(x))
    if r == "NM-CERT":
        return x.nSheets == 1 and mHi(x) < 4 * fl ** 2 / lamHi(x)
    if r == "NM-G":
        g = cfg.gUp(lamHi(x))
        mh = mHi(x)
        return x.nSheets == 1 and mh < g and (mh < THREE_QUARTERS * g or k_resolved(x))
    raise KeyError(r)


def prov_of(x: _X, name: str) -> Optional[Prov]:
    return {"isoA": x.provA, "isoLam": x.provLam, "isoM": x.provM, "eps": x.epsProv,
            "mat": x.mat["prov"] if x.mat else None, "dos": x.dos["prov"] if x.dos else None}[name]


def full_same_run(ps: list[Optional[Prov]]) -> bool:
    """Lean ``fullSameRun``: non-empty, first is fullDFPT, all equal to it (tier and run)."""
    if not ps or ps[0] is None:
        return False
    p = ps[0]
    return p.tier == "fullDFPT" and all(q is not None and q.tier == p.tier and q.run == p.run for q in ps)


def prov_ok(x: _X, r: str) -> bool:
    return full_same_run([prov_of(x, n) for n in RULE_USES[r]])


def fires(cfg: DecideConfig, x: _X, r: str) -> bool:
    return floor_r(cfg, x, r) and prov_ok(x, r)


def first_rule(cfg: DecideConfig, x: _X) -> Optional[str]:
    for r in RULES:
        if fires(cfg, x, r):
            return r
    return None


def would_fire(cfg: DecideConfig, x: _X) -> bool:
    return any(floor_r(cfg, x, r) for r in RULES)


def decide_core(cfg: DecideConfig, x: _X) -> tuple[str, Optional[str]]:
    if wf(x, cfg) and gvalid(x):
        r = first_rule(cfg, x)
        if r is not None:
            return D_REJECTED_PROVEN, r
        return (D_SCREENED_OUT if would_fire(cfg, x) else D_PASS), None
    return D_UNDECIDED, None


# ------------------------------------------------------------------------------------ tags and policy
def blocks_rejection(t: InputTag) -> bool:
    """Lean ``InputTag.blocksRejection``. A kind Lean does not know is treated as ``unknown``."""
    k = t.kind
    if k in ("fileIncomplete", "fileCompleteUnverified", "unknown"):
        return True
    if k == "dropNonpos":
        return True if t.massK is None else 0 < F(t.massK)     # None: fail-closed
    if k == "negClip":
        return False
    if k == "sigmaExcluded":
        return True if t.blocksRejection is None else bool(t.blocksRejection)   # None: fail-closed
    if k == "informational":
        return False
    return True


def blocked(tags: Sequence[InputTag]) -> bool:
    return any(blocks_rejection(t) for t in tags)


RANK = {D_REJECTED_PROVEN: 2, D_SCREENED_OUT: 1, D_UNDECIDED: 0, D_PASS: 0}


def is_removal(o: str) -> bool:
    return 1 <= RANK[o]


def apply_policy(b: bool, o: str) -> str:
    return D_UNDECIDED if (b and is_removal(o)) else o


def core_b(wfgv: bool, first: bool, would: bool, b: bool) -> str:
    """Bool skeleton of decideT (Lean ``decideB``): wf&&gvalid, firstRule.isSome, wouldFire,
    blocked. (``first`` implies ``would``; the row first & !would is unreachable.)"""
    core = (D_REJECTED_PROVEN if first else (D_SCREENED_OUT if would else D_PASS)) if wfgv else D_UNDECIDED
    return apply_policy(b, core)


# ------------------------------------------------------------------------------------ decide_t
def _inputs_used(x: DecideInputs, r: str) -> list[dict]:
    out = []
    for n in RULE_USES[r]:
        if n == "isoA":
            out.append({"name": "isoA", "A": x.iso.A, "dA": x.iso.dA, "prov": x.iso.provA.to_json()})
        elif n == "isoLam":
            out.append({"name": "isoLam", "lam": x.iso.lam, "dLam": x.iso.dLam, "prov": x.iso.provLam.to_json()})
        elif n == "isoM":
            out.append({"name": "isoM", "m": x.iso.m, "dM": x.iso.dM, "prov": x.iso.provM.to_json()})
        elif n == "eps":
            out.append({"name": "eps", "epsSys": x.eps.epsSys, "epsOmega": x.eps.epsOmega,
                        "qSpread": x.eps.qSpread, "prov": x.eps.prov.to_json()})
        elif n == "mat":
            out.append({"name": "mat", "B": x.mat.B, "dCalA": x.mat.dCalA, "prov": x.mat.prov.to_json()})
        elif n == "dos":
            out.append({"name": "dos", "w": x.dos.w, "dW": x.dos.dW, "prov": x.dos.prov.to_json()})
    return out


def decide_t(cfg: DecideConfig, x: DecideInputs) -> DecideOutcome:
    b = blocked(x.tags)
    try:
        X = exact_view(x)
    except _Missing as e:
        return DecideOutcome(outcome=D_UNDECIDED, rule="MISSING", reasons=[f"MISSING:{e}"],
                             core_outcome=D_UNDECIDED, blocked=b)
    except _Malformed as e:
        return DecideOutcome(outcome=D_UNDECIDED, rule="WELL-FORMED", reasons=[f"WF:{e}"],
                             core_outcome=D_UNDECIDED, blocked=b)
    if not wf(X, cfg):
        return DecideOutcome(outcome=D_UNDECIDED, rule="WELL-FORMED", reasons=["WF:domain"],
                             core_outcome=D_UNDECIDED, blocked=b)
    if not gvalid(X):
        return DecideOutcome(outcome=D_UNDECIDED, rule="G-VALID", reasons=gvalid_reasons(X) or ["GVALID"],
                             core_outcome=D_UNDECIDED, blocked=b)
    core, r = decide_core(cfg, X)
    out = apply_policy(b, core)
    wouldf = [q for q in RULES if floor_r(cfg, X, q)]
    reasons = []
    if core == D_SCREENED_OUT:
        reasons.append("WOULD_FIRE_LOWER_TIER:" + ",".join(wouldf))
    if out != core:
        reasons.append("BLOCKED_REJ:" + ",".join(t.name for t in x.tags if blocks_rejection(t)))
    if out != D_REJECTED_PROVEN and r is not None:      # Lean Outcome.assumptions/inputsUsed are [] then
        reasons.append(f"BLOCKED_RULE:{r}")
        r = None
    assumptions = list(COMMON_ASSUMPTIONS) + list(rule_assumptions(cfg, r)) if r else []
    return DecideOutcome(outcome=out, rule=r, assumptions=assumptions,
                         inputs_used=_inputs_used(x, r) if r else [], reasons=reasons,
                         would_fire=wouldf, core_outcome=core, blocked=b)


#: Verbatim lines of lean/TargetLogic/Decide.lean that this module transcribes (drift-checked by
#: tests/test_decide_lean_table.py; if Decide.lean changes, re-check this mirror).
LEAN_TEXTS = (
    """def decideCore (cfg : Config K) (x : Inputs K) : Outcome :=
  if wf x && gvalid x then
    match firstRule cfg x with
    | some r => .rejectedProven r
    | none => if wouldFire cfg x then .screenedOut else .pass
  else .undecided""",
    """def applyPolicy (b : Bool) (o : Outcome) : Outcome :=
  if b && o.isRemoval then .undecided else o""",
    "def decideT (cfg : Config K) (x : Inputs K) : Outcome := applyPolicy (blocked x) (decideCore cfg x)",
    "def allRules : List Rule := [.kappa, .rowmax, .dos, .iso, .nmCert, .nmG]",
    """def wf (x : Inputs K) : Bool :=
  decide (0 < x.iso.A) && decide (0 < x.iso.lam) && decide (0 < x.iso.m) &&
  decide (0 < x.val.omegaMax) &&
  (match x.val.EF with | some e => decide (0 < e) | none => true) &&
  decide (0 ≤ x.val.imagAbsMaxK) && decide (0 ≤ x.val.anharmShiftAbs) &&
  decide (1 ≤ x.val.nSheets) &&
  in0half x.iso.dA && in0half x.iso.dLam && in0half x.iso.dM && in0half x.eps.qSpread &&
  decide (1 / 10 ≤ x.eps.epsSys) && decide (0 ≤ x.eps.epsOmega) && decide (x.eps.epsOmega < 1) &&
  decide (x.iso.A ^ 2 ≤ x.iso.lam * x.iso.m / 4 * (1 + 1 / 100)) &&
  matWF x.mat && dosWF x.val.nSheets x.dos""",
    """    decide (99 / 100 ≤ ∑ i, dd.w i ∧ ∑ i, dd.w i ≤ 101 / 100) && in0half dd.dW &&
    decide (∑ i, dd.w i * (1 - dd.dW) ≤ 1)""",
    "  | .nmG => [.mIso, .certNM]",
    """def gvalid (x : Inputs K) : Bool :=
  (match x.val.EF with
    | some e => decide (0 < e) && decide (lamHi x * x.val.omegaMax ≤ e / 10)
    | none => false) &&
  decide (x.val.imagAbsMaxK ≤ 30) &&
  (!x.val.containsHLiBe || x.val.anharmonicCheck) &&
  (!decide (1 / 10 < x.val.anharmShiftAbs) || x.val.alpha2FRecomputedAnharmonic) &&
  !x.val.muStarNegativeSuspected &&
  decide (3 ≤ x.val.admittedSmearings) && decide (2 ≤ x.val.qMeshes) &&
  decide (1 / 10 ≤ x.eps.epsSys) &&
  decide (max (1 / 50) (max x.val.anharmShiftAbs x.eps.qSpread) ≤ x.eps.epsOmega) &&
  (!decide (1 < x.val.nSheets) || mbPresent x)""",
    """  | .iso => decide (x.val.nSheets = 1) && decide (Ahi x < cfg.floor) &&
      (decide (Ahi x < 3 / 4 * cfg.floor) || kResolved x)
  | .nmCert => decide (x.val.nSheets = 1) && decide (mHi x < 4 * cfg.floor ^ 2 / lamHi x)
  | .nmG => decide (x.val.nSheets = 1) && decide (mHi x < cfg.gUp (lamHi x)) &&
      (decide (mHi x < 3 / 4 * cfg.gUp (lamHi x)) || kResolved x)""",
    """def fullSameRun : List (Option Prov) → Bool
  | [] => false
  | none :: _ => false
  | some p :: rest => (p.tier == .fullDFPT) && rest.all (fun q => q == some p)""",
    """def InputTag.blocksRejection : InputTag K → Bool
  | .fileIncomplete => true
  | .fileCompleteUnverified => true
  | .dropNonpos m => decide (0 < m)
  | .negClip _ => false
  | .sigmaExcluded r _ => r
  | .informational _ => false
  | .unknown _ => true""",
)


def rule_assumptions(cfg: DecideConfig, r: str) -> tuple[str, ...]:
    if r == "NM-G" and not cfg.nmg_certified:
        return NMG_UNCERTIFIED_ASSUMPTIONS
    return RULE_ASSUMPTIONS[r]
