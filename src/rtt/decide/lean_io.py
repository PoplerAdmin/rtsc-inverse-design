"""Read the Lean export ``lean/export/decide_cases.json`` (``lean/ExportCases.lean``: decideT at
K = Q, exact rationals as strings) into rtt.decide types, so the mirror is tested case by case."""
from __future__ import annotations

from fractions import Fraction
from typing import Any

from ..records import DecideInputs, DosData, EpsData, InputTag, IsoData, MatData, Prov, Validity
from .core import DecideConfig

#: Lean export input-group names -> rtt.decide.core input names
EXPORT_INPUT_NAMES = {"A": "isoA", "lambda": "isoLam", "m": "isoM", "eps": "eps", "calA": "mat", "mat": "mat",
                      "w": "dos", "dos": "dos"}


def _q(s: Any):
    if s is None:
        return None
    return Fraction(s) if isinstance(s, str) else s


def _prov(d: dict | None) -> Prov | None:
    if d is None:
        return None
    tier = {"ML": "ml"}.get(d["tier"], d["tier"])
    return Prov(tier, d["run"])


def config_from_export(c: dict, **flags) -> DecideConfig:
    return DecideConfig(a0=Fraction(c["a0"]), Tstar=Fraction(c["Tc_target"]),
                        tbl=tuple((Fraction(l), Fraction(g)) for l, g in c["g_table"]), gInf=Fraction(c["g_inf"]),
                        gUp_source="lean/export/decide_cases.json", **flags)


def inputs_from_export(d: dict) -> DecideInputs:
    v, i, e = d["validity"], d["iso"], d["eps"]
    val = Validity(EF=_q(v["E_F"]), omegaMax=_q(v["omega_max"]), imagAbsMaxK=_q(v["imag_abs_max_K"]),
                   containsHLiBe=v["contains_H_Li_Be"], anharmonicCheck=v["anharmonic_check"],
                   anharmShiftAbs=_q(v["anharm_shift_abs"]),
                   alpha2FRecomputedAnharmonic=v["alpha2F_recomputed_anharmonic"],
                   muStarNegativeSuspected=v["mu_star_negative_suspected"], admittedSmearings=v["admitted_smearings"],
                   qMeshes=v["q_meshes"], nSheets=v["n_sheets"])
    iso = IsoData(A=_q(i["A_nom"]), lam=_q(i["lambda_nom"]), m=_q(i["m_nom"]), dA=_q(i["d_A"]),
                  dLam=_q(i["d_lambda"]), dM=_q(i["d_m"]), provA=_prov(i["prov_A"]), provLam=_prov(i["prov_lambda"]),
                  provM=_prov(i["prov_m"]))
    eps = EpsData(epsSys=_q(e["eps_sys"]), epsOmega=_q(e["eps_omega"]), qSpread=_q(e["qspread_omega_log"]),
                  prov=_prov(e["prov"]))
    mat = dos = None
    if d.get("mat") is not None:
        m = d["mat"]
        mat = MatData(calA=[[_q(a) for a in r] for r in m["calA"]], n=[_q(a) for a in m["n"]],
                      v=[_q(a) for a in m["v"]], dCalA=_q(m["d_calA"]), prov=_prov(m["prov"]))
    if d.get("dos") is not None:
        s = d["dos"]
        if s["n"] != len(s["w"]):
            raise ValueError("dos.n != len(w)")
        dos = DosData(w=[_q(a) for a in s["w"]], dW=_q(s["d_w"]), prov=_prov(s["prov"]))
    tags = [InputTag(name=t.get("name", t["kind"]), kind=t["kind"], massK=_q(t.get("massK")),
                     blocksRejection=t.get("blocksRejection"), blocksCertification=t.get("blocksCertification"))
            for t in d.get("tags", [])]
    tier = "fullDFPT"
    return DecideInputs(tier=tier, val=val, iso=iso, eps=eps, mat=mat, dos=dos, tags=tags)
