"""Composition-based metallicity prior for the S3 ranking score (doc 31 §2). Ranking only.

Contract (R2; doc 11): this prior only *reorders* candidates in the S3 rank queue. It never
removes one: metallicity is decided at S4 by DFT (H1, tetrahedron N(E_F)). A wrong prior can
delay a metal, never exclude it.

Why a composition rule and not an ML band-gap model: no band-gap model is available offline
(doc 10: the MatterGen ``dft_band_gap`` checkpoint is not local; ALIGNN is not installed and has
weight 0 for H/light elements, doc 11 §2.3). Campaign 1 (doc 30 §3, §7) showed the failure this
prior targets: charge-balanced saline hydrides (AlH3, SrH2) and closed-shell transition-metal
hydride complexes (SrPtH6, BaPtH6, NaAuH2, Ba2Pt2H8) ranked high and were then excluded by H1.

The rule is formal electron counting (textbook ionic / 18-electron counting; no fitted
parameters). Hydrogen is counted as H-. Donors have fixed oxidation states; for a hydride with
transition metals (TM), the TM formal oxidation state is ``ox = (n_H - sum donor charge) / n_TM``
and its d-count is ``d = group - ox``. A closed shell -- d0, d10, an 18-electron count
``d + 2 * (H per TM) == 18`` or a square-planar d8 with 4 H per TM (16 e) -- gives the
"insulator-likely" prior; open shells, fractional or negative TM oxidation states and excess
electrons give "metal-likely". Everything else (mixed anions, p-block frameworks such as MgB2,
H3S) is neutral 0.5. The three levels (0.15 / 0.5 / 0.85; 0.8 for all-metal compositions) are
ordinal choices, **not calibrated probabilities** (doc 31 §2 gives the campaign-1 check).
"""
from __future__ import annotations

import re
from fractions import Fraction
from typing import Any, Mapping

PRIOR_VERSION = "metallicity-prior/2 (parity + electron counting)"
P_METAL, P_INSULATOR, P_NEUTRAL, P_ALLMETAL = 0.85, 0.15, 0.5, 0.80

#: fixed formal oxidation states of electropositive donors (alkali, alkaline earth, group 3,
#: lanthanides, Al/Ga/In, group 12)
DONOR_OX: dict[str, int] = {
    **{e: 1 for e in ("Li", "Na", "K", "Rb", "Cs")},
    **{e: 2 for e in ("Be", "Mg", "Ca", "Sr", "Ba", "Zn", "Cd", "Hg")},
    **{e: 3 for e in ("Sc", "Y", "La", "Ce", "Pr", "Nd", "Pm", "Sm", "Eu", "Gd", "Tb", "Dy", "Ho",
                      "Er", "Tm", "Yb", "Lu", "Al", "Ga", "In")},
}
#: transition metals with variable valence: group number (= valence electrons in the d^n s^m count)
TM_GROUP: dict[str, int] = {
    "Ti": 4, "V": 5, "Cr": 6, "Mn": 7, "Fe": 8, "Co": 9, "Ni": 10, "Cu": 11,
    "Zr": 4, "Nb": 5, "Mo": 6, "Tc": 7, "Ru": 8, "Rh": 9, "Pd": 10, "Ag": 11,
    "Hf": 4, "Ta": 5, "W": 6, "Re": 7, "Os": 8, "Ir": 9, "Pt": 10, "Au": 11,
}
#: p-block metals: metallic as elements, no fixed hydride valence in this rule
P_METALS = frozenset({"Tl", "Sn", "Pb", "Bi", "Ge", "Sb"})
METALS = frozenset(DONOR_OX) | frozenset(TM_GROUP) | P_METALS

_TOKEN = re.compile(r"([A-Z][a-z]?)(\d*\.?\d*)")


def parse_formula(formula: str) -> dict[str, float]:
    """'Mg2IrH6' -> {'Mg': 2, 'Ir': 1, 'H': 6} (no brackets; reduced formulas only)."""
    out: dict[str, float] = {}
    for el, n in _TOKEN.findall(formula or ""):
        out[el] = out.get(el, 0.0) + (float(n) if n else 1.0)
    return out


def composition_of(features_or_metrics: Mapping[str, Any]) -> tuple[dict[str, float] | None, bool]:
    """(element counts, is_cell) from S3 features: ``composition`` (counts of the cell; parity
    usable) or ``formula`` (reduced formula; parity not usable), or a metrics dict."""
    f = features_or_metrics or {}
    if "features" in f and isinstance(f["features"], Mapping):
        f = f["features"]
    comp = f.get("composition")
    if isinstance(comp, Mapping) and comp:
        return {str(k): float(v) for k, v in comp.items()}, True
    if f.get("formula"):
        c = parse_formula(str(f["formula"]))
        return (c or None), False
    return None, False


def _res(p: float, rule: str, **detail: Any) -> dict[str, Any]:
    return {"p_metal": p, "rule": rule, "version": PRIOR_VERSION, **detail}


#: valence electrons counted for parity: s+d for groups 1-12 (= group number), group - 10 for the
#: p block, 3 for the lanthanides (4f treated as core, as in the SG15 set's La; others UNVERIFIED)
_PBLOCK = {"B": 3, "Al": 3, "Ga": 3, "In": 3, "Tl": 3, "C": 4, "Si": 4, "Ge": 4, "Sn": 4, "Pb": 4,
           "N": 5, "P": 5, "As": 5, "Sb": 5, "Bi": 5, "O": 6, "S": 6, "Se": 6, "Te": 6,
           "F": 7, "Cl": 7, "Br": 7, "I": 7}
VALENCE_E: dict[str, int] = {
    "H": 1, **{e: 1 for e in ("Li", "Na", "K", "Rb", "Cs")}, **{e: 2 for e in ("Be", "Mg", "Ca", "Sr", "Ba")},
    **{e: 3 for e in ("Sc", "Y", "La", "Ce", "Pr", "Nd", "Pm", "Sm", "Eu", "Gd", "Tb", "Dy", "Ho", "Er",
                      "Tm", "Yb", "Lu")},
    **TM_GROUP, "Zn": 12, "Cd": 12, "Hg": 12, **_PBLOCK}
P_ODD = 0.95


def valence_electrons(comp: Mapping[str, float]) -> int | None:
    """Valence-electron count of the cell (``comp`` = element counts of the cell itself: parity
    depends on the cell, so pass the primitive cell, not a reduced formula). None if an element
    is not tabulated or a count is fractional."""
    tot = 0
    for el, n in comp.items():
        if el not in VALENCE_E or float(n) != int(float(n)):
            return None
        tot += VALENCE_E[el] * int(float(n))
    return tot


def metallicity_prior(comp: Mapping[str, float] | None, *, cell: bool = True) -> dict[str, Any]:
    """{p_metal, rule, version, ...detail} for an element-count mapping. ``p_metal`` is an ordinal
    prior in [0, 1]; 0.5 = no rule applies.

    Valence-electron parity (R30 F1(c)/F8, doc 31 §2) comes first when ``cell`` is True (the
    counts are those of the primitive cell): an ODD count per cell cannot fill spin-degenerate
    bands, so a non-magnetic band structure is necessarily metallic (p = 0.95). An even count
    carries no such guarantee; the electron-counting rules below then apply. In campaign 1 parity
    separated every S4 H1 outcome of the G1 cells (R30, meas.)."""
    if not comp:
        return _res(P_NEUTRAL, "no_composition")
    ne = valence_electrons(comp) if cell else None
    if ne is not None and ne % 2 == 1:
        return _res(P_ODD, "odd_valence_electron_count", n_valence_e=ne)
    res = _counting_rule(comp)
    if ne is not None:
        res["n_valence_e"] = ne
    return res


def _counting_rule(comp: Mapping[str, float]) -> dict[str, Any]:
    c = {k: Fraction(str(v)) for k, v in comp.items() if float(v) > 0}
    els = set(c)
    nonmetals = els - METALS - {"H"}
    if "H" not in els:
        if not nonmetals:
            return _res(P_ALLMETAL, "all_metal_composition")
        return _res(P_NEUTRAL, "no_rule_nonhydride", nonmetals=sorted(nonmetals))
    if nonmetals:
        return _res(P_NEUTRAL, "no_rule_mixed_anion", nonmetals=sorted(nonmetals))
    if els & P_METALS:
        return _res(P_NEUTRAL, "no_rule_p_metal_hydride")
    nH = c["H"]
    donor_q = sum((c[e] * DONOR_OX[e] for e in els if e in DONOR_OX), Fraction(0))
    tms = [e for e in els if e in TM_GROUP]
    if not tms:
        e_excess = donor_q - nH
        if not (els - {"H"}):
            return _res(P_NEUTRAL, "no_rule_pure_hydrogen")
        if e_excess == 0:
            return _res(P_INSULATOR, "charge_balanced_saline_hydride", excess_e=0.0)
        if e_excess > 0:
            return _res(P_METAL, "excess_electrons", excess_e=float(e_excess))
        return _res(0.40, "hydrogen_beyond_hydride_count", excess_e=float(e_excess))
    n_tm = sum(c[e] for e in tms)
    group = sum(c[e] * TM_GROUP[e] for e in tms) / n_tm
    ox = (nH - donor_q) / n_tm
    L = nH / n_tm
    detail = {"tm_ox": float(ox), "H_per_tm": float(L)}
    if ox < 0:
        return _res(0.80, "negative_tm_oxidation_state", **detail)
    if ox.denominator != 1 or group.denominator != 1:
        return _res(P_METAL, "fractional_tm_oxidation_state", **detail)
    d = group - ox
    detail["d_count"] = float(d)
    if d < 0:
        return _res(P_NEUTRAL, "no_rule_h_exceeds_tm_valence", **detail)
    if d == 0:
        return _res(P_INSULATOR, "closed_shell_d0", **detail)
    if d == 10:
        return _res(P_INSULATOR, "closed_shell_d10", **detail)
    if d + 2 * L == 18:
        return _res(P_INSULATOR, "closed_shell_18_electron", **detail)
    if d == 8 and L == 4:
        return _res(P_INSULATOR, "closed_shell_square_planar_d8", **detail)
    return _res(P_METAL, "open_shell_tm", **detail)


def prior_from_metrics(metrics_or_features: Mapping[str, Any]) -> dict[str, Any]:
    comp, is_cell = composition_of(metrics_or_features)
    return metallicity_prior(comp, cell=is_cell)


__all__ = ["metallicity_prior", "prior_from_metrics", "parse_formula", "composition_of", "valence_electrons",
           "PRIOR_VERSION"]
