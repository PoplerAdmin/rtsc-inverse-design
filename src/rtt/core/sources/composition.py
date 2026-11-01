"""Chemical formulas, compositions and a small element table (no external dependencies).

Used by sources (formula parsing), descriptors (masses, light-element fractions) and validity
(magnetic / d- and f-electron elements).  Standard atomic weights (IUPAC, rounded to 3-4
significant decimals; precision is irrelevant for the heuristic uses here).
"""
from __future__ import annotations

import math
import re
from collections import OrderedDict

# Z-ordered symbols and standard atomic weights (u) for Z = 1 .. 94.
_ELEMENTS = (
    ("H", 1.008), ("He", 4.0026), ("Li", 6.94), ("Be", 9.0122), ("B", 10.81), ("C", 12.011),
    ("N", 14.007), ("O", 15.999), ("F", 18.998), ("Ne", 20.180), ("Na", 22.990), ("Mg", 24.305),
    ("Al", 26.982), ("Si", 28.085), ("P", 30.974), ("S", 32.06), ("Cl", 35.45), ("Ar", 39.948),
    ("K", 39.098), ("Ca", 40.078), ("Sc", 44.956), ("Ti", 47.867), ("V", 50.942), ("Cr", 51.996),
    ("Mn", 54.938), ("Fe", 55.845), ("Co", 58.933), ("Ni", 58.693), ("Cu", 63.546), ("Zn", 65.38),
    ("Ga", 69.723), ("Ge", 72.630), ("As", 74.922), ("Se", 78.971), ("Br", 79.904), ("Kr", 83.798),
    ("Rb", 85.468), ("Sr", 87.62), ("Y", 88.906), ("Zr", 91.224), ("Nb", 92.906), ("Mo", 95.95),
    ("Tc", 98.0), ("Ru", 101.07), ("Rh", 102.91), ("Pd", 106.42), ("Ag", 107.87), ("Cd", 112.41),
    ("In", 114.82), ("Sn", 118.71), ("Sb", 121.76), ("Te", 127.60), ("I", 126.90), ("Xe", 131.29),
    ("Cs", 132.91), ("Ba", 137.33), ("La", 138.91), ("Ce", 140.12), ("Pr", 140.91), ("Nd", 144.24),
    ("Pm", 145.0), ("Sm", 150.36), ("Eu", 151.96), ("Gd", 157.25), ("Tb", 158.93), ("Dy", 162.50),
    ("Ho", 164.93), ("Er", 167.26), ("Tm", 168.93), ("Yb", 173.05), ("Lu", 174.97), ("Hf", 178.49),
    ("Ta", 180.95), ("W", 183.84), ("Re", 186.21), ("Os", 190.23), ("Ir", 192.22), ("Pt", 195.08),
    ("Au", 196.97), ("Hg", 200.59), ("Tl", 204.38), ("Pb", 207.2), ("Bi", 208.98), ("Po", 209.0),
    ("At", 210.0), ("Rn", 222.0), ("Fr", 223.0), ("Ra", 226.0), ("Ac", 227.0), ("Th", 232.04),
    ("Pa", 231.04), ("U", 238.03), ("Np", 237.0), ("Pu", 244.0),
)
ATOMIC_NUMBER: dict[str, int] = {s: i + 1 for i, (s, _) in enumerate(_ELEMENTS)}
ATOMIC_MASS: dict[str, float] = {s: m for s, m in _ELEMENTS}

#: elements with partially filled 3d/4f shells that commonly carry local moments
MAGNETIC_ELEMENTS = frozenset({"V", "Cr", "Mn", "Fe", "Co", "Ni", "Cu",
                               "Ce", "Pr", "Nd", "Pm", "Sm", "Eu", "Gd", "Tb", "Dy", "Ho", "Er",
                               "Tm", "Yb", "U", "Np", "Pu"})
TRANSITION_METALS = frozenset(s for s, z in ATOMIC_NUMBER.items()
                              if 21 <= z <= 30 or 39 <= z <= 48 or 72 <= z <= 80)
LANTHANIDES_ACTINIDES = frozenset(s for s, z in ATOMIC_NUMBER.items()
                                  if 57 <= z <= 71 or 89 <= z <= 94)

_TOKEN = re.compile(r"([A-Z][a-z]?)|(\d*\.\d+|\d+)|([\(\)\[\]\{\}])|(\s+|·|\*)")


def is_element(symbol: str) -> bool:
    return symbol in ATOMIC_NUMBER


def parse_formula(formula: str) -> dict[str, float]:
    """Parse ``'LaH10'``, ``'Ca(BH4)2'``, ``'H3 S1'``, ``'Mg0.5B'`` into {element: amount}.

    Unknown element-like symbols are kept (so synthetic placeholders survive), callers decide.
    """
    stack: list[dict[str, float]] = [OrderedDict()]
    pos = 0
    last: dict[str, float] | str | None = None
    s = formula.strip()
    while pos < len(s):
        m = _TOKEN.match(s, pos)
        if not m:
            raise ValueError(f"cannot parse formula {formula!r} at {s[pos:]!r}")
        pos = m.end()
        sym, num, br, _ws = m.groups()
        if sym:
            stack[-1][sym] = stack[-1].get(sym, 0.0) + 1.0
            last = sym
        elif num:
            k = float(num)
            if isinstance(last, str):
                stack[-1][last] += k - 1.0
            elif isinstance(last, dict):
                for e, v in last.items():
                    stack[-1][e] = stack[-1].get(e, 0.0) + v * (k - 1.0)
            else:
                raise ValueError(f"number without element in {formula!r}")
            last = None
        elif br in "([{":
            stack.append(OrderedDict())
            last = None
        elif br in ")]}":
            if len(stack) < 2:
                raise ValueError(f"unbalanced brackets in {formula!r}")
            grp = stack.pop()
            for e, v in grp.items():
                stack[-1][e] = stack[-1].get(e, 0.0) + v
            last = grp
        else:
            last = None
    if len(stack) != 1:
        raise ValueError(f"unbalanced brackets in {formula!r}")
    out = {e: v for e, v in stack[0].items() if v != 0}
    if not out:
        raise ValueError(f"empty formula {formula!r}")
    return dict(out)


def _gcd_float(values: list[float]) -> float:
    if all(abs(v - round(v)) < 1e-8 for v in values):
        g = 0
        for v in values:
            g = math.gcd(g, int(round(v)))
        return float(g or 1)
    return 1.0


def reduced_composition(comp: dict[str, float]) -> dict[str, float]:
    vals = [v for v in comp.values() if v > 0]
    g = _gcd_float(vals) if vals else 1.0
    return {e: v / g for e, v in comp.items() if v > 0}


def _num(v: float) -> str:
    if abs(v - 1.0) < 1e-9:
        return ""
    if abs(v - round(v)) < 1e-9:
        return str(int(round(v)))
    return f"{v:.4g}"


def formula_string(comp: dict[str, float], reduce: bool = True) -> str:
    """Deterministic formula string: non-H elements by atomic number, then H ('LaH10', 'SH3')."""
    c = reduced_composition(comp) if reduce else comp
    keys = sorted((k for k in c if k != "H"), key=lambda k: (ATOMIC_NUMBER.get(k, 999), k))
    if "H" in c:
        keys.append("H")
    return "".join(f"{k}{_num(c[k])}" for k in keys)


def atomic_fractions(comp: dict[str, float]) -> dict[str, float]:
    tot = sum(v for v in comp.values() if v > 0)
    if tot <= 0:
        return {}
    return {e: v / tot for e, v in comp.items() if v > 0}


def chemsys(comp: dict[str, float]) -> str:
    return "-".join(sorted(e for e, v in comp.items() if v > 0))


def n_atoms(comp: dict[str, float]) -> float:
    return sum(v for v in comp.values() if v > 0)
