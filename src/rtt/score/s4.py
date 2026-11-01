"""S4 ranking score (architecture §5.1 ``score_s4``; doc 11 T1c and T1a).

**Ranking only (contract R2): this number orders candidates for S5; it never removes one.**

Components, each mapped to [0, 1]:

* ``s_phi``  Belli Tc_phi / T* (clipped to [0, 1]), hydrogen compounds only (the fit of Belli
  et al. 2021 was made on hydrogen-based compounds; doc 11 §3.3). Missing -> component absent.
* ``s_dos``  N(E_F) per eV per atom (both spins) / (N(E_F) + n_half), n_half = 0.3 /eV/atom
  (a saturating map; the scale is a convention, not a fit).
* ``s_hdos`` H share of N(E_F): 1 inside the doc 02 S3 window [0.4, 0.6]
  (targets.toml ``soft.S3_HDOS``), falling linearly to 0 at 0 and 1.

Weights: H compounds 0.6 s_phi + 0.25 s_dos + 0.15 s_hdos; compounds without H s_dos alone.
Missing components are dropped and the remaining weights renormalised. A candidate with no
usable component scores 0.0. The weights are conventions (doc 11 gives T1c ordering power
only for hydrides) and should be re-fitted against S5 A_coarse once the A4->5 audit has data.
"""
from __future__ import annotations

from typing import Any

N_HALF = 0.3
W_H = {"s_phi": 0.6, "s_dos": 0.25, "s_hdos": 0.15}
W_NOH = {"s_dos": 1.0}


def _get(x: Any, name: str) -> Any:
    if x is None:
        return None
    if isinstance(x, dict):
        return x.get(name)
    return getattr(x, name, None)


def _hdos_score(h: float, window: tuple[float, float]) -> float:
    lo, hi = window
    if lo <= h <= hi:
        return 1.0
    if h < lo:
        return max(0.0, h / lo) if lo > 0 else 0.0
    return max(0.0, (1.0 - h) / (1.0 - hi)) if hi < 1 else 0.0


def components(electronic: Any, phi: Any = None, *, T_star_K: float = 300.0,
               hdos_window: tuple[float, float] = (0.4, 0.6)) -> dict[str, float]:
    """Score components from an ElectronicResult (or its metrics dict) and optional Belli
    result (overrides the electronic result's phi fields)."""
    out: dict[str, float] = {}
    h_f = _get(phi, "H_f") if phi is not None else _get(electronic, "H_f")
    tc = _get(phi, "Tc_phi_K") if phi is not None else _get(electronic, "Tc_phi_K")
    hd = _get(phi, "H_DOS") if phi is not None else _get(electronic, "H_DOS")
    has_h = bool(h_f and h_f > 0)
    if has_h and tc is not None:
        out["s_phi"] = min(max(float(tc) / float(T_star_K), 0.0), 1.0)
    n = _get(electronic, "N_EF_per_eV_atom")
    if n is not None and n >= 0:
        out["s_dos"] = float(n) / (float(n) + N_HALF)
    if has_h and hd is not None:
        out["s_hdos"] = _hdos_score(float(hd), hdos_window)
    return out


def score_s4(electronic: Any, phi: Any = None, *, T_star_K: float = 300.0,
             hdos_window: tuple[float, float] = (0.4, 0.6)) -> float:
    c = components(electronic, phi, T_star_K=T_star_K, hdos_window=hdos_window)
    h_f = _get(phi, "H_f") if phi is not None else _get(electronic, "H_f")
    w = W_H if (h_f and h_f > 0) else W_NOH
    used = {k: v for k, v in w.items() if k in c}
    tot = sum(used.values())
    if tot <= 0:
        return 0.0
    return float(sum(c[k] * wk for k, wk in used.items()) / tot)
