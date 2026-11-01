"""S3(c): cheap geometric features (architecture §3.1 S3; doc 02 S2, S4, S8). Ranking only.

* ``H_f``: H atom fraction (S2).
* ``vol_per_H_A3``: cell volume / number of H (S8; lower scores better, no threshold).
* ``min_HH_A``, ``min_MH_A``: shortest H-H and (non-H)-H distances (S4).
* ``motif``: S4 bonding class from H-H geometry only:
  ``H2_molecular`` (min H-H < 0.90 Å; penalised), ``H_network`` (1.0-1.55 Å, weak-covalent),
  ``H_short_other`` (0.90-1.0 Å), ``isolated_H`` (> 1.55 Å: e.g. [MH_n] complexes or H⁻),
  ``no_H``. ELF-based classification needs DFT (S4 φ).
* ``n_H_neighbours_M``: mean number of H within 2.0 Å of each non-H atom (complex size hint).
* ``contains_HLiBe``, ``light_fraction`` (H, Li, Be, B, C, N).
"""
from __future__ import annotations

from typing import Any

import numpy as np

LIGHT = ("H", "Li", "Be", "B", "C", "N")


def _min_pair(st, a_sel, b_sel, rmax: float = 6.0) -> float | None:
    a_idx = [i for i, s in enumerate(st) if a_sel(s.specie.symbol)]
    if not a_idx:
        return None
    best = None
    for i in a_idx:
        for nb in st.get_neighbors(st[i], rmax):
            if b_sel(nb.specie.symbol) and nb.nn_distance > 1e-6:
                if best is None or nb.nn_distance < best:
                    best = nb.nn_distance
    return None if best is None else float(best)


def features(struct: Any) -> dict:
    from ..structure.canonical import to_pymatgen
    d = struct.to_json() if hasattr(struct, "to_json") else struct
    st = to_pymatgen(d.get("structure", d) if isinstance(d, dict) and "lattice" not in d else d)
    sp = [s.specie.symbol for s in st]
    n = len(sp)
    nH = sp.count("H")
    out: dict[str, Any] = {
        "n_atoms": n, "n_H": nH, "H_f": nH / n,
        "vol_per_atom_A3": st.volume / n,
        "vol_per_H_A3": (st.volume / nH) if nH else None,
        "contains_HLiBe": any(x in ("H", "Li", "Be") for x in sp),
        "light_fraction": sum(1 for x in sp if x in LIGHT) / n,
        "formula": st.composition.reduced_formula,
        "composition": {el: sp.count(el) for el in sorted(set(sp))},
    }
    from ..score.metallicity import metallicity_prior      # doc 31 §2: ranking only, recorded
    out["metallicity_prior"] = metallicity_prior(out["composition"])
    if nH:
        hh = _min_pair(st, lambda x: x == "H", lambda x: x == "H")
        mh = _min_pair(st, lambda x: x != "H", lambda x: x == "H")
        out["min_HH_A"] = hh
        out["min_MH_A"] = mh
        if hh is None or hh > 1.55:
            motif = "isolated_H"
        elif hh < 0.90:
            motif = "H2_molecular"
        elif hh < 1.0:
            motif = "H_short_other"
        else:
            motif = "H_network"
        out["motif"] = motif
        m_idx = [i for i, x in enumerate(sp) if x != "H"]
        if m_idx:
            cnt = [sum(1 for nb in st.get_neighbors(st[i], 2.0) if nb.specie.symbol == "H") for i in m_idx]
            out["n_H_neighbours_M"] = float(np.mean(cnt))
    else:
        out.update(min_HH_A=None, min_MH_A=None, motif="no_H", n_H_neighbours_M=0.0)
    return out


__all__ = ["features", "LIGHT"]
