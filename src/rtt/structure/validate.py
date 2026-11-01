"""S2 validity check (architecture §3.1 S2: ``PASS | INVALID``).

Checks (each recorded as criterion / value / threshold):

* ``finite``: lattice and coordinates finite, lattice non-degenerate, species known elements;
* ``min_distance_A``: shortest interatomic distance (periodic) >= ``d_min`` (default 0.5 Å).
  For H-H pairs the same floor applies: an H2 molecule (0.74 Å) is chemically possible and is
  only *penalised* in ranking (doc 02 S4), never invalid;
* ``volume_per_atom_A3`` inside [v_min, v_max] (defaults 1.0 and 100 Å^3/atom: catches
  collapsed or exploded MatterGen samples, not chemistry).

Thresholds come from ``tiers.toml [S2]`` when present, otherwise the defaults below with source
id ``doc21#S2-defaults`` (docs/design/21_build_front.md). INVALID is a structural statement, not
a physics screen, so it is not part of the heuristic-screen recall accounting.
"""
from __future__ import annotations

from typing import Any

import numpy as np

from ..records import INVALID, PASS, Outcome
from .canonical import to_pymatgen

DEFAULTS = {"d_min_A": 0.5, "v_min_A3": 1.0, "v_max_A3": 100.0}
DEFAULT_SOURCE = "doc21#S2-defaults"


def _param(tiers: Any, key: str) -> tuple[float, str, bool]:
    if tiers is not None:
        try:
            th = tiers.threshold(f"S2.{key}")
            return float(th.value), th.source_id, th.provisional
        except Exception:  # noqa: BLE001 - key absent: use default
            pass
    return DEFAULTS[key], DEFAULT_SOURCE, False


def min_distance(st) -> float:
    if len(st) == 1:
        # single site: nearest periodic image
        return float(min(np.linalg.norm(v) for v in _short_lattice_vectors(st.lattice.matrix)))
    dm = st.distance_matrix.copy()
    np.fill_diagonal(dm, np.inf)
    d = float(dm.min())
    # periodic self-images can be shorter than any pair distance in tiny cells
    d_self = float(min(np.linalg.norm(v) for v in _short_lattice_vectors(st.lattice.matrix)))
    return min(d, d_self)


def _short_lattice_vectors(m):
    m = np.asarray(m)
    out = []
    for i in (-1, 0, 1):
        for j in (-1, 0, 1):
            for k in (-1, 0, 1):
                if (i, j, k) != (0, 0, 0):
                    out.append(i * m[0] + j * m[1] + k * m[2])
    return out


def checks(s: Any, tiers: Any = None) -> list[dict]:
    """All checks as dicts {criterion, value, threshold, source_id, ok}."""
    out: list[dict] = []
    try:
        from pymatgen.core.periodic_table import Element
        d = s.to_json() if hasattr(s, "to_json") else s
        if isinstance(d, dict):
            lat = np.asarray(d["lattice"], float)
            fc = np.asarray(d["frac_coords"], float)
            ok = (lat.shape == (3, 3) and np.all(np.isfinite(lat)) and fc.ndim == 2 and fc.shape[1] == 3
                  and np.all(np.isfinite(fc)) and len(d["species"]) == len(fc) and len(fc) > 0
                  and abs(np.linalg.det(lat)) > 1e-3)
            for el in d["species"]:
                Element(el)
        else:
            ok = True
        st = to_pymatgen(d)
    except Exception as ex:  # noqa: BLE001
        return [{"criterion": "finite", "value": f"{type(ex).__name__}: {ex}", "threshold": "parseable",
                 "source_id": DEFAULT_SOURCE, "ok": False}]
    out.append({"criterion": "finite", "value": bool(ok), "threshold": True, "source_id": DEFAULT_SOURCE,
                "ok": bool(ok)})
    if not ok:
        return out
    dmin, src, _ = _param(tiers, "d_min_A")
    md = min_distance(st)
    out.append({"criterion": "min_distance_A", "value": round(md, 4), "threshold": dmin, "source_id": src,
                "ok": md >= dmin})
    vpa = st.volume / len(st)
    vmin, src1, _ = _param(tiers, "v_min_A3")
    vmax, src2, _ = _param(tiers, "v_max_A3")
    out.append({"criterion": "volume_per_atom_A3", "value": round(vpa, 4), "threshold": [vmin, vmax],
                "source_id": src1 if src1 == src2 else f"{src1};{src2}", "ok": vmin <= vpa <= vmax})
    return out


def validate(s: Any, tiers: Any = None) -> Outcome:
    """S2 validity: ``PASS`` or ``INVALID`` (first failing check)."""
    for c in checks(s, tiers):
        if not c["ok"]:
            return Outcome(stage="S2", outcome=INVALID, criterion=c["criterion"], value=c["value"],
                           threshold=c["threshold"], source_id=c["source_id"])
    return Outcome(stage="S2", outcome=PASS)
