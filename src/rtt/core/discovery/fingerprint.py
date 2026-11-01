"""Structure validation, fingerprints and duplicate confirmation (DISCOVERY.md section 2.4).

``fingerprint(structure, pressure)`` is a deterministic, dependency-free hash used ONLY to find
CANDIDATE duplicates.  It is invariant under the choice of unit cell (supercells, origin shifts,
site permutations) and under isotropic scaling: per site, the sorted list of (neighbour species,
distance / (V/N)^(1/3)) within 2 x (V/N)^(1/3), binned at ``DIST_BIN``, the multiset of site
environments reduced by the gcd of its multiplicities, the reduced formula and the EXACT pressure
(rounded to 1e-6 GPa; no pressure binning).  Equal fingerprints do NOT imply equal structures
(e.g. stacking polytypes differ only beyond the cutoff).

A fingerprint hit becomes a duplicate only when CONFIRMED by ``StructureMatcher.fit`` (pymatgen,
if installed and enabled) or else by ``same_structure``: equal reduced composition and equal
per-atom multisets of absolute (angstrom) pair distances per species pair, out to a radius that
covers the reduced cell of both structures (sum of the reduced lattice-vector lengths) plus 2 A,
compared within ``confirm_tol_A``.  An unconfirmed hit is a NEW material (provenance
``possible_duplicate_of``).  A false "different" costs one extra computation; a false "same" would
lose a material, so confirmation is the strict side.
"""
from __future__ import annotations

import hashlib
import itertools
import json
import math
from collections import Counter
from functools import reduce
from typing import Any

import numpy as np

from ..sources.composition import formula_string, is_element, reduced_composition

DIST_BIN = 0.02
FP_VERSION = "fp1"


class InvalidStructure(ValueError):
    pass


def validate_structure(s: Any) -> dict[str, Any]:
    """Return a clean {"lattice", "species", "frac_coords"} dict or raise InvalidStructure."""
    if not isinstance(s, dict):
        raise InvalidStructure(f"structure must be a dict, got {type(s).__name__}")
    try:
        lat = np.asarray(s["lattice"], dtype=float)
        sp = [str(x) for x in s["species"]]
        fc = np.asarray(s["frac_coords"], dtype=float)
    except Exception as exc:  # noqa: BLE001
        raise InvalidStructure(f"malformed structure: {type(exc).__name__}: {exc}") from exc
    if lat.shape != (3, 3) or not np.all(np.isfinite(lat)):
        raise InvalidStructure("lattice must be a finite 3x3 matrix")
    vol = abs(float(np.linalg.det(lat)))
    if not vol > 1e-3:
        raise InvalidStructure(f"degenerate lattice (volume {vol})")
    if not sp or fc.shape != (len(sp), 3) or not np.all(np.isfinite(fc)):
        raise InvalidStructure("species / frac_coords mismatch or non-finite coordinates")
    bad = [x for x in sp if not is_element(x)]
    if bad:
        raise InvalidStructure(f"unknown element symbols {sorted(set(bad))}")
    if len(sp) > 400:
        raise InvalidStructure(f"{len(sp)} sites: too large for the discovery loop")
    out = {"lattice": lat.tolist(), "species": sp, "frac_coords": (fc % 1.0).tolist()}
    for k in ("spacegroup", "pressure_GPa", "enthalpy_eV_per_atom"):
        if k in s:
            out[k] = s[k]
    return out


def composition_of(s: dict[str, Any]) -> dict[str, float]:
    comp: dict[str, float] = {}
    for x in s["species"]:
        comp[x] = comp.get(x, 0.0) + 1.0
    return comp


def reduced_formula(s: dict[str, Any]) -> str:
    return formula_string(reduced_composition(composition_of(s)))


def pressure_key(p: float) -> str:
    """Exact pressure (to 1e-6 GPa): materials at different pressures are never merged."""
    return f"{round(float(p), 6):.6f}"


def coarse_key(s: dict[str, Any], pressure_GPa: float) -> str:
    return f"{reduced_formula(s)}|P{pressure_key(pressure_GPa)}"


def _site_environments(s: dict[str, Any]) -> list[tuple]:
    L = np.asarray(s["lattice"], float)
    F = np.asarray(s["frac_coords"], float) % 1.0
    sp = s["species"]
    n = len(sp)
    vol = abs(float(np.linalg.det(L)))
    scale = (vol / n) ** (1.0 / 3.0)
    rc = 2.0 * scale
    # number of images needed along each axis: rc / (plane spacing)
    widths = [vol / float(np.linalg.norm(np.cross(L[(i + 1) % 3], L[(i + 2) % 3]))) for i in range(3)]
    rng = [int(math.ceil(rc / w)) + 1 for w in widths]
    shifts = np.array(list(itertools.product(*(range(-r, r + 1) for r in rng))), float)
    cart = F @ L
    images = (cart[None, :, :] + (shifts @ L)[:, None, :]).reshape(-1, 3)
    img_sp = [sp[j] for _ in range(len(shifts)) for j in range(n)]
    envs = []
    for i in range(n):
        d = np.linalg.norm(images - cart[i], axis=1)
        mask = (d > 1e-6) & (d <= rc)
        env = sorted((img_sp[k], int(round(d[k] / scale / DIST_BIN))) for k in np.nonzero(mask)[0])
        envs.append((sp[i], tuple(env)))
    return envs


def fingerprint(structure: dict[str, Any], pressure_GPa: float) -> str:
    s = validate_structure(structure)
    cnt = Counter(_site_environments(s))
    g = reduce(math.gcd, cnt.values())
    items = sorted((repr(k), v // g) for k, v in cnt.items())
    key = json.dumps([FP_VERSION, reduced_formula(s), pressure_key(pressure_GPa),
                      items], separators=(",", ":"))
    return f"{FP_VERSION}:" + hashlib.sha256(key.encode()).hexdigest()[:32]


def _reduce(L: np.ndarray) -> np.ndarray:
    """Greedy (Gauss/Lagrange-type) lattice reduction: repeatedly shorten each vector by integer
    multiples of the others.  Gives short vectors, as Niggli/LLL would."""
    B = np.array(L, float)
    for _ in range(100):
        changed = False
        for i in range(3):
            for j in range(3):
                if i == j:
                    continue
                m = round(float(B[i] @ B[j]) / float(B[j] @ B[j]))
                if m:
                    nb = B[i] - m * B[j]
                    if nb @ nb < B[i] @ B[i] - 1e-12:
                        B[i] = nb
                        changed = True
        if not changed:
            break
    return B


def _pair_spectrum(s: dict[str, Any], R: float) -> dict[tuple[str, str], np.ndarray]:
    L = np.asarray(s["lattice"], float)
    F = np.asarray(s["frac_coords"], float) % 1.0
    sp = s["species"]
    vol = abs(float(np.linalg.det(L)))
    widths = [vol / float(np.linalg.norm(np.cross(L[(i + 1) % 3], L[(i + 2) % 3]))) for i in range(3)]
    rng = [int(math.ceil(R / w)) + 1 for w in widths]
    shifts = np.array(list(itertools.product(*(range(-r, r + 1) for r in rng))), float) @ L
    cart = F @ L
    out: dict[tuple[str, str], list] = {}
    for i in range(len(sp)):
        for j in range(len(sp)):
            d = np.linalg.norm(cart[j] + shifts - cart[i], axis=1)
            d = d[(d > 1e-6) & (d <= R)]
            out.setdefault((sp[i], sp[j]), []).append(d)
    return {k: np.sort(np.concatenate(v)) for k, v in out.items()}


def same_structure(a: dict[str, Any], b: dict[str, Any], tol: float = 0.02) -> bool:
    """Absolute-angstrom confirmation that a and b are the same crystal (see module doc)."""
    a, b = validate_structure(a), validate_structure(b)
    if reduced_formula(a) != reduced_formula(b):
        return False
    na, nb = len(a["species"]), len(b["species"])
    R = 2.0 + max(float(np.linalg.norm(_reduce(np.asarray(x["lattice"])), axis=1).sum()) for x in (a, b))
    if R > 40.0:
        return False                          # too large to confirm cheaply: treat as different
    pa, pb = _pair_spectrum(a, R + 2 * tol), _pair_spectrum(b, R + 2 * tol)
    if set(pa) != set(pb):
        return False
    for k in pa:
        da, db = pa[k][pa[k] <= R], pb[k][pb[k] <= R]
        # per-atom multiplicities must agree: na*len(db) == nb*len(da) (supercells scale both)
        if len(da) * nb != len(db) * na:
            return False
        if len(da) == 0:
            continue
        g = math.gcd(na, nb)
        ra, rb = nb // g, na // g          # repeat factors to equal multiplicity
        xa, xb = np.repeat(da, ra), np.repeat(db, rb)
        if len(xa) != len(xb) or float(np.max(np.abs(xa - xb))) > tol:
            return False
    return True


def invalid_fingerprint(raw: Any) -> str:
    return "invalid:" + hashlib.sha256(json.dumps(raw, sort_keys=True, default=str).encode()).hexdigest()[:32]


class StructureMatcherBackend:
    """Optional pymatgen StructureMatcher (None when pymatgen is missing or disabled)."""

    def __init__(self, ltol: float = 0.2, stol: float = 0.3, angle_tol: float = 5.0,
                 enabled: bool = True):
        self.matcher = None
        if enabled:
            try:
                from pymatgen.analysis.structure_matcher import StructureMatcher
                self.matcher = StructureMatcher(ltol=ltol, stol=stol, angle_tol=angle_tol)
            except ImportError:
                self.matcher = None

    @property
    def available(self) -> bool:
        return self.matcher is not None

    def fit(self, a: dict[str, Any], b: dict[str, Any]) -> bool:
        if self.matcher is None:
            return False
        from pymatgen.core import Structure
        try:
            return bool(self.matcher.fit(Structure(a["lattice"], a["species"], a["frac_coords"]),
                                         Structure(b["lattice"], b["species"], b["frac_coords"])))
        except Exception:  # noqa: BLE001 - a failed comparison means "not shown equal"
            return False
