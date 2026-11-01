"""De-duplication of candidates.

``method = "structure_matcher"`` uses pymatgen's ``StructureMatcher`` within groups of equal
reduced formula and (rounded) pressure; ``"composition_spacegroup"`` groups by reduced formula +
space group + pressure; ``"auto"`` uses the matcher when pymatgen is importable and both
structures exist, else the fallback.  Space groups are computed with spglib when available
(``symprec``), otherwise the source-reported symbol is used (normalised); a candidate with no
space group is never merged by the fallback.

Duplicates are NOT dropped: the pipeline marks them DEFERRED ("duplicate of <id>; not computed
separately"), so every input still ends in exactly one status (S4).  The representative of a
group is the member with the lowest reported enthalpy (ties: first seen).
"""
from __future__ import annotations

import logging
from dataclasses import dataclass, field

from .base import Candidate
from .composition import formula_string

log = logging.getLogger(__name__)


@dataclass
class DedupeResult:
    kept: list[Candidate]
    duplicate_of: dict[str, str] = field(default_factory=dict)   # duplicate id -> kept id
    method_used: dict[str, str] = field(default_factory=dict)    # kept id -> method


def _norm_sg(sg) -> str | None:
    if sg in (None, ""):
        return None
    return str(sg).replace(" ", "").replace("_", "")


def spacegroup_of(c: Candidate, symprec: float = 0.1) -> str | None:
    if c.lattice is not None and c.species and c.frac_coords is not None:
        try:
            import spglib  # optional
            from .composition import ATOMIC_NUMBER
            nums = [ATOMIC_NUMBER.get(s, 0) for s in c.species]
            if all(nums):
                res = spglib.get_spacegroup((c.lattice, c.frac_coords, nums), symprec=symprec)
                if res:
                    return _norm_sg(res.split()[0])
        except Exception:  # noqa: BLE001 - optional dependency or degenerate cell
            pass
    return _norm_sg(c.spacegroup)


def _pressure_key(c: Candidate, tol: float) -> float:
    p = c.pressure_GPa or 0.0
    return round(p / tol) * tol if tol > 0 else p


def _better(a: Candidate, b: Candidate) -> bool:
    ea, eb = a.enthalpy_eV_per_atom, b.enthalpy_eV_per_atom
    if ea is None:
        return False
    if eb is None:
        return True
    return ea < eb


def dedupe(candidates: list[Candidate], *, method: str = "auto", ltol: float = 0.2,
           stol: float = 0.3, angle_tol: float = 5.0, symprec: float = 0.1,
           pressure_tol_GPa: float = 1.0) -> DedupeResult:
    if method not in ("auto", "structure_matcher", "composition_spacegroup", "none"):
        raise ValueError(f"unknown dedupe method {method!r}")
    if method == "none":
        return DedupeResult(list(candidates))
    matcher = None
    if method in ("auto", "structure_matcher"):
        try:
            from pymatgen.analysis.structure_matcher import StructureMatcher
            matcher = StructureMatcher(ltol=ltol, stol=stol, angle_tol=angle_tol)
        except ImportError:
            if method == "structure_matcher":
                raise
            matcher = None
    groups: dict[tuple, list[Candidate]] = {}
    passthrough: list[Candidate] = []
    for c in candidates:
        if c.properties.get("load_error") or not c.composition:
            passthrough.append(c)
            continue
        key = (formula_string(c.composition), _pressure_key(c, pressure_tol_GPa))
        groups.setdefault(key, []).append(c)
    res = DedupeResult(kept=[])
    order = {id(c): i for i, c in enumerate(candidates)}
    kept_set: list[Candidate] = list(passthrough)
    for key, members in groups.items():
        clusters: list[list[Candidate]] = []
        used = "none"
        if matcher is not None and len(members) > 1:
            with_struct = [(m, m.pymatgen_structure()) for m in members]
            if all(s is not None for _, s in with_struct):
                used = "structure_matcher"
                for m, s in with_struct:
                    for cl in clusters:
                        if matcher.fit(cl[0].pymatgen_structure(), s):
                            cl.append(m)
                            break
                    else:
                        clusters.append([m])
        if used == "none":
            if len(members) > 1:
                used = "composition_spacegroup"
                by_sg: dict[str, list[Candidate]] = {}
                for m in members:
                    sg = spacegroup_of(m, symprec)
                    if sg is None:
                        clusters.append([m])
                    else:
                        by_sg.setdefault(sg, []).append(m)
                clusters.extend(by_sg.values())
            else:
                clusters = [members]
        for cl in clusters:
            rep = cl[0]
            for m in cl[1:]:
                if _better(m, rep):
                    rep = m
            kept_set.append(rep)
            res.method_used[rep.material_id] = used
            for m in cl:
                if m is not rep:
                    res.duplicate_of[m.material_id] = rep.material_id
    res.kept = sorted(kept_set, key=lambda c: order[id(c)])
    return res
