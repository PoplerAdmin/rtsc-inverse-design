"""Crystal-structure normalisation for the QE writers.

Inputs may be a pymatgen ``Structure``, an ASE ``Atoms``, a ``CellData`` or a plain dict
``{"lattice": 3x3 (angstrom, rows = vectors), "species": [...], "frac_coords": Nx3,
"masses": {El: amu}}``. pymatgen and ASE are optional imports: we duck-type on attributes and
never import them at module load.

Atomic masses: taken from pymatgen / ASE when the structure comes from them; for plain dicts
they must be supplied explicitly (``masses``) -- we do not ship a hand-typed mass table.
"""
from __future__ import annotations

import math
import re
from dataclasses import dataclass, field
from typing import Any

import numpy as np

BOHR_TO_ANGSTROM = 0.529177210903  # CODATA 2018 (QE: BOHR_RADIUS_ANGS)


@dataclass
class CellData:
    lattice_A: np.ndarray            # (3, 3) rows are lattice vectors, angstrom
    species: list[str]               # element symbol per site
    frac_coords: np.ndarray          # (N, 3) fractional coordinates
    masses_amu: dict[str, float] = field(default_factory=dict)

    def __post_init__(self) -> None:
        self.lattice_A = np.asarray(self.lattice_A, dtype=float).reshape(3, 3)
        self.frac_coords = np.asarray(self.frac_coords, dtype=float).reshape(-1, 3)
        self.species = [str(s) for s in self.species]
        if len(self.species) != len(self.frac_coords):
            raise ValueError("species and frac_coords lengths differ")
        if abs(np.linalg.det(self.lattice_A)) < 1e-8:
            raise ValueError("singular lattice")
        for el in self.elements:
            if el not in self.masses_amu:
                raise ValueError(f"no atomic mass for element {el!r}; pass masses explicitly")

    @property
    def elements(self) -> list[str]:
        out: list[str] = []
        for s in self.species:
            if s not in out:
                out.append(s)
        return out

    @property
    def nat(self) -> int:
        return len(self.species)

    @property
    def volume_A3(self) -> float:
        return float(abs(np.linalg.det(self.lattice_A)))

    def reciprocal_lengths_2pi(self) -> np.ndarray:
        """|b_i| including the 2*pi factor (1/angstrom)."""
        rec = 2.0 * math.pi * np.linalg.inv(self.lattice_A).T
        return np.linalg.norm(rec, axis=1)

    def to_dict(self) -> dict[str, Any]:
        return {"lattice_A": self.lattice_A.tolist(), "species": list(self.species),
                "frac_coords": self.frac_coords.tolist(), "masses_amu": dict(self.masses_amu)}

    @classmethod
    def from_dict(cls, d: dict[str, Any]) -> "CellData":
        lat = d.get("lattice_A", d.get("lattice"))
        masses = d.get("masses_amu", d.get("masses", {}))
        return cls(np.asarray(lat), list(d["species"]), np.asarray(d["frac_coords"]), dict(masses))


def _strip_symbol(s: str) -> str:
    m = re.match(r"([A-Z][a-z]?)", s)
    if not m:
        raise ValueError(f"cannot read element symbol from {s!r}")
    return m.group(1)


def to_celldata(obj: Any, masses: dict[str, float] | None = None) -> CellData:
    """Convert pymatgen Structure / ASE Atoms / dict / CellData to ``CellData``.

    ``masses`` (amu) overrides masses from the source object.
    """
    masses = dict(masses or {})
    if isinstance(obj, CellData):
        m = dict(obj.masses_amu)
        m.update(masses)
        return CellData(obj.lattice_A.copy(), list(obj.species), obj.frac_coords.copy(), m)
    if isinstance(obj, dict):
        d = dict(obj)
        m = dict(d.get("masses_amu", d.get("masses", {})))
        m.update(masses)
        d["masses_amu"] = m
        d.pop("masses", None)
        return CellData.from_dict(d)
    # pymatgen Structure (duck typing: .lattice.matrix, .sites, .frac_coords)
    if hasattr(obj, "lattice") and hasattr(obj, "frac_coords") and hasattr(obj, "sites"):
        if not getattr(obj, "is_ordered", True):
            raise ValueError("disordered pymatgen Structure is not supported (order it first)")
        species, m = [], {}
        for site in obj.sites:
            sp = site.specie
            el = getattr(sp, "element", sp)
            sym = _strip_symbol(getattr(el, "symbol", str(el)))
            species.append(sym)
            if sym not in m:
                am = getattr(el, "atomic_mass", None)
                if am is not None:
                    m[sym] = float(am)
        m.update(masses)
        return CellData(np.asarray(obj.lattice.matrix), species, np.asarray(obj.frac_coords), m)
    # ASE Atoms (duck typing)
    if hasattr(obj, "get_chemical_symbols") and hasattr(obj, "get_scaled_positions"):
        syms = [_strip_symbol(s) for s in obj.get_chemical_symbols()]
        am = obj.get_masses()
        m = {}
        for s, a in zip(syms, am):
            m.setdefault(s, float(a))
        m.update(masses)
        if not all(bool(p) for p in obj.get_pbc()):
            raise ValueError("ASE Atoms must be periodic in all three directions")
        return CellData(np.asarray(obj.cell.array if hasattr(obj.cell, "array") else obj.cell),
                        syms, np.asarray(obj.get_scaled_positions(wrap=False)), m)
    raise TypeError(f"unsupported structure type {type(obj)!r}")


# ------------------------------------------------------------------------------------------
# Parse the final structure of a pw.x vc-relax / relax run
# ------------------------------------------------------------------------------------------
_FLOAT = r"[-+]?(?:\d+\.\d*|\.\d+|\d+)(?:[eEdD][-+]?\d+)?"


def _f(s: str) -> float:
    return float(s.replace("d", "e").replace("D", "e"))


def parse_pw_final_structure(text: str, template: CellData) -> CellData:
    """Read the last CELL_PARAMETERS / ATOMIC_POSITIONS block of a pw.x output.

    Prefers the block between "Begin final coordinates" and "End final coordinates";
    otherwise the last printed blocks. ``template`` provides masses and is used for the cell
    when a fixed-cell relax prints no CELL_PARAMETERS.
    """
    i0 = text.rfind("Begin final coordinates")
    block = text[i0:] if i0 >= 0 else text
    i1 = block.find("End final coordinates")
    if i0 >= 0 and i1 >= 0:
        block = block[:i1]
    lines = block.splitlines()
    lattice = None
    pos_species: list[str] = []
    pos: list[list[float]] = []
    pos_unit = None
    alat_bohr = None
    m = re.search(r"celldm\(1\)=\s*(" + _FLOAT + ")", text)
    if m:
        alat_bohr = _f(m.group(1))
    idx_cell = [i for i, l in enumerate(lines) if l.strip().startswith("CELL_PARAMETERS")]
    if idx_cell:
        i = idx_cell[-1]
        head = lines[i]
        rows = [[_f(x) for x in lines[i + k].split()[:3]] for k in (1, 2, 3)]
        lat = np.asarray(rows, dtype=float)
        ma = re.search(r"alat\s*=\s*(" + _FLOAT + ")", head)
        if "angstrom" in head.lower():
            lattice = lat
        elif "bohr" in head.lower():
            lattice = lat * BOHR_TO_ANGSTROM
        elif ma or "alat" in head.lower():
            a = _f(ma.group(1)) if ma else alat_bohr
            if a is None:
                raise ValueError("CELL_PARAMETERS in alat units but alat unknown")
            lattice = lat * a * BOHR_TO_ANGSTROM
        else:
            lattice = lat  # QE default for CELL_PARAMETERS without unit is bohr/alat; assume A
    else:
        lattice = template.lattice_A.copy()
    idx_pos = [i for i, l in enumerate(lines) if l.strip().startswith("ATOMIC_POSITIONS")]
    if not idx_pos:
        raise ValueError("no ATOMIC_POSITIONS block found in pw.x output")
    i = idx_pos[-1]
    head = lines[i].lower()
    pos_unit = ("crystal" if "crystal" in head else "angstrom" if "angstrom" in head
                else "bohr" if "bohr" in head else "alat")
    for l in lines[i + 1:i + 1 + template.nat]:
        p = l.split()
        pos_species.append(_strip_symbol(p[0]))
        pos.append([_f(x) for x in p[1:4]])
    P = np.asarray(pos, dtype=float)
    if pos_unit == "crystal":
        frac = P
    else:
        if pos_unit == "angstrom":
            cart = P
        elif pos_unit == "bohr":
            cart = P * BOHR_TO_ANGSTROM
        else:
            if alat_bohr is None:
                raise ValueError("positions in alat units but alat unknown")
            cart = P * alat_bohr * BOHR_TO_ANGSTROM
        frac = cart @ np.linalg.inv(lattice)
    return CellData(lattice, pos_species, frac, dict(template.masses_amu))
