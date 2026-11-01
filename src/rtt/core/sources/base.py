"""Candidate record shared by all ORCH stages.

A ``Candidate`` is one input material.  It carries a composition, an optional crystal structure
(either a pymatgen ``Structure`` in ``structure`` or the light-weight ``lattice`` / ``species`` /
``frac_coords`` triple, which needs no third-party package), the pressure and enthalpy reported by
its source, free-form ``properties`` (e.g. ``energy_above_hull``, ``is_metal``, ``band_gap``,
``ordering``, ``total_magnetization``, ``dos_at_fermi``, ``phonon``) and ``provenance`` (where every
number came from).  Pre-computed alpha^2F inputs may be attached as ``spectra_files`` (paths with
explicit units) or ``spectra`` (inline ``Spectrum.to_dict()`` records).

``synthetic=True`` marks invented demo data; it is propagated into every report and flag so that
synthetic numbers are never presented as data for a real compound.
"""
from __future__ import annotations

import hashlib
import math
import re
from dataclasses import dataclass, field
from typing import Any

from .composition import formula_string, parse_formula, reduced_composition

_ID_SAFE = re.compile(r"[^A-Za-z0-9._+=-]+")


def safe_id(material_id: str) -> str:
    """Filesystem-safe, injective version of a material id (used for <workdir>/<id>/).

    Ids made only of safe characters are kept as they are.  Whenever sanitising changes the id,
    a hash of the ORIGINAL id is appended, so 'X/1' -> 'X_1-<hash>' never equals 'X_1'.
    (Case-insensitive file systems are handled by a collision check in pipeline.plan.)"""
    s = _ID_SAFE.sub("_", material_id).strip("._")
    if s == material_id and s:
        return s
    h = hashlib.sha256(material_id.encode("utf-8")).hexdigest()[:10]
    return f"{s or 'material'}-{h}"


@dataclass
class Candidate:
    material_id: str
    composition: dict[str, float]
    source: str = "unknown"
    formula: str = ""
    lattice: list[list[float]] | None = None          # 3x3 rows, angstrom
    species: list[str] | None = None
    frac_coords: list[list[float]] | None = None
    spacegroup: str | None = None
    pressure_GPa: float | None = None
    enthalpy_eV_per_atom: float | None = None         # at pressure_GPa as reported by the source
    properties: dict[str, Any] = field(default_factory=dict)
    provenance: dict[str, Any] = field(default_factory=dict)
    spectra_files: list[dict[str, Any]] = field(default_factory=list)
    spectra: list[dict[str, Any]] = field(default_factory=list)
    band_coupling_files: list[dict[str, Any]] = field(default_factory=list)   # {path, table}
    synthetic: bool = False
    structure: Any = field(default=None, repr=False, compare=False)   # pymatgen Structure (opt.)

    def __post_init__(self) -> None:
        if isinstance(self.composition, str):
            self.composition = parse_formula(self.composition)
        self.composition = {str(k): float(v) for k, v in dict(self.composition).items() if v}
        if not self.formula and self.composition:
            self.formula = formula_string(self.composition)

    # ------------------------------------------------------------------ helpers
    @property
    def reduced(self) -> dict[str, float]:
        return reduced_composition(self.composition)

    @property
    def elements(self) -> list[str]:
        return sorted(self.composition)

    @property
    def n_sites(self) -> float:
        return sum(self.composition.values())

    def lattice_lengths(self) -> tuple[float, float, float] | None:
        lat = self.lattice
        if lat is None and self.structure is not None:
            try:
                return tuple(float(x) for x in self.structure.lattice.abc)  # type: ignore
            except Exception:
                return None
        if lat is None:
            return None
        return tuple(math.sqrt(sum(x * x for x in row)) for row in lat)  # type: ignore

    def volume(self) -> float | None:
        if self.lattice is None:
            return None
        a, b, c = self.lattice
        return abs(a[0] * (b[1] * c[2] - b[2] * c[1]) - a[1] * (b[0] * c[2] - b[2] * c[0])
                   + a[2] * (b[0] * c[1] - b[1] * c[0]))

    def pymatgen_structure(self):
        """Return a pymatgen Structure (building it from lattice/species/coords) or None."""
        if self.structure is not None:
            return self.structure
        if self.lattice is None or not self.species or self.frac_coords is None:
            return None
        try:
            from pymatgen.core import Structure  # optional
        except ImportError:
            return None
        return Structure(self.lattice, self.species, self.frac_coords)

    # ------------------------------------------------------------------ serialisation
    def to_dict(self) -> dict[str, Any]:
        d = {
            "material_id": self.material_id, "composition": dict(self.composition),
            "source": self.source, "formula": self.formula, "lattice": self.lattice,
            "species": self.species, "frac_coords": self.frac_coords,
            "spacegroup": self.spacegroup, "pressure_GPa": self.pressure_GPa,
            "enthalpy_eV_per_atom": self.enthalpy_eV_per_atom,
            "properties": _jsonable(self.properties), "provenance": _jsonable(self.provenance),
            "spectra_files": list(self.spectra_files), "spectra": list(self.spectra),
            "band_coupling_files": list(self.band_coupling_files),
            "synthetic": self.synthetic,
        }
        if self.structure is not None and self.lattice is None:
            try:
                d["pmg_structure"] = self.structure.as_dict()
            except Exception:
                pass
        return d

    @classmethod
    def from_dict(cls, d: dict[str, Any]) -> "Candidate":
        d = dict(d)
        pmg = d.pop("pmg_structure", None)
        c = cls(**{k: v for k, v in d.items() if k in _FIELDS})
        if pmg is not None:
            try:
                from pymatgen.core import Structure
                c.structure = Structure.from_dict(pmg)
            except ImportError:
                pass
        return c


_FIELDS = {"material_id", "composition", "source", "formula", "lattice", "species", "frac_coords",
           "spacegroup", "pressure_GPa", "enthalpy_eV_per_atom", "properties", "provenance",
           "spectra_files", "spectra", "band_coupling_files", "synthetic"}


def _jsonable(x: Any) -> Any:
    """Best-effort conversion of numpy / enum / pymatgen objects into JSON-able values."""
    if isinstance(x, dict):
        return {str(k): _jsonable(v) for k, v in x.items()}
    if isinstance(x, (list, tuple)):
        return [_jsonable(v) for v in x]
    if isinstance(x, (str, int, float, bool)) or x is None:
        return x
    if hasattr(x, "tolist"):
        return x.tolist()
    if hasattr(x, "value") and isinstance(getattr(x, "value"), (str, int, float)):
        return x.value
    if hasattr(x, "as_dict"):
        try:
            return x.as_dict()
        except Exception:
            pass
    return str(x)


def candidate_from_structure_arrays(material_id: str, lattice, species, frac_coords, *,
                                    source: str, spacegroup: str | None = None,
                                    pressure_GPa: float | None = None,
                                    enthalpy_eV_per_atom: float | None = None,
                                    provenance: dict | None = None,
                                    properties: dict | None = None) -> Candidate:
    comp: dict[str, float] = {}
    for s in species:
        comp[s] = comp.get(s, 0.0) + 1.0
    return Candidate(material_id=material_id, composition=comp, source=source,
                     lattice=[list(map(float, r)) for r in lattice], species=list(species),
                     frac_coords=[list(map(float, r)) for r in frac_coords],
                     spacegroup=spacegroup, pressure_GPa=pressure_GPa,
                     enthalpy_eV_per_atom=enthalpy_eV_per_atom,
                     provenance=dict(provenance or {}), properties=dict(properties or {}))
