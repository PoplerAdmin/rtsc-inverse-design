"""Canonical structure representation and hash (architecture §5.1 ``structure/canonical.py``).

* A *structure dict* (``sdict``) is ``{"lattice": 3x3 Å, "species": [...], "frac_coords": [...]}``,
  the form stored in the candidate record (§6.1). Extra keys (``spacegroup``, ``n_atoms``) are
  added by :func:`to_sdict`.
* :func:`canonicalize` returns the spglib standardized primitive cell (pymatgen
  ``SpacegroupAnalyzer.get_primitive_standard_structure``) at ``symprec``; if spglib fails the
  Niggli-reduced cell is used. Sites are sorted by (Z, frac coords) after wrapping into [0, 1).
* :func:`structure_sha256` hashes the canonical cell with lattice rounded to 1e-4 Å and
  fractional coordinates to 1e-5. Two cells that differ only by setting / origin / supercell
  give the same hash if spglib standardizes them identically; *equality* is decided by
  StructureMatcher in :mod:`rtt.structure.dedupe`, never by the hash alone.
"""
from __future__ import annotations

import hashlib
import json
from typing import Any

import numpy as np

SYMPREC = 1e-3          # tight: standardize without symmetrizing a distorted cell
ANGLE_TOL = 5.0
HASH_VERSION = "rtt.sha/1"


def to_pymatgen(s: Any):
    """sdict | pymatgen Structure | ase Atoms -> pymatgen Structure."""
    from pymatgen.core import Lattice, Structure
    if isinstance(s, Structure):
        return s
    if hasattr(s, "get_chemical_symbols") and hasattr(s, "cell"):
        from pymatgen.io.ase import AseAtomsAdaptor
        return AseAtomsAdaptor.get_structure(s)
    if isinstance(s, dict):
        return Structure(Lattice(np.asarray(s["lattice"], float)), list(s["species"]),
                         np.asarray(s["frac_coords"], float))
    raise TypeError(f"cannot convert {type(s).__name__} to a structure")


def spacegroup_number(s: Any, symprec: float = 0.01) -> int | None:
    from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
    try:
        return int(SpacegroupAnalyzer(to_pymatgen(s), symprec=symprec).get_space_group_number())
    except Exception:  # noqa: BLE001
        return None


def to_sdict(s: Any, *, with_symmetry: bool = True) -> dict[str, Any]:
    st = to_pymatgen(s)
    d = {"lattice": [[float(x) for x in r] for r in st.lattice.matrix],
         "species": [site.specie.symbol for site in st],
         "frac_coords": [[float(x) % 1.0 for x in site.frac_coords] for site in st],
         "n_atoms": len(st)}
    if with_symmetry:
        d["spacegroup"] = spacegroup_number(st)
    return d


def _sorted(st):
    from pymatgen.core import Structure
    wrapped = [(site.specie.Z, tuple(np.round(np.mod(site.frac_coords, 1.0), 6) % 1.0), site.specie.symbol)
               for site in st]
    wrapped.sort()
    return Structure(st.lattice, [w[2] for w in wrapped], [list(w[1]) for w in wrapped])


def canonicalize(rs: Any, symprec: float = SYMPREC):
    """Standardized primitive cell (pymatgen Structure), deterministic site order.

    ``rs`` may be a RawStructure-like object with a ``structure`` attribute, an sdict, a pymatgen
    Structure or ase Atoms."""
    from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
    obj = getattr(rs, "structure", rs)
    st = to_pymatgen(obj)
    try:
        prim = SpacegroupAnalyzer(st, symprec=symprec, angle_tolerance=ANGLE_TOL) \
            .get_primitive_standard_structure()
        if prim is None or len(prim) == 0 or abs(prim.volume * len(st) / len(prim) - st.volume) > 1e-3 * st.volume:
            raise ValueError("spglib primitive inconsistent with input volume")
    except Exception:  # noqa: BLE001 - degenerate cells: fall back to Niggli
        prim = st.get_reduced_structure("niggli")
    return _sorted(prim)


def structure_sha256(s: Any) -> str:
    st = to_pymatgen(s)
    payload = {"v": HASH_VERSION,
               "lattice": np.round(st.lattice.matrix, 4).tolist(),
               "species": [site.specie.symbol for site in st],
               "frac": (np.round(np.mod(st.frac_coords, 1.0), 5) % 1.0).tolist()}
    txt = json.dumps(payload, sort_keys=True, separators=(",", ":"))
    return hashlib.sha256(txt.encode()).hexdigest()


def reduced_formula(s: Any) -> str:
    return to_pymatgen(s).composition.reduced_formula


def chemsys_of(s: Any) -> str:
    return "-".join(sorted({site.specie.symbol for site in to_pymatgen(s)}))


#: tolerance (Å) at which the S1-S3 records label the space group (``spacegroup_number``)
DFT_SYMPREC = 0.01


def symmetrized_for_dft(s: Any, symprec: float = DFT_SYMPREC) -> tuple[dict[str, Any], dict[str, Any]]:
    """(sdict, info) for DFT input.

    Campaign 1 found that MatterGen cells carry ~1e-3-1e-2 Å noise: spglib labels them (e.g.
    Pm-3m at 0.01 Å, the tolerance of ``spacegroup_number``) but they are P1 at QE's tolerance,
    so pw.x/ph.x ran with "No symmetry found" (14 instead of 4 irreducible q on 3³ for a
    perovskite, no k-point reduction). If the cell has more symmetry at ``symprec`` than at
    1e-5 Å, return spglib's idealized standardized primitive cell; otherwise return the cell
    unchanged (exact cells keep their setting, so finished QE runs stay adoptable). The
    symmetry is then preserved by the DFT relaxation (a constraint, recorded in ``info``)."""
    import spglib
    d = s if isinstance(s, dict) else to_sdict(s, with_symmetry=False)
    els = sorted(set(d["species"]))
    cell = (np.asarray(d["lattice"], float), np.asarray(d["frac_coords"], float),
            [els.index(x) + 1 for x in d["species"]])
    try:
        exact = spglib.get_symmetry_dataset(cell, symprec=1e-5)
        loose = spglib.get_symmetry_dataset(cell, symprec=symprec)
    except Exception:  # noqa: BLE001
        exact = loose = None
    n_exact = exact.number if exact is not None else 1
    n_loose = loose.number if loose is not None else 1
    info = {"symprec_A": symprec, "spacegroup_exact": n_exact, "spacegroup_symprec": n_loose,
            "symmetrized": False}
    if loose is None or len(loose.rotations) <= (len(exact.rotations) if exact is not None else 1):
        return {k: d[k] for k in ("lattice", "species", "frac_coords")}, info
    lat, pos, nums = spglib.standardize_cell(cell, to_primitive=True, no_idealize=False, symprec=symprec)
    v_in = abs(np.linalg.det(cell[0])) / len(d["species"])
    if abs(abs(np.linalg.det(lat)) / len(nums) - v_in) > 0.05 * v_in:          # defensive; never seen
        return ({k: d[k] for k in ("lattice", "species", "frac_coords")},
                {**info, "note": "standardize changed the volume per atom; not symmetrized"})
    out = {"lattice": [[float(x) for x in r] for r in lat], "species": [els[n - 1] for n in nums],
           "frac_coords": [[float(x) % 1.0 for x in r] for r in pos]}
    info.update(symmetrized=True, n_atoms_in=len(d["species"]), n_atoms_out=len(nums))
    return out, info
