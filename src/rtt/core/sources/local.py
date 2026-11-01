"""Local directory source: CIF, POSCAR/CONTCAR/*.vasp, extended XYZ and AIRSS .res files.

Reader preference (``backend``): ``"auto"`` tries pymatgen, then ASE, then the built-in
dependency-free readers in ``fileformats``.  A file that cannot be read still yields a
``Candidate`` (with ``properties['load_error']``) so that the run summary accounts for it (S4);
the pipeline turns it into an ERROR verdict.

Optional manifest (``manifest.toml`` or ``manifest.json`` in the directory, or an explicit path)
attaches metadata per file stem or material id::

    ["LaH10_Fm-3m"]                     # file stem (or material id)
    material_id = "LaH10-Fm-3m-200GPa"  # optional override
    pressure_GPa = 200.0
    enthalpy_eV_per_atom = -3.21
    properties = { is_metal = true, energy_above_hull = 0.0 }
    spectra_files = [ { path = "a2F/LaH10.a2F.dos3", units = "Ry", label = "smearing 3" } ]
    band_coupling = "bands/LaH10_bc.toml"   # optional rtt.core.dfpt.BandCoupling file(s), for
                                         # [certify] rejection_model = "multiband"

Paths in ``spectra_files`` and ``band_coupling`` are resolved relative to the manifest's directory.
"""
from __future__ import annotations

import json
import logging
import os
import tomllib
from pathlib import Path
from typing import Any, Iterable

from .base import Candidate, candidate_from_structure_arrays
from .fileformats import RawStructure, read_cif, read_extxyz, read_poscar, read_res

log = logging.getLogger(__name__)

STRUCTURE_PATTERNS = ("*.cif", "*.vasp", "*.poscar", "POSCAR*", "CONTCAR*", "*.xyz", "*.extxyz",
                      "*.res")


def _kind(path: Path) -> str | None:
    n = path.name.lower()
    if n.endswith(".cif"):
        return "cif"
    if n.endswith((".xyz", ".extxyz")):
        return "extxyz"
    if n.endswith(".res"):
        return "res"
    if n.endswith((".vasp", ".poscar")) or n.startswith(("poscar", "contcar")):
        return "poscar"
    return None


def load_manifest(path: Path) -> dict[str, dict[str, Any]]:
    if path.suffix == ".toml":
        with open(path, "rb") as fh:
            data = tomllib.load(fh)
    else:
        data = json.loads(path.read_text())
    base = path.parent
    out = {}
    for key, meta in data.items():
        meta = dict(meta)
        files = []
        for sf in meta.get("spectra_files", []):
            sf = dict(sf)
            p = Path(sf["path"])
            sf["path"] = str(p if p.is_absolute() else (base / p).resolve())
            files.append(sf)
        if files:
            meta["spectra_files"] = files
        bcs = meta.get("band_coupling")
        if bcs is not None:
            bcs = bcs if isinstance(bcs, list) else [bcs]
            norm = []
            for bc in bcs:
                bc = {"path": bc} if isinstance(bc, str) else dict(bc)
                bp = Path(bc["path"])
                bc["path"] = str(bp if bp.is_absolute() else (base / bp).resolve())
                norm.append(bc)
            meta["band_coupling"] = norm
        out[str(key)] = meta
    return out


def _read_builtin(path: Path, kind: str) -> list[RawStructure]:
    text = path.read_text(errors="replace")
    if kind == "cif":
        return read_cif(text)
    if kind == "extxyz":
        return read_extxyz(text)
    if kind == "res":
        return [read_res(text)]
    return [read_poscar(text)]


def _read_pymatgen(path: Path, kind: str):
    from pymatgen.core import Structure  # noqa: F401  (ImportError -> next backend)
    if kind == "cif":
        from pymatgen.io.cif import CifParser
        structs = CifParser(str(path)).parse_structures(primitive=False)
        return [(s, {}) for s in structs]
    if kind == "res":
        raise ImportError("use builtin .res reader (keeps AIRSS TITL metadata)")
    if kind == "extxyz":
        raise ImportError("use ASE / builtin for extxyz (keeps per-frame info)")
    from pymatgen.core import Structure as S
    return [(S.from_file(str(path)), {})]


def _read_ase(path: Path, kind: str) -> list[RawStructure]:
    if kind == "res":
        raise ImportError("use builtin .res reader")
    import ase.io  # noqa: F401
    frames = ase.io.read(str(path), index=":")
    out = []
    for at in frames:
        info = {k: v for k, v in at.info.items() if isinstance(v, (int, float, str))}
        try:
            e = at.get_potential_energy()
            info.setdefault("energy", float(e))
        except Exception:
            pass
        out.append(RawStructure(at.cell.tolist(), at.get_chemical_symbols(),
                                at.get_scaled_positions(wrap=True).tolist(), info))
    return out


def _spacegroup_from_info(info: dict) -> str | None:
    sg = info.get("spacegroup")
    return str(sg).strip() if sg not in (None, "") else None


def read_structure_file(path: Path, backend: str = "auto") -> tuple[list[Any], str]:
    """Return ([RawStructure | (pymatgen Structure, info)], backend_used)."""
    kind = _kind(path)
    if kind is None:
        raise ValueError(f"unrecognised structure file type: {path.name}")
    order = {"auto": ("pymatgen", "ase", "builtin"), "pymatgen": ("pymatgen",),
             "ase": ("ase",), "builtin": ("builtin",)}[backend]
    last_exc: Exception | None = None
    for b in order:
        try:
            if b == "pymatgen":
                return _read_pymatgen(path, kind), b
            if b == "ase":
                return _read_ase(path, kind), b
            return _read_builtin(path, kind), b
        except ImportError as exc:
            last_exc = exc
            continue
    raise ImportError(f"no reader backend available for {path} ({last_exc})")


def _energy_per_atom(info: dict, natoms: int) -> tuple[float | None, str | None]:
    """Enthalpy (preferred) or energy per atom from file metadata, with the key used."""
    for key in ("enthalpy_per_atom", "enthalpy_eV_per_atom"):
        if key in info:
            return float(info[key]), key
    for key in ("enthalpy_eV", "enthalpy", "energy", "free_energy"):
        if key in info:
            return float(info[key]) / natoms, key + " (per cell / n_atoms)"
    return None, None


def load_local(path: str | os.PathLike, *, recursive: bool = False,
               patterns: Iterable[str] = STRUCTURE_PATTERNS, pressure_GPa: float | None = None,
               backend: str = "auto", manifest: str | os.PathLike | None = None,
               id_prefix: str = "") -> list[Candidate]:
    root = Path(path)
    if not root.is_dir():
        raise FileNotFoundError(f"local source directory not found: {root}")
    man: dict[str, dict] = {}
    if manifest is not None:
        man = load_manifest(Path(manifest))
    else:
        for name in ("manifest.toml", "manifest.json"):
            if (root / name).is_file():
                man = load_manifest(root / name)
                break
    files: set[Path] = set()
    for pat in patterns:
        files.update(root.rglob(pat) if recursive else root.glob(pat))
    out: list[Candidate] = []
    for f in sorted(p for p in files if p.is_file()):
        stem = f.stem if f.suffix else f.name
        rel = str(f.relative_to(root))
        try:
            items, used = read_structure_file(f, backend)
        except Exception as exc:  # noqa: BLE001 - recorded, becomes ERROR (S4)
            mid = id_prefix + stem
            meta = man.get(stem, man.get(mid, {}))
            out.append(Candidate(material_id=meta.get("material_id", mid), composition={},
                                 source="local", properties={"load_error": f"{type(exc).__name__}: {exc}"},
                                 provenance={"file": str(f.resolve())}))
            continue
        for k, item in enumerate(items):
            mid = id_prefix + (stem if len(items) == 1 else f"{stem}#{k}")
            meta = dict(man.get(mid[len(id_prefix):], man.get(stem, {})) if len(items) == 1
                        else man.get(mid[len(id_prefix):], {}))
            prov = {"source": "local", "file": str(f.resolve()), "relpath": rel, "frame": k,
                    "reader": used}
            if isinstance(item, RawStructure):
                info = dict(item.info)
                e_pa, e_key = _energy_per_atom(info, len(item.species))
                p = info.get("pressure_GPa", info.get("pressure", pressure_GPa))
                c = candidate_from_structure_arrays(
                    meta.get("material_id", mid), item.lattice, item.species, item.frac_coords,
                    source="local", spacegroup=_spacegroup_from_info(info),
                    pressure_GPa=None if p is None else float(p), enthalpy_eV_per_atom=e_pa,
                    provenance=prov)
                if e_key:
                    c.provenance["enthalpy_from"] = e_key
                if info.get("partial_occupancy"):
                    c.properties["partial_occupancy"] = True
                c.provenance["file_info"] = {k2: v for k2, v in info.items()
                                             if isinstance(v, (int, float, str, bool))}
            else:
                s, _info = item
                comp = {str(el): float(v) for el, v in s.composition.get_el_amt_dict().items()}
                c = Candidate(material_id=meta.get("material_id", mid), composition=comp,
                              source="local", lattice=s.lattice.matrix.tolist(),
                              species=[str(sp.specie) if hasattr(sp, "specie") else str(sp)
                                       for sp in s],
                              frac_coords=s.frac_coords.tolist(), pressure_GPa=pressure_GPa,
                              provenance=prov)
                if not s.is_ordered:
                    c.properties["partial_occupancy"] = True
                    c.species = None
                    c.structure = s
            _apply_meta(c, meta)
            out.append(c)
    return out


def _apply_meta(c: Candidate, meta: dict[str, Any]) -> None:
    if not meta:
        return
    for key in ("pressure_GPa", "enthalpy_eV_per_atom", "spacegroup"):
        if key in meta:
            setattr(c, key, meta[key])
    c.properties.update(meta.get("properties", {}))
    c.spectra_files.extend(meta.get("spectra_files", []))
    c.band_coupling_files.extend(meta.get("band_coupling", []))
    if meta.get("synthetic"):
        c.synthetic = True
    c.provenance["manifest"] = {k: v for k, v in meta.items()
                                if k not in ("spectra_files", "band_coupling")}
