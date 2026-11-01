"""Database baseline lane D (architecture §3.3): known compounds from the offline Alex-MP
reference set (doc 10 §3.3; 845,997 entries, MP2020-corrected PBE(+U) energies).

Query = all entries of the chemical system and its subsystems (read in the MatterGen venv by
``mgvenv_tool.py refquery``), kept when

* the entry has >= 2 elements, and contains H when the query system contains H;
* its DFT E_above_hull on the hull of the same entries is <= ``ehull_max`` (0.1 eV/atom).

**Deviation (doc 21):** the reference set stores no band gap (entries carry energy, structure,
POTCARs and ids only), so ``gap_max`` cannot be applied here. Every lane-D structure is flagged
``GAP_UNKNOWN``; metallicity is decided at S4 (H1), which is where §3.1 puts it anyway.

Ordering: H fraction descending, then number of atoms ascending, then E_hull.
"""
from __future__ import annotations

from typing import Any

from ..records import Origin, RawStructure, StructureData
from .external import refquery
from .lanes import norm_chemsys


def _hf(species: list[str]) -> float:
    return sum(1 for s in species if s == "H") / max(1, len(species))


def query_alexmp(chemsys: str, gap_max: float | None = 0.0, ehull_max: float = 0.1, *,
                 machine: Any = None, cache_dir: Any = None, max_atoms: int | None = None) -> list[RawStructure]:
    cs = norm_chemsys(chemsys)
    want_h = "H" in cs.split("-")
    data = refquery(cs, machine=machine, cache_dir=cache_dir)
    rows = []
    for r in data["rows"]:
        els = r["chemsys"].split("-")
        if len(els) < 2 or (want_h and "H" not in els):
            continue
        eh = r.get("e_above_hull_dft_eV")
        if eh is None or eh > ehull_max:
            continue
        if max_atoms is not None and r["n_atoms"] > max_atoms:
            continue
        rows.append(r)
    rows.sort(key=lambda r: (-_hf(r["structure"]["species"]), r["n_atoms"], r["e_above_hull_dft_eV"]))
    out = []
    for r in rows:
        sd = StructureData(lattice=r["structure"]["lattice"], species=r["structure"]["species"],
                           frac_coords=[[x % 1.0 for x in fc] for fc in r["structure"]["frac_coords"]])
        origin = Origin(lane="D", generator="alexmp_reference_MP2020correction", conditioning=None,
                        batch_id=f"D-{cs}", flags=["KNOWN_COMPOUND", "GAP_UNKNOWN"])
        out.append(RawStructure(structure=sd, origin=origin, pressure_GPa=0.0,
                                source_ref=f"alexmp:{r['entry_id']}",
                                meta={"entry_id": r["entry_id"], "formula": r["formula"],
                                      "e_above_hull_dft_eV": r["e_above_hull_dft_eV"],
                                      "dft_energy_per_atom_eV": r["energy_per_atom_eV"],
                                      "query_chemsys": cs, "ehull_max": ehull_max,
                                      "gap_max_requested": gap_max, "gap_filter_applied": False,
                                      "complete_hull": data.get("complete_hull")}))
    return out


def query_many(systems: list[str], **kw: Any) -> list[RawStructure]:
    """Union over systems; an entry reached from several queries is returned once."""
    seen, out = set(), []
    for cs in systems:
        for rs in query_alexmp(cs, **kw):
            if rs.source_ref in seen:
                continue
            seen.add(rs.source_ref)
            out.append(rs)
    return out


__all__ = ["query_alexmp", "query_many"]
