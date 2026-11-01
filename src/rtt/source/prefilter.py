"""S1 composition / size prefilter (architecture §3.1 S1: ``PASS | SCREENED_OUT | DEFERRED``).

R30 F1(b) / doc 31 §7: rules 2 and 3 are MatterGen-domain limits. They remove only generated
candidates (lanes G*). For lanes D and C, rule 2 gives ``DEFERRED(cost_class)`` and rule 3 is
not applied (PASS with a tag).

Removes (SCREENED_OUT, heuristic and reversible, R2) when

1. ``no_pseudo``: an element has no SG15 ONCV pseudopotential (``<El>_ONCV_PBE-*.upf``, scalar
   relativistic) in ``tiers [S1] pseudo_dir`` -- the DFT stages could not run it;
2. ``n_atoms_primitive``: the primitive cell has more than ``tiers [S1] n_atoms_max`` atoms (20;
   doc 10 §1, MatterGen training <= 20 atoms);
3. ``outside_mattergen_domain``: an element outside MatterGen's domain (noble gases, Tc, Pm,
   Z > 83).

Also computes (recorded, not used to remove): cost class (§8.1), H_f, element set.

Cost class (§8.1): C1 = <= 5 atoms with cubic or hexagonal lattice; C2 = 6-10 atoms cubic;
C3 = everything else up to ``n_atoms_max`` (6-10 atoms of moderate symmetry, and larger cells,
which §8.1 does not price; they are C3 here so the scheduler prices them as expensive).
"""
from __future__ import annotations

import re
from functools import lru_cache
from pathlib import Path
from typing import Any

from ..records import DEFERRED, PASS, SCREENED_OUT, Outcome
from .external import path_of
from .lanes import MATTERGEN_ELEMENTS

DEFAULT_N_ATOMS_MAX = (20, "doc10#1 (MatterGen training <=20 atoms)")


@lru_cache(maxsize=8)
def pseudo_elements(pseudo_dir: str) -> frozenset[str]:
    p = Path(pseudo_dir).expanduser()
    out = set()
    if p.is_dir():
        for f in p.iterdir():
            m = re.match(r"^([A-Z][a-z]?)_ONCV_PBE-[0-9.]+\.upf$", f.name)
            if m:
                out.add(m.group(1))
    return frozenset(out)


def _tier(tiers: Any, key: str, default: tuple[Any, str]) -> tuple[Any, str, bool]:
    if tiers is not None:
        try:
            th = tiers.threshold(key)
            return th.value, th.source_id, th.provisional
        except Exception:  # noqa: BLE001
            pass
    return default[0], default[1], False


def primitive_info(rs: Any) -> dict:
    """n_atoms of the primitive cell, crystal system, H_f, elements, cost class."""
    from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
    from ..structure.canonical import to_pymatgen
    sd = rs.structure.to_json() if hasattr(rs, "structure") else rs
    st = to_pymatgen(sd)
    try:
        sga = SpacegroupAnalyzer(st, symprec=0.01)
        prim = sga.find_primitive() or st
        system = sga.get_crystal_system()
        sg = sga.get_space_group_number()
    except Exception:  # noqa: BLE001
        prim, system, sg = st, None, None
    n = len(prim)
    if n <= 5 and system in ("cubic", "hexagonal"):
        cc = "C1"
    elif 6 <= n <= 10 and system == "cubic":
        cc = "C2"
    elif n <= 5:
        cc = "C2"      # <= 5 atoms but low symmetry: priced like a cubic 6-10-atom cell
    else:
        cc = "C3"
    species = [s.specie.symbol for s in st]
    return {"n_atoms_primitive": n, "crystal_system": system, "spacegroup": sg, "cost_class": cc,
            "H_f": sum(1 for x in species if x == "H") / len(species),
            "elements": sorted(set(species))}


def prefilter(rs: Any, tiers: Any = None, *, info: dict | None = None) -> Outcome:
    """S1 outcome for one RawStructure. ``info`` (from :func:`primitive_info`) may be passed in."""
    info = info or primitive_info(rs)
    els = info["elements"]
    pdir = None
    if tiers is not None:
        try:
            pdir = tiers.get("S1.pseudo_dir", None)
        except Exception:  # noqa: BLE001
            pdir = None
    pdir = str(Path(pdir).expanduser()) if pdir else str(path_of("pseudo_dir"))
    have = pseudo_elements(pdir)
    missing = [e for e in els if e not in have]
    if missing:
        return Outcome(stage="S1", outcome=SCREENED_OUT, criterion="no_pseudo", value=missing,
                       threshold=f"SG15 file in {Path(pdir).name}", source_id="doc12#3.1:S1;tiers:S1.pseudo_dir")
    lane = getattr(getattr(rs, "origin", None), "lane", None)
    generated = lane is None or str(lane).startswith("G")
    nmax, src, prov = _tier(tiers, "S1.n_atoms_max", DEFAULT_N_ATOMS_MAX)
    if info["n_atoms_primitive"] > int(nmax):
        if not generated:      # R30 F1(b): a generator-domain limit, not physics: park, do not remove
            return Outcome(stage="S1", outcome=DEFERRED, criterion="cost_class", value=info["n_atoms_primitive"],
                           threshold=int(nmax), source_id=f"{src}; doc31#7 (generator limit, lane {lane})",
                           provisional=prov, tags=["GENERATOR_DOMAIN_LIMIT_NOT_APPLIED"],
                           reason=(f"n_atoms_primitive {info['n_atoms_primitive']} > {int(nmax)}: beyond the MatterGen "
                                   f"training size, which does not apply to lane {lane}; deferred as too costly "
                                   "for the current DFT lanes (cost class above C3), not removed"))
        return Outcome(stage="S1", outcome=SCREENED_OUT, criterion="n_atoms_primitive",
                       value=info["n_atoms_primitive"], threshold=int(nmax), source_id=src, provisional=prov)
    bad = [e for e in els if e not in MATTERGEN_ELEMENTS]
    if bad:
        if not generated:      # same reasoning: MatterGen's element set does not bind lanes D and C
            return Outcome(stage="S1", outcome=PASS, tags=["OUTSIDE_MATTERGEN_DOMAIN_NOT_APPLIED"])
        return Outcome(stage="S1", outcome=SCREENED_OUT, criterion="outside_mattergen_domain", value=bad,
                       threshold="mattergen SELECTED_ATOMIC_NUMBERS", source_id="doc10#1;doc12#3.1:S1")
    return Outcome(stage="S1", outcome=PASS)


__all__ = ["prefilter", "primitive_info", "pseudo_elements"]
