"""Elemental references for the S4 heuristic hull (architecture §3.1 S4; doc 22 §4).

Campaign 1 found that ``S4Electronic`` never passed ``e_form_eV`` to ``decide_s4``, so every S4
outcome carried ``H7_HEUR_NOT_CHECKED``. This module supplies the missing piece: the energy
per atom of each element at the **same** profile as the candidate (s4: 60 Ry, SG15), computed
with :func:`rtt.dft.refhull.phase_energy` (vc-relax + scf, cached by structure, pseudopotential
md5 and profile hash).

Reference structures:

* H: H2 in a 10 A box (:func:`rtt.dft.refhull.h2_molecule`), fixed-cell relax at Gamma.
* every other element: the lowest-energy-per-atom elemental entry of the offline Alex-MP
  reference set (``rtt.source.external.refquery``) with at most ``max_atoms`` atoms. Choosing a
  higher-energy polymorph can only raise mu_el and lower E_form, so the S4 quantity stays a
  lower bound on E_hull (the removal direction stays conservative); the chosen entry id is
  recorded.

E_form = (E - sum n_el mu_el)/N (``formation_energy_per_atom_eV``). It is a lower bound on E_hull
and not ZPE-corrected; tiers ``S4.ehull_heur_max_eV`` already carries the ZPE allowance.
"""
from __future__ import annotations

import contextlib
from pathlib import Path
from typing import Any, Callable

from .refhull import RY_EV, h2_molecule, phase_energy  # noqa: F401  (RY_EV re-exported)


@contextlib.contextmanager
def _lock(path: Path):
    try:
        import fcntl
    except ImportError:          # non-POSIX: no lock (single worker)
        yield
        return
    with open(path, "w") as fh:
        fcntl.flock(fh, fcntl.LOCK_EX)
        try:
            yield
        finally:
            fcntl.flock(fh, fcntl.LOCK_UN)


def elemental_structure(el: str, *, machine: Any = None, max_atoms: int = 8,
                        query: Callable[..., dict] | None = None) -> tuple[dict, str]:
    """(structure dict, provenance id) of the reference phase of ``el``."""
    if el == "H":
        return h2_molecule(), "H2_box10A"
    if query is None:
        from ..source.external import refquery as query
    data = query(el, machine=machine)
    rows = [r for r in data["rows"] if r["chemsys"] == el and r["n_atoms"] <= max_atoms]
    if not rows:
        raise LookupError(f"no elemental Alex-MP entry for {el} with <= {max_atoms} atoms")
    best = min(rows, key=lambda r: (r["energy_per_atom_eV"], r["n_atoms"]))
    st = best["structure"]
    return ({"lattice": st["lattice"], "species": st["species"],
             "frac_coords": [[x % 1.0 for x in fc] for fc in st["frac_coords"]]},
            f"alexmp:{best['entry_id']}")


def element_refs(elements: list[str], profile: Any, workroot: str | Path, *, launcher: Any = None,
                 cache_dir: str | Path | None = None, machine: Any = None,
                 energy_fn: Callable[..., Any] | None = None,
                 query: Callable[..., dict] | None = None) -> tuple[dict[str, float], dict]:
    """mu_el (eV/atom) per element at ``profile`` and an info dict {el: {ref, e_per_atom_eV} |
    {error}}. A failed element is reported, not raised; the caller must then not compute
    E_form (H7-heur stays NOT CHECKED)."""
    energy_fn = energy_fn or phase_energy
    root = Path(workroot).expanduser()
    mu: dict[str, float] = {}
    info: dict[str, dict] = {}
    root.mkdir(parents=True, exist_ok=True)
    for el in sorted(set(elements)):
        try:
            st, ref = elemental_structure(el, machine=machine, query=query)
            with _lock(root / f"elem_{el}.lock"):       # parallel S4 workers share the cache
                pe = energy_fn(f"elem_{el}", st, profile, root / f"elem_{el}", launcher=launcher,
                               cache_dir=cache_dir, molecule=(el == "H"))
            mu[el] = float(pe.e_per_atom_eV)
            info[el] = {"ref": ref, "e_per_atom_eV": mu[el], "profile": getattr(profile, "name", None)}
        except Exception as ex:  # noqa: BLE001 - reported; no E_form then
            info[el] = {"error": f"{type(ex).__name__}: {ex}"[:300]}
    return mu, info


def formation_energy(e_total_Ry: float | None, composition: dict[str, int],
                     mu: dict[str, float]) -> float | None:
    if e_total_Ry is None or not composition or any(el not in mu for el in composition):
        return None
    n = sum(composition.values())
    return (e_total_Ry * RY_EV - sum(c * mu[el] for el, c in composition.items())) / n


__all__ = ["elemental_structure", "element_refs", "formation_energy"]
