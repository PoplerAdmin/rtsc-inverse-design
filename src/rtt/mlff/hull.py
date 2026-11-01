"""S3(b): ML convex hulls (architecture §3.1 S3; doc 10 §3.3).

Two definitions, as recommended by doc 10 §3.3:

* **MatterSim-only hull** (``ehull_self``; ranks and screens): every Alex-MP reference entry of
  the chemical system and its subsystems is relaxed with the *same* MatterSim model at the *same*
  pressure; candidate enthalpies are compared with that hull. Systematic MatterSim-vs-DFT offsets
  cancel to first order. No MP2020 corrections. Ported from ``src/ml/ehull_self.py``.
* **DFT-reference hull** (``ehull_dftref``; second opinion, recorded only): MatterSim energy +
  MP2020 correction (mattergen ``VasprunLike``) against the DFT reference energies. 0 GPa only.
  Ported from ``src/ml/ehull_ref.py``; runs in the MatterGen venv through ``mgvenv_tool.py``.

Reference cells with more than ``max_ref_atoms`` atoms are skipped (cost). Skipping a phase can
only *lower* the hull's coverage, i.e. make E_hull smaller, which is the lenient direction for
the only removal this hull feeds (``SCREENED_OUT(H7-ml)``); the count is recorded.
"""
from __future__ import annotations

import json
import os
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Sequence

from ..source.external import path_of, refquery, run_mgvenv_tool


@dataclass
class SelfHull:
    chemsys: str
    pressure_GPa: float
    model: str
    entries: list[dict]                        # {entry_id, formula, enthalpy_per_atom_eV, converged}
    n_ref_total: int = 0
    n_ref_skipped_size: int = 0
    n_ref_unconverged: int = 0
    elements_covered: list[str] = field(default_factory=list)
    _pd: Any = None

    @property
    def complete(self) -> bool:
        return set(self.chemsys.split("-")) <= set(self.elements_covered)

    def phase_diagram(self):
        if self._pd is None:
            from pymatgen.analysis.phase_diagram import PDEntry, PhaseDiagram
            from pymatgen.core import Composition
            ents = []
            for e in self.entries:
                if not e["converged"]:
                    continue
                c = Composition(e["formula"])
                ents.append(PDEntry(c, e["enthalpy_per_atom_eV"] * c.num_atoms, name=e["entry_id"]))
            self._pd = PhaseDiagram(ents)
        return self._pd


def _cache_file(chemsys: str, pressure_GPa: float, model: str, machine: Any = None) -> Path:
    d = path_of("front_cache", machine) / "mlhull" / f"{model}_P{pressure_GPa:g}"
    d.mkdir(parents=True, exist_ok=True)
    return d / f"{chemsys}.json"


def build_self_hull(chemsys: str, *, pressure_GPa: float = 0.0, model: str = "MatterSim-v1.0.0-5M",
                    device: str | None = None, machine: Any = None, max_ref_atoms: int = 64,
                    relax_fn: Any = None) -> SelfHull:
    """Relax the reference entries of ``chemsys`` (cached per chemsys / P / model)."""
    cs = "-".join(sorted(set(chemsys.split("-"))))
    f = _cache_file(cs, pressure_GPa, model, machine)
    if f.exists():
        d = json.loads(f.read_text())
        return SelfHull(**d)
    data = refquery(cs, machine=machine)
    rows = [r for r in data["rows"] if r["n_atoms"] <= max_ref_atoms]
    skipped = len(data["rows"]) - len(rows)
    if relax_fn is None:
        from .relax import relax_batch as relax_fn  # noqa: N813
    res = relax_fn([r["structure"] for r in rows], model=model, pressure_GPa=pressure_GPa, device=device) if rows else []
    entries = []
    for r, x in zip(rows, res):
        if x.error or x.enthalpy_per_atom_eV != x.enthalpy_per_atom_eV:   # NaN
            continue
        entries.append({"entry_id": r["entry_id"], "formula": r["formula"],
                        "enthalpy_per_atom_eV": x.enthalpy_per_atom_eV, "converged": bool(x.converged)})
    covered = sorted({e["formula"] for e in entries if e["converged"] and _is_element(e["formula"])})
    h = SelfHull(chemsys=cs, pressure_GPa=pressure_GPa, model=model, entries=entries,
                 n_ref_total=len(data["rows"]), n_ref_skipped_size=skipped,
                 n_ref_unconverged=sum(1 for e in entries if not e["converged"]),
                 elements_covered=[_element_of(x) for x in covered])
    tmp = f.with_suffix(".tmp")
    tmp.write_text(json.dumps({k: v for k, v in h.__dict__.items() if k != "_pd"}))
    os.replace(tmp, f)
    return h


def _is_element(formula: str) -> bool:
    from pymatgen.core import Composition
    return len(Composition(formula).elements) == 1


def _element_of(formula: str) -> str:
    from pymatgen.core import Composition
    return Composition(formula).elements[0].symbol


def ehull_self(struct: Any, enthalpy_per_atom_eV: float, chemsys: str | None = None, *,
               hull: SelfHull | None = None, **kw: Any) -> tuple[float | None, dict]:
    """(E_hull in eV/atom on the MatterSim-only hull, info). None if the hull is incomplete."""
    from pymatgen.analysis.phase_diagram import PDEntry
    from ..structure.canonical import to_pymatgen
    d = struct.to_json() if hasattr(struct, "to_json") else struct
    st = to_pymatgen(d)
    cs = chemsys or "-".join(sorted({s.specie.symbol for s in st}))
    hull = hull or build_self_hull(cs, **kw)
    info = {"chemsys": hull.chemsys, "n_ref": len(hull.entries), "n_ref_skipped_size": hull.n_ref_skipped_size,
            "pressure_GPa": hull.pressure_GPa, "model": hull.model}
    if not set(st.composition.chemical_system.split("-")) <= set(hull.elements_covered):
        info["note"] = "incomplete hull (missing elemental reference)"
        return None, info
    c = st.composition
    e = PDEntry(c, enthalpy_per_atom_eV * c.num_atoms, name="cand")
    decomp, eh = hull.phase_diagram().get_decomp_and_e_above_hull(e, allow_negative=True)
    info["decomposition"] = {k.composition.reduced_formula: round(v, 3) for k, v in decomp.items()}
    return round(float(eh), 5), info


def ehull_dftref(items: Sequence[tuple[str, Any, float]], *, machine: Any = None,
                 workdir: str | Path | None = None) -> dict[str, dict]:
    """Second-opinion hull for many structures at once.

    ``items``: (key, structure, MatterSim total energy eV) at 0 GPa. Returns key -> row
    ({e_above_hull_eV, mp2020_correction_eV_per_atom, note?})."""
    if not items:
        return {}
    wd = Path(workdir or path_of("front_cache", machine) / "dftref_tmp")
    wd.mkdir(parents=True, exist_ok=True)
    tag = f"{os.getpid()}_{abs(hash(tuple(k for k, _, _ in items))) % 10**8}"
    inp, out = wd / f"in_{tag}.json", wd / f"out_{tag}.json"
    payload = []
    for k, s, e in items:
        d = s.to_json() if hasattr(s, "to_json") else s
        payload.append({"key": k, "structure": {x: d[x] for x in ("lattice", "species", "frac_coords")},
                        "energy_eV": float(e)})
    inp.write_text(json.dumps({"items": payload}))
    try:
        run_mgvenv_tool(["ehull-dftref", "--lmdb", str(path_of("alexmp_lmdb", machine)),
                         "--inp", str(inp), "--out", str(out)], machine=machine)
        rows = json.loads(out.read_text())["rows"]
    finally:
        inp.unlink(missing_ok=True)
        out.unlink(missing_ok=True)
    return {r["key"]: r for r in rows}


__all__ = ["SelfHull", "build_self_hull", "ehull_self", "ehull_dftref"]
