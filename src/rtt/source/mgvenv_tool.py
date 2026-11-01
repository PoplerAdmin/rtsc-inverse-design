"""Helper that runs INSIDE the MatterGen venv (~/mattergen/.venv, Python 3.10).

It is executed as a standalone script by :mod:`rtt.source.database` and :mod:`rtt.mlff.hull`
through ``subprocess``; it must not import ``rtt`` (the main venv is Python 3.12 and has no
``lmdb``; the MatterGen venv has ``lmdb``, ``pymatgen`` and ``mattergen``).

Ported from the research scripts ``src/ml/ehull_ref.py`` (reader, MP2020 correction via
mattergen's ``VasprunLike``) -- logic copied, not imported.

Subcommands (all read/write JSON):

``refquery --lmdb L --chemsys A-B-C --out OUT``
    Every reference entry of the chemical system and of all its subsystems, with its structure
    (lattice / species / frac_coords), DFT energy (MP2020-corrected, as stored) and its DFT
    E_above_hull on the hull of those entries.

``ehull-dftref --lmdb L --inp IN --out OUT``
    IN: {"items": [{"key": ..., "structure": {...}, "energy_eV": float}, ...]} -- MatterSim
    total energies at 0 GPa. Applies MP2020 corrections (VasprunLike) and returns
    E_above_hull against the DFT reference hull ("second opinion" hull, doc 10 §3.3).

Security note: the LMDB values are pickles from Microsoft's official MatterGen data release
(doc 10 §3.3). Never point ``--lmdb`` at any other file.
"""
from __future__ import annotations

import argparse
import itertools
import json
import pickle
import sys


def _open(path):
    import lmdb
    return lmdb.open(str(path), subdir=False, readonly=True, lock=False, readahead=False)


def _get(env, key, default=None):
    with env.begin() as txn:
        v = txn.get(key.encode())
    return default if v is None else pickle.loads(v)


def entries_by_chemsys(env, cs):
    from pymatgen.entries.computed_entries import ComputedStructureEntry
    out = []
    for rf in _get(env, f"{cs}.reduced_formulas", []) or []:
        for n in range(_get(env, f"{cs}.{rf}.length", 0)):
            out.append(ComputedStructureEntry.from_dict(_get(env, f"{cs}.{rf}.{n}")))
    return out


def ref_entries(env, elements):
    out = []
    els = sorted(set(elements))
    for k in range(1, len(els) + 1):
        for sub in itertools.combinations(els, k):
            out.extend(entries_by_chemsys(env, "-".join(sub)))
    return out


def _has_all_elements(ents, els):
    have = {e.composition.elements[0].symbol for e in ents if e.composition.is_element}
    return all(x in have for x in els)


def _sdict(s):
    return {"lattice": [list(map(float, r)) for r in s.lattice.matrix],
            "species": [site.specie.symbol for site in s],
            "frac_coords": [list(map(float, site.frac_coords)) for site in s]}


def cmd_refquery(a):
    from pymatgen.analysis.phase_diagram import PhaseDiagram
    env = _open(a.lmdb)
    els = sorted(set(a.chemsys.split("-")))
    ents = ref_entries(env, els)
    complete = bool(ents) and _has_all_elements(ents, els)
    pd = PhaseDiagram(ents) if complete else None
    rows = []
    for e in ents:
        r = {"entry_id": str(e.entry_id), "formula": e.composition.reduced_formula,
             "chemsys": "-".join(sorted({x.symbol for x in e.composition.elements})),
             "n_atoms": int(e.structure.num_sites), "energy_eV": float(e.energy),
             "energy_per_atom_eV": float(e.energy) / e.structure.num_sites,
             "structure": _sdict(e.structure)}
        if pd is not None:
            r["e_above_hull_dft_eV"] = round(float(pd.get_e_above_hull(e, allow_negative=True)), 5)
        else:
            r["e_above_hull_dft_eV"] = None
        rows.append(r)
    json.dump({"chemsys": "-".join(els), "complete_hull": complete, "n": len(rows), "rows": rows},
              open(a.out, "w"))


def cmd_ehull_dftref(a):
    from pymatgen.analysis.phase_diagram import PhaseDiagram
    from pymatgen.core import Lattice, Structure
    from pymatgen.entries.compatibility import MaterialsProject2020Compatibility
    from mattergen.evaluation.utils.vasprunlike import VasprunLike
    env = _open(a.lmdb)
    items = json.load(open(a.inp))["items"]
    compat = MaterialsProject2020Compatibility()
    cache, rows = {}, []
    for it in items:
        sd = it["structure"]
        s = Structure(Lattice(sd["lattice"]), sd["species"], sd["frac_coords"])
        els = tuple(sorted({e.symbol for e in s.composition.elements}))
        if els not in cache:
            ents = ref_entries(env, els)
            cache[els] = PhaseDiagram(ents) if ents and _has_all_elements(ents, els) else None
        pd = cache[els]
        r = {"key": it["key"], "formula": s.composition.reduced_formula}
        try:
            entry = VasprunLike(s, float(it["energy_eV"])).get_computed_entry(
                inc_structure=False, energy_correction_scheme=compat)
        except Exception as ex:  # noqa: BLE001
            r.update(e_above_hull_eV=None, note=f"correction failed: {type(ex).__name__}: {ex}")
            rows.append(r)
            continue
        r["mp2020_correction_eV_per_atom"] = entry.correction / s.num_sites
        if pd is None:
            r.update(e_above_hull_eV=None, note="no complete reference hull for this chemical system")
        else:
            eh = pd.get_e_above_hull(entry, allow_negative=True)
            r["e_above_hull_eV"] = None if eh is None else round(float(eh), 5)
        rows.append(r)
    json.dump({"rows": rows}, open(a.out, "w"))


def main(argv=None):
    ap = argparse.ArgumentParser()
    sub = ap.add_subparsers(dest="cmd", required=True)
    q = sub.add_parser("refquery")
    q.add_argument("--lmdb", required=True)
    q.add_argument("--chemsys", required=True)
    q.add_argument("--out", required=True)
    h = sub.add_parser("ehull-dftref")
    h.add_argument("--lmdb", required=True)
    h.add_argument("--inp", required=True)
    h.add_argument("--out", required=True)
    a = ap.parse_args(argv)
    {"refquery": cmd_refquery, "ehull-dftref": cmd_ehull_dftref}[a.cmd](a)


if __name__ == "__main__":
    sys.exit(main())
