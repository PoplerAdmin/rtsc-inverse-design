"""Energy above hull of MatterSim-relaxed structures vs the MatterGen Alex-MP reference hull.

Runs in the MatterGen venv (~/mattergen/.venv: pymatgen 2024.10.29, lmdb, mattergen evaluation code).
This reproduces what `mattergen-evaluate` does for its E_hull metric (MatterSim energy + MP2020
corrections via mattergen's VasprunLike, compared to DFT-PBE(+U)+MP2020 reference entries of the
845,997-entry Alexandria+MP reference set), but without loading the whole set into RAM:
the reference LMDB (3.85 GB uncompressed) is memory-mapped and queried per chemical system.

Prepare once (streamed decompression; the mattergen deserializer would read 1.7+ GB into RAM):
  ~/mattergen/.venv/bin/python -c "from mattergen.evaluation.reference.presets import \
      get_reference_dataset_path as g; print(g('reference_MP2020correction.gz'))"
  zcat <that path> > ~/rtsc_targeted_work/T3-mlaudit/ref/reference_MP2020correction.lmdb

Usage:
  ~/mattergen/.venv/bin/python src/ml/ehull_ref.py --lmdb REF.lmdb --relaxed W/relax_X/relaxed.json \
      --extxyz W/relax_X/relaxed.extxyz --out W/relax_X/ehull.json

Caveats (documented in docs/design/10_ml_capabilities.md):
  * mixes MatterSim (ML surrogate of PBE/PBE+U) energies with DFT reference energies: MatterSim
    errors (tens of meV/atom) go straight into E_hull;
  * 0 GPa only (reference entries are ambient-pressure DFT energies, not enthalpies);
  * the reference set contains entries pickled by Microsoft; unpickling a downloaded file is only
    acceptable because it is the official MatterGen data release.
"""
from __future__ import annotations

import argparse
import itertools
import json
from pathlib import Path

import pickle

import lmdb
from ase.io import read
from pymatgen.analysis.phase_diagram import PhaseDiagram
from pymatgen.entries.compatibility import MaterialsProject2020Compatibility
from pymatgen.entries.computed_entries import ComputedStructureEntry
from pymatgen.io.ase import AseAtomsAdaptor

from mattergen.evaluation.utils.vasprunlike import VasprunLike


class RefLMDB:
    """Minimal read-only reader of the MatterGen reference LMDB (key layout documented in
    mattergen/evaluation/reference/reference_dataset_serializer.py). Written because
    LMDBBackedReferenceDatasetImpl opens the same environment twice and fails with
    'environment already open in this process' on lmdb 3.0 (observed 2026-10-28)."""

    def __init__(self, path):
        self.env = lmdb.open(str(path), subdir=False, readonly=True, lock=False, readahead=False)

    def _get(self, key, default=None):
        with self.env.begin() as txn:
            v = txn.get(key.encode())
        return default if v is None else pickle.loads(v)

    def entries_by_chemsys(self, cs):
        out = []
        for rf in self._get(f"{cs}.reduced_formulas", []) or []:
            for n in range(self._get(f"{cs}.{rf}.length", 0)):
                out.append(ComputedStructureEntry.from_dict(self._get(f"{cs}.{rf}.{n}")))
        return out


def ref_entries(ref, elements):
    out = []
    for k in range(1, len(elements) + 1):
        for sub in itertools.combinations(sorted(elements), k):
            out.extend(ref.entries_by_chemsys("-".join(sorted(sub))))
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--lmdb", required=True)
    ap.add_argument("--relaxed", help="relaxed.json from mattersim_screen.py relax")
    ap.add_argument("--extxyz")
    ap.add_argument("--out")
    ap.add_argument("--export-ref", nargs=2, metavar=("CHEMSYS", "OUT_EXTXYZ"),
                    help="only write all reference structures of CHEMSYS (e.g. H-Ir-Mg) and its "
                         "subsystems to OUT_EXTXYZ (info: ref_id, dft_energy_corrected_eV), then exit; "
                         "used to build a MatterSim-only hull with ehull_self.py")
    a = ap.parse_args()
    ref = RefLMDB(Path(a.lmdb))
    if a.export_ref:
        from ase.io import write
        ents = ref_entries(ref, a.export_ref[0].split("-"))
        frames = []
        for i, e in enumerate(ents):
            at = AseAtomsAdaptor.get_atoms(e.structure)
            at.info = {"ref_id": str(e.entry_id), "ref_formula": e.composition.reduced_formula,
                       "dft_energy_corrected_eV": float(e.energy)}
            frames.append(at)
        write(a.export_ref[1], frames)
        print(f"wrote {len(frames)} reference structures, max natoms {max(len(f) for f in frames)}")
        return
    data = json.load(open(a.relaxed))
    if data["meta"].get("pressure_GPa", 0.0) != 0.0:
        raise SystemExit("reference hull is a 0 GPa DFT hull; refusing to compare enthalpies at P>0")
    frames = read(a.extxyz, ":")
    compat = MaterialsProject2020Compatibility()
    cache, rows = {}, []
    for row, at in zip(data["rows"], frames):
        s = AseAtomsAdaptor.get_structure(at)
        els = tuple(sorted({e.symbol for e in s.composition.elements}))
        if els not in cache:
            ents = ref_entries(ref, els)
            cache[els] = (PhaseDiagram(ents) if ents and all(
                any(e.composition.is_element and e.composition.elements[0].symbol == x for e in ents)
                for x in els) else None, len(ents))
        pd, nref = cache[els]
        entry = VasprunLike(s, row["energy_eV"]).get_computed_entry(
            inc_structure=False, energy_correction_scheme=compat)
        r = dict(index=row["index"], formula=s.composition.reduced_formula, chemsys="-".join(els),
                 n_ref_entries=nref, mp2020_correction_eV_per_atom=entry.correction / s.num_sites)
        if pd is None:
            r["e_above_hull_eV_per_atom"] = None
            r["note"] = "no complete reference hull for this chemical system"
        else:
            decomp, eh = pd.get_decomp_and_e_above_hull(entry, allow_negative=True)
            r["e_above_hull_eV_per_atom"] = round(float(eh), 4)
            r["decomposition"] = {k.composition.reduced_formula: round(v, 3) for k, v in decomp.items()}
        rows.append(r)
    stable = [r for r in rows if r["e_above_hull_eV_per_atom"] is not None]
    summ = dict(n=len(rows), n_with_hull=len(stable),
                n_below_0=sum(r["e_above_hull_eV_per_atom"] <= 0.0 for r in stable),
                n_below_0p1=sum(r["e_above_hull_eV_per_atom"] <= 0.1 for r in stable),
                ref_entries_per_chemsys={"-".join(k): v[1] for k, v in cache.items()})
    json.dump(dict(summary=summ, rows=rows), open(a.out, "w"), indent=1)
    print(json.dumps(summ, indent=1))
    for r in sorted(stable, key=lambda r: r["e_above_hull_eV_per_atom"])[:15]:
        print(r["index"], r["formula"], r["e_above_hull_eV_per_atom"], r.get("decomposition"))


if __name__ == "__main__":
    main()
