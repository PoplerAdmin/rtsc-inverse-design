"""MatterSim-only (self-consistent) convex hull: every energy, reference and candidate, from the same
MatterSim model, so systematic MatterSim-vs-DFT offsets cancel to first order. No MP2020 corrections.

Runs in the main venv (~/.rtt/envs/py).  Inputs are relaxed.json/.extxyz pairs written by
`mattersim_screen.py relax`: one for the reference structures (exported with
`ehull_ref.py --export-ref`) and one for the candidates.  Works at any pressure (uses the enthalpy
per atom H = E + PV that mattersim_screen.py stores), provided both sets were relaxed at the same P.

  $PY src/ml/ehull_self.py --ref W/relax_ref_HIrMg/relaxed.json --cand W/relax_MgIrH64/relaxed.json \
      --out W/relax_MgIrH64/ehull_self.json
"""
from __future__ import annotations

import argparse
import json

from pymatgen.analysis.phase_diagram import PDEntry, PhaseDiagram
from pymatgen.core import Composition


def entries(path, tag):
    d = json.load(open(path))
    out = []
    for r in d["rows"]:
        c = Composition(r["formula"])
        e = PDEntry(c, r["enthalpy_per_atom_eV"] * c.num_atoms, name=f"{tag}:{r['index']}")
        out.append((r, e))
    return d["meta"], out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--ref", required=True)
    ap.add_argument("--cand", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--include-candidates-in-hull", action="store_true",
                    help="also let candidates define the hull (E_hull can then not be negative)")
    a = ap.parse_args()
    mref, ref = entries(a.ref, "ref")
    mcand, cand = entries(a.cand, "cand")
    if abs(mref["pressure_GPa"] - mcand["pressure_GPa"]) > 1e-9 or mref["model"] != mcand["model"]:
        raise SystemExit("reference and candidates relaxed at different P or with different models")
    hull_entries = [e for r, e in ref if r["converged"]]
    if a.include_candidates_in_hull:
        hull_entries += [e for r, e in cand if r["converged"]]
    pd = PhaseDiagram(hull_entries)
    rows = []
    for r, e in cand:
        try:
            decomp, eh = pd.get_decomp_and_e_above_hull(e, allow_negative=True)
        except ValueError as ex:  # element outside the reference system
            rows.append(dict(index=r["index"], formula=e.composition.reduced_formula, error=str(ex)))
            continue
        rows.append(dict(index=r["index"], formula=e.composition.reduced_formula,
                         converged=r["converged"], e_above_hull_self_eV_per_atom=round(float(eh), 4),
                         decomposition={k.composition.reduced_formula: round(v, 3) for k, v in decomp.items()}))
    ok = [x for x in rows if "e_above_hull_self_eV_per_atom" in x]
    summ = dict(model=mref["model"], pressure_GPa=mref["pressure_GPa"], n_ref_converged=len(hull_entries),
                n=len(rows), n_below_0=sum(x["e_above_hull_self_eV_per_atom"] <= 0 for x in ok),
                n_below_0p05=sum(x["e_above_hull_self_eV_per_atom"] <= 0.05 for x in ok),
                n_below_0p1=sum(x["e_above_hull_self_eV_per_atom"] <= 0.1 for x in ok))
    json.dump(dict(summary=summ, rows=rows), open(a.out, "w"), indent=1)
    print(json.dumps(summ, indent=1))
    for x in sorted(ok, key=lambda x: x["e_above_hull_self_eV_per_atom"])[:12]:
        print(x["index"], x["formula"], x["e_above_hull_self_eV_per_atom"], x["decomposition"])


if __name__ == "__main__":
    main()
