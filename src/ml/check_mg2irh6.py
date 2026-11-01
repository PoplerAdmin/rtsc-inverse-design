"""MatterSim cross-check on Fm-3m Mg2IrH6 (Dolui et al., PRL 132, 166001 (2024), arXiv:2310.07562:
PBE a = 6.655 A, harmonic-stable, lambda ~ 2.5, Tc 160-175 K).

1. MatterSim full relax (cell + positions) of the DFT reference structure from the Alex-MP set.
2. Phonons at the MatterSim lattice and at the PBE lattice (positions relaxed, cell fixed).
Writes JSON to stdout. Runs in ~/.rtt/envs/py on the GPU.
  $PY src/ml/check_mg2irh6.py W/ref/ref_H-Ir-Mg.extxyz
"""
import json
import os
import sys

os.environ.setdefault("TORCH_DISABLE_NATIVE_JIT", "1")
import numpy as np
from ase.filters import FrechetCellFilter
from ase.io import read
from ase.optimize import FIRE
from mattersim.forcefield.potential import MatterSimCalculator
from phonopy import Phonopy
from phonopy.structure.atoms import PhonopyAtoms
from pymatgen.io.ase import AseAtomsAdaptor
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer

THZ_TO_CM1 = 33.35640951981521
calc = MatterSimCalculator(load_path="MatterSim-v1.0.0-5M", device="cuda")


def conv_a(at):
    s = SpacegroupAnalyzer(AseAtomsAdaptor.get_structure(at), symprec=0.05).get_conventional_standard_structure()
    return round(float(s.lattice.a), 4), SpacegroupAnalyzer(s, symprec=0.05).get_space_group_symbol()


def phon(at, n=3, mesh=16):
    ph = Phonopy(PhonopyAtoms(symbols=at.get_chemical_symbols(), cell=at.cell.array,
                              scaled_positions=at.get_scaled_positions()),
                 supercell_matrix=np.eye(3, dtype=int) * n, primitive_matrix="auto")
    ph.generate_displacements(distance=0.01)
    from ase import Atoms
    F = []
    for s in ph.supercells_with_displacements:
        a = Atoms(s.symbols, cell=s.cell, scaled_positions=s.scaled_positions, pbc=True); a.calc = calc
        f = a.get_forces(); F.append(f - f.mean(0))
    ph.forces = np.array(F); ph.produce_force_constants(); ph.symmetrize_force_constants()
    out = {}
    for lab, q in {"G": [0, 0, 0], "X": [0.5, 0, 0.5], "L": [0.5, 0.5, 0.5], "W": [0.5, 0.25, 0.75]}.items():
        out[lab] = np.round(np.sort(ph.get_frequencies(q)) * THZ_TO_CM1, 1)[:6].tolist()
    ph.run_mesh([mesh] * 3, is_gamma_center=True)
    fr = ph.mesh.frequencies * THZ_TO_CM1
    return dict(lowest6_cm1_at=out, min_mesh_cm1=round(float(fr.min()), 1), max_cm1=round(float(fr.max()), 1))


frames = read(sys.argv[1], ":")
ref = [f for f in frames if sorted(f.get_chemical_symbols()) == sorted(["Mg"] * 2 * (len(f) // 9) + ["Ir"] * (len(f) // 9) + ["H"] * 6 * (len(f) // 9))]
res = []
for f in ref:
    at = f.copy(); at.calc = calc
    a_dft, sg = conv_a(at)
    FIRE(FrechetCellFilter(at), logfile=None).run(fmax=0.005, steps=2000)
    a_ms, sg2 = conv_a(at)
    r = dict(ref_id=f.info["ref_id"], sg_ref=sg, a_ref_dft=a_dft, sg_mattersim=sg2, a_mattersim=a_ms,
             phonons_mattersim_lattice=phon(at))
    # PBE lattice of Dolui et al.: scale to a = 6.655 A, relax positions only
    b = at.copy(); b.calc = calc
    b.set_cell(at.cell.array * (6.655 / a_ms), scale_atoms=True)
    FIRE(b, logfile=None).run(fmax=0.002, steps=2000)
    r["a_pbe_fixed"] = conv_a(b)[0]
    r["pressure_at_pbe_lattice_GPa"] = round(float(-np.mean(b.get_stress()[:3]) * 160.21766208), 2)
    r["phonons_pbe_lattice"] = phon(b)
    res.append(r)
print(json.dumps(res, indent=1))
