"""MatterSim-5M (CPU) Gamma-point phonons of Fm-3m Mg2IrH6 at the QE/PBE geometry, with irreps,
to identify WHICH Gamma mode MatterSim makes imaginary. Gamma is exact in the 1x1x1 cell
(finite displacement, 0.01 A). Two geometries: (a) QE-relaxed a and x_H unchanged,
(b) a fixed, x_H relaxed with MatterSim.  Usage: python mattersim_gamma_irreps.py a_A x_H
"""
import json
import os
import sys

os.environ.setdefault("TORCH_DISABLE_NATIVE_JIT", "1")
import numpy as np
from ase import Atoms
from ase.constraints import FixSymmetry
from ase.optimize import BFGS
from mattersim.forcefield.potential import MatterSimCalculator
from phonopy import Phonopy
from phonopy.structure.atoms import PhonopyAtoms

THZ_TO_CM1 = 33.35640951981521
a, x = float(sys.argv[1]), float(sys.argv[2])
calc = MatterSimCalculator(load_path="MatterSim-v1.0.0-5M", device="cpu")
cell = a / 2 * np.array([[-1, 0, 1], [0, 1, 1], [-1, 1, 0]])
cart = np.array([[0, 0, 0], [.25, .25, .25], [-.25, -.25, -.25], [x, 0, 0], [-x, 0, 0],
                 [0, x, 0], [0, -x, 0], [0, 0, x], [0, 0, -x]]) * a
at0 = Atoms("IrMg2H6", positions=cart, cell=cell, pbc=True)


def gamma(at):
    ph = Phonopy(PhonopyAtoms(symbols=at.get_chemical_symbols(), cell=at.cell.array,
                              scaled_positions=at.get_scaled_positions()),
                 supercell_matrix=np.eye(3, dtype=int), primitive_matrix="P")
    ph.generate_displacements(distance=0.01)
    F = []
    for s in ph.supercells_with_displacements:
        b = Atoms(s.symbols, cell=s.cell, scaled_positions=s.scaled_positions, pbc=True)
        b.calc = calc
        f = b.get_forces()
        F.append(f - f.mean(0))
    ph.forces = np.array(F)
    ph.produce_force_constants()
    ph.symmetrize_force_constants()
    fr = np.sort(ph.get_frequencies([0, 0, 0])) * THZ_TO_CM1
    out = {"freqs_cm1": [round(float(v), 1) for v in fr], "irreps_table": None}
    try:
        import contextlib
        import io
        ph.run_irreps([0, 0, 0], degeneracy_tolerance=1e-3)
        buf = io.StringIO()
        with contextlib.redirect_stdout(buf):
            ph.irreps.show()
        out["irreps_table"] = buf.getvalue()
    except Exception as e:  # irreps are a convenience; frequencies are the result
        out["irreps_table"] = f"irreps failed: {e!r}"
    return out


at = at0.copy(); at.calc = calc
res = {"a_A": a, "x_H_in": x, "max_force_at_QE_geometry_eV_A": float(np.abs(at.get_forces()).max()),
       "pressure_at_QE_geometry_GPa": float(-np.mean(at.get_stress()[:3]) * 160.21766208),
       "gamma_QE_geometry": gamma(at)}
b = at0.copy(); b.calc = calc
b.set_constraint(FixSymmetry(b))
BFGS(b, logfile=None).run(fmax=1e-4, steps=500)
b.set_constraint()
res["x_H_mattersim_relaxed"] = float(b.positions[3, 0] / a)
res["gamma_mattersim_relaxed_x"] = gamma(b)
print(json.dumps(res, indent=1))
