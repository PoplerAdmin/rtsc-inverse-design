"""MatterSim GPU screen for MatterGen output: batched relaxation, energies, finite-displacement phonons.

Runs in the main venv (~/.rtt/envs/py: mattersim 1.2.5, torch 2.14+cu130, phonopy 4.8.6).
HEURISTIC ONLY: MatterSim numbers may rank or defer candidates; they never decide anything.

Subcommands
-----------
relax    : relax all structures of one or more extxyz files at pressure P (BatchRelaxer, FIRE +
           FrechetCellFilter), single-point the final frames, write relaxed.extxyz + relaxed.json.
phonons  : phonopy finite displacements with MatterSim forces for selected relaxed structures.

Examples
--------
  export OMP_NUM_THREADS=2 OPENBLAS_NUM_THREADS=2 MKL_NUM_THREADS=2 TORCH_DISABLE_NATIVE_JIT=1
  PY=~/.rtt/envs/py/bin/python
  nice -n 19 $PY src/ml/mattersim_screen.py relax --out W/relax_LiMgH W/gen_LiMgH64/generated_crystals.extxyz
  nice -n 19 $PY src/ml/mattersim_screen.py phonons --relaxed W/relax_LiMgH/relaxed.extxyz --indices 3 7 --out W/ph_LiMgH
"""
from __future__ import annotations

import argparse
import json
import math
import os
import time
from pathlib import Path

# torch 2.14 dispatches the bmm in MatterSim's force backward to a Triton kernel that must be
# compiled with a C compiler (absent in this WSL) -> RuntimeError on CUDA. Must be set before
# torch is imported (it is imported lazily below).
os.environ.setdefault("TORCH_DISABLE_NATIVE_JIT", "1")

import numpy as np

EV_A3_TO_GPA = 160.21766208
THZ_TO_CM1 = 33.35640951981521
CM1_TO_K = 1.4387768775039338


def _gpu_peak_mib():
    import torch
    return round(torch.cuda.max_memory_allocated() / 2**20, 1) if torch.cuda.is_available() else None


def load_potential(model: str, device: str):
    from mattersim.forcefield.potential import Potential
    return Potential.from_checkpoint(load_path=model, device=device, load_training_state=False)


def cmd_relax(a):
    import torch
    from ase.io import read, write
    from ase.filters import FrechetCellFilter
    from mattersim.applications.batch_relax import BatchRelaxer
    from mattersim.forcefield.potential import MatterSimCalculator

    out = Path(a.out); out.mkdir(parents=True, exist_ok=True)
    atoms_list, src = [], []
    for f in a.inputs:
        for i, at in enumerate(read(f, ":")):
            atoms_list.append(at); src.append(f"{f}#{i}")
    t0 = time.time()
    pot = load_potential(a.model, a.device)
    t_load = time.time() - t0
    if a.pressure_GPa != 0.0:
        # FrechetCellFilter takes scalar_pressure in eV/A^3; BatchRelaxer instantiates the filter
        # with defaults only, so wrap it.
        p = a.pressure_GPa / EV_A3_TO_GPA

        class PFilter(FrechetCellFilter):
            def __init__(self, atoms):
                super().__init__(atoms, scalar_pressure=p)
        filt = PFilter
    else:
        filt = "FRECHETCELLFILTER"
    relaxer = BatchRelaxer(pot, optimizer="FIRE", filter=filt, fmax=a.fmax,
                           max_natoms_per_batch=a.max_natoms_per_batch, max_n_steps=a.max_steps)
    torch.cuda.reset_peak_memory_stats() if torch.cuda.is_available() else None
    t1 = time.time()
    traj = relaxer.relax(atoms_list)
    t_relax = time.time() - t1
    finals = [traj[i][-1].copy() for i in range(len(atoms_list))]
    nsteps = [len(traj[i]) for i in range(len(atoms_list))]
    # single-point re-evaluation of the final frames with the ASE calculator
    calc = MatterSimCalculator(potential=pot)
    rows = []
    t2 = time.time()
    for i, at in enumerate(finals):
        at.info.clear(); at.calc = calc
        e = float(at.get_potential_energy())
        fmax = float(np.abs(at.get_forces()).max())
        s = at.get_stress(voigt=True) * EV_A3_TO_GPA  # GPa (ASE sign: positive = tensile)
        p_int = -float(np.mean(s[:3]))
        n = len(at)
        rows.append(dict(index=i, source=src[i], formula=at.get_chemical_formula(mode="metal"),
                         n_atoms=n, energy_eV=e, energy_per_atom_eV=e / n,
                         enthalpy_per_atom_eV=(e + a.pressure_GPa / EV_A3_TO_GPA * at.get_volume()) / n,
                         volume_per_atom_A3=at.get_volume() / n, max_force_eV_A=fmax,
                         internal_pressure_GPa=p_int, n_steps=nsteps[i],
                         converged=bool(nsteps[i] < a.max_steps)))
        at.info["mattersim_energy_eV"] = e
        at.info["index"] = i
        at.calc = None
    t_sp = time.time() - t2
    write(out / "relaxed.extxyz", finals)
    meta = dict(model=a.model, device=a.device, pressure_GPa=a.pressure_GPa, fmax=a.fmax,
                max_steps=a.max_steps, max_natoms_per_batch=a.max_natoms_per_batch,
                n_structures=len(atoms_list), n_atoms_total=int(sum(len(x) for x in atoms_list)),
                t_model_load_s=round(t_load, 2), t_batch_relax_s=round(t_relax, 2),
                s_per_structure_relax=round(t_relax / max(1, len(atoms_list)), 3),
                t_single_points_s=round(t_sp, 2), total_ionic_steps=int(sum(nsteps)),
                n_converged=int(sum(r["converged"] for r in rows)), gpu_peak_alloc_MiB=_gpu_peak_mib())
    json.dump(dict(meta=meta, rows=rows), open(out / "relaxed.json", "w"), indent=1)
    print(json.dumps(meta, indent=1))


def _auto_supercell(cell, min_len, max_atoms, nat):
    lengths = np.linalg.norm(cell, axis=1)
    # use perpendicular widths, not lengths, for skewed cells
    vol = abs(np.linalg.det(cell))
    widths = [vol / np.linalg.norm(np.cross(cell[(i + 1) % 3], cell[(i + 2) % 3])) for i in range(3)]
    rep = [max(1, math.ceil(min_len / w)) for w in widths]
    while nat * np.prod(rep) > max_atoms and max(rep) > 1:
        j = int(np.argmax([r * w for r, w in zip(rep, widths)]))  # shrink largest dimension
        rep[j] -= 1
        rep = [max(1, r) for r in rep]
    return np.diag(rep), lengths


def cmd_phonons(a):
    import torch
    from ase.io import read
    from mattersim.forcefield.potential import MatterSimCalculator
    from phonopy import Phonopy
    from phonopy.structure.atoms import PhonopyAtoms

    out = Path(a.out); out.mkdir(parents=True, exist_ok=True)
    frames = read(a.relaxed, ":")
    pot = load_potential(a.model, a.device)
    calc = MatterSimCalculator(potential=pot)
    results = []
    for idx in a.indices:
        at = frames[idx]
        t0 = time.time()
        sc, _ = _auto_supercell(at.cell.array, a.min_supercell_length, a.max_supercell_atoms, len(at))
        unit = PhonopyAtoms(symbols=at.get_chemical_symbols(), cell=at.cell.array,
                            scaled_positions=at.get_scaled_positions())
        ph = Phonopy(unit, supercell_matrix=sc, primitive_matrix="auto")
        ph.generate_displacements(distance=a.displacement)
        scs = ph.supercells_with_displacements
        forces = []
        from ase import Atoms
        t1 = time.time()
        for s in scs:
            sa = Atoms(s.symbols, cell=s.cell, scaled_positions=s.scaled_positions, pbc=True)
            sa.calc = calc
            f = sa.get_forces()
            forces.append(f - f.mean(axis=0))  # remove drift
        t_forces = time.time() - t1
        ph.forces = np.array(forces)
        ph.produce_force_constants(); ph.symmetrize_force_constants()
        ph.run_mesh([a.mesh] * 3, is_gamma_center=True, with_eigenvectors=True, is_mesh_symmetry=False)
        mesh = ph.mesh
        freqs = np.array(mesh.frequencies) * THZ_TO_CM1           # (nq, 3N)
        w = np.array(mesh.weights, float); w /= w.sum()
        fmin = float(freqs.min()); qmin = mesh.qpoints[int(np.unravel_index(freqs.argmin(), freqs.shape)[0])]
        imag_frac = float((w[:, None] * (freqs < -a.imag_tol_cm1)).sum() / freqs.shape[1])
        # H weight of modes above 1000 cm^-1 (|eigvec|^2 on H sites; phonopy eigvecs are
        # mass-weighted-normalised, so this is a heuristic participation, not a DOS projection)
        ev = np.array(mesh.eigenvectors)                         # (nq, 3N, 3N): [q, comp, band]
        prim = ph.primitive
        syms = prim.symbols
        H = [i for i, s in enumerate(syms) if s == "H"]
        hfrac = None
        pos = freqs > 1.0
        lnw = np.where(pos, np.log(np.clip(freqs, 1e-9, None)), 0.0)
        # log-average of phonon DOS (not alpha^2F) over positive modes, weight 1/omega like omega_log
        wq = w[:, None] * pos / np.where(pos, freqs, 1.0)
        omega_log_dos_K = float(math.exp((wq * lnw).sum() / wq.sum()) * CM1_TO_K) if wq.sum() > 0 else None
        if H:
            comps = np.abs(ev) ** 2                               # (nq, 3N, 3N)
            hidx = [3 * i + k for i in H for k in range(3)]
            hw = comps[:, hidx, :].sum(axis=1)                    # (nq, bands)
            hfrac = float((w[:, None] * hw * (freqs > 1000.0)).sum() / (w[:, None] * hw).sum())
        r = dict(index=idx, formula=at.get_chemical_formula(mode="metal"), n_atoms=len(at),
                 supercell=sc.diagonal().tolist(), n_supercell_atoms=int(len(at) * np.prod(sc.diagonal())),
                 n_displacements=len(scs), t_forces_s=round(t_forces, 2),
                 t_total_s=round(time.time() - t0, 2), min_freq_cm1=round(fmin, 1),
                 min_freq_qpoint=np.round(qmin, 3).tolist(), imag_mode_fraction=round(imag_frac, 4),
                 max_freq_cm1=round(float(freqs.max()), 1),
                 H_weight_fraction_above_1000cm1=hfrac, omega_log_phDOS_K=omega_log_dos_K,
                 dynamically_stable_ml=bool(fmin > -a.imag_tol_cm1), mesh=a.mesh,
                 gpu_peak_alloc_MiB=_gpu_peak_mib())
        print(json.dumps(r))
        results.append(r)
    json.dump(dict(model=a.model, displacement_A=a.displacement, imag_tol_cm1=a.imag_tol_cm1,
                   results=results), open(out / "phonons.json", "w"), indent=1)


def main():
    ap = argparse.ArgumentParser()
    sub = ap.add_subparsers(dest="cmd", required=True)
    r = sub.add_parser("relax")
    r.add_argument("inputs", nargs="+")
    r.add_argument("--out", required=True)
    r.add_argument("--model", default="MatterSim-v1.0.0-5M")
    r.add_argument("--device", default="cuda")
    r.add_argument("--pressure-GPa", dest="pressure_GPa", type=float, default=0.0)
    r.add_argument("--fmax", type=float, default=0.02)
    r.add_argument("--max-steps", dest="max_steps", type=int, default=1000)
    r.add_argument("--max-natoms-per-batch", dest="max_natoms_per_batch", type=int, default=2048)
    p = sub.add_parser("phonons")
    p.add_argument("--relaxed", required=True)
    p.add_argument("--indices", type=int, nargs="+", required=True)
    p.add_argument("--out", required=True)
    p.add_argument("--model", default="MatterSim-v1.0.0-5M")
    p.add_argument("--device", default="cuda")
    p.add_argument("--displacement", type=float, default=0.01)
    p.add_argument("--min-supercell-length", dest="min_supercell_length", type=float, default=9.0)
    p.add_argument("--max-supercell-atoms", dest="max_supercell_atoms", type=int, default=250)
    p.add_argument("--mesh", type=int, default=16)
    p.add_argument("--imag-tol-cm1", dest="imag_tol_cm1", type=float, default=15.0)
    a = ap.parse_args()
    {"relax": cmd_relax, "phonons": cmd_phonons}[a.cmd](a)


if __name__ == "__main__":
    main()
