"""S3(d): MatterSim harmonic phonons -- **ranking only** (architecture §3.1 S3, §4 R2; doc 10 §3.4).

Ported from ``src/ml/mattersim_screen.py phonons``: phonopy finite displacements (0.01 Å),
supercell with perpendicular widths >= 9 Å and <= 250 atoms, Γ-centred mesh.

**Never a removal criterion.** MatterSim gives a false dynamical instability for Fm-3m Mg2IrH6
(doc 10 §3.4: Γ -255, X -287 cm⁻¹, while DFT/SSCHA find it stable). An imaginary MatterSim
mode only sets ``needs_phonon_check``; :func:`rtt.mlff.screen.decide_s3` cannot use it.

Summary quantities (``PhononSummary``):

* ``min_freq_cm1`` (negative = imaginary), ``imag_mode_fraction`` (q-weighted share of modes below
  -``imag_tol_cm1``), ``omega_max_cm1``;
* ``omega2_light_K``: sqrt of the second moment of the phonon spectrum, each mode weighted by its
  |eigenvector|² share on light atoms (H, Li, Be, B, C, N); positive modes only. A *phonon*-DOS
  moment, not an α²F moment: it is the ω₂^proxy ranking feature of §3.1 notes (λ_needed =
  2·floor/ω₂^proxy), nothing more;
* ``zpe_ml_meV_atom``: harmonic ZPE ½Σħω per atom over positive modes (ranking only; the H7 ZPE
  correction uses DFT at S7);
* ``omega_log_phDOS_K``: log-average of the phonon DOS with 1/ω weight (as in the research script).
"""
from __future__ import annotations

import os

os.environ.setdefault("TORCH_DISABLE_NATIVE_JIT", "1")

import math  # noqa: E402
import time  # noqa: E402
from dataclasses import dataclass, field  # noqa: E402
from typing import Any  # noqa: E402

import numpy as np  # noqa: E402

THZ_TO_CM1 = 33.35640951981521
CM1_TO_K = 1.4387768775039338
CM1_TO_MEV = 0.12398419843320026
LIGHT = ("H", "Li", "Be", "B", "C", "N")


@dataclass
class PhononSummary:
    min_freq_cm1: float
    omega_max_cm1: float
    omega2_light_K: float | None
    zpe_ml_meV_atom: float
    imag_mode_fraction: float
    omega_log_phDOS_K: float | None
    light_weight_fraction: float | None
    min_freq_qpoint: list = field(default_factory=list)
    supercell: list = field(default_factory=list)
    n_supercell_atoms: int = 0
    n_displacements: int = 0
    mesh: int = 0
    stable_ml: bool = True
    needs_phonon_check: bool = False
    wall_s: float = 0.0
    model: str = ""
    imag_tol_cm1: float = 15.0

    def to_json(self) -> dict:
        return dict(self.__dict__)


def auto_supercell(cell, min_len: float = 9.0, max_atoms: int = 250, nat: int = 1) -> list[int]:
    cell = np.asarray(cell, float)
    vol = abs(np.linalg.det(cell))
    widths = [vol / np.linalg.norm(np.cross(cell[(i + 1) % 3], cell[(i + 2) % 3])) for i in range(3)]
    rep = [max(1, math.ceil(min_len / w)) for w in widths]
    while nat * np.prod(rep) > max_atoms and max(rep) > 1:
        # shrink the longest *reducible* direction (campaign 1: argmax over all directions could
        # pick one already at 1 and loop forever, e.g. a long c axis with rep [2, 2, 1])
        j = int(np.argmax([r * w if r > 1 else -np.inf for r, w in zip(rep, widths)]))
        rep[j] -= 1
    return rep


def summarize(freqs_cm1: np.ndarray, weights: np.ndarray, light_w: np.ndarray | None, n_atoms: int,
              imag_tol_cm1: float = 15.0) -> dict:
    """Pure function over mesh data: freqs (nq, nb) cm⁻¹, q weights (nq,), light weight (nq, nb)."""
    f = np.asarray(freqs_cm1, float)
    w = np.asarray(weights, float)
    w = w / w.sum()
    W = np.broadcast_to(w[:, None], f.shape)
    pos = f > 1.0
    imag_frac = float((W * (f < -imag_tol_cm1)).sum() / f.shape[1])
    zpe = float(0.5 * (W * np.where(pos, f, 0.0)).sum() * CM1_TO_MEV / n_atoms)
    wq = W * pos / np.where(pos, f, 1.0)
    lnw = np.where(pos, np.log(np.clip(f, 1e-9, None)), 0.0)
    olog = float(math.exp((wq * lnw).sum() / wq.sum()) * CM1_TO_K) if wq.sum() > 0 else None
    o2, lfrac = None, None
    if light_w is not None:
        lw = np.asarray(light_w, float)
        den = (W * lw * pos).sum()
        if den > 0:
            o2 = float(math.sqrt((W * lw * pos * f ** 2).sum() / den) * CM1_TO_K)
        lfrac = float((W * lw).sum() / W.sum())
    fmin = float(f.min())
    return {"min_freq_cm1": round(fmin, 2), "omega_max_cm1": round(float(f.max()), 2),
            "omega2_light_K": None if o2 is None else round(o2, 2), "zpe_ml_meV_atom": round(zpe, 3),
            "imag_mode_fraction": round(imag_frac, 5), "omega_log_phDOS_K": None if olog is None else round(olog, 2),
            "light_weight_fraction": None if lfrac is None else round(lfrac, 4),
            "stable_ml": fmin > -imag_tol_cm1, "needs_phonon_check": fmin <= -imag_tol_cm1}


def ml_phonons(struct: Any, *, calculator: Any = None, model: str = "MatterSim-v1.0.0-5M",
               device: str | None = None, displacement: float = 0.01, min_supercell_length: float = 9.0,
               max_supercell_atoms: int = 250, mesh: int | None = None,
               imag_tol_cm1: float = 15.0) -> PhononSummary:
    from ase import Atoms
    from phonopy import Phonopy
    from phonopy.structure.atoms import PhonopyAtoms
    from .relax import load_potential, to_atoms

    t0 = time.time()
    at = to_atoms(struct)
    if calculator is None:
        from mattersim.forcefield.potential import MatterSimCalculator
        calculator = MatterSimCalculator(potential=load_potential(model, device))
    sc = auto_supercell(at.cell.array, min_supercell_length, max_supercell_atoms, len(at))
    unit = PhonopyAtoms(symbols=at.get_chemical_symbols(), cell=at.cell.array,
                        scaled_positions=at.get_scaled_positions())
    ph = Phonopy(unit, supercell_matrix=np.diag(sc), primitive_matrix="auto")
    ph.generate_displacements(distance=displacement)
    scs = ph.supercells_with_displacements
    forces = []
    for s in scs:
        sa = Atoms(s.symbols, cell=s.cell, scaled_positions=s.scaled_positions, pbc=True)
        sa.calc = calculator
        f = sa.get_forces()
        forces.append(f - f.mean(axis=0))
    ph.forces = np.array(forces)
    ph.produce_force_constants()
    ph.symmetrize_force_constants()
    nprim = len(ph.primitive.symbols)
    m = int(mesh or (16 if nprim <= 10 else 10))
    ph.run_mesh([m] * 3, is_gamma_center=True, with_eigenvectors=True, is_mesh_symmetry=False)
    md = ph.mesh
    freqs = np.array(md.frequencies) * THZ_TO_CM1
    ev = np.array(md.eigenvectors)                    # (nq, 3N, 3N) [q, component, band]
    syms = list(ph.primitive.symbols)
    lidx = [3 * i + k for i, s in enumerate(syms) if s in LIGHT for k in range(3)]
    light_w = (np.abs(ev[:, lidx, :]) ** 2).sum(axis=1) if lidx else None
    s = summarize(freqs, np.array(md.weights, float), light_w, nprim, imag_tol_cm1)
    qmin = md.qpoints[int(np.unravel_index(freqs.argmin(), freqs.shape)[0])]
    from .relax import release_gpu_cache
    release_gpu_cache()
    return PhononSummary(**s, min_freq_qpoint=np.round(qmin, 3).tolist(), supercell=list(map(int, sc)),
                         n_supercell_atoms=int(len(at) * int(np.prod(sc))), n_displacements=len(scs), mesh=m,
                         wall_s=round(time.time() - t0, 2), model=model, imag_tol_cm1=imag_tol_cm1)


__all__ = ["PhononSummary", "ml_phonons", "summarize", "auto_supercell"]
