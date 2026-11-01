"""S3(a): batched MatterSim relaxation (architecture §3.1 S3; §5.1 ``mlff/relax.py``).

Ported from ``src/ml/mattersim_screen.py relax`` (doc 10 §3.1, §3.3): ``BatchRelaxer`` with FIRE
+ FrechetCellFilter, fmax 0.02 eV/Å, then a single-point re-evaluation of every final frame with
the ASE calculator. 0.6-0.7 s/structure on the GPU at batch 64 (doc 10, meas.). Sequential
single-structure relaxation is ~25x slower; always batch.

``TORCH_DISABLE_NATIVE_JIT=1`` is set before torch is imported (doc 10 §3.1: without it every
CUDA force call fails in the main venv). Heuristic only: nothing here decides anything.
"""
from __future__ import annotations

import os

os.environ.setdefault("TORCH_DISABLE_NATIVE_JIT", "1")

import time  # noqa: E402
from dataclasses import dataclass, field  # noqa: E402
from typing import Any, Sequence  # noqa: E402

import numpy as np  # noqa: E402

EV_A3_TO_GPA = 160.21766208
DEFAULT_MODEL = "MatterSim-v1.0.0-5M"

_POT_CACHE: dict[tuple[str, str], Any] = {}


@dataclass
class RelaxResult:
    index: int
    structure: dict                    # relaxed sdict (lattice, species, frac_coords)
    formula: str
    n_atoms: int
    energy_eV: float
    energy_per_atom_eV: float
    enthalpy_per_atom_eV: float
    volume_per_atom_A3: float
    max_force_eV_A: float
    internal_pressure_GPa: float
    n_steps: int
    converged: bool
    pressure_GPa: float
    model: str
    error: str | None = None
    timing: dict = field(default_factory=dict)

    def to_json(self) -> dict:
        return dict(self.__dict__)


def resolve_device(device: str | None = None) -> str:
    if device:
        return device
    try:
        import torch
        return "cuda" if torch.cuda.is_available() else "cpu"
    except Exception:  # noqa: BLE001
        return "cpu"


def load_potential(model: str = DEFAULT_MODEL, device: str | None = None):
    from mattersim.forcefield.potential import Potential
    dev = resolve_device(device)
    key = (model, dev)
    if key not in _POT_CACHE:
        if dev == "cpu":
            try:
                import torch
                torch.set_num_threads(2)
            except Exception:  # noqa: BLE001
                pass
        _POT_CACHE[key] = Potential.from_checkpoint(load_path=model, device=dev, load_training_state=False)
    return _POT_CACHE[key]


def to_atoms(s: Any):
    from ase import Atoms
    if hasattr(s, "get_chemical_symbols"):
        return s.copy()
    d = s.to_json() if hasattr(s, "to_json") else s
    if hasattr(d, "lattice") and hasattr(d, "sites"):          # pymatgen Structure
        from pymatgen.io.ase import AseAtomsAdaptor
        return AseAtomsAdaptor.get_atoms(d)
    d = d.get("structure", d) if isinstance(d, dict) and "structure" in d and "lattice" not in d else d
    return Atoms(symbols=list(d["species"]), cell=np.asarray(d["lattice"], float),
                 scaled_positions=np.asarray(d["frac_coords"], float), pbc=True)


def atoms_to_sdict(at) -> dict:
    return {"lattice": [[float(x) for x in r] for r in at.cell.array],
            "species": list(at.get_chemical_symbols()),
            "frac_coords": [[float(x) % 1.0 for x in r] for r in at.get_scaled_positions(wrap=True)]}


def relax_batch(structs: Sequence[Any], model: str = DEFAULT_MODEL, fmax: float = 0.02, *,
                pressure_GPa: float | Sequence[float] = 0.0, device: str | None = None,
                max_steps: int = 1000, max_natoms_per_batch: int = 2048) -> list[RelaxResult]:
    """Relax ``structs`` (sdicts / StructureData / pymatgen / ase) at ``pressure_GPa``.

    ``pressure_GPa`` may be a list (one per structure); structures are then grouped by pressure
    (BatchRelaxer takes one cell filter per call). Enthalpy per atom H = (E + PV)/N.
    """
    if isinstance(pressure_GPa, (int, float)):
        pressures = [float(pressure_GPa)] * len(structs)
    else:
        pressures = [float(p) for p in pressure_GPa]
        if len(pressures) != len(structs):
            raise ValueError("one pressure per structure")
    results: list[RelaxResult | None] = [None] * len(structs)
    for p in sorted(set(pressures)):
        idx = [i for i, q in enumerate(pressures) if q == p]
        for i, r in zip(idx, _relax_group([structs[i] for i in idx], model, fmax, p, device, max_steps,
                                           max_natoms_per_batch)):
            r.index = i
            results[i] = r
    return [r for r in results if r is not None]


def _relax_group(structs, model, fmax, p, device, max_steps, max_natoms_per_batch) -> list[RelaxResult]:
    from ase.filters import FrechetCellFilter
    from mattersim.applications.batch_relax import BatchRelaxer
    from mattersim.forcefield.potential import MatterSimCalculator

    atoms_list = [to_atoms(s) for s in structs]
    t0 = time.time()
    pot = load_potential(model, device)
    t_load = time.time() - t0
    if p != 0.0:
        sp = p / EV_A3_TO_GPA

        class PFilter(FrechetCellFilter):
            def __init__(self, atoms):
                super().__init__(atoms, scalar_pressure=sp)
        filt: Any = PFilter
    else:
        filt = "FRECHETCELLFILTER"
    relaxer = BatchRelaxer(pot, optimizer="FIRE", filter=filt, fmax=fmax,
                           max_natoms_per_batch=max_natoms_per_batch, max_n_steps=max_steps)
    t1 = time.time()
    traj = relaxer.relax([a.copy() for a in atoms_list])
    t_relax = time.time() - t1
    calc = MatterSimCalculator(potential=pot)
    out = []
    for i in range(len(atoms_list)):
        frames = traj.get(i) if isinstance(traj, dict) else traj[i]
        at = frames[-1].copy()
        nsteps = len(frames)
        at.info.clear()
        at.calc = calc
        try:
            e = float(at.get_potential_energy())
            f = float(np.abs(at.get_forces()).max())
            s = at.get_stress(voigt=True) * EV_A3_TO_GPA
            n = len(at)
            out.append(RelaxResult(
                index=i, structure=atoms_to_sdict(at), formula=at.get_chemical_formula(mode="metal"),
                n_atoms=n, energy_eV=e, energy_per_atom_eV=e / n,
                enthalpy_per_atom_eV=(e + p / EV_A3_TO_GPA * at.get_volume()) / n,
                volume_per_atom_A3=at.get_volume() / n, max_force_eV_A=f,
                internal_pressure_GPa=-float(np.mean(s[:3])), n_steps=nsteps,
                converged=bool(nsteps < max_steps), pressure_GPa=p, model=model,
                timing={"batch_relax_s": round(t_relax, 2), "model_load_s": round(t_load, 2),
                        "n_in_batch": len(atoms_list)}))
        except Exception as ex:  # noqa: BLE001
            out.append(RelaxResult(index=i, structure=atoms_to_sdict(at), formula=at.get_chemical_formula(),
                                   n_atoms=len(at), energy_eV=float("nan"), energy_per_atom_eV=float("nan"),
                                   enthalpy_per_atom_eV=float("nan"), volume_per_atom_A3=float("nan"),
                                   max_force_eV_A=float("nan"), internal_pressure_GPa=float("nan"),
                                   n_steps=nsteps, converged=False, pressure_GPa=p, model=model,
                                   error=f"{type(ex).__name__}: {ex}"))
        at.calc = None
    release_gpu_cache()
    return out


def release_gpu_cache() -> None:
    """Give cached CUDA blocks back after each batch (the GPU is shared with other agents)."""
    try:
        import torch
        if torch.cuda.is_available():
            torch.cuda.empty_cache()
    except Exception:  # noqa: BLE001
        pass


__all__ = ["RelaxResult", "relax_batch", "load_potential", "to_atoms", "atoms_to_sdict", "resolve_device",
           "DEFAULT_MODEL"]
