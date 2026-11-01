"""Own-DFT reference hull (architecture S7a; also the elemental references of the S4 heuristic).

Every phase is relaxed and evaluated with the *same* profile (cutoff, k density, smearing, SG15
pseudopotentials) as the candidate it is compared with; energies from different profiles are
never mixed (the cache key contains the profile hash and the pseudopotential md5s).

* :func:`phase_energy` runs relax (``vc-relax`` for crystals, fixed-cell ``relax`` at Gamma for
  molecules in a box, e.g. :func:`h2_molecule`) + scf through :class:`rtt.dft.workflow.Chain`,
  and caches the result as JSON under ``cache_dir`` by
  ``sha256(structure | pseudo md5s | profile sha256)``.
* :func:`ref_hull` evaluates a set of reference phases for one chemical system and returns a
  :class:`RefHull`; :meth:`RefHull.e_above_hull` uses pymatgen's ``PhaseDiagram`` with optional
  harmonic ZPE added per phase (H7: ZPE-corrected E_hull for H compounds).
* :func:`zpe` = harmonic ZPE (meV/atom) from a matdyn phonon DOS.

Choosing the reference phases (which competing compounds exist) is the caller's job (database
lane D / Alex-MP; architecture §3.1 S7a: alpha-B instead of beta-B).
"""
from __future__ import annotations

import hashlib
import json
from dataclasses import asdict, dataclass, field
from pathlib import Path
from typing import Any

from .coarse import zpe_meV_per_atom
from .profiles import PROFILES, Profile, pseudo_md5s, sg15_pseudos
from .workflow import Chain, ChainStep, RunPolicy, default_launcher

RY_EV = 13.605693122994


def zpe(phonon_dos: str | Path, nat: int) -> float | None:
    """Harmonic zero-point energy in meV/atom from a matdyn ``phonon.dos`` file."""
    return zpe_meV_per_atom(Path(phonon_dos), nat)


def h2_molecule(box_A: float = 10.0, bond_A: float = 0.7414) -> dict:
    """H2 in a cubic box (the H reference of formation energies); fixed-cell relax at Gamma."""
    b = bond_A / box_A
    return {"lattice": [[box_A, 0, 0], [0, box_A, 0], [0, 0, box_A]], "species": ["H", "H"],
            "frac_coords": [[0.5 - b / 2, 0.5, 0.5], [0.5 + b / 2, 0.5, 0.5]],
            "masses_amu": {"H": 1.00794}, "molecule": True}


@dataclass
class PhaseEnergy:
    name: str
    composition: dict
    nat: int
    e_total_Ry: float
    e_per_atom_eV: float
    profile: str
    cache_key: str
    zpe_meV_atom: float | None = None
    molecule: bool = False
    workdir: str | None = None
    wall_s: float = 0.0

    def to_dict(self) -> dict[str, Any]:
        return asdict(self)


@dataclass
class RefHull:
    chemsys: list[str]
    profile: str
    phases: list[PhaseEnergy] = field(default_factory=list)
    errors: dict = field(default_factory=dict)

    def element_refs(self) -> dict[str, float]:
        """Lowest energy per atom (eV) of the single-element phases."""
        out: dict[str, float] = {}
        for p in self.phases:
            if len(p.composition) == 1:
                el = next(iter(p.composition))
                out[el] = min(out.get(el, float("inf")), p.e_per_atom_eV)
        return out

    def e_above_hull(self, composition: dict[str, int], e_per_atom_eV: float,
                     zpe_meV_atom: float | None = None, with_zpe: bool = False) -> float:
        """E_hull (eV/atom) of a candidate against these phases (pymatgen PhaseDiagram).
        ``with_zpe`` adds every phase's ZPE (a phase without ZPE then raises)."""
        from pymatgen.analysis.phase_diagram import PDEntry, PhaseDiagram
        from pymatgen.core import Composition

        def e_tot(c: dict, epa: float, z: float | None, nat: int) -> float:
            if with_zpe:
                if z is None:
                    raise ValueError("ZPE missing for a phase while with_zpe=True")
                epa = epa + z / 1000.0
            return epa * nat

        entries = [PDEntry(Composition(p.composition), e_tot(p.composition, p.e_per_atom_eV,
                                                             p.zpe_meV_atom, p.nat), name=p.name)
                   for p in self.phases]
        n = sum(composition.values())
        cand = PDEntry(Composition(composition), e_tot(composition, e_per_atom_eV, zpe_meV_atom, n),
                       name="candidate")
        missing = set(composition) - set(self.element_refs())
        if missing:
            raise ValueError(f"no elemental reference for {sorted(missing)}")
        pd = PhaseDiagram(entries + [cand])
        return float(pd.get_e_above_hull(cand))

    def to_dict(self) -> dict[str, Any]:
        return {"chemsys": self.chemsys, "profile": self.profile,
                "phases": [p.to_dict() for p in self.phases], "errors": self.errors}


def _canon(cell: Any) -> str:
    d = cell.to_dict()
    return json.dumps({"lattice_A": [[round(x, 8) for x in r] for r in d["lattice_A"]],
                       "species": d["species"],
                       "frac_coords": [[round(x, 8) for x in r] for r in d["frac_coords"]]},
                      sort_keys=True)


def phase_energy(name: str, struct: Any, profile: Profile, workdir: str | Path, *,
                 launcher: Any = None, cache_dir: str | Path | None = None,
                 policy: RunPolicy | None = None, molecule: bool | None = None) -> PhaseEnergy:
    from rtt.core.dfpt import qe_inputs as qi
    from rtt.core.dfpt.grids import grid_from_spacing
    from rtt.core.dfpt.structure import CellData, parse_pw_final_structure, to_celldata
    from .electronic import parse_total_energy_Ry
    from .workflow import _masses_for

    is_mol = bool(molecule if molecule is not None else
                  (isinstance(struct, dict) and struct.get("molecule")))
    s = {k: v for k, v in struct.items() if k != "molecule"} if isinstance(struct, dict) else struct
    cell = to_celldata(s, _masses_for(s))
    pseudos = sg15_pseudos(cell.elements, pattern=profile.pseudo_pattern)
    key = hashlib.sha256("|".join([_canon(cell), json.dumps(pseudo_md5s(pseudos, cell.elements),
                                                            sort_keys=True),
                                   profile.sha256(), str(is_mol)]).encode()).hexdigest()
    cdir = Path(cache_dir).expanduser() if cache_dir else None
    if cdir is not None and (cdir / f"{key}.json").is_file():
        return PhaseEnergy(**json.loads((cdir / f"{key}.json").read_text()))
    w = Path(workdir).expanduser()
    launcher = launcher or default_launcher(2)
    common = dict(prefix="pwscf", outdir="./tmp/",
                  pseudo_dir=str(Path(pseudos.directory).resolve()) + "/",
                  pseudo_files={el: pseudos.files[el] for el in cell.elements},
                  ecutwfc=float(profile.ecutwfc_Ry), ecutrho=profile.ecutrho_Ry,
                  smearing=profile.smearing, degauss=float(profile.degauss_Ry),
                  conv_thr=profile.conv_thr, mixing_beta=profile.mixing_beta)

    def kg(c: Any) -> tuple[int, int, int]:
        return (1, 1, 1) if is_mol else grid_from_spacing(c.reciprocal_lengths_2pi(), profile.k_spacing)

    def relaxed() -> Any:
        return CellData.from_dict(json.loads((w / "relaxed_structure.json").read_text()))

    def after(wd: Path) -> None:
        c = parse_pw_final_structure((wd / "relax.out").read_text(errors="replace"), cell)
        (wd / "relaxed_structure.json").write_text(json.dumps(c.to_dict(), indent=1))

    def chk(name: str, needle: str, extra: tuple[str, ...] = ()):
        def f(wd: Path) -> tuple[bool, str]:
            o = wd / f"{name}.out"
            if not o.exists():
                return False, "no output"
            t = o.read_text(errors="replace")
            if "JOB DONE" not in t or needle not in t:
                return False, f"{needle!r} / JOB DONE missing"
            for e in extra:
                if not (wd / e).exists():
                    return False, f"{e} missing"
            return True, ""
        return f

    calc = "relax" if is_mol else "vc-relax"
    steps = [
        ChainStep("relax", "pw.x",
                  lambda: qi.write_pw_input(calculation=calc, cell=cell, kgrid=kg(cell),
                                            press_kbar=0.0,
                                            press_conv_thr_kbar=profile.press_conv_thr_kbar,
                                            etot_conv_thr=profile.etot_conv_thr,
                                            forc_conv_thr=profile.forc_conv_thr,
                                            nstep=profile.relax_nstep,
                                            cell_dofree=profile.cell_dofree, **common),
                  chk("relax", "End final coordinates", ("relaxed_structure.json",)), after=after),
        ChainStep("scf", "pw.x",
                  lambda: qi.write_pw_input(calculation="scf", cell=relaxed(), kgrid=kg(relaxed()),
                                            **common),
                  chk("scf", "convergence has been achieved"), upstream=("relax",)),
    ]
    chain = Chain(w, steps, launcher, material_id=name, profile=profile.name, policy=policy)
    ok, err = chain.run()
    if not ok:
        raise RuntimeError(f"reference phase {name}: {err}")
    e = parse_total_energy_Ry((w / "scf.out").read_text(errors="replace"))
    if e is None:
        raise RuntimeError(f"reference phase {name}: no total energy in scf.out")
    comp: dict[str, int] = {}
    for sp in cell.species:
        comp[sp] = comp.get(sp, 0) + 1
    pe = PhaseEnergy(name=name, composition=comp, nat=cell.nat, e_total_Ry=e,
                     e_per_atom_eV=e * RY_EV / cell.nat, profile=profile.name, cache_key=key,
                     molecule=is_mol, workdir=str(w),
                     wall_s=float(sum(t.get("wall_s", 0.0) for t in chain.timings)))
    if cdir is not None:
        cdir.mkdir(parents=True, exist_ok=True)
        (cdir / f"{key}.json").write_text(json.dumps(pe.to_dict(), indent=1))
    return pe


def ref_hull(chemsys: list[str] | str, profile: Profile | None, phases: dict[str, Any],
             workroot: str | Path, *, launcher: Any = None, cache_dir: str | Path | None = None,
             policy: RunPolicy | None = None, zpe_by_phase: dict[str, float] | None = None
             ) -> RefHull:
    """Own-DFT reference energies for ``chemsys`` (e.g. "Mg-Ir-H"). ``phases`` maps a name to
    a structure (core-CellData-compatible; dicts may set ``"molecule": True``). Failed
    phases are recorded in ``RefHull.errors`` and left out. Fewer phases can only raise the hull
    surface, so the computed E_hull is then a lower bound on the true one (a removal for
    E_hull > window stays conservative), but the evidence statement must list the missing
    phases."""
    prof = profile or PROFILES["ref_hull"]
    els = chemsys.split("-") if isinstance(chemsys, str) else list(chemsys)
    rh = RefHull(chemsys=sorted(els), profile=prof.name)
    root = Path(workroot).expanduser()
    for name, st in phases.items():
        try:
            pe = phase_energy(name, st, prof, root / name, launcher=launcher, cache_dir=cache_dir,
                              policy=policy)
            if zpe_by_phase and name in zpe_by_phase:
                pe.zpe_meV_atom = zpe_by_phase[name]
            rh.phases.append(pe)
        except Exception as exc:
            rh.errors[name] = str(exc)
    return rh
