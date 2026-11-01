"""Input-file writers for pw.x (vc-relax, scf), ph.x (interpolated e-ph), q2r.x, matdyn.x.

Templates follow QE's own electron-phonon example (``PHonon/examples/example03/run_example``,
tag ``qe-6.7.0``, fetched from GitHub): a dense-k scf with ``la2F=.true.`` (writes
``<prefix>.a2Fsave``), then the coarse scf, then ``ph.x`` with ``electron_phonon='interpolated'``,
``el_ph_sigma``, ``el_ph_nsigma``, ``trans=.true.``, ``ldisp=.true.``, ``nq1..3``; ``q2r.x``
with ``la2F=.true.``; ``matdyn.x`` with ``dos=.true., la2F=.true.``, ``nk1..3``, ``ndos``.
Namelist variable names were checked against the NAMELIST statements in ``q2r.f90`` and
``matdyn.f90`` (both define ``el_ph_nsigma``; matdyn ``DeltaE`` is in cm^-1).
"""
from __future__ import annotations

from typing import Any, Iterable

from .structure import CellData


def fortran_value(v: Any) -> str:
    if isinstance(v, bool):
        return ".true." if v else ".false."
    if isinstance(v, int):
        return str(v)
    if isinstance(v, float):
        r = repr(v)
        if "e" in r or "E" in r:
            m, e = r.lower().split("e")
            if "." not in m:
                m += ".0"
            return f"{m}d{int(e)}"
        return r
    if isinstance(v, str):
        if "'" in v:
            raise ValueError(f"quote in namelist string {v!r}")
        return f"'{v}'"
    raise TypeError(f"unsupported namelist value {v!r}")


def namelist(name: str, values: dict[str, Any]) -> str:
    lines = [f" &{name}"]
    for k, v in values.items():
        if v is None:
            continue
        lines.append(f"    {k} = {fortran_value(v)}")
    lines.append(" /")
    return "\n".join(lines) + "\n"


def _species_card(cell: CellData, pseudo_files: dict[str, str]) -> str:
    out = ["ATOMIC_SPECIES"]
    for el in cell.elements:
        out.append(f" {el:<3s} {cell.masses_amu[el]:.8f} {pseudo_files[el]}")
    return "\n".join(out) + "\n"


def _cell_cards(cell: CellData) -> str:
    out = ["CELL_PARAMETERS angstrom"]
    for row in cell.lattice_A:
        out.append("  " + " ".join(f"{x:20.14f}" for x in row))
    out.append("ATOMIC_POSITIONS crystal")
    for s, p in zip(cell.species, cell.frac_coords):
        out.append(f" {s:<3s} " + " ".join(f"{x:20.14f}" for x in p))
    return "\n".join(out) + "\n"


def _kpoints(grid: Iterable[int]) -> str:
    g = list(grid)
    return f"K_POINTS automatic\n {g[0]} {g[1]} {g[2]} 0 0 0\n"


def write_pw_input(*, calculation: str, cell: CellData, prefix: str, outdir: str,
                   pseudo_dir: str, pseudo_files: dict[str, str], kgrid: Iterable[int],
                   ecutwfc: float, ecutrho: float | None = None,
                   smearing: str = "mv", degauss: float = 0.02, la2F: bool = False,
                   conv_thr: float = 1e-10, mixing_beta: float = 0.5,
                   press_kbar: float | None = None, press_conv_thr_kbar: float = 0.5,
                   etot_conv_thr: float = 1e-5, forc_conv_thr: float = 1e-4,
                   nstep: int = 100, cell_dofree: str = "all",
                   extra_control: dict | None = None, extra_system: dict | None = None,
                   extra_electrons: dict | None = None, kpoints_card: str | None = None) -> str:
    """pw.x input for calculation in {'scf', 'nscf', 'vc-relax', 'relax'}. Non-spin-polarised
    metal. ``kpoints_card`` replaces the automatic K_POINTS card (e.g. an explicit nscf list)."""
    if calculation not in ("scf", "nscf", "vc-relax", "relax"):
        raise ValueError(f"unsupported calculation {calculation!r}")
    control = {"calculation": calculation, "prefix": prefix, "outdir": outdir,
               "pseudo_dir": pseudo_dir, "tprnfor": True, "tstress": True,
               "verbosity": "high" if calculation in ("scf", "nscf") else None}
    if calculation not in ("scf", "nscf"):
        control.update({"etot_conv_thr": float(etot_conv_thr),
                        "forc_conv_thr": float(forc_conv_thr), "nstep": int(nstep)})
    control.update(extra_control or {})
    system = {"ibrav": 0, "nat": cell.nat, "ntyp": len(cell.elements),
              "ecutwfc": float(ecutwfc),
              "ecutrho": float(ecutrho) if ecutrho else None,
              "occupations": "smearing", "smearing": smearing, "degauss": float(degauss)}
    if la2F:
        system["la2F"] = True
    system.update(extra_system or {})
    electrons = {"conv_thr": float(conv_thr), "mixing_beta": float(mixing_beta)}
    electrons.update(extra_electrons or {})
    s = namelist("control", control) + namelist("system", system) + namelist("electrons", electrons)
    if calculation in ("relax", "vc-relax"):
        s += namelist("ions", {"ion_dynamics": "bfgs"})
    if calculation == "vc-relax":
        s += namelist("cell", {"cell_dynamics": "bfgs",
                               "press": float(press_kbar or 0.0),
                               "press_conv_thr": float(press_conv_thr_kbar),
                               "cell_dofree": cell_dofree})
    s += _species_card(cell, pseudo_files) + _cell_cards(cell)
    s += kpoints_card if kpoints_card is not None else _kpoints(kgrid)
    return s


def write_ph_input(*, prefix: str, outdir: str, cell: CellData, qgrid: Iterable[int],
                   el_ph_sigma: float, el_ph_nsigma: int, tr2_ph: float = 1e-14,
                   fildyn: str = "dyn", fildvscf: str | None = "dvscf",
                   recover: bool = False, alpha_mix: float | None = None,
                   title: str = "electron-phonon (interpolated)",
                   extra: dict | None = None) -> str:
    q = list(qgrid)
    vals: dict[str, Any] = {"prefix": prefix, "outdir": outdir, "tr2_ph": float(tr2_ph),
                            "fildyn": fildyn, "fildvscf": fildvscf,
                            "electron_phonon": "interpolated",
                            "el_ph_sigma": float(el_ph_sigma), "el_ph_nsigma": int(el_ph_nsigma),
                            "trans": True, "ldisp": True,
                            "nq1": int(q[0]), "nq2": int(q[1]), "nq3": int(q[2]),
                            "recover": True if recover else None,
                            "alpha_mix(1)": float(alpha_mix) if alpha_mix else None}
    for i, el in enumerate(cell.elements, start=1):
        vals[f"amass({i})"] = float(cell.masses_amu[el])
    vals.update(extra or {})
    return title[:75] + "\n" + namelist("inputph", vals)


def write_q2r_input(*, fildyn: str = "dyn", flfrc: str = "ifc.fc", zasr: str = "crystal",
                    el_ph_nsigma: int = 10) -> str:
    return namelist("input", {"fildyn": fildyn, "flfrc": flfrc, "zasr": zasr, "la2F": True,
                              "el_ph_nsigma": int(el_ph_nsigma)})


def write_matdyn_dos_input(*, cell: CellData, flfrc: str = "ifc.fc",
                           kgrid: Iterable[int] = (16, 16, 16), asr: str = "crystal",
                           ndos: int | None = 200, deltaE_cm1: float | None = None,
                           fldos: str = "phonon.dos", el_ph_nsigma: int = 10) -> str:
    """matdyn.x with dos=.true., la2F=.true.: writes a2F.dos1..N and 'lambda' in its cwd."""
    k = list(kgrid)
    vals: dict[str, Any] = {"asr": asr, "flfrc": flfrc, "flfrq": "freq.dos",
                            "dos": True, "la2F": True, "fldos": fldos,
                            "nk1": int(k[0]), "nk2": int(k[1]), "nk3": int(k[2]),
                            "el_ph_nsigma": int(el_ph_nsigma)}
    if deltaE_cm1 is not None:
        vals["deltaE"] = float(deltaE_cm1)
    elif ndos is not None:
        vals["ndos"] = int(ndos)
    for i, el in enumerate(cell.elements, start=1):
        vals[f"amass({i})"] = float(cell.masses_amu[el])
    return namelist("input", vals)
