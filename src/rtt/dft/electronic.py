"""S4 dft-electronic stage: relax, metallicity, magnetism, N(E_F), E_F, n_sheets, Belli phi
(architecture §3.1 S4; doc 02 H1, H2, H7; doc 01 G-VALID item 1 for E_F).

Chain (profile ``s4``; workdir ``<w>``; QE scratch ``<w>/tmp``)::

    relax    pw.x vc-relax, non-spin-polarised, press = 10 * pressure_GPa kbar
    scf      pw.x scf on the relaxed cell (k from k_spacing)
    elf      pp.x plot_num=8 -> elf.cube                          (Belli phi)
    nscf     pw.x nscf, smearing, denser k mesh (nscf_k_spacing)  (for the ILDOS)
    ildos    pp.x plot_num=10, E_F +- window -> ildos.cube     (only if H present; H_DOS)
    nscf_tet pw.x nscf, occupations='tetrahedra_opt', same mesh -> eigenvalues (XML)
    dos      dos.x (tetrahedra)                                 -> pwscf.dos: N(E_F)
    spin     pw.x scf, nspin=2, starting_magnetization = 0.5 on every species (ferromagnetic
             start), separate prefix, on the relaxed cell        -> moments and E_spin

Derived quantities:

* ``N_EF_per_eV_atom``: **tetrahedron** DOS at E_Fermi (dos.x, both spins), per atom. The
  Gaussian estimates from the eigenvalues at sigma and 2 sigma are kept as cross-checks
  (``N_EF_gauss_*``; flag ``NEF_MESH_SENSITIVE`` when they differ by > 25 %).
* ``n_sheets``: number of bands with min_k e < E_Fermi < max_k e on the nscf mesh (per spin,
  summed); a band with several pockets counts once (lower bound on the number of sheets).
  ``crossing_bands`` gives each band's depth, height, width and DOS fraction at E_F;
  ``sheet_w_min`` is the smallest DOS fraction (min over sigma and 2 sigma Gaussians), the
  w_min input of N1-DOS / the S5 multiband rule (doc 01 §5.1).
* ``gap_eV``: 0 if any band crosses; else min(e > E_Fermi) - max(e <= E_Fermi) on the mesh.
* ``E_F_eV`` (G-VALID item 1, literal): min over crossing bands of min(depth, width/2).
* moments: |total| and absolute magnetization of the spin scf; ``e_total_spin_Ry`` vs the NM
  ``e_total_Ry`` (same cell, k mesh, cutoff and smearing) decides which is the ground state.

Decision (:func:`decide_s4`, revised after review RC1):

* **H1** (doc 02: "N(E_F) > 0 at the PBE level"): metallic iff the largest of the N(E_F)
  estimates exceeds ``nef_min`` (numerical floor, provisional). Otherwise EXCLUDED only if the
  gap exceeds gap_eV_max + ``marginal_gap`` (mesh resolution + smearing); a gap inside the
  marginal band gives ``DEFERRED(H1-marginal)`` (doc 02: hybrid check for marginal gaps).
* **H2** (doc 02: non-magnetic *ground state*): the moment per atom is max(|M_tot|, |M|_abs)/nat
  (so ferri-/antiferromagnetic solutions count). Above the per-atom tolerance it excludes only
  if the spin solution is lower in energy than the NM one by more than ``mag_energy_tol``;
  a higher-energy (metastable) magnetic solution passes with a flag; missing energies or a moment
  within 2x tolerance give ``DEFERRED``.
* **H7-heur**: formation energy vs elemental references (a lower bound on E_hull) above the
  window (+ an optional ZPE allowance) gives SCREENED_OUT, but ``DEFERRED`` when the run is
  flagged ``RELAX_PRESSURE_MISMATCH`` (cutoff unconverged).
* Every outcome decided with a provisional threshold carries ``provisional=True`` and a
  ``PROVISIONAL_THRESHOLD(<check>)`` tag; checks not made carry ``<check>_NOT_CHECKED``.
"""
from __future__ import annotations

import json
import math
import re
import xml.etree.ElementTree as ET
from dataclasses import asdict, dataclass, field
from pathlib import Path
from typing import Any

import numpy as np

from .belli import BelliResult, phi_and_hdos
from .profiles import PROFILES, Profile, pseudo_md5s, sg15_pseudos, th_parts
from .workflow import Chain, ChainStep, DFPTRun, RunPolicy, default_launcher

HA_EV = 27.211386245988
RY_EV = HA_EV / 2.0
GPA_TO_KBAR = 10.0


# ---------------------------------------------------------------------- band analysis
@dataclass
class Bands:
    eig_eV: np.ndarray            # (nspin, nk, nbnd)
    weights: np.ndarray           # (nk,) normalised to sum 1
    fermi_eV: float
    nat: int
    lsda: bool
    nelec: float | None = None

    @property
    def spin_degeneracy(self) -> float:
        return 1.0 if self.lsda else 2.0


def parse_pw_xml(path: str | Path) -> Bands:
    """Eigenvalues, k weights and Fermi energy from a pw.x data-file-schema XML (QE >= 6.2)."""
    root = ET.parse(str(path)).getroot()
    out = root.find("output")
    bs = out.find("band_structure")
    lsda = (bs.findtext("lsda") or "false").strip().lower() == "true"
    ef = bs.findtext("fermi_energy")
    if ef is None:
        hl = bs.findtext("highestOccupiedLevel")
        if hl is None:
            raise ValueError(f"{path}: no fermi_energy")
        ef = hl
    nat = int(out.find("atomic_structure").get("nat"))
    ne = bs.findtext("nelec")
    eigs, wts = [], []
    for ks in bs.findall("ks_energies"):
        kp = ks.find("k_point")
        wts.append(float(kp.get("weight")))
        eigs.append([float(x) for x in ks.findtext("eigenvalues").split()])
    e = np.asarray(eigs) * HA_EV
    if lsda:
        nup = int(bs.findtext("nbnd_up"))
        e = np.stack([e[:, :nup], e[:, nup:]])
    else:
        e = e[None, :, :]
    w = np.asarray(wts, float)
    return Bands(e, w / w.sum(), float(ef) * HA_EV, nat, lsda,
                 float(ne) if ne is not None else None)


@dataclass
class BandAnalysis:
    gap_eV: float
    n_sheets: int
    n_sheets_by_spin: list[int]
    E_F_eV: float | None              # G-VALID definition
    E_F_hole_eV: float | None
    N_EF_per_eV_cell: float           # Gaussian, sigma
    N_EF_per_eV_atom: float
    fermi_level_eV: float
    vbm_eV: float | None
    cbm_eV: float | None
    crossing_bands: list = field(default_factory=list)
    N_EF_2sigma_per_eV_atom: float | None = None
    sheet_w_min: float | None = None


def _gauss(x: np.ndarray, s: float) -> np.ndarray:
    return np.exp(-0.5 * (x / s) ** 2) / (s * math.sqrt(2 * math.pi))


def analyse_bands(b: Bands, dos_sigma_eV: float = 0.136) -> BandAnalysis:
    ef = b.fermi_eV
    sheets, efs, holes = [], [], []
    crossing: list[dict] = []
    g1 = _gauss(b.eig_eV - ef, dos_sigma_eV)
    g2 = _gauss(b.eig_eV - ef, 2 * dos_sigma_eV)
    per_band1 = b.spin_degeneracy * np.einsum("skn,k->sn", g1, b.weights)
    per_band2 = b.spin_degeneracy * np.einsum("skn,k->sn", g2, b.weights)
    n_cell, n2 = float(per_band1.sum()), float(per_band2.sum())
    for s in range(b.eig_eV.shape[0]):
        e = b.eig_eV[s]
        lo, hi = e.min(axis=0), e.max(axis=0)
        cross = (lo < ef) & (hi > ef)
        sheets.append(int(cross.sum()))
        for n in np.nonzero(cross)[0]:
            efs.append(min(ef - lo[n], 0.5 * (hi[n] - lo[n])))
            holes.append(min(ef - lo[n], hi[n] - ef))
            f1 = float(per_band1[s, n] / n_cell) if n_cell > 0 else 0.0
            f2 = float(per_band2[s, n] / n2) if n2 > 0 else 0.0
            crossing.append({"spin": s, "band": int(n), "depth_eV": float(ef - lo[n]),
                             "height_eV": float(hi[n] - ef), "width_eV": float(hi[n] - lo[n]),
                             "dos_fraction": f1, "dos_fraction_2sigma": f2})
    allv = b.eig_eV.reshape(-1)
    below, above = allv[allv <= ef], allv[allv > ef]
    vbm = float(below.max()) if below.size else None
    cbm = float(above.min()) if above.size else None
    n_sh = int(sum(sheets))
    gap = 0.0 if n_sh > 0 or vbm is None or cbm is None else max(cbm - vbm, 0.0)
    w_min = (min(min(c["dos_fraction"], c["dos_fraction_2sigma"]) for c in crossing)
             if crossing else None)
    return BandAnalysis(gap_eV=float(gap), n_sheets=n_sh, n_sheets_by_spin=sheets,
                        E_F_eV=float(min(efs)) if efs else None,
                        E_F_hole_eV=float(min(holes)) if holes else None,
                        N_EF_per_eV_cell=n_cell, N_EF_per_eV_atom=n_cell / b.nat,
                        fermi_level_eV=float(ef), vbm_eV=vbm, cbm_eV=cbm,
                        crossing_bands=crossing, N_EF_2sigma_per_eV_atom=n2 / b.nat,
                        sheet_w_min=w_min)


_DOS_HEAD = re.compile(r"EFermi\s*=\s*([-+]?\d+\.\d+)")


def parse_dosx(path: str | Path) -> dict[str, float]:
    """dos.x output -> N(E_F) per cell (both spins) by linear interpolation at the header's
    EFermi, and the integrated DOS there (which must equal nelec for a correct file)."""
    lines = Path(path).read_text(errors="replace").splitlines()
    m = _DOS_HEAD.search(lines[0])
    if not m:
        raise ValueError(f"{path}: no EFermi in header")
    ef = float(m.group(1))
    spin = "dosup" in lines[0]
    rows = np.array([[float(x) for x in ln.split()] for ln in lines[1:] if ln.strip()
                     and not ln.lstrip().startswith("#")])
    E = rows[:, 0]
    dos = rows[:, 1] + rows[:, 2] if spin else rows[:, 1]
    idos = rows[:, 3] if spin else rows[:, 2]
    return {"EFermi_eV": ef, "N_EF_per_eV_cell": float(np.interp(ef, E, dos)),
            "int_dos_at_EF": float(np.interp(ef, E, idos))}


_TOTMAG = re.compile(r"total magnetization\s+=\s+([-+]?\d+\.\d+)\s+Bohr mag/cell")
_ABSMAG = re.compile(r"absolute magnetization\s+=\s+([-+]?\d+\.\d+)\s+Bohr mag/cell")
_ETOT = re.compile(r"^!\s+total energy\s+=\s+([-+]?\d+\.\d+)\s+Ry", re.M)
_EF_OUT = re.compile(r"the Fermi energy is\s+([-+]?\d+\.\d+)\s+ev")
_PRESS = re.compile(r"P=\s*([-+]?\d+\.\d+)")
PULAY_TOL_KBAR = 10.0
NEF_REL_TOL = 0.25      # N(E_F) at sigma vs 2 sigma: larger relative change = mesh too coarse


def parse_pressure_kbar(text: str) -> float | None:
    p = _PRESS.findall(text)
    return float(p[-1]) if p else None


def parse_magnetization(text: str) -> tuple[float | None, float | None]:
    t = _TOTMAG.findall(text)
    a = _ABSMAG.findall(text)
    return (abs(float(t[-1])) if t else None, float(a[-1]) if a else None)


def parse_total_energy_Ry(text: str) -> float | None:
    e = _ETOT.findall(text)
    return float(e[-1]) if e else None


# ---------------------------------------------------------------------- result
@dataclass
class ElectronicResult:
    run: DFPTRun
    gap_eV: float | None = None
    moment_muB: float | None = None              # |total magnetization| per cell
    abs_moment_muB: float | None = None          # absolute magnetization per cell
    N_EF_per_eV_atom: float | None = None        # tetrahedron (dos.x) when available
    N_EF_per_eV_cell: float | None = None
    N_EF_gauss_per_eV_atom: float | None = None
    N_EF_2sigma_per_eV_atom: float | None = None
    int_dos_at_EF: float | None = None
    nelec: float | None = None
    E_F_eV: float | None = None
    E_F_hole_eV: float | None = None
    fermi_level_eV: float | None = None
    n_sheets: int | None = None
    n_sheets_by_spin: list[int] = field(default_factory=list)
    crossing_bands: list[dict] = field(default_factory=list)
    sheet_w_min: float | None = None
    phi: float | None = None
    H_f: float | None = None
    H_DOS: float | None = None
    Tc_phi_K: float | None = None
    e_total_Ry: float | None = None
    e_total_spin_Ry: float | None = None
    scf_pressure_kbar: float | None = None      # fresh-basis scf on the relaxed cell
    nat: int | None = None
    composition: dict = field(default_factory=dict)
    relaxed_cell: dict | None = None
    flags: list[str] = field(default_factory=list)
    error: str = ""

    @property
    def ok(self) -> bool:
        return self.run.ok and not self.error

    @property
    def e_per_atom_eV(self) -> float | None:
        return None if self.e_total_Ry is None or not self.nat else self.e_total_Ry * RY_EV / self.nat

    @property
    def moment_per_atom_muB(self) -> float | None:
        vals = [v for v in (self.moment_muB, self.abs_moment_muB) if v is not None]
        if not vals or not self.nat:
            return None
        return max(abs(v) for v in vals) / self.nat

    @property
    def magnetic_energy_meV_atom(self) -> float | None:
        """E_spin - E_NM per atom (negative: the magnetic solution is lower)."""
        if self.e_total_Ry is None or self.e_total_spin_Ry is None or not self.nat:
            return None
        return (self.e_total_spin_Ry - self.e_total_Ry) * RY_EV * 1000.0 / self.nat

    def metrics(self) -> dict[str, Any]:
        """S4 metrics with the §6.1 key names."""
        return {"gap_eV": self.gap_eV, "moment_muB": self.moment_muB,
                "abs_moment_muB": self.abs_moment_muB,
                "moment_per_atom_muB": self.moment_per_atom_muB,
                "magnetic_energy_meV_atom": self.magnetic_energy_meV_atom,
                "N_EF_per_eV_atom": self.N_EF_per_eV_atom,
                "N_EF_gauss_per_eV_atom": self.N_EF_gauss_per_eV_atom,
                "N_EF_2sigma_per_eV_atom": self.N_EF_2sigma_per_eV_atom,
                "int_dos_at_EF": self.int_dos_at_EF, "nelec": self.nelec, "H_DOS": self.H_DOS,
                "E_F_eV": self.E_F_eV, "E_F_hole_eV": self.E_F_hole_eV,
                "fermi_level_eV": self.fermi_level_eV, "n_sheets": self.n_sheets,
                "sheet_w_min": self.sheet_w_min, "crossing_bands": self.crossing_bands,
                "phi": self.phi, "H_f": self.H_f, "Tc_phi_K": self.Tc_phi_K,
                "e_total_Ry": self.e_total_Ry, "e_total_spin_Ry": self.e_total_spin_Ry,
                "e_per_atom_eV": self.e_per_atom_eV,
                "scf_pressure_kbar": self.scf_pressure_kbar, "flags": list(self.flags),
                "wall_s": self.run.wall_s, "core_h": self.run.core_h}

    def to_dict(self) -> dict[str, Any]:
        d = asdict(self)
        d["run"] = self.run.to_dict()
        d["metrics"] = self.metrics()
        return d


# ---------------------------------------------------------------------- chain
def _pp_input(filplot: str, plot_num: int, fileout: str, extra: dict | None = None) -> str:
    from rtt.core.dfpt.qe_inputs import namelist
    inp = {"prefix": "pwscf", "outdir": "./tmp/", "filplot": filplot, "plot_num": plot_num}
    inp.update(extra or {})
    return namelist("inputpp", inp) + namelist("plot", {"iflag": 3, "output_format": 6,
                                                        "fileout": fileout})


DOS_WINDOW_EV = 1.0     # dos.x only around E_F: the full range (SG15 semicore ~ -100 eV) took 8 min for Al


def _dos_input(ef_eV: float | None = None) -> str:
    from rtt.core.dfpt.qe_inputs import namelist
    v: dict[str, Any] = {"prefix": "pwscf", "outdir": "./tmp/", "fildos": "pwscf.dos",
                         "DeltaE": 0.005}
    if ef_eV is not None:
        v["Emin"] = round(ef_eV - DOS_WINDOW_EV, 6)
        v["Emax"] = round(ef_eV + DOS_WINDOW_EV, 6)
    return namelist("dos", v)


def _check_out(name: str, *needles: str, files: tuple[str, ...] = ()):
    def chk(w: Path) -> tuple[bool, str]:
        out = w / f"{name}.out"
        if not out.exists():
            return False, "no output"
        text = out.read_text(errors="replace")
        if name.startswith(("relax", "scf", "nscf", "spin")) and "JOB DONE" not in text:
            return False, "no 'JOB DONE'"
        for n in needles:
            if n not in text:
                return False, f"{n!r} not in output"
        for f in files:
            if not (w / f).exists():
                return False, f"{f} missing"
        return True, ""
    return chk


S4_ARTIFACTS = ("relax.in", "relax.out", "scf.in", "scf.out", "nscf.in", "nscf.out",
                "nscf_tet.in", "nscf_tet.out", "dos.in", "pwscf.dos", "spin.in", "spin.out",
                "elf.in", "ildos.in", "relaxed_structure.json", "chain_state.json",
                "relax.nstep_retry", "relax.nstep1.out", "relax.fallback",
                "timing.jsonl", "electronic_result.json", "tmp/pwscf.save/data-file-schema.xml")


def run_electronic(struct: Any, workdir: str | Path, *, profile: Profile | None = None,
                   launcher: Any = None, material_id: str = "material",
                   pressure_GPa: float = 0.0, policy: RunPolicy | None = None
                   ) -> ElectronicResult:
    """Run the S4 chain and derive gap, moments, N(E_F), E_F, n_sheets and Belli phi."""
    from rtt.core.dfpt import qe_inputs as qi
    from rtt.core.dfpt.grids import grid_from_spacing
    from rtt.core.dfpt.structure import CellData, parse_pw_final_structure, to_celldata
    from .workflow import _masses_for

    prof = profile or PROFILES["s4"]
    launcher = launcher or default_launcher(2)
    w = Path(workdir).expanduser()
    w.mkdir(parents=True, exist_ok=True)
    cell0 = to_celldata(struct, _masses_for(struct))
    pseudos = sg15_pseudos(cell0.elements, pattern=prof.pseudo_pattern)
    has_h = "H" in cell0.elements
    common = dict(prefix="pwscf", outdir="./tmp/",
                  pseudo_dir=str(Path(pseudos.directory).resolve()) + "/",
                  pseudo_files={el: pseudos.files[el] for el in cell0.elements},
                  ecutwfc=float(prof.ecutwfc_Ry), ecutrho=prof.ecutrho_Ry,
                  smearing=prof.smearing, degauss=float(prof.degauss_Ry),
                  conv_thr=prof.conv_thr, mixing_beta=prof.mixing_beta)

    def relaxed() -> Any:
        if not prof.do_relax:
            return cell0
        p = w / "relaxed_structure.json"
        if not p.exists():
            raise FileNotFoundError("relaxed_structure.json (run relax first)")
        return CellData.from_dict(json.loads(p.read_text()))

    def kgrid(cell: Any, spacing: float) -> tuple[int, int, int]:
        return grid_from_spacing(cell.reciprocal_lengths_2pi(), spacing)

    def after_relax(wd: Path) -> None:
        c = parse_pw_final_structure((wd / "relax.out").read_text(errors="replace"), cell0)
        (wd / "relaxed_structure.json").write_text(json.dumps(c.to_dict(), indent=1))

    def nscf_fermi() -> float:
        t = (w / "nscf.out").read_text(errors="replace")
        m = _EF_OUT.findall(t)
        if not m:
            raise ValueError("no Fermi energy in nscf.out")
        return float(m[-1])

    nscf_k = lambda: kgrid(relaxed(), prof.nscf_k_spacing or prof.k_spacing / 2)  # noqa: E731
    steps: list[ChainStep] = []
    if prof.do_relax:
        steps.append(ChainStep(
            "relax", "pw.x",
            lambda: qi.write_pw_input(calculation="vc-relax", cell=cell0,
                                      kgrid=kgrid(cell0, prof.k_spacing),
                                      press_kbar=pressure_GPa * GPA_TO_KBAR,
                                      press_conv_thr_kbar=prof.press_conv_thr_kbar,
                                      etot_conv_thr=prof.etot_conv_thr,
                                      forc_conv_thr=prof.forc_conv_thr, nstep=prof.relax_nstep,
                                      cell_dofree=prof.cell_dofree, **common),
            _check_out("relax", "End final coordinates", files=("relaxed_structure.json",)),
            after=after_relax))
    up_r = ("relax",) if prof.do_relax else ()
    steps.append(ChainStep(
        "scf", "pw.x",
        lambda: qi.write_pw_input(calculation="scf", cell=relaxed(),
                                  kgrid=kgrid(relaxed(), prof.k_spacing), **common),
        _check_out("scf", "convergence has been achieved"), upstream=up_r))
    steps.append(ChainStep("elf", "pp.x", lambda: _pp_input("elf.dat", 8, "elf.cube"),
                           _check_out("elf", files=("elf.cube",)), upstream=("scf",),
                           parallel=False))
    steps.append(ChainStep(
        "nscf", "pw.x",
        lambda: qi.write_pw_input(calculation="nscf", cell=relaxed(), kgrid=nscf_k(), **common),
        _check_out("nscf", "the Fermi energy is"), upstream=("elf",)))
    last = "nscf"
    if has_h:
        steps.append(ChainStep(
            "ildos", "pp.x",
            lambda: _pp_input("ildos.dat", 10, "ildos.cube",
                              {"emin": round(nscf_fermi() - prof.ildos_window_eV, 6),
                               "emax": round(nscf_fermi() + prof.ildos_window_eV, 6)}),
            _check_out("ildos", files=("ildos.cube",)), upstream=("nscf",), parallel=False))
        last = "ildos"
    tet_sys = {"occupations": "tetrahedra_opt", "smearing": None, "degauss": None}
    steps.append(ChainStep(
        "nscf_tet", "pw.x",
        lambda: qi.write_pw_input(calculation="nscf", cell=relaxed(), kgrid=nscf_k(),
                                  extra_system=tet_sys, **common),
        _check_out("nscf_tet", "the Fermi energy is"), upstream=(last,)))
    def tet_fermi() -> float:
        m = _EF_OUT.findall((w / "nscf_tet.out").read_text(errors="replace"))
        if not m:
            raise ValueError("no Fermi energy in nscf_tet.out")
        return float(m[-1])

    steps.append(ChainStep("dos", "dos.x", lambda: _dos_input(tet_fermi()),
                           _check_out("dos", files=("pwscf.dos",)),
                           upstream=("nscf_tet",), parallel=False))
    spin_sys: dict[str, Any] = {"nspin": 2}
    for i, _el in enumerate(cell0.elements, start=1):
        spin_sys[f"starting_magnetization({i})"] = float(prof.starting_magnetization)
    spin_common = dict(common, prefix="pwscf_spin")
    steps.append(ChainStep(
        "spin", "pw.x",
        lambda: qi.write_pw_input(calculation="scf", cell=relaxed(),
                                  kgrid=kgrid(relaxed(), prof.k_spacing),
                                  extra_system=spin_sys, **spin_common),
        _check_out("spin", "convergence has been achieved"), upstream=up_r))

    run = DFPTRun(material_id=material_id, profile=prof.to_dict(), workdir=str(w),
                  status="failed", cell0=cell0.to_dict(), pseudo_dir=str(pseudos.directory),
                  pseudo_files=dict(pseudos.files),
                  pseudo_md5=pseudo_md5s(pseudos, cell0.elements))
    import time as _t
    run.started = _t.time()
    chain = Chain(w, steps, launcher, material_id=material_id, profile=prof.name, policy=policy)
    ok, err = chain.run()
    run.steps = chain.timings
    run.finished = _t.time()
    res = ElectronicResult(run=run, nat=cell0.nat)
    comp: dict[str, int] = {}
    for s in cell0.species:
        comp[s] = comp.get(s, 0) + 1
    res.composition = comp
    if not ok:
        run.error = res.error = err
        _write(run, res)
        return res
    run.status = "done"
    try:
        rc = relaxed()
        res.relaxed_cell = run.relaxed_cell = rc.to_dict()
        scf_text = (w / "scf.out").read_text(errors="replace")
        res.e_total_Ry = parse_total_energy_Ry(scf_text)
        res.scf_pressure_kbar = parse_pressure_kbar(scf_text)
        target = pressure_GPa * GPA_TO_KBAR
        if (prof.do_relax and res.scf_pressure_kbar is not None
                and abs(res.scf_pressure_kbar - target) > max(PULAY_TOL_KBAR, 0.05 * abs(target))):
            # vc-relax keeps the initial plane-wave basis; a fresh-basis scf far from the
            # target pressure means the cutoff is unconverged for these pseudopotentials
            # (e.g. SG15 Al-1.2 with 11 valence electrons at 25 Ry: P = -6800 kbar).
            res.flags.append(f"RELAX_PRESSURE_MISMATCH(P={res.scf_pressure_kbar:.1f}kbar,"
                             f"target={target:.1f})")
        stext = (w / "spin.out").read_text(errors="replace")
        res.moment_muB, res.abs_moment_muB = parse_magnetization(stext)
        res.e_total_spin_Ry = parse_total_energy_Ry(stext)
        bands = parse_pw_xml(_xml_path(w))
        ba = analyse_bands(bands, prof.dos_sigma_eV)
        res.nelec = bands.nelec
        res.gap_eV, res.n_sheets, res.n_sheets_by_spin = ba.gap_eV, ba.n_sheets, ba.n_sheets_by_spin
        res.E_F_eV, res.E_F_hole_eV, res.fermi_level_eV = ba.E_F_eV, ba.E_F_hole_eV, ba.fermi_level_eV
        res.crossing_bands, res.sheet_w_min = ba.crossing_bands, ba.sheet_w_min
        res.N_EF_gauss_per_eV_atom = ba.N_EF_per_eV_atom
        res.N_EF_2sigma_per_eV_atom = ba.N_EF_2sigma_per_eV_atom
        if ba.N_EF_per_eV_atom > 0 and abs(ba.N_EF_2sigma_per_eV_atom / ba.N_EF_per_eV_atom - 1) > NEF_REL_TOL:
            res.flags.append(f"NEF_MESH_SENSITIVE(sigma:{ba.N_EF_per_eV_atom:.3f},"
                             f"2sigma:{ba.N_EF_2sigma_per_eV_atom:.3f})")
        d = parse_dosx(w / "pwscf.dos")
        res.N_EF_per_eV_cell = d["N_EF_per_eV_cell"]
        res.N_EF_per_eV_atom = d["N_EF_per_eV_cell"] / cell0.nat
        # dos.x runs on a window around E_F, so its integral starts at E_F - window, not at the
        # band bottom: record it but do not compare with nelec (checked once on Al with the full
        # range: 11.0 = nelec, i.e. dos.x reports the spin-summed DOS).
        res.int_dos_at_EF = d["int_dos_at_EF"]
        if abs(d["EFermi_eV"] - ba.fermi_level_eV) > 0.05:
            res.flags.append(f"DOS_EF_MISMATCH(dos.x={d['EFermi_eV']:.3f},xml={ba.fermi_level_eV:.3f})")
        if 0 < ba.gap_eV < 0.05:
            res.flags.append("GAP_SMALL_MESH_LIMITED")
        if (res.moment_muB is not None and res.abs_moment_muB is not None
                and res.moment_muB < 0.05 < res.abs_moment_muB):
            res.flags.append("ANTIFERRO_OR_FERRI_SUSPECTED")
        br: BelliResult = phi_and_hdos(w)
        res.phi, res.H_f, res.H_DOS, res.Tc_phi_K = br.phi, br.H_f, br.H_DOS, br.Tc_phi_K
        from rtt.core.dfpt.parsers import parse_qe_version
        run.qe_version = parse_qe_version(scf_text)
    except Exception as exc:
        run.status = "failed"
        run.error = res.error = f"post-processing failed: {exc}"
    _write(run, res)
    return res


def _xml_path(w: Path) -> Path:
    for p in (w / "tmp" / "pwscf.save" / "data-file-schema.xml", w / "tmp" / "pwscf.xml"):
        if p.is_file():
            return p
    raise FileNotFoundError("pw.x XML data file not found under tmp/")


def _write(run: DFPTRun, res: ElectronicResult) -> None:
    Path(run.workdir, "electronic_result.json").write_text(
        json.dumps(res.to_dict(), indent=1, default=str))


def electronic_artifacts(res: ElectronicResult, *, ledger: Any = None,
                         tier_run_id: int | None = None) -> list:
    """sha256 :class:`ArtifactRef` of every S4 file the decision rests on (review RC1 N1:
    H1/H2 removals must be re-verifiable). Registered in the ledger when one is given. Call
    before any final pruning of ``tmp/``."""
    from rtt.records import ArtifactRef
    from .prune import sha256_file
    w = Path(res.run.workdir)
    out = []
    for name in S4_ARTIFACTS:
        p = w / name
        if not p.is_file():
            continue
        if ledger is not None and tier_run_id is not None:
            from rtt.ledger.artifacts import register
            out.append(register(ledger, tier_run_id, p))
        else:
            out.append(ArtifactRef(path=str(p), sha256=sha256_file(p), bytes=p.stat().st_size))
    return out


# ---------------------------------------------------------------------- heuristic hull + decide
def formation_energy_per_atom_eV(e_total_Ry: float, composition: dict[str, int],
                                 element_refs_eV_per_atom: dict[str, float]) -> float:
    """(E - sum_el n_el * mu_el) / N with mu_el the reference energy per atom of the element at
    the *same* profile (cutoff, pseudopotentials). A lower bound on E_hull, because the convex
    hull lies at or below the elemental tie line."""
    n = sum(composition.values())
    ref = sum(cnt * element_refs_eV_per_atom[el] for el, cnt in composition.items())
    return (e_total_Ry * RY_EV - ref) / n


_th = th_parts          # backwards-compatible name

#: numerical floor for "N(E_F) > 0" (states/eV/atom, both spins). Provisional.
NEF_MIN = {"value": 0.01, "provisional": True,
           "source": "doc02#5a:H1 (N(E_F)>0); numerical floor doc22#RC1-F2"}
#: band above gap_eV_max inside which a gap is mesh/smearing-limited (eV): ~ degauss 0.02 Ry.
MARGINAL_GAP = {"value": 0.3, "provisional": True,
                "source": "doc02#5a:H1 (hybrid check if marginal); doc22#RC1-F3 (smearing 0.27 eV)"}
#: E_spin must be this far below E_NM (meV/atom) for a magnetic ground state.
MAG_ENERGY_TOL = {"value": 1.0, "provisional": True, "source": "doc22#RC1-F4"}


def decide_s4(res: ElectronicResult, *, gap_max: Any, moment_tol: Any,
              ehull_heur_max: Any = None, e_form_eV: float | None = None,
              nef_min: Any = NEF_MIN, marginal_gap: Any = MARGINAL_GAP,
              mag_energy_tol: Any = MAG_ENERGY_TOL, zpe_allowance_eV: Any = None,
              enable_h1: bool = True, enable_h2: bool = True, enable_h7: bool = True,
              tier_run_id: int | None = None) -> Any:
    """Pure S4 decision -> :class:`rtt.records.Outcome` (stage S4); rules in the module
    docstring. ``moment_tol`` is a tolerance **per atom** (muB/atom). Thresholds may be
    :class:`rtt.config.Threshold`, inline TOML tables (``{gap_eV_max=..., source=...}``) or
    numbers. Order: ERROR -> H1 -> H2 -> H7-heur -> PASS."""
    from rtt.records import (DEFERRED, ERROR, EXCLUDED_ASSUMPTION, PASS, SCREENED_OUT,
                             Outcome)
    if not res.ok:
        return Outcome(stage="S4", outcome=ERROR, reason=res.error or res.run.error or "failed",
                       tier_run_id=tier_run_id)
    tags = list(res.flags)
    prov_used: list[str] = []

    def out(outcome: str, **kw: Any) -> Any:
        t = tags + [f"PROVISIONAL_THRESHOLD({c})" for c in dict.fromkeys(prov_used)]
        return Outcome(stage="S4", outcome=outcome, tier_run_id=tier_run_id,
                       provisional=bool(prov_used) or kw.pop("provisional", False), tags=t, **kw)

    # ------------------------------------------------------------------ H1
    g, gsrc, gprov = th_parts(gap_max)
    nmin, nsrc, nprov = th_parts(nef_min)
    mg, mgsrc, mgprov = th_parts(marginal_gap)
    # Campaign 1 fix: the tetrahedron DOS decides when present. The Gaussian values (sigma and
    # 2 sigma = 0.27 eV) are a fallback only: their tails reach across a gap (NaAuH2, PBE gap
    # 2.0 eV, tetrahedron N(E_F) = 0, 2-sigma Gaussian 0.017 > nef_min) and made an insulator PASS H1.
    if res.N_EF_per_eV_atom is not None:
        nef = res.N_EF_per_eV_atom
    else:
        nefs = [v for v in (res.N_EF_gauss_per_eV_atom, res.N_EF_2sigma_per_eV_atom) if v is not None]
        nef = max(nefs) if nefs else None
    if not enable_h1:
        tags.append("H1_NOT_CHECKED")
    elif ((nef is not None and nmin is not None and nef > float(nmin))
          or (res.n_sheets or 0) > 0 or (nef is None and res.gap_eV == 0.0)):
        # states at E_F, or a band crossing E_Fermi (gap exactly 0 = a band crosses)
        if nprov and nef is not None and nef > float(nmin or 0) and not (res.n_sheets or 0) > 0:
            prov_used.append("H1")
    elif res.gap_eV is not None and g is not None:
        band = float(mg) if mg is not None else 0.0
        if mgprov or gprov:
            prov_used.append("H1")
        if res.gap_eV > float(g) + band:
            return out(EXCLUDED_ASSUMPTION, criterion="H1", value=res.gap_eV,
                       threshold=float(g) + band, source_id=gsrc or "doc02#5a:H1",
                       reason=f"no states at E_F (N(E_F)={nef}) and PBE gap {res.gap_eV:.3f} eV "
                              f"> gap_eV_max + marginal band ({mgsrc})")
        return out(DEFERRED, criterion="H1-marginal", value=res.gap_eV,
                   threshold=float(g) + band, source_id=mgsrc or "doc02#5a:H1",
                   reason="gap within the mesh/smearing-limited band: hybrid-functional / denser "
                          "mesh check required (doc 02 H1)")
    else:
        tags.append("H1_NOT_CHECKED")

    # ------------------------------------------------------------------ H2
    m, msrc, mprov = th_parts(moment_tol)
    et, etsrc, etprov = th_parts(mag_energy_tol)
    mpa = res.moment_per_atom_muB
    if not enable_h2 or m is None or mpa is None:
        tags.append("H2_NOT_CHECKED")
    else:
        if mprov:
            prov_used.append("H2")
        tol = float(m)
        if mpa > tol:
            de = res.magnetic_energy_meV_atom
            kind = ("ferri/antiferromagnetic" if (res.moment_muB or 0) / max(res.nat or 1, 1) <= tol
                    else "ferromagnetic")
            if de is None:
                return out(DEFERRED, criterion="H2-energy-unknown", value=mpa, threshold=tol,
                           source_id=msrc or "doc02#5a:H2",
                           reason=f"{kind} solution ({mpa:.3f} muB/atom) but E_spin or E_NM "
                                  "missing: ground state not established")
            if etprov:
                prov_used.append("H2")
            if de < -float(et):
                if mpa > 2 * tol:
                    return out(EXCLUDED_ASSUMPTION, criterion="H2", value=mpa, threshold=tol,
                               source_id=msrc or "doc02#5a:H2",
                               reason=f"{kind} ground state: {mpa:.3f} muB/atom, "
                                      f"E_spin - E_NM = {de:.2f} meV/atom")
                return out(DEFERRED, criterion="H2-marginal", value=mpa, threshold=tol,
                           source_id=msrc or "doc02#5a:H2",
                           reason=f"magnetic ground state with a moment within 2x tolerance "
                                  f"({mpa:.3f} muB/atom)")
            tags.append(f"H2_METASTABLE_MAGNETIC_SOLUTION(dE={de:.2f}meV/atom,"
                        f"m={mpa:.3f}muB/atom)")

    # ------------------------------------------------------------------ H7-heur
    h, hsrc, hprov = th_parts(ehull_heur_max)
    z, zsrc, zprov = th_parts(zpe_allowance_eV)
    if enable_h7 and e_form_eV is not None and h is not None:
        lim = float(h) + (float(z) if z is not None else 0.0)
        if hprov or zprov:
            prov_used.append("H7")
        if e_form_eV > lim:
            if any(f.startswith("RELAX_PRESSURE_MISMATCH") for f in res.flags):
                return out(DEFERRED, criterion="H7-heur-unconverged", value=e_form_eV,
                           threshold=lim, source_id=hsrc or "doc02#5a:H7",
                           reason="formation energy above the window but the cutoff is not "
                                  "converged (RELAX_PRESSURE_MISMATCH)")
            return out(SCREENED_OUT, criterion="H7-heur", value=e_form_eV, threshold=lim,
                       source_id=hsrc or "doc02#5a:H7",
                       reason="formation energy vs elemental references (a lower bound on "
                              "E_hull) above the heuristic window")
    else:
        tags.append("H7_HEUR_NOT_CHECKED")
    return out(PASS)
