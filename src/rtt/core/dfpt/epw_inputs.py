"""EPW (anisotropic Migdal-Eliashberg) input writer and helpers.

Template: QE's own EPW example ``EPW/examples/mgb2/epw/epw.in`` (tag qe-6.7.0 = EPW 5.3), which is
the reference anisotropic run (ephwrite, eliashberg, laniso, limag, iverbosity=2). Key names are
validated against the ``NAMELIST / inputepw /`` statement of ``EPW/src/epw_readin.f90`` at tag
qe-6.7.0 (list below); other EPW versions renamed some variables, so pass ``version=None`` to skip
the check deliberately.

Workflow EPW expects (from the example tree)::

    pw.x scf  -> ph.x (ldisp, fildvscf, coarse nq)  -> collect_ph_for_epw()   [= EPW/bin/pp.py]
    pw.x scf  -> pw.x nscf on the FULL coarse k grid (nk1..3, explicit list, see
                 nscf_kpoints_card) with nbnd >= nbndsub + skipped bands -> epw.x

Material-specific Wannierisation input (nbndsub, proj, windows, bands_skipped) cannot be chosen
generically and MUST be given by the caller; the writer refuses to guess it.
Outputs used here: <prefix>.lambda_FS (per-state lambda_nk), <prefix>.a2f (omega_max via the
grid end), see band_coupling.py.
"""
from __future__ import annotations

import shutil
from pathlib import Path
from typing import Any, Sequence

from .qe_inputs import fortran_value

#: NAMELIST / inputepw / of EPW 5.3 (EPW/src/epw_readin.f90, tag qe-6.7.0), copied from source
EPW53_KEYS = frozenset("""
amass outdir prefix iverbosity fildvscf elph nq1 nq2 nq3 nk1 nk2 nk3 nbndsub filukk epbread
epbwrite epwread epwwrite etf_mem kmaps eig_read wepexst epexst vme degaussw fsthick nsmear
delta_smear dvscf_dir ngaussw epmatkqread selecqread wannierize dis_win_max dis_win_min
dis_froz_min dis_froz_max num_iter proj bands_skipped wdata iprint write_wfn wmin wmax nw
eps_acustic a2f nest_fn plselfen elecselfen phonselfen use_ws nc rand_q rand_nq rand_k rand_nk
specfun_pl nqf1 nqf2 nqf3 nkf1 nkf2 nkf3 mp_mesh_k mp_mesh_q filqf filkf ephwrite band_plot
fermi_plot degaussq delta_qsmear nqsmear nqstep nswfc nswc nswi pwc wsfc wscut system_2d
broyden_beta broyden_ndim nstemp temps conv_thr_raxis conv_thr_iaxis conv_thr_racon gap_edge
nsiter muc lreal limag lpade lacon liso laniso lpolar npade lscreen scr_typ fermi_diff smear_rpa
cumulant bnd_cum lifc asr_typ lunif kerwrite kerread imag_read eliashberg ep_coupling fila2f
max_memlt efermi_read fermi_energy specfun_el specfun_ph wmin_specfun wmax_specfun nw_specfun
delta_approx scattering int_mob scissor ncarrier carrier iterative_bte scattering_serta
scattering_0rta longrange shortrange scatread restart restart_step restart_filq prtgkk nel meff
epsiheg lphase omegamin omegamax omegastep n_r lindabs mob_maxiter auto_projections scdm_proj
scdm_entanglement scdm_mu scdm_sigma assume_metal wannier_plot wannier_plot_list reduce_unk
wannier_plot_supercell wannier_plot_scale wannier_plot_radius fixsym epw_no_t_rev epw_tr
epw_nosym epw_noinv epw_crysym
""".split())


def _base(key: str) -> str:
    return key.split("(")[0].strip()


def write_epw_aniso_input(*, prefix: str, masses: Sequence[float], nk: Sequence[int],
                          nq: Sequence[int], nkf: Sequence[int], nqf: Sequence[int],
                          nbndsub: int, proj: Sequence[str], dvscf_dir: str = "./save",
                          outdir: str = "./", bands_skipped: str | None = None,
                          dis_win_max: float | None = None, dis_win_min: float | None = None,
                          dis_froz_max: float | None = None, dis_froz_min: float | None = None,
                          num_iter: int = 500, fsthick_eV: float = 0.4, degaussw_eV: float = 0.10,
                          degaussq_meV: float = 0.5, nqstep: int = 500,
                          eps_acustic_cm1: float = 2.0, muc: float = 0.16,
                          temps_K: Sequence[float] = (15.0,), wscut_eV: float = 1.0,
                          nsiter: int = 500, conv_thr_iaxis: float = 1e-4,
                          lpade: bool = True, mp_mesh_k: bool = True,
                          extra: dict[str, Any] | None = None, version: str | None = "5.3") -> str:
    """&inputepw for: Wannierisation, ephwrite, anisotropic imaginary-axis Eliashberg, iverbosity=2
    (writes .lambda_FS, .lambda_aniso, .a2f). Units follow EPW: fsthick/degaussw/wscut in eV,
    degaussq in meV, eps_acustic in cm^-1, temps in K."""
    nk, nq, nkf, nqf = (list(map(int, g)) for g in (nk, nq, nkf, nqf))
    for name, g in (("nk", nk), ("nq", nq), ("nkf", nkf), ("nqf", nqf)):
        if len(g) != 3 or min(g) < 1:
            raise ValueError(f"{name} must be three positive integers")
    for i in range(3):
        if nkf[i] % nqf[i]:
            # EPW/src/supercond.f90 eliashberg_init: 'eliashberg requires nkf1,nkf2,nkf3 to be
            # multiple of nqf1,nqf2,nqf3 when fila2f is not used'
            raise ValueError("nkf must be a multiple of nqf for eliashberg (EPW eliashberg_init)")
    if not proj or nbndsub < 1:
        raise ValueError("Wannier projections and nbndsub are material specific and required")
    v: dict[str, Any] = {"prefix": prefix, "outdir": outdir}
    for i, m in enumerate(masses, start=1):
        v[f"amass({i})"] = float(m)
    v.update({"elph": True, "kmaps": False, "epbwrite": True, "epbread": False,
              "epwwrite": True, "epwread": False, "etf_mem": 1,
              "nbndsub": int(nbndsub), "bands_skipped": bands_skipped,
              "wannierize": True, "num_iter": int(num_iter),
              "dis_win_max": dis_win_max, "dis_win_min": dis_win_min,
              "dis_froz_max": dis_froz_max, "dis_froz_min": dis_froz_min})
    for i, p in enumerate(proj, start=1):
        v[f"proj({i})"] = str(p)
    v.update({"iverbosity": 2, "eps_acustic": float(eps_acustic_cm1), "ephwrite": True,
              "fsthick": float(fsthick_eV), "degaussw": float(degaussw_eV), "nsmear": 1,
              "degaussq": float(degaussq_meV), "nqstep": int(nqstep),
              "eliashberg": True, "laniso": True, "limag": True, "lpade": bool(lpade),
              "conv_thr_iaxis": float(conv_thr_iaxis), "wscut": float(wscut_eV),
              "nstemp": len(temps_K), "nsiter": int(nsiter), "muc": float(muc),
              "dvscf_dir": dvscf_dir,
              "nk1": nk[0], "nk2": nk[1], "nk3": nk[2], "nq1": nq[0], "nq2": nq[1], "nq3": nq[2],
              "mp_mesh_k": bool(mp_mesh_k),
              "nkf1": nkf[0], "nkf2": nkf[1], "nkf3": nkf[2],
              "nqf1": nqf[0], "nqf2": nqf[1], "nqf3": nqf[2]})
    for i, t in enumerate(temps_K, start=1):
        v[f"temps({i})"] = float(t)
    v.update(extra or {})
    if version == "5.3":
        bad = sorted({_base(k) for k in v} - EPW53_KEYS)
        if bad:
            raise ValueError(f"not in the EPW 5.3 inputepw namelist: {bad}")
    elif version is not None:
        raise ValueError("only version='5.3' key lists are bundled; pass version=None to skip")
    lines = ["--", "&inputepw"]
    for k, val in v.items():
        if val is None:
            continue
        lines.append(f"  {k} = {fortran_value(val)}")
    lines.append(" /")
    return "\n".join(lines) + "\n"


def nscf_kpoints_card(grid: Sequence[int]) -> str:
    """Full unshifted grid in crystal coordinates, weights 1/N, z fastest -- the list produced by
    EPW/Wannier90's kmesh.pl and used in the EPW examples' nscf.in."""
    n1, n2, n3 = (int(x) for x in grid)
    N = n1 * n2 * n3
    out = ["K_POINTS crystal", str(N)]
    for i in range(n1):
        for j in range(n2):
            for k in range(n3):
                out.append(f"  {i / n1:.8f} {j / n2:.8f} {k / n3:.8f} {1.0 / N:.6e}")
    return "\n".join(out) + "\n"


def collect_ph_for_epw(phdir: str | Path, prefix: str, dest: str = "save") -> Path:
    """Non-interactive port of EPW/bin/pp.py (qe-6.7.0): copy dyn files, dvscf and phsave of a
    ph.x run (outdir = phdir) into ``phdir/dest`` with EPW's names. Handles XML/non-XML dyn files,
    sequential (``prefix.dvscf``) and parallel (``prefix.dvscf1``) runs, PAW, and .fc files."""
    d = Path(phdir)
    ph0 = d / "_ph0"
    ctrl = (ph0 / f"{prefix}.phsave" / "control_ph.xml").read_text().splitlines()
    i = next(n for n, l in enumerate(ctrl) if "NUMBER_OF_Q_POINTS" in l)
    nqpt = int(ctrl[i + 1].strip())
    xml = (d / f"{prefix}.dyn1.xml").is_file()
    if not xml and not (d / f"{prefix}.dyn1").is_file():
        raise FileNotFoundError("no dyn1 file found; cannot tell if xml format was used")
    schema = (d / f"{prefix}.save" / "data-file-schema.xml")
    paw = schema.is_file() and "<paw>true</paw>" in schema.read_text().replace(" ", "")
    seq = (ph0 / f"{prefix}.dvscf").is_file()
    suf = "" if seq else "1"
    out = d / dest
    out.mkdir(exist_ok=True)
    if xml:
        shutil.copy(d / f"{prefix}.dyn0", d / f"{prefix}.dyn0.xml")
    for iq in range(1, nqpt + 1):
        ext = ".xml" if xml else ""
        shutil.copy(d / f"{prefix}.dyn{iq}{ext}", out / f"{prefix}.dyn_q{iq}{ext}")
        src = ph0 if iq == 1 else ph0 / f"{prefix}.q_{iq}"
        shutil.copy(src / f"{prefix}.dvscf{suf}", out / f"{prefix}.dvscf_q{iq}")
        if paw:
            shutil.copy(src / f"{prefix}.dvscf_paw{suf}", out / f"{prefix}.dvscf_paw_q{iq}")
        if iq == 1:
            if (out / f"{prefix}.phsave").exists():
                shutil.rmtree(out / f"{prefix}.phsave")
            shutil.copytree(ph0 / f"{prefix}.phsave", out / f"{prefix}.phsave")
            for fc, tgt in ((f"{prefix}.fc.xml", "ifc.q2r.xml"), (f"{prefix}.fc", "ifc.q2r")):
                if (d / fc).is_file():
                    shutil.copy(d / fc, out / tgt)
                    break
    return out
