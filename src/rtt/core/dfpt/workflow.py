"""Quantum ESPRESSO interpolated electron-phonon workflow for ONE material.

Steps (all run in ``workdir``; QE scratch in ``workdir/tmp``)::

    relax    pw.x  vc-relax at press = 10 * pressure_GPa kbar          (optional)
    scf_fit  pw.x  scf on the DENSE k grid with la2F=.true. -> tmp/<prefix>.a2Fsave
    scf      pw.x  scf on the coarse k grid (overwrites <prefix>.save, keeps a2Fsave)
    ph       ph.x  electron_phonon='interpolated', ldisp, trans, el_ph_sigma/nsigma
    q2r      q2r.x la2F=.true.     -> ifc.fc, elph_dir/a2Fmatdyn.61..
    matdyn   matdyn.x dos=.true., la2F=.true. -> a2F.dos1..N, lambda

The dense-before-coarse order is the one of QE's own example (PHonon/examples/example03).

Resumability: a step is skipped when (a) its recorded token (hash of its input text and of the
tokens of the steps it depends on) matches and (b) its output products pass the completion
check (``JOB DONE`` plus the files it must produce). If the state file is missing but the
outputs pass and the existing ``<step>.in`` is byte-identical to what would be written now,
the step is adopted. Rerunning a step invalidates everything downstream (token chain).
A ph.x step interrupted with the same token is restarted with ``recover=.true.``.
"""
from __future__ import annotations

import hashlib
import json
import re
import shutil
import subprocess
import time
from dataclasses import asdict, dataclass, field
from pathlib import Path
from typing import Any

import numpy as np

from . import qe_inputs as qi
from .checks import Tolerances, check_qe_spectrum, spread_report
from .grids import EPhGrids, select_grids
from .launcher import Launcher, RunRecord, run_command
from .parsers import (QEA2FDos, mode_negative_fractions, parse_qe_a2f_dos, parse_qe_lambda_file,
                      parse_qe_version, qe_a2f_dos_to_spectrum)
from .pseudo import PseudoTable
from .structure import CellData, parse_pw_final_structure, to_celldata
from .units import GPA_TO_KBAR, to_kelvin_factor

STEPS = ("relax", "scf_fit", "scf", "ph", "q2r", "matdyn")
CODES = {"relax": "pw.x", "scf_fit": "pw.x", "scf": "pw.x", "ph": "ph.x", "q2r": "q2r.x",
         "matdyn": "matdyn.x"}
UPSTREAM = {"relax": (), "scf_fit": ("relax",), "scf": ("scf_fit",), "ph": ("scf",),
            "q2r": ("ph",), "matdyn": ("q2r",)}


class WorkflowError(RuntimeError):
    def __init__(self, msg: str, record: RunRecord | None = None, result: Any = None):
        super().__init__(msg)
        self.record = record
        self.result = result


@dataclass
class EPhSettings:
    pressure_GPa: float = 0.0
    do_relax: bool = True
    prefix: str = "pwscf"
    # plane waves / electrons
    ecutwfc_Ry: float | None = None          # None: max over SSSP/UPF suggestions
    ecutrho_Ry: float | None = None
    smearing: str = "mv"
    degauss_Ry: float = 0.02
    conv_thr: float = 1e-10
    mixing_beta: float = 0.5
    # relax
    press_conv_thr_kbar: float = 0.5
    etot_conv_thr: float = 1e-5
    forc_conv_thr: float = 1e-4
    relax_nstep: int = 100
    cell_dofree: str = "all"
    # grids
    k_spacing: float = 0.25                  # 1/angstrom incl. 2pi
    q_spacing: float = 0.8
    fit_factor: int = 2
    k_fit_spacing: float | None = None
    matdyn_spacing: float | None = None
    k_scf: tuple[int, int, int] | None = None
    k_fit: tuple[int, int, int] | None = None
    q: tuple[int, int, int] | None = None
    k_matdyn: tuple[int, int, int] | None = None
    # e-ph
    el_ph_sigma_Ry: float = 0.005            # broadenings are isig * el_ph_sigma, isig=1..nsigma
    el_ph_nsigma: int = 10
    use_smearings: list[int] | None = None   # None: every smearing is an alternative spectrum
    # "mode_positive": admit only the contiguous block of smearings, starting at the smallest one
    # whose per-mode alpha^2F is non-negative (see admit_smearings); "all": no admission rule
    smearing_admission: str = "all"
    tr2_ph: float = 1e-14
    alpha_mix: float | None = None
    fildvscf: str | None = "dvscf"
    asr: str = "crystal"
    ndos: int | None = 200
    deltaE_cm1: float | None = None
    # running
    timeouts_s: dict[str, float] = field(default_factory=dict)  # per step; missing = no limit
    adopt_existing: bool = True
    extra_pw_system: dict = field(default_factory=dict)
    extra_pw_control: dict = field(default_factory=dict)
    extra_ph: dict = field(default_factory=dict)

    @classmethod
    def from_config(cls, cfg: dict[str, Any] | None) -> "EPhSettings":
        cfg = dict(cfg or {})
        if "press_kbar" in cfg:
            cfg["pressure_GPa"] = float(cfg.pop("press_kbar")) / GPA_TO_KBAR
        for k in ("k_scf", "k_fit", "q", "k_matdyn"):
            if cfg.get(k) is not None:
                cfg[k] = tuple(int(x) for x in cfg[k])
        known = set(cls.__dataclass_fields__)
        unknown = set(cfg) - known
        if unknown:
            raise ValueError(f"unknown dfpt settings: {sorted(unknown)}")
        return cls(**cfg)

    @property
    def press_kbar(self) -> float:
        return float(self.pressure_GPa) * GPA_TO_KBAR

    def to_dict(self) -> dict[str, Any]:
        return asdict(self)


@dataclass
class EPhResult:
    material_id: str
    spectra: list                    # list[rtt.core.spectrum.Spectrum], one per smearing kept
    checks: list                     # list[ConsistencyCheck]
    spread: dict
    warnings: list[str]
    flags: list[str]
    provenance: dict
    lambda_entries: list = field(default_factory=list)
    errors: list[str] = field(default_factory=list)   # fatal: some alternative unusable

    @property
    def checks_ok(self) -> bool:
        return all(c.ok for c in self.checks)

    @property
    def ok(self) -> bool:
        """False if any requested smearing could not be turned into a Spectrum. Callers must then
        NOT use the remaining spectra alone (S5: a verdict needs every alternative) -> ERROR."""
        return not self.errors

    def to_dict(self) -> dict[str, Any]:
        return {"material_id": self.material_id,
                "spectra": [s.to_dict() for s in self.spectra],
                "checks": [c.to_dict() for c in self.checks],
                "checks_ok": self.checks_ok, "ok": self.ok, "errors": list(self.errors),
                "spread": self.spread, "warnings": list(self.warnings), "flags": list(self.flags),
                "provenance": self.provenance,
                "lambda_entries": [asdict(e) for e in self.lambda_entries]}


def _sha(s: str) -> str:
    return hashlib.sha256(s.encode()).hexdigest()


def _json_default(o: Any) -> Any:
    if isinstance(o, (np.floating, np.integer)):
        return o.item()
    if isinstance(o, np.ndarray):
        return o.tolist()
    if isinstance(o, Path):
        return str(o)
    raise TypeError(type(o))


class EPhWorkflow:
    def __init__(self, material_id: str, structure: Any, workdir: str | Path,
                 settings: EPhSettings, pseudos: PseudoTable, launcher: Launcher | None = None,
                 masses: dict[str, float] | None = None):
        self.material_id = material_id
        self.cell0: CellData = to_celldata(structure, masses)
        self.workdir = Path(workdir)
        self.settings = settings
        self.pseudos = pseudos
        self.launcher = launcher or Launcher()
        self.state_path = self.workdir / "dfpt_state.json"
        self.state: dict[str, Any] = {"steps": {}}
        if self.state_path.exists():
            try:
                self.state = json.loads(self.state_path.read_text())
            except json.JSONDecodeError:
                self.state = {"steps": {}}
        self.state.setdefault("steps", {})
        if settings.use_smearings is not None:
            bad = [i for i in settings.use_smearings if not 1 <= i <= settings.el_ph_nsigma]
            if bad:
                raise ValueError(f"use_smearings {bad} outside 1..{settings.el_ph_nsigma}")
        if settings.smearing_admission not in SMEARING_ADMISSION:
            raise ValueError(f"smearing_admission must be one of {SMEARING_ADMISSION}")
        if settings.smearing_admission != "all" and settings.use_smearings is not None:
            raise ValueError("use_smearings and smearing_admission are mutually exclusive")
        self._grids: EPhGrids | None = None

    # ------------------------------------------------------------------ helpers
    @property
    def outdir(self) -> Path:
        return self.workdir / "tmp"

    def _save_state(self) -> None:
        self.workdir.mkdir(parents=True, exist_ok=True)
        tmp = self.state_path.with_suffix(".json.tmp")
        tmp.write_text(json.dumps(self.state, indent=1, default=_json_default))
        tmp.replace(self.state_path)

    def steps(self) -> tuple[str, ...]:
        return STEPS if self.settings.do_relax else STEPS[1:]

    def _upstream(self, step: str) -> tuple[str, ...]:
        ups = UPSTREAM[step]
        return tuple(u for u in ups if u in self.steps())

    def working_cell(self) -> CellData:
        """Relaxed cell if the relax step is complete, else (no relax) the input cell."""
        if not self.settings.do_relax:
            return self.cell0
        p = self.workdir / "relaxed_structure.json"
        if not p.exists():
            raise WorkflowError("relaxed structure not available yet (run the relax step)")
        return CellData.from_dict(json.loads(p.read_text()))

    def grids(self, cell: CellData | None = None) -> EPhGrids:
        cell = cell or (self.working_cell() if self._relax_complete_or_skipped() else self.cell0)
        s = self.settings
        return select_grids(cell.reciprocal_lengths_2pi(), k_spacing=s.k_spacing,
                            q_spacing=s.q_spacing, fit_factor=s.fit_factor,
                            k_fit_spacing=s.k_fit_spacing, matdyn_spacing=s.matdyn_spacing,
                            k_scf=s.k_scf, k_fit=s.k_fit, q=s.q, k_matdyn=s.k_matdyn)

    def _relax_complete_or_skipped(self) -> bool:
        return (not self.settings.do_relax) or (self.workdir / "relaxed_structure.json").exists()

    def _cutoffs(self, cell: CellData) -> tuple[float, float | None]:
        w, r = self.pseudos.cutoffs(cell.elements)
        ecw = self.settings.ecutwfc_Ry or w
        ecr = self.settings.ecutrho_Ry or r
        if not ecw:
            raise WorkflowError("ecutwfc unknown: set dfpt.ecutwfc_Ry or provide SSSP cutoffs")
        return float(ecw), (float(ecr) if ecr else None)

    # ------------------------------------------------------------------ inputs
    def input_text(self, step: str, recover: bool = False) -> str:
        s = self.settings
        cell = self.cell0 if step == "relax" else self.working_cell()
        ecw, ecr = self._cutoffs(cell)
        files = {el: self.pseudos.files[el] for el in cell.elements}
        common = dict(cell=cell, prefix=s.prefix, outdir="./tmp/",
                      pseudo_dir=str(Path(self.pseudos.directory).resolve()) + "/",
                      pseudo_files=files, ecutwfc=ecw, ecutrho=ecr, smearing=s.smearing,
                      degauss=s.degauss_Ry, conv_thr=s.conv_thr, mixing_beta=s.mixing_beta,
                      extra_system=s.extra_pw_system, extra_control=s.extra_pw_control)
        if step == "relax":
            g = self.grids(self.cell0)
            return qi.write_pw_input(calculation="vc-relax", kgrid=g.k_scf, press_kbar=s.press_kbar,
                                     press_conv_thr_kbar=s.press_conv_thr_kbar,
                                     etot_conv_thr=s.etot_conv_thr, forc_conv_thr=s.forc_conv_thr,
                                     nstep=s.relax_nstep, cell_dofree=s.cell_dofree, **common)
        g = self.grids(cell)
        if step == "scf_fit":
            return qi.write_pw_input(calculation="scf", kgrid=g.k_fit, la2F=True, **common)
        if step == "scf":
            return qi.write_pw_input(calculation="scf", kgrid=g.k_scf, la2F=False, **common)
        if step == "ph":
            return qi.write_ph_input(prefix=s.prefix, outdir="./tmp/", cell=cell, qgrid=g.q,
                                     el_ph_sigma=s.el_ph_sigma_Ry, el_ph_nsigma=s.el_ph_nsigma,
                                     tr2_ph=s.tr2_ph, fildyn="dyn", fildvscf=s.fildvscf,
                                     recover=recover, alpha_mix=s.alpha_mix, extra=s.extra_ph,
                                     title=f"e-ph {self.material_id}")
        if step == "q2r":
            return qi.write_q2r_input(fildyn="dyn", flfrc="ifc.fc", zasr=s.asr,
                                      el_ph_nsigma=s.el_ph_nsigma)
        if step == "matdyn":
            return qi.write_matdyn_dos_input(cell=cell, flfrc="ifc.fc", kgrid=g.k_matdyn,
                                             asr=s.asr, ndos=s.ndos, deltaE_cm1=s.deltaE_cm1,
                                             el_ph_nsigma=s.el_ph_nsigma)
        raise ValueError(step)

    def token(self, step: str) -> str:
        parts = [self.input_text(step)]
        for u in self._upstream(step):
            parts.append(self.state["steps"].get(u, {}).get("token", "MISSING"))
        return _sha("\n#--#\n".join(parts))

    # ------------------------------------------------------------------ completion checks
    def products_ok(self, step: str) -> tuple[bool, str]:
        w = self.workdir
        out = w / f"{step}.out"
        if not out.exists():
            return False, "no output"
        text = out.read_text(errors="replace")
        if "JOB DONE" not in text:
            return False, "no 'JOB DONE' in output"
        s = self.settings
        if step == "relax":
            if "End final coordinates" not in text:
                return False, "vc-relax did not reach 'End final coordinates' (not converged)"
            if not (w / "relaxed_structure.json").exists():
                return False, "relaxed_structure.json missing"
            return True, ""
        if step == "scf_fit":
            ok = (self.outdir / f"{s.prefix}.a2Fsave").exists()
            return ok, "" if ok else f"tmp/{s.prefix}.a2Fsave missing"
        if step == "scf":
            if "convergence has been achieved" not in text:
                return False, "scf not converged"
            return True, ""
        if step == "ph":
            dyn0 = w / "dyn0"
            if not dyn0.exists():
                return False, "dyn0 missing"
            try:
                nqirr = int(dyn0.read_text().split("\n")[1].split()[0])
            except (IndexError, ValueError):
                return False, "dyn0 unreadable"
            for iq in range(1, nqirr + 1):
                if not (w / f"dyn{iq}").exists():
                    return False, f"dyn{iq} missing"
                for isig in range(1, s.el_ph_nsigma + 1):
                    if not (w / "elph_dir" / f"a2Fq2r.{50 + isig}.{iq}").exists():
                        return False, f"elph_dir/a2Fq2r.{50 + isig}.{iq} missing"
            return True, ""
        if step == "q2r":
            if not (w / "ifc.fc").exists():
                return False, "ifc.fc missing"
            for isig in range(1, s.el_ph_nsigma + 1):
                # q2r writes elph_dir/a2Fmatdyn.<60+isig>; matdyn.f90 a2Fdos opens the same
                if not (w / "elph_dir" / f"a2Fmatdyn.{60 + isig}").exists():
                    return False, f"elph_dir/a2Fmatdyn.{60 + isig} missing"
            return True, ""
        if step == "matdyn":
            for isig in range(1, s.el_ph_nsigma + 1):
                if not (w / f"a2F.dos{isig}").exists():
                    return False, f"a2F.dos{isig} missing"
            lf = w / "lambda"
            if not lf.exists() or len(parse_qe_lambda_file(lf)) != s.el_ph_nsigma:
                return False, "lambda file missing or incomplete"
            return True, ""
        raise ValueError(step)

    def is_complete(self, step: str) -> tuple[bool, str]:
        try:
            tok = self.token(step)
        except WorkflowError as e:
            return False, str(e)
        rec = self.state["steps"].get(step, {})
        ok, why = self.products_ok(step)
        if not ok:
            return False, why
        if rec.get("token") == tok and rec.get("status") == "done":
            return True, "done (state)"
        if self.settings.adopt_existing and not rec:
            # adopt only when there is NO record at all (an invalidated record blocks it), the
            # input is identical, upstream is done and our output is newer than upstream's
            inp = self.workdir / f"{step}.in"
            out_t = (self.workdir / f"{step}.out").stat().st_mtime
            newer = all((self.workdir / f"{u}.out").exists()
                        and (self.workdir / f"{u}.out").stat().st_mtime <= out_t
                        for u in self._upstream(step))
            if inp.exists() and inp.read_text() == self.input_text(step) and newer \
                    and all(self.state["steps"].get(u, {}).get("status") == "done"
                            for u in self._upstream(step)):
                self.state["steps"][step] = {"token": tok, "status": "done", "adopted": True,
                                             "time": time.time()}
                self._save_state()
                return True, "adopted existing outputs"
        return False, "inputs or upstream changed since last run"

    def status(self) -> dict[str, str]:
        out = {}
        blocked = False
        for st in self.steps():
            if blocked:
                out[st] = "pending (upstream incomplete)"
                continue
            ok, why = self.is_complete(st)
            out[st] = "done" if ok else f"pending: {why}"
            blocked = not ok
        return out

    # ------------------------------------------------------------------ running
    def plan(self) -> list[dict[str, Any]]:
        """Commands that ``run`` would execute (inputs of later steps may not exist yet)."""
        out = []
        for st in self.steps():
            ok, why = self.is_complete(st)
            out.append({"step": st, "complete": ok, "reason": why,
                        "command": self.launcher.command(CODES[st], f"{st}.in")})
        return out

    def _after_success(self, step: str) -> None:
        if step == "relax":
            text = (self.workdir / "relax.out").read_text(errors="replace")
            cell = parse_pw_final_structure(text, self.cell0)
            (self.workdir / "relaxed_structure.json").write_text(
                json.dumps(cell.to_dict(), indent=1))
            info = parse_relax_summary(text, self.settings.press_kbar)
            self.state["relax_summary"] = info

    def run_step(self, step: str) -> RunRecord:
        self.workdir.mkdir(parents=True, exist_ok=True)
        (self.workdir / "tmp").mkdir(exist_ok=True)
        tok = self.token(step)
        rec = self.state["steps"].get(step, {})
        recover = (step == "ph" and rec.get("token") == tok and rec.get("status") != "done"
                   and (self.outdir / "_ph0").exists())
        text = self.input_text(step, recover=recover)
        (self.workdir / f"{step}.in").write_text(text)
        if step == "relax" and (self.workdir / "relaxed_structure.json").exists():
            (self.workdir / "relaxed_structure.json").unlink()
        self.state["steps"][step] = {"token": tok, "status": "running", "time": time.time(),
                                     "recover": recover}
        # invalidate downstream records (an 'invalidated' record also blocks adoption)
        seen = False
        for st in self.steps():
            if seen:
                self.state["steps"][st] = {"status": "invalidated", "by": step,
                                           "time": time.time()}
            if st == step:
                seen = True
        self._save_state()
        cmd = self.launcher.command(CODES[step], f"{step}.in")
        r = run_command(step, cmd, self.workdir, env=self.launcher.environment(),
                        timeout_s=self.settings.timeouts_s.get(step))
        status = "failed"
        if r.ok:
            try:
                self._after_success(step)
                ok, why = self.products_ok(step)
            except Exception as exc:  # parse problems become step failures
                ok, why = False, f"post-processing failed: {exc}"
            if ok:
                status = "done"
            else:
                r.error_message = (r.error_message + " || " if r.error_message else "") + why
        self.state["steps"][step].update({"status": status, "record": r.to_dict()})
        self._save_state()
        return r

    def run(self, until: str | None = None, collect: bool = True) -> EPhResult | None:
        """Run every incomplete step in order (stop after ``until``). Raises WorkflowError
        with the captured diagnosis on the first failing step."""
        if until is not None and until not in self.steps():
            raise ValueError(f"unknown step {until!r}")
        prov_ps = self.pseudos.validate(self.cell0.elements)
        self.state["pseudopotentials"] = prov_ps
        self._save_state()
        for st in self.steps():
            ok, why = self.is_complete(st)
            if not ok:
                r = self.run_step(st)
                if self.state["steps"][st]["status"] != "done":
                    raise WorkflowError(f"[{self.material_id}] step {st} failed "
                                        f"(rc={r.returncode}, timeout={r.timed_out}): "
                                        f"{r.error_message}", r)
            if st == until:
                return None
        if not collect:
            return None
        res = self.collect()
        if not res.ok:
            raise WorkflowError(f"[{self.material_id}] QE finished but some spectra are unusable: "
                                + "; ".join(res.errors), None, res)
        return res

    # ------------------------------------------------------------------ results
    def collect(self, tol: Tolerances | None = None) -> EPhResult:
        s = self.settings
        cell = self.working_cell()
        g = self.grids(cell)
        version = None
        mo = self.workdir / "matdyn.out"
        if mo.exists():
            version = parse_qe_version(mo.read_text(errors="replace"))
        prov = {"material_id": self.material_id, "synthetic": False,
                "code": "Quantum ESPRESSO", "code_version": version,
                "workdir": str(self.workdir.resolve()), "pressure_GPa": s.pressure_GPa,
                "relaxed": s.do_relax, "grids": g.to_dict(),
                "el_ph_sigma_Ry": s.el_ph_sigma_Ry, "el_ph_nsigma": s.el_ph_nsigma,
                "degauss_Ry": s.degauss_Ry, "smearing": s.smearing,
                "ecut_Ry": list(self._cutoffs(cell)), "ndos": s.ndos,
                "deltaE_cm1": s.deltaE_cm1,
                "pseudopotentials": self.state.get("pseudopotentials"),
                "relax_summary": self.state.get("relax_summary"),
                "cell": cell.to_dict()}
        ro = self.workdir / "relax.out"
        if s.do_relax and ro.exists():
            self.state["relax_summary"] = parse_relax_summary(ro.read_text(errors="replace"),
                                                              s.press_kbar)
            prov["relax_summary"] = self.state["relax_summary"]
        res = collect_qe_eph(self.workdir, self.material_id, el_ph_nsigma=s.el_ph_nsigma,
                             use_smearings=s.use_smearings, tol=tol, base_provenance=prov,
                             smearing_admission=s.smearing_admission)
        for wmsg in (self.state.get("relax_summary") or {}).get("warnings", []):
            res.warnings.append("relax: " + wmsg)
        ps_warn = (self.state.get("pseudopotentials") or {}).get("_warnings", [])
        res.warnings.extend("pseudopotentials: " + x for x in ps_warn)
        if not res.checks_ok:
            res.warnings.append("WARNING: consistency checks failed: "
                                + ", ".join(f for f in res.flags if f.startswith("DFPT_CHECK_FAIL")))
        return res


SMEARING_ADMISSION = ("all", "mode_positive")
#: per-mode negative mass tolerated as print/rounding noise (same scale as atoms_to_spectrum's)
MODE_NEGATIVE_TOL_REL = 1e-6


def admit_smearings(tabs: dict[int, QEA2FDos], tol_rel: float = MODE_NEGATIVE_TOL_REL
                    ) -> tuple[list[int], dict[str, Any], str | None]:
    """Rule 'mode_positive': smearing i passes iff every mode's alpha^2F_nu is non-negative at
    omega > 0 (within tol_rel of the total mass; see parsers.mode_negative_fractions). The admitted
    set is the contiguous block from the smallest passing smearing upwards (as in protocol P3).

    Failing smearings below that block are excluded: they are too narrow for the k grid, which is
    the only way per-mode negativity arises from a monotone ladder. A failure ABOVE a passing
    smearing has no such explanation, so it is an error, as is no smearing passing at all.
    Returns (admitted, report, error); error is None on success.
    """
    order = sorted(tabs)
    neg = {i: mode_negative_fractions(tabs[i]) for i in order}
    passing = [i for i in order if max(neg[i], default=0.0) <= tol_rel]
    report: dict[str, Any] = {
        "rule": "mode_positive", "tol_rel": tol_rel,
        "max_mode_negative_fraction": {str(i): max(neg[i], default=0.0) for i in order},
        "lambda_footer": {str(i): tabs[i].lambda_reported for i in order}}
    if not passing:
        report.update(admitted=[], excluded=order)
        return [], report, ("smearing admission: no smearing has non-negative per-mode alpha^2F "
                            f"(max negative fraction per smearing {report['max_mode_negative_fraction']})")
    block = [i for i in order if i >= passing[0]]
    holes = [i for i in block if i not in passing]
    excluded = [i for i in order if i < passing[0]]
    report.update(admitted=block, excluded=excluded)
    if holes:
        return [], report, (f"smearing admission: smearings {holes} fail per-mode positivity above "
                            f"the passing smearing {passing[0]} (non-contiguous; not explained by "
                            "under-resolved smearing)")
    return block, report, None


def _exclusion_tag(report: dict[str, Any], tabs: dict[int, QEA2FDos]) -> str | None:
    """Assumption carried by every admitted spectrum when smearings were excluded:
    'A-SIGMA-EXCLUDED(sigma1:lam=<l>:A=<a>K,...)' with the excluded spectra's own (signed) lambda
    and A. rtt.core.certify_policy.blocking_tags compares them with the admitted alternatives' envelope
    boxes and blocks the proof direction(s) they fall outside of."""
    ex = report.get("excluded") or []
    if not ex:
        return None
    f = to_kelvin_factor("Ry")
    parts = []
    for i in ex:
        t = tabs[i]
        pos = t.E_Ry > 0
        lam = float(2.0 * np.sum(t.a2f_total[pos] / t.E_Ry[pos]) * t.delta_Ry)
        A = float(np.sum(t.a2f_total[pos]) * t.delta_Ry * f)
        parts.append(f"sigma{i}:lam={lam!r}:A={A!r}K")
    return f"A-SIGMA-EXCLUDED({','.join(parts)})"


def collect_qe_eph(workdir: str | Path, material_id: str, *, el_ph_nsigma: int | None = None,
                   use_smearings: list[int] | None = None, tol: Tolerances | None = None,
                   base_provenance: dict | None = None, imaginary_tol_cm1: float = 1.0,
                   smearing_admission: str = "all", **spectrum_kw: Any) -> EPhResult:
    """Parse a finished QE matdyn directory (a2F.dos*, lambda, freq.dos) into an EPhResult.

    Usable on runs made outside this workflow (``base_provenance`` then documents them).
    """
    w = Path(workdir)
    lam_entries = parse_qe_lambda_file(w / "lambda") if (w / "lambda").exists() else []
    if el_ph_nsigma is None:
        el_ph_nsigma = len(sorted(w.glob("a2F.dos*")))
    idx = list(use_smearings) if use_smearings is not None else list(range(1, el_ph_nsigma + 1))
    warnings: list[str] = []
    flags: list[str] = []
    if use_smearings is not None:
        warnings.append(f"only smearing indices {idx} of 1..{el_ph_nsigma} kept (explicit config)")
    spectra, checks = [], []
    base = dict(base_provenance or {})
    errors: list[str] = []
    excl_tag = None
    if smearing_admission not in SMEARING_ADMISSION:
        raise ValueError(f"smearing_admission must be one of {SMEARING_ADMISSION}")
    if smearing_admission == "mode_positive":
        if use_smearings is not None:
            raise ValueError("use_smearings and smearing_admission are mutually exclusive")
        tabs = {}
        for isig in idx:
            try:
                tabs[isig] = parse_qe_a2f_dos(w / f"a2F.dos{isig}")
            except (ValueError, OSError):
                pass        # reported as DFPT_SPECTRUM_ERROR by the main loop below
        admitted, report, err = admit_smearings(tabs)
        base["smearing_admission"] = report
        if err:
            errors.append(err)
            flags.append("SMEARING_ADMISSION_FAIL")
            idx = []
        else:
            unparsed = [i for i in idx if i not in tabs]
            idx = sorted(set(admitted) | set(unparsed))
            if report["excluded"]:
                warnings.append(
                    f"only smearing indices {idx} of 1..{el_ph_nsigma} kept (smearing_admission="
                    f"mode_positive: {['sigma%d' % i for i in report['excluded']]} have negative "
                    "per-mode alpha^2F, i.e. are not resolved by the k grid)")
                excl_tag = _exclusion_tag(report, tabs)
    fmin = None
    fq = w / "freq.dos"
    if fq.exists():
        try:
            fmin = min_frequency_cm1(fq)
        except Exception as exc:  # pragma: no cover - diagnostic only
            warnings.append(f"could not parse freq.dos: {exc}")
    for isig in idx:
        entry = lam_entries[isig - 1] if len(lam_entries) >= isig else None
        b = entry.broadening_Ry if entry else None
        prov = dict(base)
        prov["qe_lambda_file"] = asdict(entry) if entry else None
        if excl_tag:
            prov["assumptions"] = list(prov.get("assumptions", [])) + [excl_tag]
        try:
            tab = parse_qe_a2f_dos(w / f"a2F.dos{isig}")
            sp = qe_a2f_dos_to_spectrum(tab, label=f"{material_id}:qe:sigma{isig}"
                                        + (f"({b:g}Ry)" if b is not None else ""),
                                        broadening_Ry=b, extra_provenance=prov, **spectrum_kw)
        except (ValueError, OSError) as exc:
            errors.append(f"sigma{isig}: {exc}")
            flags.append(f"DFPT_SPECTRUM_ERROR:sigma{isig}")
            continue
        for wmsg in sp.provenance.get("warnings", []):
            warnings.append(f"sigma{isig}: {wmsg}")
            if "omega <= 0" in wmsg and "NEGATIVE_FREQ_A2F_WEIGHT" not in flags:
                flags.append("NEGATIVE_FREQ_A2F_WEIGHT")
        cs = check_qe_spectrum(sp, tab, entry, tol)
        checks.extend(cs)
        for c in cs:
            if not c.ok:
                flags.append(f"DFPT_CHECK_FAIL:{c.name}:sigma{isig}")
        spectra.append(sp)
    if fmin is not None:
        for sp in spectra:
            sp.provenance["min_freq_matdyn_grid_cm1"] = fmin
        if fmin < -abs(imaginary_tol_cm1):
            flags.append(f"IMAGINARY_MODES(min={fmin:.2f}cm-1)")
            warnings.append(f"WARNING: imaginary phonon frequencies on the matdyn grid "
                            f"(min {fmin:.3f} cm^-1)")
            if errors:
                errors.append(f"note: imaginary modes present (min {fmin:.2f} cm^-1); matdyn's "
                              "a2F is unreliable/NaN/zero for dynamically unstable structures")
    prov = dict(base)
    prov["n_alternatives"] = len(spectra)
    if errors:
        warnings.append("ERROR: some smearings gave no usable spectrum; the material must not be "
                        "classified from the remaining ones (ARCH S5)")
    return EPhResult(material_id, spectra, checks, spread_report(spectra), warnings, flags,
                     prov, lam_entries, errors)


def min_frequency_cm1(path: str | Path) -> float:
    """Minimum frequency in a matdyn ``flfrq`` file (&plot header, then q line / freq lines)."""
    vals: list[float] = []
    lines = Path(path).read_text().splitlines()
    m = re.search(r"nbnd=\s*(\d+),\s*nks=\s*(\d+)", lines[0])
    if not m:
        raise ValueError("not a matdyn flfrq file")
    nbnd = int(m.group(1))
    i = 1
    while i < len(lines):
        i += 1  # q-point line
        got: list[float] = []
        while len(got) < nbnd and i < len(lines):
            got.extend(float(x) for x in lines[i].split())
            i += 1
        vals.extend(got)
    return float(min(vals))


def parse_relax_summary(text: str, target_kbar: float | None = None) -> dict[str, Any]:
    info: dict[str, Any] = {}
    info["bfgs_converged"] = "bfgs converged" in text
    info["max_steps_reached"] = "maximum number of steps has been reached" in text
    info["warnings"] = []
    ps = re.findall(r"P=\s*([-+]?\d+\.\d+)", text)
    if ps:
        info["final_scf_pressure_kbar"] = float(ps[-1])
    i_conv = text.find("bfgs converged")
    if i_conv >= 0:
        pc = re.findall(r"P=\s*([-+]?\d+\.\d+)", text[:i_conv])
        if pc:
            info["relaxed_pressure_kbar"] = float(pc[-1])
    if info["max_steps_reached"]:
        info["warnings"].append("vc-relax reached the maximum number of steps")
    if target_kbar is not None and "final_scf_pressure_kbar" in info:
        dp = info["final_scf_pressure_kbar"] - target_kbar
        info["final_scf_pressure_minus_target_kbar"] = dp
        if abs(dp) > max(2.0, 0.05 * abs(target_kbar)):
            info["warnings"].append(
                f"WARNING: final scf after vc-relax gives P={info['final_scf_pressure_kbar']} kbar "
                f"vs target {target_kbar} kbar (Pulay stress: increase ecutwfc or rerun vc-relax)")
    en = re.findall(r"!\s+total energy\s+=\s+([-+]?\d+\.\d+)\s+Ry", text)
    if en:
        info["final_total_energy_Ry"] = float(en[-1])
    h = re.findall(r"Final enthalpy\s+=\s+([-+]?\d+\.\d+)\s+Ry", text)
    if h:
        info["final_enthalpy_Ry"] = float(h[-1])
    v = re.findall(r"unit-cell volume\s+=\s+([\d.]+)\s+\(a\.u\.\)\^3", text)
    if v:
        info["final_volume_bohr3"] = float(v[-1])
    return info


def build_fortify_shim(dest_dir: str | Path, cc: str = "gcc") -> Path:
    """Compile contrib/snprintf_chk_shim.c into ``dest_dir`` (for LD_PRELOAD).

    Only needed for the Ubuntu 24.04 ``quantum-espresso`` 6.7 package, whose md5 helper trips
    glibc's _FORTIFY_SOURCE check at start-up ("buffer overflow detected"). See the C file.
    """
    src = Path(__file__).with_name("contrib") / "snprintf_chk_shim.c"
    dest = Path(dest_dir)
    dest.mkdir(parents=True, exist_ok=True)
    so = dest / "libsnprintf_chk_shim.so"
    if not so.exists() or so.stat().st_mtime < src.stat().st_mtime:
        if shutil.which(cc) is None:
            raise RuntimeError(f"C compiler {cc!r} not found")
        subprocess.run([cc, "-O2", "-shared", "-fPIC", "-o", str(so), str(src)], check=True)
    return so


def workflow_from_config(material_id: str, structure: Any, workdir: str | Path,
                         cfg: dict[str, Any], masses: dict[str, float] | None = None) -> EPhWorkflow:
    """Build an EPhWorkflow from a config mapping (e.g. the ``[dfpt]`` table of config.toml)::

        [dfpt]            -> EPhSettings fields (pressure_GPa or press_kbar, k_spacing, q, ...)
        [dfpt.pseudo]     -> PseudoTable.from_config (dir, sssp_json, files, cutoffs_Ry)
        [dfpt.launcher]   -> Launcher.from_config (kind, nprocs, npool, mpi_args, bin_dir, env)
        [dfpt.shim]       -> optional {build = true, dir = "..."}: LD_PRELOAD fortify shim
    """
    cfg = dict(cfg)
    pseudo_cfg = cfg.pop("pseudo")
    launcher_cfg = dict(cfg.pop("launcher", {}) or {})
    shim_cfg = cfg.pop("shim", None)
    if shim_cfg and shim_cfg.get("build"):
        so = build_fortify_shim(shim_cfg.get("dir", Path(workdir) / ".shim"))
        env = dict(launcher_cfg.get("env", {}))
        env["LD_PRELOAD"] = str(so)
        launcher_cfg["env"] = env
    return EPhWorkflow(material_id, structure, workdir, EPhSettings.from_config(cfg),
                       PseudoTable.from_config(pseudo_cfg), Launcher.from_config(launcher_cfg),
                       masses)
