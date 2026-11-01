"""Profile runner around the core QE workflow (architecture §5.1 ``rtt.dft.workflow``).

``run_profile(struct, profile, workdir, launcher)`` runs one profile for one structure:

* DFPT profiles (``coarse``, ``full``) drive the **frozen core**
  :class:`rtt.core.dfpt.workflow.EPhWorkflow` (relax -> scf_fit -> scf -> ph -> q2r ->
  matdyn) step by step, adding what the architecture asks of the wrapper only:
  - adoption of finished steps: the core token chain (``dfpt_state.json``) and its
    adopt-if-identical-input rule, plus a whole-run marker ``profile_done.json`` that stores the
    output sha256s, so a run whose ``tmp/`` was pruned is adopted without re-running;
  - per-step timing (wall, children user+sys CPU, core-hours) appended to ``timing.jsonl``;
  - per-q scratch pruning during ph.x (hash-then-delete, :mod:`rtt.dft.prune`) and a final
    prune after success;
  - a disk guard before every step and a pause/resume watcher during ph.x;
  - a machine-wide MPI-rank budget: before a parallel step, wait until the number of running
    QE processes plus this step's ranks is within ``max_machine_ranks`` (other agents run QE).
* ``electronic`` (``s4``) and ``hull`` (``ref_hull``) profiles use the small :class:`Chain`
  step runner below, with the same adoption, timing and guard rules
  (:mod:`rtt.dft.electronic`, :mod:`rtt.dft.refhull`).

Nothing here raises on a QE failure: the returned :class:`DFPTRun` has ``status = "failed"``
and the core diagnosis in ``error``; the caller records ``ERROR`` with a retry count.
"""
from __future__ import annotations

import contextlib
import hashlib
import json
import os
import resource
import subprocess
import time
from dataclasses import asdict, dataclass, field
from pathlib import Path
from typing import Any, Callable

from .profiles import Profile, pseudo_md5s, qe_bin, sg15_pseudos
from .prune import (DiskSpaceLow, ScratchWatcher, disk_guard, prune_final, sha256_file)

QE_PROCS = ("pw.x", "ph.x", "q2r.x", "matdyn.x", "pp.x", "dos.x", "epw.x", "projwfc.x")
DONE_MARKER = "profile_done.json"
TIMING_LOG = "timing.jsonl"


# ---------------------------------------------------------------------- launcher / ranks
def default_launcher(ranks: int = 2, *, npool: int | None = None, nice: int = 19) -> Any:
    """mpirun launcher for QE 7.6 (``~/micromamba/envs/qe/bin``), ``nice 19``, one OpenMP
    thread per rank (CONTEXT.md compute rules; architecture §7.3)."""
    from rtt.core.dfpt.launcher import Launcher
    b = qe_bin()
    return Launcher(kind="mpirun" if ranks > 1 else "plain", nprocs=int(ranks),
                    npool=npool if npool is not None else int(ranks),
                    mpirun=str(b / "mpirun"), mpi_args=["--bind-to", "none"], bin_dir=str(b),
                    nice=nice, omp_threads=1,
                    env={"OPENBLAS_NUM_THREADS": "1", "MKL_NUM_THREADS": "1"})


def running_qe_ranks() -> int:
    """Number of running QE executable processes on this machine (all users)."""
    n = 0
    for code in QE_PROCS:
        try:
            out = subprocess.run(["pgrep", "-c", "-x", code], capture_output=True, text=True)
            n += int(out.stdout.strip() or 0)
        except (OSError, ValueError):
            pass
    return n


def wait_for_ranks(need: int, max_total: int, *, poll_s: float = 30.0,
                   timeout_s: float | None = None, counter: Callable[[], int] = running_qe_ranks
                   ) -> float:
    """Block until ``counter() + need <= max_total``. Returns seconds waited."""
    t0 = time.monotonic()
    while counter() + need > max_total:
        if timeout_s is not None and time.monotonic() - t0 > timeout_s:
            raise TimeoutError(f"machine rank budget {max_total} still full after {timeout_s} s")
        time.sleep(poll_s)
    return time.monotonic() - t0


@dataclass
class RunPolicy:
    """Machine-side knobs (machine.toml [disk] and the rank budget)."""
    min_free_gb: float | None = None          # disk guard; None = off
    disk_path: str | None = None              # file system to watch (default: workdir)
    central_manifest_dir: str | None = None   # second copy of prune manifests
    max_machine_ranks: int | None = 12        # None = do not wait
    rank_poll_s: float = 30.0
    rank_wait_timeout_s: float | None = None
    watcher_interval_s: float = 60.0

    @classmethod
    def from_machine(cls, machine: Any | None) -> "RunPolicy":
        if machine is None:
            return cls()
        raw = getattr(machine, "raw", machine)
        d = raw.get("disk", {}) or {}
        return cls(min_free_gb=d.get("min_free_GB"), disk_path=d.get("watch_path"),
                   central_manifest_dir=d.get("manifest_dir"),
                   max_machine_ranks=(raw.get("ranks", {}) or {}).get("machine_max", 12))


# ---------------------------------------------------------------------- records
@dataclass
class StepTiming:
    step: str
    state: str                    # "ran" | "adopted" | "done (state)" | "failed"
    wall_s: float = 0.0
    cpu_s: float = 0.0            # user + sys of child processes (all ranks)
    ranks: int = 1
    waited_for_ranks_s: float = 0.0
    paused_s: float = 0.0
    pruned_bytes: int = 0
    started: float | None = None
    error: str = ""

    @property
    def core_h(self) -> float:
        return self.cpu_s / 3600.0


@dataclass
class DFPTRun:
    """Result of one profile run (also used for the electronic and hull chains)."""
    material_id: str
    profile: dict
    workdir: str
    status: str                                  # "done" | "failed"
    q: list[int] | None = None
    grids: dict | None = None
    settings: dict | None = None                 # core EPhSettings.to_dict()
    cell0: dict | None = None                    # input structure (core CellData dict)
    relaxed_cell: dict | None = None
    pseudo_dir: str | None = None
    pseudo_files: dict | None = None
    pseudo_md5: dict | None = None
    qe_version: str | None = None
    steps: list[dict] = field(default_factory=list)
    error: str = ""
    adopted: bool = False
    pruned_bytes: int = 0
    started: float | None = None
    finished: float | None = None

    @property
    def profile_name(self) -> str:
        return self.profile["name"]

    @property
    def spectrum_tier(self) -> str | None:
        return self.profile.get("spectrum_tier")

    @property
    def ok(self) -> bool:
        return self.status == "done"

    @property
    def wall_s(self) -> float:
        return float(sum(s.get("wall_s", 0.0) for s in self.steps))

    @property
    def core_h(self) -> float:
        return float(sum(s.get("cpu_s", 0.0) for s in self.steps)) / 3600.0

    def to_dict(self) -> dict[str, Any]:
        d = asdict(self)
        d["wall_s"] = self.wall_s
        d["core_h"] = self.core_h
        return d

    @classmethod
    def from_dict(cls, d: dict[str, Any]) -> "DFPTRun":
        keys = set(cls.__dataclass_fields__)
        return cls(**{k: v for k, v in d.items() if k in keys})


def _log_timing(workdir: Path, material_id: str, profile: str, t: StepTiming) -> None:
    row = asdict(t)
    row.update(material_id=material_id, profile=profile, core_h=t.core_h)
    with open(Path(workdir) / TIMING_LOG, "a") as f:
        f.write(json.dumps(row) + "\n")


def _children_cpu() -> float:
    ru = resource.getrusage(resource.RUSAGE_CHILDREN)
    return ru.ru_utime + ru.ru_stime


def _guard(workdir: Path, policy: RunPolicy) -> None:
    if policy.min_free_gb is not None:
        disk_guard(policy.disk_path or workdir, policy.min_free_gb)


def _watcher(workdir: Path, policy: RunPolicy, prune: bool, with_dvscf: bool) -> ScratchWatcher:
    return ScratchWatcher(workdir, interval_s=policy.watcher_interval_s, prune=prune,
                          with_dvscf=with_dvscf,
                          central_dir=Path(policy.central_manifest_dir).expanduser()
                          if policy.central_manifest_dir else None,
                          disk_path=Path(policy.disk_path) if policy.disk_path else workdir,
                          min_free_gb=policy.min_free_gb)


# ---------------------------------------------------------------------- DFPT profiles
def _dfpt_outputs(workdir: Path, nsigma: int) -> list[Path]:
    w = Path(workdir)
    names = [f"a2F.dos{i}" for i in range(1, nsigma + 1)] + ["lambda", "freq.dos", "phonon.dos",
                                                              "ifc.fc", "dyn0"]
    names += sorted(p.name for p in w.glob("dyn[0-9]*") if p.name != "dyn0")
    names += [f"{s}.in" for s in ("relax", "scf_fit", "scf", "ph", "q2r", "matdyn")]
    names += [f"{s}.out" for s in ("relax", "scf_fit", "scf", "ph", "q2r", "matdyn")]
    names += ["relaxed_structure.json", "dfpt_state.json"]
    return [w / n for n in names if (w / n).is_file()]


def output_hashes(workdir: Path, nsigma: int) -> dict[str, str]:
    return {p.name: sha256_file(p) for p in _dfpt_outputs(workdir, nsigma)}


def _marker_valid(workdir: Path, wf: Any, nsigma: int) -> dict | None:
    m = Path(workdir) / DONE_MARKER
    if not m.is_file():
        return None
    try:
        d = json.loads(m.read_text())
        if d.get("matdyn_token") != wf.token("matdyn"):
            return None
        hashes = d.get("output_sha256", {})
        if not hashes:
            return None
        for name, h in hashes.items():
            if name in ("dfpt_state.json", DONE_MARKER):
                continue
            p = Path(workdir) / name
            if not p.is_file() or sha256_file(p) != h:
                return None
        return d
    except Exception:
        return None


def build_eph_workflow(material_id: str, struct: Any, workdir: Path, profile: Profile,
                       launcher: Any, *, q: Any = None, pressure_GPa: float = 0.0,
                       do_relax: bool | None = None) -> Any:
    from rtt.core.dfpt.structure import to_celldata
    from rtt.core.dfpt.workflow import EPhWorkflow
    cell = to_celldata(struct, _masses_for(struct))
    pseudos = sg15_pseudos(cell.elements, pattern=profile.pseudo_pattern)
    settings = profile.eph_settings(q=q, pressure_GPa=pressure_GPa, do_relax=do_relax)
    return EPhWorkflow(material_id, cell, workdir, settings, pseudos, launcher)


def _masses_for(struct: Any) -> dict[str, float] | None:
    """Standard atomic masses (ASE table) for plain-dict structures without masses."""
    if isinstance(struct, dict) and not (struct.get("masses_amu") or struct.get("masses")):
        from ase.data import atomic_masses, atomic_numbers
        return {el: float(atomic_masses[atomic_numbers[el]]) for el in set(struct["species"])}
    return None


def run_profile(struct: Any, profile: Profile, workdir: str | Path, launcher: Any = None, *,
                material_id: str = "material", q: Any = None, pressure_GPa: float = 0.0,
                do_relax: bool | None = None, policy: RunPolicy | None = None) -> DFPTRun:
    """Run one DFPT profile (``coarse``/``full``) for ``struct`` in ``workdir``.

    For electronic and hull profiles use :func:`rtt.dft.electronic.run_electronic` and
    :func:`rtt.dft.refhull.ref_hull` (this function dispatches to the former for ``s4``)."""
    if profile.kind == "electronic":
        from .electronic import run_electronic
        return run_electronic(struct, workdir, profile=profile, launcher=launcher,
                              material_id=material_id, pressure_GPa=pressure_GPa,
                              policy=policy).run
    if profile.kind != "dfpt":
        raise ValueError(f"run_profile: profile {profile.name} ({profile.kind}) is not runnable here")
    policy = policy or RunPolicy()
    launcher = launcher or default_launcher(2)
    w = Path(workdir).expanduser()
    w.mkdir(parents=True, exist_ok=True)
    wf = build_eph_workflow(material_id, struct, w, profile, launcher, q=q,
                            pressure_GPa=pressure_GPa, do_relax=do_relax)
    s = wf.settings
    run = DFPTRun(material_id=material_id, profile=profile.to_dict(), workdir=str(w),
                  status="failed", q=list(s.q) if s.q else None, settings=s.to_dict(),
                  cell0=wf.cell0.to_dict(), pseudo_dir=str(wf.pseudos.directory),
                  pseudo_files=dict(wf.pseudos.files),
                  pseudo_md5=pseudo_md5s(wf.pseudos, wf.cell0.elements), started=time.time())
    try:
        marker = _marker_valid(w, wf, s.el_ph_nsigma)
    except Exception:
        marker = None
    if marker is not None:
        run.status, run.adopted = "done", True
        run.steps = marker.get("steps", [])
        _finish(run, wf)
        return run
    central = Path(policy.central_manifest_dir).expanduser() if policy.central_manifest_dir else None
    for st in wf.steps():
        ok, why = wf.is_complete(st)
        if ok:
            t = StepTiming(st, "adopted" if "adopted" in why else "done (state)")
            run.steps.append(asdict(t))
            _log_timing(w, material_id, profile.name, t)
            continue
        parallel = st in ("relax", "scf_fit", "scf", "ph")
        ranks = int(launcher.nprocs) if parallel else 1
        t = StepTiming(st, "ran", ranks=ranks, started=time.time())
        try:
            _guard(w, policy)
            if policy.max_machine_ranks is not None:
                t.waited_for_ranks_s = wait_for_ranks(ranks, policy.max_machine_ranks,
                                                      poll_s=policy.rank_poll_s,
                                                      timeout_s=policy.rank_wait_timeout_s)
        except (DiskSpaceLow, TimeoutError) as exc:
            t.state, t.error = "failed", str(exc)
            run.steps.append(asdict(t))
            _log_timing(w, material_id, profile.name, t)
            run.error = f"step {st} not started: {exc}"
            _finish(run, wf)
            return run
        if st == "ph":
            try:        # not on a core recover restart (same token interrupted): dyn files are reused
                same = wf.state["steps"].get("ph", {}).get("token") == wf.token("ph")
            except Exception:  # noqa: BLE001
                same = False
            if not same:
                _clear_stale_ph_outputs(w)
        cpu0, t0 = _children_cpu(), time.monotonic()
        use_watch = st == "ph" and (profile.prune_per_q or policy.min_free_gb is not None)
        watch = (_watcher(w, policy, profile.prune_per_q, profile.prune_dvscf_per_q)
                 if use_watch else contextlib.nullcontext())
        err = ""
        with watch as wobj:
            try:
                wf.run(until=st, collect=False)
            except Exception as exc:          # core WorkflowError carries the diagnosis
                err = str(exc)
        t.wall_s = time.monotonic() - t0
        t.cpu_s = _children_cpu() - cpu0
        if isinstance(wobj, ScratchWatcher):
            t.pruned_bytes, t.paused_s = wobj.pruned_bytes, wobj.paused_s
            run.pruned_bytes += wobj.pruned_bytes
        rec = wf.state["steps"].get(st, {})
        if err or rec.get("status") != "done":
            t.state, t.error = "failed", err or f"status {rec.get('status')}"
            run.steps.append(asdict(t))
            _log_timing(w, material_id, profile.name, t)
            run.error = t.error
            _finish(run, wf)
            return run
        run.steps.append(asdict(t))
        _log_timing(w, material_id, profile.name, t)
    run.status = "done"
    hashes = output_hashes(w, s.el_ph_nsigma)
    run.pruned_bytes += prune_final(w, profile.prune_final, central_dir=central)
    (w / DONE_MARKER).write_text(json.dumps({
        "matdyn_token": wf.token("matdyn"), "output_sha256": hashes, "steps": run.steps,
        "profile_sha256": profile.sha256(), "finished": time.time()}, indent=1))
    _finish(run, wf)
    return run


def _clear_stale_ph_outputs(w: Path) -> list[str]:
    """Before a (re)run of ph.x, move dyn files of an earlier run aside (campaign 1). The per-q
    ScratchWatcher treats q as finished when its dyn file is complete and the next one exists;
    stale dyn0..dynN of a killed run with different input made it prune the *current* q's
    wavefunctions ("file pwscf.wfc not found"). Files are moved to ``stale_ph_<time>/``, not
    deleted."""
    stale = sorted(p for p in w.glob("dyn*") if p.is_file() and p.name[3:].isdigit())
    if not stale:
        return []
    dest = w / f"stale_ph_{int(time.time())}"
    dest.mkdir(exist_ok=True)
    for p in stale:
        p.rename(dest / p.name)
    return [p.name for p in stale]


def _finish(run: DFPTRun, wf: Any) -> None:
    run.finished = time.time()
    try:
        run.grids = wf.grids().to_dict()
    except Exception:
        pass
    try:
        run.relaxed_cell = wf.working_cell().to_dict()
    except Exception:
        run.relaxed_cell = None
    from rtt.core.dfpt.parsers import parse_qe_version
    for name in ("matdyn.out", "ph.out", "scf.out"):
        p = Path(run.workdir) / name
        if p.is_file():
            v = parse_qe_version(p.read_text(errors="replace"))
            if v:
                run.qe_version = v
                break
    (Path(run.workdir) / "dfpt_run.json").write_text(json.dumps(run.to_dict(), indent=1,
                                                                default=str))


def run_full(struct: Any, profile: Profile, workroot: str | Path, launcher: Any = None, *,
             material_id: str = "material", pressure_GPa: float = 0.0,
             policy: RunPolicy | None = None) -> list[DFPTRun]:
    """Run ``profile`` on each of ``profile.q_meshes`` (coarser first), each in
    ``workroot/q<n1>x<n2>x<n3>``. The first mesh relaxes (at the profile cutoff and k; R12 #14
    for S6); later meshes reuse that relaxed cell without relaxing again, so the q-spread
    (delta^q, eps_omega) compares the same structure. S6: ``full`` (4^3, 6^3); S5: ``coarse``
    (3^3, 4^3; doc 13 T2-SCREEN, :func:`run_coarse`)."""
    meshes = list(profile.q_meshes) or [profile.q]
    root = Path(workroot).expanduser()
    out: list[DFPTRun] = []
    cell = struct
    for i, qm in enumerate(meshes):
        sub = root / ("q" + "x".join(str(x) for x in qm))
        r = run_profile(cell, profile, sub, launcher, material_id=material_id, q=qm,
                        pressure_GPa=pressure_GPa, do_relax=(i == 0) and profile.do_relax,
                        policy=policy)
        out.append(r)
        if not r.ok:
            break
        if i == 0 and r.relaxed_cell:
            cell = r.relaxed_cell
    return out


def run_coarse(struct: Any, profile: Profile, workroot: str | Path, launcher: Any = None, *,
               material_id: str = "material", pressure_GPa: float = 0.0,
               policy: RunPolicy | None = None) -> list[DFPTRun]:
    """S5: the ``coarse`` profile on q 3^3 and q 4^3 (``profile.q_meshes``), same relaxed cell.
    Feed the last run's metrics to :func:`rtt.dft.coarse.decide_coarse` as ``metrics`` and the
    first's as ``metrics_coarser`` (doc 13 §5 T2-SCREEN)."""
    return run_full(struct, profile, workroot, launcher, material_id=material_id,
                    pressure_GPa=pressure_GPa, policy=policy)

# ---------------------------------------------------------------------- generic step chain
@dataclass
class ChainStep:
    name: str
    code: str                                  # "pw.x", "pp.x", ...
    text: Callable[[], str]                    # input text (called lazily: may need upstream)
    check: Callable[[Path], tuple[bool, str]]  # products check
    upstream: tuple[str, ...] = ()
    parallel: bool = True
    after: Callable[[Path], None] | None = None


#: pw.x failures that a more robust solver usually fixes (campaign 1: La hydrides in nscf,
#: YH3 vc-relax). Matched against the step's error message / output.
PW_FALLBACK_ERRORS = ("too many bands are not converged", "dE0s is positive",
                      "S matrix not positive definite", "problems computing cholesky")


def pw_fallback_text(text: str) -> str:
    """Robust variant of a pw.x input: conjugate-gradient diagonalisation, and damped
    ion/cell dynamics instead of BFGS for relaxations. Physics (cutoff, k, smearing, pressure,
    convergence thresholds) is unchanged."""
    out = []
    for line in text.splitlines():
        low = line.strip().lower()
        if low == "&electrons":
            out += [line, "    diagonalization = 'cg'"]
            continue
        if low.startswith("ion_dynamics"):
            line = "    ion_dynamics = 'damp'"
        elif low.startswith("cell_dynamics"):
            line = "    cell_dynamics = 'damp-w'"
        out.append(line)
    return "\n".join(out) + ("\n" if text.endswith("\n") else "")


#: pw.x relax that stopped at ``nstep`` without converging (campaign 1: BaIrH6 S4 vc-relax). QE
#: then prints no "End final coordinates" block, so the products check fails. Matched in the output.
NSTEP_EXHAUSTED = ("maximum number of steps has been reached",)
NSTEP_RETRY_FACTOR = 2
BOHR_A = 0.529177210903


def _last_block(lines: list[str], head: str, n: int) -> tuple[str, list[str]] | None:
    idx = [i for i, ln in enumerate(lines) if ln.strip().startswith(head)]
    if not idx:
        return None
    i = idx[-1]
    rows = lines[i + 1: i + 1 + n]
    return (lines[i], rows) if len(rows) == n else None


def pw_continue_text(inp: str, out: str, factor: int = NSTEP_RETRY_FACTOR) -> tuple[str, dict]:
    """Input for one continuation of a relax/vc-relax that hit ``nstep``: the same input with
    ``nstep`` x ``factor`` and the LAST cell and positions printed in ``out`` (QE prints them
    after every ionic step). Physics settings are unchanged. Returns (text, info)."""
    import re
    info: dict[str, Any] = {}
    m = re.search(r"(?im)^(\s*nstep\s*=\s*)(\d+)", inp)
    if m:
        old = int(m.group(2))
        new = old * int(factor)
        txt = inp[:m.start(2)] + str(new) + inp[m.end(2):]
    else:                                    # QE default nstep for relaxations is 50
        old, new = 50, 50 * int(factor)
        txt = re.sub(r"(?im)^(\s*&control\s*)$", lambda mm: mm.group(1) + f"\n    nstep = {new}", inp, count=1)
    info.update(nstep_old=old, nstep_new=new)
    nat_m = re.search(r"(?im)^\s*nat\s*=\s*(\d+)", inp)
    lines_out = out.splitlines()
    lines_in = txt.splitlines()
    geo = []
    cell = _last_block(lines_out, "CELL_PARAMETERS", 3)
    ci = [i for i, ln in enumerate(lines_in) if ln.strip().upper().startswith("CELL_PARAMETERS")]
    if cell is not None and ci:
        head, rows = cell
        try:
            vals = [[float(x.replace("d", "e").replace("D", "e")) for x in r.split()[:3]] for r in rows]
            if any(len(v) != 3 for v in vals):
                raise ValueError("short cell row")
            hl = head.lower()
            ma = re.search(r"alat\s*=\s*([-+0-9.eEdD]+)", head)
            if "angstrom" in hl:
                f = 1.0
            elif "bohr" in hl:
                f = BOHR_A
            elif ma:
                f = float(ma.group(1).replace("d", "e").replace("D", "e")) * BOHR_A
            else:
                raise ValueError(f"cell units unknown: {head.strip()}")
            i = ci[-1]
            lines_in[i:i + 4] = ["CELL_PARAMETERS angstrom"] + [
                "  " + " ".join(f"{x * f:20.14f}" for x in v) for v in vals]
            geo.append("cell")
        except ValueError as exc:
            info["cell_note"] = str(exc)
    if nat_m:
        nat = int(nat_m.group(1))
        pos = _last_block(lines_out, "ATOMIC_POSITIONS", nat)
        pi = [i for i, ln in enumerate(lines_in) if ln.strip().upper().startswith("ATOMIC_POSITIONS")]
        if pos is not None and pi and all(len(r.split()) >= 4 for r in pos[1]):
            i = pi[-1]
            lines_in[i:i + 1 + nat] = [pos[0].strip()] + [" " + r.strip() for r in pos[1]]
            geo.append("positions")
    info["geometry"] = ("last printed " + " and ".join(geo)) if geo else "initial (no ionic step printed)"
    return "\n".join(lines_in) + ("\n" if txt.endswith("\n") else ""), info


def _sha(s: str) -> str:
    return hashlib.sha256(s.encode()).hexdigest()


class Chain:
    """Token-chained, resumable runner for small pw.x/pp.x chains (same rules as the core
    workflow: a step is done when its token -- hash of its input and upstream tokens -- matches
    and its products check passes; without a record it is adopted if the existing ``.in`` is
    byte-identical and upstream is done)."""

    def __init__(self, workdir: Path, steps: list[ChainStep], launcher: Any, *,
                 material_id: str, profile: str, policy: RunPolicy | None = None,
                 state_name: str = "chain_state.json"):
        self.w = Path(workdir)
        self.w.mkdir(parents=True, exist_ok=True)
        (self.w / "tmp").mkdir(exist_ok=True)
        self.steps = steps
        self.by = {s.name: s for s in steps}
        self.launcher = launcher
        self.material_id = material_id
        self.profile = profile
        self.policy = policy or RunPolicy()
        self.state_path = self.w / state_name
        try:
            self.state = json.loads(self.state_path.read_text()) if self.state_path.exists() else {}
        except json.JSONDecodeError:
            self.state = {}
        self.state.setdefault("steps", {})
        self.timings: list[dict] = []
        self.retries: list[dict] = []          # nstep continuations made by this run (recorded)

    def _save(self) -> None:
        tmp = self.state_path.with_suffix(".tmp")
        tmp.write_text(json.dumps(self.state, indent=1))
        tmp.replace(self.state_path)

    def _fallback_marker(self, name: str) -> Path:
        return self.w / f"{name}.fallback"

    def _nstep_marker(self, name: str) -> Path:
        return self.w / f"{name}.nstep_retry"

    def text(self, name: str) -> str:
        """Step input; the continuation input once an nstep-retry marker exists (it stores the
        exact text), then the robust pw.x variant once a fallback marker exists (both sticky, so
        a finished retry stays adoptable)."""
        st = self.by[name]
        nm = self._nstep_marker(name)
        if st.code == "pw.x" and nm.exists():
            try:
                t = json.loads(nm.read_text())["text"]
            except (json.JSONDecodeError, KeyError):
                t = st.text()
        else:
            t = st.text()
        return pw_fallback_text(t) if st.code == "pw.x" and self._fallback_marker(name).exists() else t

    def _try_nstep_retry(self, st: ChainStep, msg: str) -> bool:
        """Campaign 1 (BaIrH6): a relax that hit ``nstep`` gets ONE continuation from its last
        geometry with nstep x 2, recorded in ``<step>.nstep_retry`` (JSON: reason, old/new
        nstep, geometry source, the first attempt's output sha256 and a copy of it)."""
        if st.code != "pw.x" or self._nstep_marker(st.name).exists():
            return False
        out = self.w / f"{st.name}.out"
        blob = out.read_text(errors="replace") if out.exists() else ""
        if not any(p in blob for p in NSTEP_EXHAUSTED):
            return False
        inp = st.text()                     # base input (fallback is applied on top by text())
        if "relax" not in inp:
            return False
        new, info = pw_continue_text(inp, blob)
        keep = self.w / f"{st.name}.nstep1.out"
        keep.write_text(blob)
        self._nstep_marker(st.name).write_text(json.dumps({
            "reason": "nstep exhausted (maximum number of steps has been reached)", **info,
            "first_output": keep.name, "first_output_sha256": _sha(blob), "time": time.time(),
            "text": new}, indent=1))
        self.retries.append({"step": st.name, "kind": "nstep_retry", **info})
        return True

    def token(self, name: str) -> str:
        st = self.by[name]
        parts = [self.text(name)] + [self.state["steps"].get(u, {}).get("token", "MISSING")
                               for u in st.upstream]
        return _sha("\n#--#\n".join(parts))

    def is_complete(self, name: str) -> tuple[bool, str]:
        st = self.by[name]
        try:
            tok = self.token(name)
        except Exception as exc:
            return False, f"input not available: {exc}"
        ok, why = st.check(self.w)
        if not ok:
            return False, why
        rec = self.state["steps"].get(name, {})
        if rec.get("token") == tok and rec.get("status") == "done":
            return True, "done (state)"
        inp = self.w / f"{name}.in"
        if not rec and inp.exists() and inp.read_text() == self.text(name) and all(
                self.state["steps"].get(u, {}).get("status") == "done" for u in st.upstream):
            self.state["steps"][name] = {"token": tok, "status": "done", "adopted": True,
                                         "time": time.time()}
            self._save()
            return True, "adopted existing outputs"
        return False, "inputs or upstream changed since last run"

    def run(self) -> tuple[bool, str]:
        for st in self.steps:
            ok, msg = self._run_step(st)
            if not ok and self._try_nstep_retry(st, msg):
                ok, msg = self._run_step(st)
            if not ok and st.code == "pw.x" and not self._fallback_marker(st.name).exists():
                out = self.w / f"{st.name}.out"
                blob = msg + (out.read_text(errors="replace")[-20000:] if out.exists() else "")
                hit = next((p for p in PW_FALLBACK_ERRORS if p in blob), None)
                if hit:          # one retry with the robust solver; recorded by the marker file
                    self._fallback_marker(st.name).write_text(json.dumps({"error": hit, "time": time.time()}))
                    ok, msg = self._run_step(st)
            if not ok:
                return False, msg
        return True, ""

    def _run_step(self, st: ChainStep) -> tuple[bool, str]:
        from rtt.core.dfpt.launcher import run_command
        for st in (st,):
            ok, why = self.is_complete(st.name)
            if ok:
                t = StepTiming(st.name, "adopted" if "adopted" in why else "done (state)")
                self.timings.append(asdict(t))
                _log_timing(self.w, self.material_id, self.profile, t)
                continue
            ranks = int(self.launcher.nprocs) if st.parallel else 1
            t = StepTiming(st.name, "ran", ranks=ranks, started=time.time())
            try:
                _guard(self.w, self.policy)
                if self.policy.max_machine_ranks is not None:
                    t.waited_for_ranks_s = wait_for_ranks(
                        ranks, self.policy.max_machine_ranks, poll_s=self.policy.rank_poll_s,
                        timeout_s=self.policy.rank_wait_timeout_s)
            except (DiskSpaceLow, TimeoutError) as exc:
                t.state, t.error = "failed", str(exc)
                self.timings.append(asdict(t))
                _log_timing(self.w, self.material_id, self.profile, t)
                return False, f"step {st.name} not started: {exc}"
            tok = self.token(st.name)
            (self.w / f"{st.name}.in").write_text(self.text(st.name))
            self.state["steps"][st.name] = {"token": tok, "status": "running", "time": time.time()}
            seen = False
            for other in self.steps:            # invalidate downstream
                if seen:
                    self.state["steps"][other.name] = {"status": "invalidated", "by": st.name}
                seen = seen or other.name == st.name
            self._save()
            if st.parallel:
                cmd = self.launcher.command(st.code, f"{st.name}.in")
            else:
                cmd = ["nice", "-n", str(self.launcher.nice or 19),
                       self.launcher.exe(st.code), "-in", f"{st.name}.in"]
            cpu0, t0 = _children_cpu(), time.monotonic()
            r = run_command(st.name, cmd, self.w, env=self.launcher.environment(),
                            timeout_s=None)
            t.wall_s, t.cpu_s = time.monotonic() - t0, _children_cpu() - cpu0
            status, msg = "failed", r.error_message
            if r.ok:
                try:
                    if st.after:
                        st.after(self.w)
                    ok2, why2 = st.check(self.w)
                except Exception as exc:
                    ok2, why2 = False, f"post-processing failed: {exc}"
                if ok2:
                    status = "done"
                else:
                    msg = (msg + " || " if msg else "") + why2
            self.state["steps"][st.name].update({"status": status, "record": r.to_dict()})
            self._save()
            t.state = "ran" if status == "done" else "failed"
            t.error = "" if status == "done" else msg
            self.timings.append(asdict(t))
            _log_timing(self.w, self.material_id, self.profile, t)
            if status != "done":
                return False, f"step {st.name} failed (rc={r.returncode}): {msg}"
        return True, ""
