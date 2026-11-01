"""Spectra and artifact references from a finished DFPT run (architecture §5.1, R4).

``spectra_from_run(run)`` turns a :class:`rtt.dft.workflow.DFPTRun` into core
:class:`Spectrum` objects -- one per *admitted* smearing -- using only the core parsers and
admission rule (``EPhWorkflow.collect`` -> ``collect_qe_eph`` -> ``admit_smearings`` with
``smearing_admission = "mode_positive"``, ``qe_a2f_dos_to_spectrum``), and returns the sha256
:class:`rtt.records.ArtifactRef` of every file the spectra and the run record depend on.

Each spectrum's provenance gains ``source_sha256`` (its a2F.dos file), ``profile``,
``spectrum_tier``, ``spectrum_digest`` and, when a ledger is given, ``source_artifact_id``.
If any requested smearing is unusable, or admission fails, :func:`spectra_from_run` raises
:class:`CollectError` (core rule: a verdict needs every alternative; the stage is ERROR).
"""
from __future__ import annotations

import hashlib
import json
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

from .prune import sha256_file
from .workflow import DFPTRun


class CollectError(RuntimeError):
    def __init__(self, msg: str, collected: "Collected | None" = None):
        super().__init__(msg)
        self.collected = collected


def spectrum_digest(sp: Any) -> str:
    """sha256 of the exact numerical content (positions, weights, envelope; floats via repr)."""
    payload = json.dumps({"omega_K": [repr(x) for x in sp.omega_K],
                          "weight_K": [repr(x) for x in sp.weight_K],
                          "rel_envelope": repr(sp.rel_envelope)}, separators=(",", ":"))
    return hashlib.sha256(payload.encode()).hexdigest()


@dataclass
class Collected:
    spectra: list
    artifacts: list                 # list[rtt.records.ArtifactRef]
    eph: Any                        # core EPhResult (warnings, flags, checks, provenance)
    admission: dict | None
    errors: list[str] = field(default_factory=list)
    flags: list[str] = field(default_factory=list)
    warnings: list[str] = field(default_factory=list)

    @property
    def ok(self) -> bool:
        return not self.errors and bool(self.spectra)


ARTIFACT_GLOBS = ("a2F.dos*", "lambda", "freq.dos", "phonon.dos", "dyn[0-9]*", "*.in", "*.out",
                  "relaxed_structure.json", "dfpt_state.json", "dfpt_run.json", "timing.jsonl",
                  "profile_done.json", "PRUNED_MANIFEST.tsv")


def run_artifact_paths(workdir: Path) -> list[Path]:
    w = Path(workdir)
    seen: dict[str, Path] = {}
    for g in ARTIFACT_GLOBS:
        for p in sorted(w.glob(g)):
            if p.is_file():
                seen[p.name] = p
    return [seen[k] for k in sorted(seen)]


def _workflow_for(run: DFPTRun) -> Any:
    from rtt.core.dfpt.launcher import Launcher
    from rtt.core.dfpt.pseudo import PseudoTable
    from rtt.core.dfpt.structure import CellData
    from rtt.core.dfpt.workflow import EPhSettings, EPhWorkflow
    settings = EPhSettings.from_config(dict(run.settings or {}))
    pseudos = PseudoTable(run.pseudo_dir or "", dict(run.pseudo_files or {}))
    return EPhWorkflow(run.material_id, CellData.from_dict(run.cell0), run.workdir, settings,
                       pseudos, Launcher())


def collect_run(run: DFPTRun, *, ledger: Any = None, tier_run_id: int | None = None) -> Collected:
    """Parse, admit and hash. Never raises on bad QE output (see :attr:`Collected.errors`)."""
    from rtt.records import ArtifactRef
    if run.profile.get("kind") != "dfpt":
        raise ValueError("collect_run needs a DFPT run")
    if not run.ok:
        return Collected([], [], None, None, errors=[f"run not done: {run.error}"])
    w = Path(run.workdir)
    # RC1 F17: hash the files the parser reads BEFORE parsing and again after; the spectra carry
    # the pre-parse hash, and any change in between is an error (the bytes are not bound).
    parsed_inputs = [p for p in run_artifact_paths(w)
                     if p.name.startswith("a2F.dos") or p.name in ("lambda", "freq.dos")]
    pre = {str(p.resolve()): sha256_file(p) for p in parsed_inputs}
    wf = _workflow_for(run)
    eph = wf.collect()
    errors: list[str] = list(getattr(eph, "errors", []))
    for path, h in pre.items():
        now = sha256_file(path) if Path(path).is_file() else None
        if now != h:
            errors.append(f"{Path(path).name} changed while it was parsed (sha256 {h[:12]} -> "
                          f"{(now or 'missing')[:12]}): spectra not bound to their source")
    arts: list = []
    by_path: dict[str, Any] = {}
    for p in run_artifact_paths(w):
        if ledger is not None and tier_run_id is not None:
            from rtt.ledger.artifacts import register
            ref = register(ledger, tier_run_id, p)
        else:
            ref = ArtifactRef(path=str(p), sha256=sha256_file(p), bytes=p.stat().st_size)
        arts.append(ref)
        by_path[str(p.resolve())] = ref
    for sp in eph.spectra:
        src = Path(sp.provenance.get("source_file", ""))
        key = str(src.resolve()) if str(src) else ""
        sp.provenance["source_sha256"] = pre.get(key) or (
            sha256_file(src) if src.is_file() else None)
        ref = by_path.get(key)
        if ref is not None and ref.id is not None and ref.sha256 == sp.provenance["source_sha256"]:
            sp.provenance["source_artifact_id"] = ref.id
        sp.provenance["profile"] = run.profile_name
        sp.provenance["spectrum_tier"] = run.spectrum_tier
        sp.provenance["q"] = run.q
        sp.provenance["spectrum_digest"] = spectrum_digest(sp)
    return Collected(spectra=list(eph.spectra), artifacts=arts, eph=eph,
                     admission=eph.provenance.get("smearing_admission"),
                     errors=errors, flags=list(eph.flags), warnings=list(eph.warnings))


def spectra_from_run(run: DFPTRun, *, ledger: Any = None, tier_run_id: int | None = None
                     ) -> tuple[list, list]:
    """(spectra, artifact refs) per architecture §5.1. Raises :class:`CollectError` when any
    alternative is unusable or admission fails."""
    c = collect_run(run, ledger=ledger, tier_run_id=tier_run_id)
    if not c.ok:
        raise CollectError("; ".join(c.errors) or "no spectra", c)
    return c.spectra, c.artifacts
