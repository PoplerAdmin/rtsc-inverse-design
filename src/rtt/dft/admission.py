"""S5 per-mesh admission failure -> DEFERRED, not ERROR (doc 31 §4; campaign 1 doc 30 §5.2, §7).

The core mode-positive smearing admission (``rtt.core.dfpt.workflow.admit_smearings``,
``MODE_NEGATIVE_TOL_REL = 1e-6``) is a soundness rule and is sha-pinned; it is NOT changed here.
When it rejects every smearing of one coarse mesh (NaIrH3 q 3^3: max negative per-mode a2F
fraction 3e-6..3e-3 > 1e-6), ``S5Coarse.run`` reports ``ok = False`` and the stage became ERROR,
discarding a finished, informative q 4^3 result and inviting a pointless identical re-run.

Doc 13 §5 (T2-SCREEN) needs admission on BOTH meshes for a removal. A mesh that fails admission
therefore only blocks the removal: the right outcome is ``DEFERRED(admission_failed)`` with the
failing mesh, the core admission report and the admitted mesh's numbers in the record.
Real failures (a pw.x/ph.x crash, a missing or changed file, any collect error that is not the
admission rule) stay ERROR: :func:`admission_deferral` returns None for them.

Call site (S5Coarse.decide, first line; the payload keys are the ones ``S5Coarse.run`` sets)::

    if (o := admission_deferral(res)) is not None: return o
"""
from __future__ import annotations

from typing import Any, Sequence

from ..records import DEFERRED, Outcome

ADMISSION_PREFIX = "smearing admission:"          # core admit_smearings error text
ADMISSION_FLAG = "SMEARING_ADMISSION_FAIL"         # core collect_qe_eph flag
CRITERION = "admission_failed"
SOURCE_ID = "doc13#5 T2-SCREEN (removal needs admission on both meshes); doc31#4"


def _qtag(q: Any) -> str:
    try:
        return "q" + "x".join(str(int(x)) for x in q)
    except Exception:  # noqa: BLE001
        return f"q{q}"


def _strip_mesh_prefix(piece: str) -> str:
    """'q[3, 3, 3]: smearing admission: ...' -> 'smearing admission: ...'."""
    s = piece.strip()
    if s.startswith("q") and "]: " in s[:40]:
        s = s.split("]: ", 1)[1]
    return s.strip()


def only_admission_errors(errors: Sequence[str]) -> bool:
    """True if every error text is the core admission message (and there is at least one)."""
    errs = [e for e in (errors or []) if str(e).strip()]
    return bool(errs) and all(_strip_mesh_prefix(str(e)).startswith(ADMISSION_PREFIX) for e in errs)


def mesh_status(run: Any, col_or_metrics: Any, error_text: str | None = None) -> str:
    """'ok' | 'admission_failed' | 'error' for one coarse mesh.

    ``col_or_metrics`` is a :class:`rtt.dft.collect.Collected` (exact: its ``errors`` list) or a
    :class:`rtt.dft.coarse.CoarseMetrics` (its ``admission_ok`` / ``flags``; then the stage's
    error text, ``error_text``, must consist of admission messages only)."""
    if run is None or not getattr(run, "ok", False) or col_or_metrics is None:
        return "error"
    c = col_or_metrics
    if hasattr(c, "errors"):                                   # Collected
        if c.ok:
            return "ok"
        return "admission_failed" if only_admission_errors(c.errors) else "error"
    if getattr(c, "admission_ok", True):                       # CoarseMetrics
        return "ok"
    if ADMISSION_FLAG not in (getattr(c, "flags", None) or []):
        return "error"
    if error_text is not None:
        pieces = [p for p in str(error_text).split("; ") if p.strip()]
        if not only_admission_errors(pieces):
            return "error"
    return "admission_failed"


def _mesh_summary(q: Any, status: str, c: Any) -> dict[str, Any]:
    d: dict[str, Any] = {"q": list(q) if q is not None else None, "status": status}
    rep = getattr(c, "admission", None)
    if isinstance(rep, dict):
        d["max_mode_negative_fraction"] = rep.get("max_mode_negative_fraction")
        d["tol_rel"] = rep.get("tol_rel")
    if hasattr(c, "to_metrics") and status == "ok":
        m = c.to_metrics()
        d.update(A_mid_K=m.get("A_mid_K"), A_lo_K=m.get("A_lo_K"), admitted_sigmas=m.get("admitted_sigmas"))
    return d


def s5_admission_outcome(meshes: Sequence[Any], runs: Sequence[Any], cols: Sequence[Any], *,
                         tier_run_id: int | None = None, error_text: str | None = None) -> Outcome | None:
    """DEFERRED(admission_failed) when every mesh either succeeded or failed ONLY the core
    admission rule, and at least one failed it; None otherwise (caller keeps its ERROR/decision).

    The outcome records, per mesh, the status, the core max negative per-mode fraction and
    tolerance (when the Collected report is available) and A_mid/A_lo of the admitted mesh(es).
    It never removes and never passes: a mesh without admission cannot support either."""
    if not runs or len(runs) != len(cols):
        return None
    meshes = list(meshes) if meshes else [getattr(c, "q", None) for c in cols]
    st = [mesh_status(r, c, error_text) for r, c in zip(runs, cols)]
    if "error" in st or "admission_failed" not in st:
        return None
    detail = [_mesh_summary(q, s, c) for q, s, c in zip(meshes, st, cols)]
    failed = [_qtag(q) for q, s in zip(meshes, st) if s == "admission_failed"]
    admitted = [d for d in detail if d["status"] == "ok"]
    info = "; ".join(f"{_qtag(d['q'])} admitted: A_mid {d.get('A_mid_K')} K, A_lo {d.get('A_lo_K')} K"
                     for d in admitted) or "no mesh admitted"
    tags = [f"ADMISSION_FAILED:{t}" for t in failed] + ["NO_REMOVAL_WITHOUT_ADMISSION"]
    if admitted:
        tags.append("SINGLE_MESH_ADMITTED" if len(admitted) == 1 else "PARTIAL_ADMISSION")
    return Outcome(stage="S5", outcome=DEFERRED, criterion=CRITERION, value={"meshes": detail},
                   source_id=SOURCE_ID, tier_run_id=tier_run_id, tags=tags,
                   reason=(f"mode-positive admission failed on {', '.join(failed)} (core rule, tol unchanged); "
                           f"removal needs admission on both meshes (doc 13 §5), so no screen: {info}"))


def admission_deferral(res: Any) -> Outcome | None:
    """:func:`s5_admission_outcome` from an S5 ``StageResult`` (payload keys of S5Coarse.run:
    ``runs``, ``metrics_by_mesh`` and, if present, ``collected_by_mesh``). Returns None when the
    stage succeeded, raised before collecting, or failed for any reason other than admission."""
    if res is None or getattr(res, "ok", False):
        return None
    pay = getattr(res, "payload", None) or {}
    runs = pay.get("runs")
    cols = pay.get("collected_by_mesh") or pay.get("metrics_by_mesh")
    if not runs or not cols:
        return None
    meshes = [getattr(r, "q", None) for r in runs]
    return s5_admission_outcome(meshes, runs, cols, tier_run_id=getattr(res, "tier_run_id", None),
                                error_text=getattr(res, "error", None))


__all__ = ["admission_deferral", "s5_admission_outcome", "mesh_status", "only_admission_errors"]
