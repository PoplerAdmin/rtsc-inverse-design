"""Artifact hashing and prune events (architecture §6.2, R4).

``register(ledger, tier_run_id, path)`` hashes a file and inserts an ``artifacts`` row.
``verify(ledger, artifact_id)`` re-hashes a retained file and records a ``verified``,
``missing`` or ``mismatch`` event. Pruning never deletes rows: ``mark_pruned`` records an event.
"""
from __future__ import annotations

import hashlib
import os
from pathlib import Path
from typing import TYPE_CHECKING

from ..records import ArtifactRef

if TYPE_CHECKING:
    from .store import Ledger


def file_sha256(path: str | os.PathLike) -> tuple[str, int]:
    h = hashlib.sha256()
    n = 0
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
            n += len(chunk)
    return h.hexdigest(), n


def register(ledger: "Ledger", tier_run_id: int, path: str | os.PathLike, retained: bool = True) -> ArtifactRef:
    p = Path(path)
    sha, n = file_sha256(p)
    aid = ledger.insert_artifact(tier_run_id, str(p), sha, n, retained)
    return ArtifactRef(path=str(p), sha256=sha, bytes=n, id=aid, retained=retained)


def verify(ledger: "Ledger", artifact_id: int) -> str:
    row = ledger.q1("SELECT path, sha256 FROM artifacts WHERE id = ?", (artifact_id,))
    if row is None:
        raise KeyError(artifact_id)
    if ledger.artifact_pruned(artifact_id):
        return "pruned"
    path, sha = row
    if not Path(path).exists():
        ledger.insert_artifact_event(artifact_id, "missing", None)
        return "missing"
    got, _ = file_sha256(path)
    ev = "verified" if got == sha else "mismatch"
    ledger.insert_artifact_event(artifact_id, ev, None if ev == "verified" else got)
    return ev


def mark_pruned(ledger: "Ledger", artifact_id: int, note: str = "") -> None:
    ledger.insert_artifact_event(artifact_id, "pruned", note)
