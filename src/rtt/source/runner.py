"""S1 runners: one MatterGen lane job, or one database query, into the ledger."""
from __future__ import annotations

from typing import Any

from .database import query_many
from .ingest import ingest_many
from .lanes import LaneJob
from .mattergen import generate


def run_lane_job(ledger: Any, job: LaneJob, *, targets: Any, tiers: Any, machine: Any = None,
                 workdir: Any = None) -> list[str]:
    """Generate one lane job with MatterGen (GPU) and ingest the structures (S1 outcomes)."""
    raws = generate(job, workdir=workdir, machine=machine)
    return ingest_many(ledger, raws, targets=targets, tiers=tiers)


def run_database_lane(ledger: Any, systems: list[str], *, targets: Any, tiers: Any, machine: Any = None,
                      ehull_max: float = 0.1, limit: int | None = None) -> list[str]:
    raws = query_many(systems, ehull_max=ehull_max, machine=machine)
    if limit is not None:
        raws = raws[:limit]
    return ingest_many(ledger, raws, targets=targets, tiers=tiers)
