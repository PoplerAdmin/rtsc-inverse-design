"""S1 ingestion: RawStructure -> canonical structure -> candidate row -> S1 outcome.

Shared helper :func:`write_outcome` is used by S1-S3 to append outcome rows with the right
targets/tiers hashes (a control's outcome is recorded against its surrogate target, §3.3; the
surrogate targets text is registered in ``inputs`` first so the R7 trigger accepts it).
"""
from __future__ import annotations

from typing import Any

from ..records import Origin, Outcome, RawStructure, StructureData, new_candidate_id
from ..structure.canonical import canonicalize, spacegroup_number, structure_sha256
from .prefilter import prefilter, primitive_info

ORIGIN_S1 = "rtt.source.prefilter"


def targets_for(origin: Origin, targets: Any) -> Any:
    if origin.lane == "C" and origin.T_star_surrogate_K and hasattr(targets, "with_T_star"):
        return targets.with_T_star(origin.T_star_surrogate_K)
    return targets


def write_outcome(ledger: Any, cid: str, o: Outcome, *, origin_obj: Origin, targets: Any, tiers: Any,
                  origin: str) -> int:
    t = targets_for(origin_obj, targets)
    ledger.register_input("targets", t.sha256, t.text)
    ledger.register_input("tiers", tiers.sha256, tiers.text)
    return ledger.insert_outcome(cid, o, targets_sha256=t.sha256, tiers_sha256=tiers.sha256, origin=origin)


def candidate_pressure(rs_or_origin: Any, default: float = 0.0) -> float:
    o = rs_or_origin.origin if hasattr(rs_or_origin, "origin") else rs_or_origin
    if o.pressure_override_GPa is not None:
        return float(o.pressure_override_GPa)
    return float(getattr(rs_or_origin, "pressure_GPa", default) or default)


def ingest_raw(ledger: Any, rs: RawStructure, *, targets: Any, tiers: Any, context: Any = None) -> str:
    """Insert one RawStructure as a candidate (canonical primitive cell) and record S1.

    Idempotent: a candidate id (structure hash + lane + pressure) already in the ledger is
    returned without new rows."""
    if context is not None and targets is None:
        targets, tiers = context.targets, context.tiers
    info = primitive_info(rs)
    st = canonicalize(rs.structure.to_json())
    sha = structure_sha256(st)
    p = candidate_pressure(rs)
    cid = new_candidate_id(sha, rs.origin.lane, p)
    if ledger.has_candidate(cid):
        return cid
    sd = StructureData(lattice=[[float(x) for x in r] for r in st.lattice.matrix],
                       species=[s.specie.symbol for s in st],
                       frac_coords=[[float(x) for x in s.frac_coords] for s in st],
                       spacegroup=spacegroup_number(st))
    origin = Origin.from_json(rs.origin.to_json())
    if rs.source_ref:
        # Origin has no source field; keep the provenance (file#index / LMDB id) as a flag
        origin.flags = list(dict.fromkeys([*origin.flags, f"source_ref={rs.source_ref}"]))
    ledger.insert_candidate(cid, sha, sd, origin, pressure_GPa=p, cost_class=info["cost_class"])
    o = prefilter(rs, tiers, info=info)
    write_outcome(ledger, cid, o, origin_obj=origin, targets=targets, tiers=tiers, origin=ORIGIN_S1)
    return cid


def ingest_many(ledger: Any, raws: list[RawStructure], *, targets: Any, tiers: Any) -> list[str]:
    return [ingest_raw(ledger, rs, targets=targets, tiers=tiers) for rs in raws]


__all__ = ["ingest_raw", "ingest_many", "write_outcome", "targets_for", "candidate_pressure", "ORIGIN_S1"]
