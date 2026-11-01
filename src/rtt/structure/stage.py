"""S2 runner: validity -> duplicate -> novelty, written to the ledger.

Emits only ``PASS | INVALID | DUPLICATE`` (§3.1). Novelty against Alex-MP is a metric (stored
in the PASS outcome's ``value``), never a removal. Candidates are processed only if their current
S1 outcome is PASS and they have no S2 outcome yet.

Controls (lane C) are never DUPLICATE and never representatives: a generated cell that matches a
control is evaluated at the real T* (the control only at its surrogate T*_s). It PASSes with the
link recorded (``value.same_structure_as_control``, tag ``SAME_AS_CONTROL``) -- linked, not
dropped (review RC1 F16).
"""
from __future__ import annotations

from typing import Any, Iterable

from ..records import DUPLICATE, PASS, Outcome
from ..source.ingest import write_outcome
from .dedupe import DedupeIndex, is_novel
from .validate import validate

ORIGIN_S2 = "rtt.structure.stage"


def _novelty(rec: Any, refs_cache: dict, machine: Any = None) -> dict:
    from ..source.external import refquery
    cs = "-".join(sorted(set(rec.structure.species)))
    try:
        if cs not in refs_cache:
            refs_cache[cs] = [(r["entry_id"], r["structure"]) for r in refquery(cs, machine=machine)["rows"]]
        novel, rid = is_novel(rec.structure.to_json(), refs_cache[cs])
        return {"novel": novel, "alexmp_match": rid}
    except Exception as ex:  # noqa: BLE001 - novelty is informational
        return {"novel": None, "novelty_error": f"{type(ex).__name__}: {ex}"[:300]}


def _spacegroup(sd: dict) -> int:
    try:
        from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
        from .canonical import to_pymatgen
        return int(SpacegroupAnalyzer(to_pymatgen(sd), symprec=0.01).get_space_group_number())
    except Exception:  # noqa: BLE001 - ordering only
        return 0


def representative_key(ledger: Any, cid: str, lane: str) -> tuple:
    """Processing order = choice of the duplicate representative (R30 F3, doc 31 §7): controls,
    then lane D (known database entries stay attributed to the database; a generated copy becomes
    its DUPLICATE tagged GENERATED_REDISCOVERY), then generated lanes; within a lane group the
    smallest cell, then the highest space-group number (a symmetry proxy), then the id."""
    group = 0 if lane == "C" else 1 if lane == "D" else 2
    if ledger.current_outcome(cid, "S2") is not None:      # already done: order irrelevant
        return (9, 0, 0, cid)
    rec = ledger.candidate_record(cid)
    sd = rec.structure.to_json()
    return (group, len(sd.get("species", [])), -_spacegroup(sd), cid)


def run_s2(ledger: Any, cids: Iterable[str] | None = None, *, targets: Any, tiers: Any,
           novelty: bool = True, machine: Any = None) -> dict[str, Outcome]:
    """Run S2 for ``cids`` (default: every candidate without S2). Returns cid -> Outcome."""
    # index of earlier S2-PASS structures (the representatives)
    idx = DedupeIndex()
    ctrl_idx = DedupeIndex()
    for cid in ledger.candidate_ids():
        o = ledger.current_outcome(cid, "S2")
        if o and o["outcome"] == PASS:
            rec = ledger.candidate_record(cid)
            if rec.origin.lane == "C":
                # N6/F16: a generated copy of a control is still evaluated at the real T* (linked)
                ctrl_idx.add(cid, rec.structure.to_json(), rec.pressure_GPa, sha=rec.structure_sha256)
                continue
            idx.add(cid, rec.structure.to_json(), rec.pressure_GPa, sha=rec.structure_sha256)
    todo = list(cids) if cids is not None else ledger.candidate_ids()
    lanes = {c: ledger.q1("SELECT lane FROM candidates WHERE id=?", (c,))[0] for c in todo}
    todo.sort(key=lambda c: representative_key(ledger, c, lanes[c]))   # R30 F3: C, then D, small, symmetric
    refs_cache: dict = {}
    out: dict[str, Outcome] = {}
    for cid in todo:
        s1 = ledger.current_outcome(cid, "S1")
        if not s1 or s1["outcome"] != PASS or ledger.current_outcome(cid, "S2"):
            continue
        rec = ledger.candidate_record(cid)
        sd = rec.structure.to_json()
        o = validate(sd, tiers)
        if o.outcome == PASS:
            # controls are injected on purpose and must traverse every stage (§3.3): never a
            # duplicate themselves, and never a representative (N6: other lanes run at the real T*)
            dup = None if rec.origin.lane == "C" else idx.find(sd, rec.pressure_GPa, sha=rec.structure_sha256)
            if dup is not None:
                ledger.insert_duplicate(cid, dup, "structure_matcher")
                rep_lane = ledger.q1("SELECT lane FROM candidates WHERE id=?", (dup,))[0]
                dtags = []
                if rec.origin.lane.startswith("G") and rep_lane == "D":
                    dtags = ["GENERATED_REDISCOVERY", f"REDISCOVERS:{dup}"]      # R30 F3: novel = false
                elif rec.origin.lane == "D" and rep_lane.startswith("G"):
                    dtags = ["KNOWN_ENTRY_OF_GENERATED", f"GENERATED_COPY:{dup}"]
                o = Outcome(stage="S2", outcome=DUPLICATE, criterion="duplicate_of", value=dup,
                            threshold="StructureMatcher(ltol=0.2,stol=0.1,angle_tol=5,scale=False)+bonding_guard",
                            source_id="doc12#3.1:S2", reason=f"duplicate of {dup} (lane {rep_lane})", tags=dtags)
            else:
                meta = _novelty(rec, refs_cache, machine) if novelty else {}
                tags = []
                if rec.origin.lane == "C":
                    ctrl_idx.add(cid, sd, rec.pressure_GPa, sha=rec.structure_sha256)
                else:
                    same = ctrl_idx.find(sd, rec.pressure_GPa, sha=rec.structure_sha256)
                    if same is not None:
                        meta = {**meta, "same_structure_as_control": same}
                        tags.append("SAME_AS_CONTROL")
                    idx.add(cid, sd, rec.pressure_GPa, sha=rec.structure_sha256)
                o = Outcome(stage="S2", outcome=PASS, value=meta or None, tags=tags)
        write_outcome(ledger, cid, o, origin_obj=rec.origin, targets=targets, tiers=tiers, origin=ORIGIN_S2)
        out[cid] = o
    return out


__all__ = ["run_s2", "ORIGIN_S2"]
