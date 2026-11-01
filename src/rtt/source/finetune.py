"""A-adapter support (architecture §3.3 lane GA; doc 10 §4.2). **Stub.**

* :func:`export_labels` writes the ledger's ``labels`` rows for one property as a MatterGen
  fine-tune CSV (columns ``material_id, cif, chemical_system, <prop>``; the format of
  doc 10 §4.2). Per §3.3, our own S4/S5 labels
  are a *validation* set for the adapter, not its training set.
* :func:`train_adapter` is not implemented here: the adapter is being trained by the label agent
  (``src/ml/labels/run_finetune_A.sh``, property ``rtt_log10A``). Lane GA stays refused by
  :class:`~rtt.source.lanes.LaneSpec` until the acceptance test of §3.3 has been run.
"""
from __future__ import annotations

import csv
import json
from pathlib import Path
from typing import Any


def _cif(sd: dict, name: str) -> str:
    from pymatgen.io.cif import CifWriter
    from ..structure.canonical import to_pymatgen
    return str(CifWriter(to_pymatgen(sd)))


def export_labels(ledger: Any, prop: str = "rtt_A_K", out: str | Path | None = None, *,
                  stage: str | None = None) -> Path:
    conn = ledger.conn if hasattr(ledger, "conn") else ledger
    q = ("SELECT l.candidate_id, l.value, s.json FROM labels l JOIN candidates c ON c.id = l.candidate_id "
         "JOIN structures s ON s.sha256 = c.structure_sha256 WHERE l.name = ?")
    args: list[Any] = [prop]
    if stage:
        q += " AND l.stage = ?"
        args.append(stage)
    rows = conn.execute(q, args).fetchall()
    out = Path(out or f"labels_{prop}.csv")
    out.parent.mkdir(parents=True, exist_ok=True)
    with open(out, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["material_id", "cif", "chemical_system", prop])
        for cid, val, js in rows:
            d = json.loads(js)
            sd = d.get("structure", d)
            w.writerow([cid, _cif(sd, cid), "-".join(sorted(set(sd["species"]))), float(val)])
    return out


def train_adapter(*_: Any, **__: Any) -> None:
    raise NotImplementedError(
        "adapter training is owned by the label agent (src/ml/labels/run_finetune_A.sh; doc 10 §4.2 "
        "template). Lane GA is refused until the §3.3 acceptance test passes.")


__all__ = ["export_labels", "train_adapter"]
