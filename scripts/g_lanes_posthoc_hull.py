"""Doc 34 §2: post-hoc MatterSim-only E_hull for H-containing S3 cells whose in-run ML-self hull
failed (CUDA OOM while relaxing reference entries under the 0.4 x 24 GB cap that let S3 share the
GPU with MatterGen). Uses the same cached hull function as the run (entry-level cache,
campaign2_lane_p.make_cached_hull) and the S3-relaxed structure + enthalpy stored in the ledger;
the ledger outcome is NOT changed (a missing hull cannot remove; this is for ranking/reporting).

    RTT_GPU_FRACTION=0.6 ~/.rtt/envs/py/bin/python scripts/g_lanes_posthoc_hull.py \
        --ledger ~/rtsc_targeted_work/G-lanes/ledger.sqlite --work ~/rtsc_targeted_work/G-lanes
Writes <work>/posthoc_hull.json  {cid: {e_hull_ml_self_eV, hull_self}}.
"""
from __future__ import annotations

import argparse
import json
import os
import sqlite3
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
sys.path.insert(0, str(Path(__file__).resolve().parent))


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--ledger", required=True)
    ap.add_argument("--work", required=True)
    a = ap.parse_args()
    work = Path(a.work).expanduser()
    import torch
    if torch.cuda.is_available():
        torch.cuda.set_per_process_memory_fraction(float(os.environ.get("RTT_GPU_FRACTION", "0.6")))
    from campaign2_lane_p import make_cached_hull
    from rtt.mlff.hull import ehull_self
    hull_fn = make_cached_hull(work, lambda m: print(m, flush=True))
    c = sqlite3.connect(f"file:{Path(a.ledger).expanduser()}?mode=ro", uri=True)
    out_f = work / "posthoc_hull.json"
    out = json.loads(out_f.read_text()) if out_f.exists() else {}
    for cid, mj in c.execute("SELECT t.candidate_id, t.metrics_json FROM tier_runs t JOIN outcomes o ON "
                             "o.tier_run_id = t.id WHERE t.stage='S3' AND t.state='done'"):
        m = json.loads(mj or "{}") or {}
        rs = m.get("relaxed_structure")
        if cid in out or not rs or "H" not in rs["species"] or m.get("e_hull_ml_self_eV") is not None:
            continue
        if not (m.get("hull_self") or {}).get("error"):
            continue
        cs = "-".join(sorted(set(rs["species"])))
        try:
            h = hull_fn(cs, pressure_GPa=0.0, device="cuda")
            eh, info = ehull_self(rs, m["relax"]["enthalpy_per_atom_eV"], cs, hull=h)
        except Exception as ex:  # noqa: BLE001
            eh, info = None, {"error": f"{type(ex).__name__}: {ex}"[:300]}
        out[cid] = {"e_hull_ml_self_eV": eh, "hull_self": info, "in_run_error": m["hull_self"]["error"][:120]}
        print(cid, cs, eh, flush=True)
        out_f.write_text(json.dumps(out, indent=1))
    out_f.write_text(json.dumps(out, indent=1))
    return 0


if __name__ == "__main__":
    sys.exit(main())
