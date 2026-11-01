"""Controls through S1-S3 ahead of the generated lanes (campaign 1), same ledger."""
import os, sys, json
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from rtt.config import load_machine, load_targets, load_tiers
from rtt.ledger.store import Ledger
from rtt.source.controls import inject_controls
from rtt.structure.stage import run_s2
from rtt.mlff.screen import run_s3
t, ti, m = load_targets(), load_tiers(), load_machine()
L = Ledger(Path(sys.argv[1]).expanduser(), run_id="campaign1")
L.register_configs(t, ti, m)
cc = inject_controls(L, targets=t, tiers=ti)
print(run_s2(L, cc, targets=t, tiers=ti, novelty=True, machine=m))
res = run_s3(L, cc, targets=t, tiers=ti, machine=m, device="cuda", workdir=Path(sys.argv[2]).expanduser() / "S3", alignn=False, rank_cut=0)
for c in cc:
    r = L.candidate_record(c)
    print(c, r.formula, r.origin.lane, getattr(r.origin, "name", None), {s: L.current_outcome(c, s)["outcome"] for s in ("S1","S2","S3") if L.current_outcome(c, s)})
