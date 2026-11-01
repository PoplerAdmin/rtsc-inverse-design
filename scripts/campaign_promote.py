"""Campaign 1 generated-lane quota: promote named DEFERRED(rank_cut) S3 candidates to PASS with tag
CAMPAIGN_GENERATED_QUOTA (a recorded deviation from the global rank cut; see doc 30)."""
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from rtt.config import load_targets, load_tiers
from rtt.ledger.store import Ledger
from rtt.mlff.screen import ORIGIN_S3
from rtt.records import Outcome
from rtt.source.ingest import write_outcome
t, ti = load_targets(), load_tiers()
L = Ledger(Path(sys.argv[1]).expanduser())
for cid in sys.argv[2:]:
    cur = L.current_outcome(cid, "S3")
    assert cur and cur["outcome"] == "DEFERRED" and cur.get("criterion") == "rank_cut", (cid, cur)
    rec = L.candidate_record(cid)
    o = Outcome(stage="S3", outcome="PASS", tags=["PROMOTED_FROM_RANK_CUT", "CAMPAIGN_GENERATED_QUOTA"],
                tier_run_id=cur.get("tier_run_id"))
    write_outcome(L, cid, o, origin_obj=rec.origin, targets=t, tiers=ti, origin=ORIGIN_S3)
    print(cid, rec.formula, "promoted")
