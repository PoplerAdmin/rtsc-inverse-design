"""Doc 34: verify the four campaign-2 MatterGen checkpoints (space_group, dft_band_gap,
dft_mag_density, ml_bulk_modulus): sha256 via rtt.source.mattergen.checkpoint_sha256, size,
and that the Lightning checkpoint deserialises in the MatterGen venv with the expected
property_embeddings_adapt (= rtt.source.lanes.CHECKPOINT_PROPERTIES). Prints JSON.

    ~/.rtt/envs/py/bin/python scripts/g_lanes_ckpt_check.py > OUT.json
"""
import json, subprocess, sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from rtt.source.external import path_of
from rtt.source.lanes import CHECKPOINT_PROPERTIES
from rtt.source.mattergen import checkpoint_sha256

NAMES = ("space_group", "dft_band_gap", "dft_mag_density", "ml_bulk_modulus")
PROBE = r'''
import json, sys, torch, yaml
ck, cfg = sys.argv[1], sys.argv[2]
d = torch.load(ck, map_location="cpu", weights_only=False)
sd = d.get("state_dict", {})
c = yaml.safe_load(open(cfg))
def find(o, key):
    if isinstance(o, dict):
        for k, v in o.items():
            if k == key: return v
            r = find(v, key)
            if r is not None: return r
    if isinstance(o, list):
        for v in o:
            r = find(v, key)
            if r is not None: return r
    return None
pea = find(c, "property_embeddings_adapt")
print(json.dumps({"n_tensors": len(sd), "n_params": int(sum(v.numel() for v in sd.values())),
                  "epoch": d.get("epoch"), "global_step": d.get("global_step"),
                  "property_embeddings_adapt": sorted(pea.keys()) if isinstance(pea, dict) else pea,
                  "adapter_keys": sorted({k.split(".")[3] for k in sd if "property_embeddings_adapt" in k and len(k.split("."))>3})}))
'''
out = {}
hub = path_of("hf_hub") / "models--microsoft--mattergen" / "snapshots"
py = str(path_of("mattergen_venv") / "bin" / "python")
for n in NAMES:
    cks = sorted(hub.glob(f"*/checkpoints/{n}/checkpoints/last.ckpt"))
    r = {"sha256": checkpoint_sha256(n), "expected_properties": sorted(CHECKPOINT_PROPERTIES[n])}
    if cks:
        ck = cks[0]
        r["snapshot"] = ck.parts[ck.parts.index("snapshots") + 1]
        r["size_bytes"] = ck.resolve().stat().st_size
        cfg = ck.parent.parent / "config.yaml"
        cp = subprocess.run([py, "-I", "-c", PROBE, str(ck), str(cfg)], capture_output=True, text=True, timeout=600)
        try:
            r["load"] = json.loads(cp.stdout.strip().splitlines()[-1])
            r["loads"] = True
        except Exception:
            r["loads"] = False
            r["stderr_tail"] = cp.stderr[-800:]
        pea = (r.get("load") or {}).get("property_embeddings_adapt")
        r["properties_match"] = pea == r["expected_properties"]
    out[n] = r
print(json.dumps(out, indent=1))
