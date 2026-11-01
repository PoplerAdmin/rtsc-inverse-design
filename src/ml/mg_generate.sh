#!/usr/bin/env bash
# Run mattergen-generate with timing + peak-GPU-memory logging.
# Usage: mg_generate.sh OUTDIR PRETRAINED_NAME BATCH NBATCH SEED [PROPS_DICT] [GUIDANCE]
#   e.g. mg_generate.sh ~/rtsc_targeted_work/T3-mlaudit/gen_uncond mattergen_base 64 1 1
#        mg_generate.sh ~/rtsc_targeted_work/T3-mlaudit/gen_LiMgH chemical_system 64 1 1 "{'chemical_system':'Li-Mg-H'}" 2.0
# Runs in the MatterGen venv (~/mattergen/.venv, python 3.10, torch 2.2.1+cu121).
# Checkpoints are fetched from huggingface.co/microsoft/mattergen on first use (~0.5 GB each).
set -euo pipefail
OUT=$1; NAME=$2; BS=$3; NB=$4; SEED=$5; PROPS=${6:-}; GUID=${7:-0.0}
mkdir -p "$OUT"
export OMP_NUM_THREADS=2 OPENBLAS_NUM_THREADS=2 MKL_NUM_THREADS=2
# sample GPU memory every 1 s in the background
( while true; do nvidia-smi --query-gpu=memory.used --format=csv,noheader,nounits; sleep 1; done ) > "$OUT/gpumem.log" 2>/dev/null &
MON=$!
ARGS=("$OUT" "--pretrained_name=$NAME" "--batch_size=$BS" "--num_batches=$NB" "--seed=$SEED" "--record_trajectories=False")
if [ -n "$PROPS" ]; then ARGS+=("--properties_to_condition_on=$PROPS" "--diffusion_guidance_factor=$GUID"); fi
t0=$(date +%s.%N)
nice -n 19 ~/mattergen/.venv/bin/mattergen-generate "${ARGS[@]}" > "$OUT/mattergen.log" 2>&1 || rc=$?
t1=$(date +%s.%N)
kill $MON 2>/dev/null || true
peak=$(sort -n "$OUT/gpumem.log" | tail -1)
base=$(head -1 "$OUT/gpumem.log")
python3 - "$OUT" "$NAME" "$BS" "$NB" "$PROPS" "$GUID" "$t0" "$t1" "${rc:-0}" "$peak" "$base" <<'EOF'
import json, sys
out, name, bs, nb, props, guid, t0, t1, rc, peak, base = sys.argv[1:]
n = int(bs) * int(nb)
wall = float(t1) - float(t0)
json.dump({"pretrained_name": name, "batch_size": int(bs), "num_batches": int(nb), "n_structures": n,
           "properties": props, "guidance": float(guid), "wall_s": round(wall, 1),
           "s_per_structure": round(wall / n, 2), "returncode": int(rc),
           "gpu_mem_peak_MiB": int(peak), "gpu_mem_before_MiB": int(base)},
          open(f"{out}/timing.json", "w"), indent=1)
print(open(f"{out}/timing.json").read())
EOF
