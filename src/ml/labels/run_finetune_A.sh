#!/usr/bin/env bash
# D1-labels: fine-tune mattergen_base on rtt_log10A = log10(A/K), A = int alpha^2F d omega.
# Input: CSV splits written by build_dataset.py into $W/ds/csv.  Output run dir: $W/ft_log10A.
# Micro-batch 16, accumulate 8 => effective batch 128 (micro-batch 32 = 20 GB does not fit beside other GPU users) (more optimizer steps than the default 512 for a
# ~8k-row set); lr 5e-6 (MatterGen default); 200 epochs; val every 5 epochs; best+last checkpoints.
set -eux
W=${W:-$HOME/rtsc_targeted_work/D1-labels}
S=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)   # src/ml
MGPY=~/mattergen/.venv/bin/python
EP=${EP:-200}
export OMP_NUM_THREADS=2 OPENBLAS_NUM_THREADS=2 MKL_NUM_THREADS=2
export MG_EXTRA_PROPS=rtt_log10A
cd "$W"
rm -rf "$W/cache/log10A"
nice -n 19 $MGPY "$S/mg_newprop.py" csv-to-dataset --csv-folder "$W/ds/csv" --dataset-name log10A --cache-folder "$W/cache"
DENS=$(python3 -c "import json;print(json.load(open('$W/ds/summary.json'))['average_density_atoms_per_A3'])")
export OUTPUT_DIR=$W/ft_log10A${TAG:-}
mkdir -p "$OUTPUT_DIR"
( while true; do nvidia-smi --query-gpu=memory.used --format=csv,noheader,nounits; sleep 10; done ) > "$OUTPUT_DIR/gpumem.log" 2>/dev/null &
MON=$!
t0=$(date +%s)
cd ~/mattergen
nice -n 19 $MGPY "$S/mg_newprop.py" finetune \
  adapter.pretrained_name=mattergen_base data_module=mp_20 \
  data_module.root_dir="$W/cache/log10A" \
  data_module.average_density=$DENS \
  "hydra.searchpath=[file://$S/mg_conf]" \
  "+lightning_module/diffusion_module/model/property_embeddings@adapter.adapter.property_embeddings_adapt.rtt_log10A=rtt_log10A" \
  '~trainer.logger' 'data_module.properties=["rtt_log10A"]' \
  data_module.batch_size.train=${MB:-16} trainer.accumulate_grad_batches=${ACC:-8} \
  trainer.max_epochs=$EP trainer.check_val_every_n_epoch=5
rc=$?
kill $MON || true
echo "FT_DONE rc=$rc wall_s=$(( $(date +%s) - t0 )) peak_MiB=$(sort -n "$OUTPUT_DIR/gpumem.log" | tail -1)"
