#!/usr/bin/env bash
# T3-mlaudit: measure MatterGen adapter fine-tuning cost on the RTX 4090.
# Uses an EXISTING property id (dft_band_gap) on a 2000/200-structure subset of MP-20 so that no
# file inside ~/mattergen has to be modified; the step time is independent of which scalar is used.
# Dataset prep (done once, CPU, ~9 s):
#   python3 -c "import zipfile; zipfile.ZipFile('dl/mp_20.zip').extractall('.')"   # mp_20.zip from GitHub LFS
#   (take the first 2000 train / 200 val rows into mp20_small/{train,val}.csv)
#   ~/mattergen/.venv/bin/csv-to-dataset --csv-folder mp20_small --dataset-name mp20_small --cache-folder cache
set -x
FT=${FT:-$HOME/rtsc_targeted_work/T3-mlaudit/ft}
export OMP_NUM_THREADS=2 OPENBLAS_NUM_THREADS=2 MKL_NUM_THREADS=2
export OUTPUT_DIR=$FT/out_bandgap_small_acc${ACC:-8}_ep${EP:-3}
cd ~/mattergen
P=dft_band_gap
( while true; do nvidia-smi --query-gpu=memory.used --format=csv,noheader,nounits; sleep 2; done ) > "$FT/ft_gpumem.log" 2>/dev/null &
MON=$!
t0=$(date +%s)
nice -n 19 ~/mattergen/.venv/bin/mattergen-finetune \
  adapter.pretrained_name=mattergen_base \
  data_module=mp_20 \
  data_module.root_dir="$FT/cache/mp20_small" \
  "+lightning_module/diffusion_module/model/property_embeddings@adapter.adapter.property_embeddings_adapt.$P=$P" \
  '~trainer.logger' \
  "data_module.properties=[\"$P\"]" \
  trainer.max_epochs=${EP:-3} \
  trainer.accumulate_grad_batches=${ACC:-8} \
  trainer.check_val_every_n_epoch=${EP:-3}
rc=$?
t1=$(date +%s)
kill $MON
echo "FT_DONE rc=$rc wall_s=$((t1-t0)) peak_MiB=$(sort -n "$FT/ft_gpumem.log" | tail -1)"
