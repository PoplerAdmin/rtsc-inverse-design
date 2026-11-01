#!/usr/bin/env bash
# D1-labels: sample the rtt_log10A fine-tuned model unconditionally and at several A targets.
# 64 structures per setting, guidance gamma=2 (as in docs/design/10), seed 1, best checkpoint.
set -ux
W=${W:-$HOME/rtsc_targeted_work/D1-labels}
S=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)   # src/ml
RUN=${RUN:-$W/ft_log10A}
G=${G:-2.0}
CK=${CK:-best}
export OMP_NUM_THREADS=2 OPENBLAS_NUM_THREADS=2 MKL_NUM_THREADS=2 MG_EXTRA_PROPS=rtt_log10A
cd ~/mattergen
gen() {  # out props
  local out=$1; shift
  [ -f "$out/generated_crystals.extxyz" ] && return 0
  mkdir -p "$out"; t0=$(date +%s)
  nice -n 19 ~/mattergen/.venv/bin/python "$S/mg_newprop.py" generate "$out" --model_path="$RUN" \
    --checkpoint_epoch=$CK --batch_size=64 --num_batches=1 --seed=1 --record_trajectories=False "$@" \
    > "$out/mattergen.log" 2>&1
  echo "$out rc=$? wall_s=$(( $(date +%s) - t0 ))"
}
gen "$W/gen/uncond"
for A in ${TARGETS:-30 100 300 800 1500}; do
  L=$(python3 -c "import math;print(round(math.log10($A),4))")
  gen "$W/gen/A$A" "--properties_to_condition_on={'rtt_log10A':$L}" --diffusion_guidance_factor=$G
done
echo GEN_ALL_DONE
