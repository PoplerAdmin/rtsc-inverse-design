#!/usr/bin/env bash
# T3-mlaudit smoke test, stage 2: MatterSim relax (GPU) + E_hull vs Alex-MP DFT reference.
# Launch: setsid nohup bash run_smoke_screen.sh > run_screen.log 2>&1 < /dev/null &
set -x
export OMP_NUM_THREADS=2 OPENBLAS_NUM_THREADS=2 MKL_NUM_THREADS=2
# torch 2.14 routes bmm (used in MatterSim force backward) to a Triton kernel that needs a C
# compiler, which WSL does not have -> RuntimeError on GPU. This env var restores the ATen path.
export TORCH_DISABLE_NATIVE_JIT=1
W=${W:-$HOME/rtsc_targeted_work/T3-mlaudit}
S=$(cd "$(dirname "$0")" && pwd)
PY=~/.rtt/envs/py/bin/python
MG=~/mattergen/.venv/bin/python
REF=$W/ref/reference_MP2020correction.lmdb
for tag in uncond64 LiMgH64 MgIrH64; do
  nice -n 19 $PY "$S/mattersim_screen.py" relax --out "$W/relax_$tag" "$W/gen_$tag/generated_crystals.extxyz"
  nice -n 19 $MG "$S/ehull_ref.py" --lmdb "$REF" --relaxed "$W/relax_$tag/relaxed.json" \
      --extxyz "$W/relax_$tag/relaxed.extxyz" --out "$W/relax_$tag/ehull.json"
done
echo ALLDONE
