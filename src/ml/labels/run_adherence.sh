#!/usr/bin/env bash
# D1-labels: proxies for the generated sets (composition kNN + MatterSim relax/phonon proxy).
set -ux
W=${W:-$HOME/rtsc_targeted_work/D1-labels}
S=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)   # src/ml
PY=~/.rtt/envs/py/bin/python
export OMP_NUM_THREADS=2 OPENBLAS_NUM_THREADS=2 MKL_NUM_THREADS=2 TORCH_DISABLE_NATIVE_JIT=1
SETS="uncond A30 A100 A300 A800 A1500"
for s in $SETS; do
  d=$W/gen/$s
  [ -f $d/relax/relaxed.json ] || nice -n 19 $PY $S/mattersim_screen.py relax --out $d/relax $d/generated_crystals.extxyz
  N=$(python3 -c "import json;print(len(json.load(open('$d/relax/relaxed.json'))['rows']))")
  [ -f $d/ph/phonons.json ] || nice -n 19 $PY $S/mattersim_screen.py phonons --relaxed $d/relax/relaxed.extxyz \
      --indices $(seq 0 $((N-1))) --out $d/ph --mesh 12 > $d/ph.log 2>&1
done
nice -n 19 $PY -I $S/labels/gen_proxies.py $W/ds/manifest.tsv $W/gen/proxies_comp.json $W/gen/uncond \
  $W/gen/A30 $W/gen/A100 $W/gen/A300 $W/gen/A800 $W/gen/A1500 ~/rtsc_targeted_work/T3-mlaudit/gen_uncond64
nice -n 19 $PY -I $S/labels/gen_phonon_proxy.py $W/calib/proxy_calib.json $W/gen/proxies_phonon.json $W/gen/uncond \
  $W/gen/A30 $W/gen/A100 $W/gen/A300 $W/gen/A800 $W/gen/A1500
echo ADH_DONE
