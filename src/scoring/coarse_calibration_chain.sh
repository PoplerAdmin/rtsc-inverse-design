#!/usr/bin/env bash
# K1 calibration chain: one QE job at a time, 6 ranks, nice 19.
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
CD=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)/coarse_dfpt.py
cd ~/rtsc_targeted_work/K1-calibrate
PY=~/.rtt/envs/py/bin/python
declare -A SIG=( [Al]=0.005 [Pb]=0.01 [Nb]=0.01 )
declare -A E23=( [Al]=40 [Pb]=34 [Nb]=80 )
run() { # tag M ecut ks kf nq
  local tag=$1 M=$2 E=$3 KS=$4 KF=$5 NQ=$6
  local wd=runs/${M}_${tag}_q${NQ}
  if [ -f $wd/coarse_run.json ] && grep -q '"ok": true' $wd/coarse_run.json; then echo "skip $wd"; return; fi
  local ea=""; [ "$E" != ref ] && ea="--ecutwfc $E"
  echo "$(date +%T) start $wd"
  $PY $CD --pw-template ref/$M/scf_fit.in --workdir $wd $ea --kscf $KS --kfit $KF --nq $NQ \
     --sigma0 ${SIG[$M]} --nsigma 10 --ranks 6 --ph-timeout 7200 > $wd.log 2>&1
  echo "$(date +%T) end $wd rc=$?"
  rm -rf $wd/tmp/_ph0/*.dvscf* $wd/tmp/_ph0/*/  2>/dev/null   # keep disk small; spectra/dyn kept
}
mkdir -p runs
for M in Al Nb Pb; do for q in 2 3 4; do run P1 $M ref 12 24 $q; done; done
for M in Al Nb Pb; do for q in 3 4; do run P2 $M ${E23[$M]} 12 24 $q; done; done
for M in Al Nb Pb; do for q in 2 4; do run P3 $M ref 8 16 $q; done; done
echo CHAIN_DONE
