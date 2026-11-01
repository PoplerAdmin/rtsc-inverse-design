#!/bin/bash
# Waits for the q2 phonon-only run, then runs the coarse e-ph chain (q 3^3) with coarse_dfpt.py.
W=~/rtsc_targeted_work/P1-control
S=$(cd "$(dirname "${BASH_SOURCE[0]}")/../../../.." && pwd)   # repo root
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
until grep -q "JOB DONE" $W/phon_q2_e80_k8/ph.out 2>/dev/null; do sleep 30; done
python3 $S/runs/controls/mg2irh6/scripts/make_template.py 6.652064 0.2586234 80 $W/template_e80.in
date; echo "start eph q3"
python3 $S/src/scoring/coarse_dfpt.py --pw-template $W/template_e80.in --workdir $W/eph_q3_k8_kf24 \
   --kscf 8 --kfit 24 --nq 3 --sigma0 0.005 --nsigma 10 --degauss 0.02 --matdyn-nk 32 --ranks 4
date; echo "end eph q3"
