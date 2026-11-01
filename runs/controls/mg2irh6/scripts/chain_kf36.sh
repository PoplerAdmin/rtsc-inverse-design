#!/bin/bash
# Denser electron k_fit (36^3) for the SAME q 3^3 phonons/dvscf: re-run only the e-ph sum
# (ph.x trans=.false. reads dyn* and _ph0/*.dvscf*), then q2r/matdyn. Tests the smearing
# convergence of lambda without redoing the linear response.
set -e
W=~/rtsc_targeted_work/P1-control
B=~/micromamba/envs/qe/bin
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
MPI="nice -n 19 $B/mpirun -np 4 --bind-to none"
SRC=$W/eph_q3_k8_kf24
D=$W/eph_q3_k8_kf36
F=$W/densefit36
mkdir -p $D $F
date; echo "dense scf k36"
sed 's/^ *24 24 24 0 0 0/ 36 36 36 0 0 0/' $SRC/scf_fit.in > $F/scf_fit.in
grep -q "36 36 36" $F/scf_fit.in
( cd $F && /usr/bin/time -v $MPI $B/pw.x -nk 4 -in scf_fit.in > scf_fit.out 2> scf_fit.err )
date; echo "copy coarse run"
cp $SRC/*.in $SRC/dyn[0-9]* $D/
rm -rf $D/tmp; cp -r $SRC/tmp $D/tmp
cp $F/tmp/pwscf.a2Fsave $D/tmp/pwscf.a2Fsave
sed 's/trans = .true./trans = .false./' $SRC/ph.in > $D/ph.in
grep -q "trans = .false." $D/ph.in
cd $D
date; echo "ph elphsum only"
/usr/bin/time -v $MPI $B/ph.x -nk 4 -in ph.in > ph.out 2> ph.err
date; echo "q2r/matdyn"
nice -n 19 $B/q2r.x -in q2r.in > q2r.out 2> q2r.err
nice -n 19 $B/matdyn.x -in matdyn.in > matdyn.out 2> matdyn.err
date; echo "end kf36"
