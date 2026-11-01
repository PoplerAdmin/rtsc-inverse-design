#!/usr/bin/env bash
# Tier-2 cost/accuracy probe: "lambda from coarse DFPT" for Im-3m H3S at 200 GPa.
#
# Runs the same QE chain as the rtt.core DFPT workflow (scf_fit -> scf -> ph.x -> q2r -> matdyn)
# but with deliberately coarse settings, on 2 MPI ranks at nice 19, and records wall time per step.
# Reference (earlier benchmark run, 8 ranks): ecut 90/360 Ry, k_scf 18^3, k_fit 36^3, q 6^3 (16 irr. q),
# lambda(sigma=0.01 Ry) = 2.40, omega_log = 1244 K, ph.x wall 1 h 58 min.
#
# Usage: coarse_dfpt_probe.sh <workdir> <ecutwfc> <kscf> <kfit> <nq> [ph_timeout_s]
# Example: setsid nohup coarse_dfpt_probe.sh ~/rtsc_targeted_work/T4-surrogates/h3s_c1 60 8 16 2 3000 > log 2>&1 &
set -euo pipefail
WD=$1; ECUT=$2; KS=$3; KF=$4; NQ=$5; PHTO=${6:-3000}
QE=~/micromamba/envs/qe/bin
PSEUDO=~/pseudo/sg15/
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1   # 2 ranks x 1 thread = 2 threads
RUN="nice -n 19 $QE/mpirun -np 2 --bind-to none"
mkdir -p "$WD"; cd "$WD"
ECUTRHO=$((ECUT*4))

cell() { cat <<EOF
ATOMIC_SPECIES
 S   32.06500000 S_ONCV_PBE-1.2.upf
 H   1.00794000 H_ONCV_PBE-1.2.upf
CELL_PARAMETERS angstrom
     -1.49412404400000     1.49412404400000     1.49412404400000
      1.49412404400000    -1.49412404400000     1.49412404400000
      1.49412404400000     1.49412404400000    -1.49412404400000
ATOMIC_POSITIONS crystal
 S       0.0 0.0 0.0
 H       0.0 0.5 0.5
 H       0.5 0.0 0.5
 H       0.5 0.5 0.0
EOF
}
pwin() { # $1 = k, $2 = extra system lines
cat <<EOF
 &control
    calculation = 'scf'
    prefix = 'pwscf'
    outdir = './tmp/'
    pseudo_dir = '$PSEUDO'
 /
 &system
    ibrav = 0, nat = 4, ntyp = 2
    ecutwfc = $ECUT, ecutrho = $ECUTRHO
    occupations = 'smearing', smearing = 'mv', degauss = 0.02
    $2
 /
 &electrons
    conv_thr = 1.0d-10
    mixing_beta = 0.5
 /
$(cell)
K_POINTS automatic
 $1 $1 $1 0 0 0
EOF
}
pwin "$KF" "la2F = .true." > scf_fit.in
pwin "$KS" "" > scf.in
cat > ph.in <<EOF
e-ph coarse probe
 &inputph
    prefix = 'pwscf'
    outdir = './tmp/'
    tr2_ph = 1.0d-14
    fildyn = 'dyn'
    fildvscf = 'dvscf'
    electron_phonon = 'interpolated'
    el_ph_sigma = 0.01
    el_ph_nsigma = 10
    ldisp = .true.
    nq1 = $NQ, nq2 = $NQ, nq3 = $NQ
    amass(1) = 32.065
    amass(2) = 1.00794
 /
EOF
cat > q2r.in <<EOF
 &input
    fildyn = 'dyn', zasr = 'crystal', flfrc = 'ifc.fc', la2F = .true.
 /
EOF
cat > matdyn.in <<EOF
 &input
    asr = 'crystal', flfrc = 'ifc.fc', flfrq = 'freq.dos', dos = .true., la2F = .true.,
    fldos = 'phonon.dos', nk1 = 16, nk2 = 16, nk3 = 16, el_ph_nsigma = 10, ndos = 200,
    amass(1) = 32.065, amass(2) = 1.00794
 /
EOF
step() { local name=$1; shift; local t0=$(date +%s.%N)
  "$@" > "$name.out" 2> "$name.err" || echo "STEP $name FAILED rc=$?" >> timings.txt
  echo "$name $(echo "$(date +%s.%N) - $t0" | bc) s" >> timings.txt; }
: > timings.txt
echo "settings ecut=$ECUT/$ECUTRHO kscf=$KS kfit=$KF nq=$NQ ranks=2" >> timings.txt
step scf_fit $RUN $QE/pw.x -nk 2 -in scf_fit.in
step scf     $RUN $QE/pw.x -nk 2 -in scf.in
step ph      timeout "$PHTO" $RUN $QE/ph.x -nk 2 -in ph.in
step q2r     nice -n 19 $QE/q2r.x -in q2r.in
step matdyn  nice -n 19 $QE/matdyn.x -in matdyn.in
grep -A12 "Broadening" lambda >> timings.txt 2>/dev/null || true
echo DONE >> timings.txt
