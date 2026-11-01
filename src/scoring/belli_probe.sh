#!/usr/bin/env bash
# Tier-1 physics proxy cost probe: Belli networking value inputs for Im-3m H3S (200 GPa cell).
# scf (k 12^3) -> pp.x ELF cube -> nscf (k 16^3) -> pp.x ILDOS at E_F (H Voronoi fraction of N(E_F)) -> phi_elf.py
# 2 MPI ranks, nice 19.  Usage: belli_probe.sh <workdir> [ecutwfc]
set -euo pipefail
WD=$1; ECUT=${2:-60}; ECUTRHO=$((ECUT*4))
QE=~/micromamba/envs/qe/bin; PSEUDO=~/pseudo/sg15/
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
RUN="nice -n 19 $QE/mpirun -np 2 --bind-to none"
HERE=$(cd "$(dirname "$0")" && pwd)
mkdir -p "$WD"; cd "$WD"
pw() { cat <<EOF
 &control
    calculation = '$1', prefix = 'pwscf', outdir = './tmp/', pseudo_dir = '$PSEUDO'
 /
 &system
    ibrav = 0, nat = 4, ntyp = 2, ecutwfc = $ECUT, ecutrho = $ECUTRHO
    occupations = 'smearing', smearing = 'mv', degauss = 0.02
    $3
 /
 &electrons
    conv_thr = 1.0d-9
 /
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
K_POINTS automatic
 $2 $2 $2 0 0 0
EOF
}
pw scf 12 "" > scf.in
pw nscf 16 "nosym = .false." > nscf.in
cat > elf.in <<EOF
 &inputpp
    prefix = 'pwscf', outdir = './tmp/', filplot = 'elf.dat', plot_num = 8
 /
 &plot
    iflag = 3, output_format = 6, fileout = 'elf.cube'
 /
EOF
step() { local name=$1; shift; local t0=$(date +%s.%N)
  "$@" > "$name.out" 2> "$name.err" || echo "STEP $name FAILED rc=$?" >> timings.txt
  echo "$name $(echo "$(date +%s.%N) - $t0" | bc) s" >> timings.txt; }
: > timings.txt
step scf     $RUN $QE/pw.x -nk 2 -in scf.in
step elf     nice -n 19 $QE/pp.x -in elf.in
step nscf    $RUN $QE/pw.x -nk 2 -in nscf.in
# SG15 pseudopotentials carry no atomic wavefunctions, so projwfc.x cannot run ("Cannot project on
# zero atomic wavefunctions").  Substitute: ILDOS in E_F +/- 0.136 eV (0.01 Ry), split by H Voronoi cells.
EF=$(grep "the Fermi energy is" nscf.out | tail -1 | awk '{print $5}')
cat > ildos.in <<EOF2
 &inputpp
    prefix = 'pwscf', outdir = './tmp/', filplot = 'ildos.dat', plot_num = 10,
    emin = $(echo "$EF - 0.136" | bc), emax = $(echo "$EF + 0.136" | bc)
 /
 &plot
    iflag = 3, output_format = 6, fileout = 'ildos.cube'
 /
EOF2
step ildos   nice -n 19 $QE/pp.x -in ildos.in
step phi     nice -n 19 ~/.rtt/envs/py/bin/python "$HERE/phi_elf.py" elf.cube --ildos-cube ildos.cube
cat phi.out >> timings.txt
echo DONE >> timings.txt
