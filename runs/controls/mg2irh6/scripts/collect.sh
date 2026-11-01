#!/bin/bash
# Copy small text inputs/outputs of the Mg2IrH6 control into the repo and write a sha256
# manifest of everything (incl. big outputs that stay in scratch). Paths are written relative
# to the scratch root so no home path ends up in the repo.
set -e
W=~/rtsc_targeted_work/P1-control
R=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)   # runs/controls/mg2irh6
cd $W
for d in ecut kconv relax_e80 relax_e100 relax_e80_from_ms2 phon_q2_e80_k8 phon_G_e80_k12 phon_G_e80_k16 eph_q3_k8_kf24 eph_q3_k6_kf36 eph_q3_k8_kf36_rerun_INVALID eph_q3_k8_kf36_INVALID_noncommensurate densefit36; do
  [ -d $d ] || continue
  mkdir -p $R/$d
  for f in $d/*.in $d/*.json $d/dyn[0-9]* $d/lambda $d/a2F.dos* $d/phonon.dos $d/freq.mesh.gp; do
    [ -f "$f" ] || continue
    sz=$(stat -c %s "$f")
    [ "$sz" -lt 300000 ] && sed "s#$HOME#~#g" "$f" > $R/$f
  done
  # pw/ph stdout: keep only a short tail-free summary (they can be MBs)
  for f in $d/*.out; do
    [ -f "$f" ] || continue
    grep -E "^!|P=|Fermi energy|number of k points|freq \(|lambda\(|WALL|JOB DONE|CELL_PARAMETERS|^(Ir|Mg|H) " "$f" \
      | sed "s#$HOME#~#g" > $R/$f.summary || true
  done
done
cp $W/mattersim_gamma.json $R/ 2>/dev/null || true
find . -type f -print0 | sort -z | xargs -0 nice -n 19 sha256sum > $R/sha256_manifest.txt   # includes tmp/ (wfc, dvscf, a2Fsave)
du -sh $W/*/ 2>/dev/null | sed "s#$HOME#~#g" > $R/scratch_sizes.txt
echo collected
