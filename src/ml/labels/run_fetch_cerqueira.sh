#!/usr/bin/env bash
# Stream all 17 batches of Materials Cloud 3kbt5-r3n56 and keep the label files (see fetch_cerqueira.py).
# ~50 GB download, ~0.4 MB kept per material.  Long: launch with setsid nohup.
set -u
export OMP_NUM_THREADS=2 OPENBLAS_NUM_THREADS=2 MKL_NUM_THREADS=2
S=$(cd "$(dirname "$0")" && pwd)
OUT=${OUT:-$HOME/rtsc_targeted_work/D1-labels/raw/mc_3kbt5/extract}
mkdir -p "$OUT"
for b in ${BATCHES:-a b c d e f g h i j k l m n o p q}; do
  if [ -f "$OUT/batch-$b.stream_hash.json" ]; then echo "skip $b"; continue; fi
  echo "start batch-$b $(date -Is)"
  for try in 1 2 3; do
    curl -sSL --retry 5 "https://archive.materialscloud.org/records/3kbt5-r3n56/files/batch-$b.tar.bz2?download=1" \
      | nice -n 19 python3 -I "$S/fetch_cerqueira.py" "$OUT" "batch-$b" && break
    echo "retry $b ($try)"; rm -rf "$OUT/batch-$b"
  done
  echo "done batch-$b $(date -Is)"
done
echo FETCH_ALL_DONE
