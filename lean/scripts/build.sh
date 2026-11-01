#!/usr/bin/env bash
# Sync lean/ sources into an ext4 build dir (whose .lake/packages holds a prebuilt Mathlib
# v4.35.0) and build with 2 threads at nice 19. Usage: build.sh <src lean dir> <build dir> [targets...]
set -euo pipefail
SRC=${1:?src}; BLD=${2:?build dir}; shift 2
export PATH=$HOME/.elan/bin:$PATH
export LEAN_NUM_THREADS=2 OMP_NUM_THREADS=2 OPENBLAS_NUM_THREADS=2 MKL_NUM_THREADS=2
mkdir -p "$BLD"
rsync -a --delete --exclude .lake --exclude logs "$SRC"/ "$BLD"/
cd "$BLD"
lake --version
# lake 5.0 has no -j flag: pin to 2 CPUs instead (LEAN_NUM_THREADS=2 as well)
taskset -c 14,15 nice -n 19 lake build "$@"
