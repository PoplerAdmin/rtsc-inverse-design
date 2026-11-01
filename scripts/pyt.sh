#!/bin/bash
# Run pytest with the compute rules of docs/CONTEXT.md (nice 19, <= 2 threads).
export OMP_NUM_THREADS=2 OPENBLAS_NUM_THREADS=2 MKL_NUM_THREADS=2
cd "$(dirname "$0")/.."
exec nice -n 19 ~/.rtt/envs/py/bin/python -m pytest "$@"
