#!/usr/bin/env bash
# Axiom audit + source scan. Usage: check.sh <build dir> <out dir>
set -euo pipefail
BLD=${1:?}; OUT=${2:?}
export PATH=$HOME/.elan/bin:$PATH LEAN_NUM_THREADS=2
mkdir -p "$OUT"; cd "$BLD"
for f in Check CheckPipeline CheckMultiBand CheckTarget; do
  taskset -c 14,15 nice -n 19 lake env lean $f.lean > "$OUT/${f}_axioms.txt" 2>&1
  n=$(grep -c "depends on axioms" "$OUT/${f}_axioms.txt" || true)
  bad=$(grep "depends on axioms" "$OUT/${f}_axioms.txt" \
        | sed -e 's/.*\[//' -e 's/\].*//' | tr ',' '\n' | sed 's/ //g' \
        | grep -v -x -E 'propext|Classical.choice|Quot.sound|' | sort -u || true)
  nodep=$(grep -c "does not depend on any axioms" "$OUT/${f}_axioms.txt" || true)
  err=$(grep -c -i "error\|sorry" "$OUT/${f}_axioms.txt" || true)
  echo "$f: $n theorems with axioms, $nodep axiom-free, disallowed axioms: [${bad}], error/sorry lines: $err"
done
echo "source scan (sorry/admit/axiom/native_decide/implemented_by/extern):"
grep -n -E '\bsorry\b|\badmit\b|^\s*axiom\b|native_decide|implemented_by|@\[extern' \
  TcR1/*.lean PipelineLogic/*.lean MultiBand/*.lean TargetLogic/*.lean *.lean || echo "  none"
