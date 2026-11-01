#!/bin/bash
# Whole rtt test suite in parts (docs/CONTEXT.md: nice 19, <= 2 threads, never one big run).
cd "$(dirname "$0")/.."
parts=(tests/test_core_pins.py tests/test_decide_core.py tests/test_decide_lean_table.py
       tests/test_decide_lean_cases.py tests/test_ledger_triggers.py tests/test_engine.py tests/test_certify_e2e.py
       tests/front tests/back tests/review
       "tests/core/test_certify_types.py tests/core/test_dfpt_smearing_admission.py tests/core/test_dfpt_soundness.py tests/core/test_dfpt_parsers.py"
       tests/core/test_certify.py tests/core/test_certify_multiband.py
       tests/core/test_certify_multiband_review.py tests/core/test_correspondence.py
       tests/core/audit_certify/test_forgery.py tests/core/oracle/test_oracle_vs_certify.py)
for p in "${parts[@]}"; do
  echo "== $p"
  bash scripts/pyt.sh -q -p no:cacheprovider $p 2>&1 | grep -E "^(FAILED|ERROR)|passed|failed|error" | tail -15
done
echo "== DONE"
