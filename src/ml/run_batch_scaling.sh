#!/usr/bin/env bash
# T3-mlaudit: does a larger batch improve MatterGen throughput? (conditional, CFG gamma=2)
set -x
W=${W:-$HOME/rtsc_targeted_work/T3-mlaudit}
S=$(cd "$(dirname "$0")" && pwd)
rm -rf "$W/gen_MgIrH128_bs128"
bash "$S/mg_generate.sh" "$W/gen_MgIrH128_bs128" chemical_system 128 1 3 "{'chemical_system':'Mg-Ir-H'}" 2.0
echo BSDONE
