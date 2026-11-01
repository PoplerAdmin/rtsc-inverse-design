#!/usr/bin/env bash
# T3-mlaudit smoke test: download two checkpoints, then 3 generation runs of 64 structures each.
# Launch: setsid nohup bash run_smoke_generation.sh > run_gen.log 2>&1 < /dev/null &
set -x
W=${W:-$HOME/rtsc_targeted_work/T3-mlaudit}
S=$(cd "$(dirname "$0")" && pwd)
cd ~/mattergen
time ~/mattergen/.venv/bin/python -c "
from huggingface_hub import snapshot_download as s
for n in ('mattergen_base', 'chemical_system'):
    s('microsoft/mattergen', allow_patterns=f'checkpoints/{n}/*')
"
bash "$S/mg_generate.sh" "$W/gen_uncond64" mattergen_base 64 1 1
bash "$S/mg_generate.sh" "$W/gen_LiMgH64" chemical_system 64 1 1 "{'chemical_system':'Li-Mg-H'}" 2.0
bash "$S/mg_generate.sh" "$W/gen_MgIrH64" chemical_system 64 1 2 "{'chemical_system':'Mg-Ir-H'}" 2.0
echo ALLDONE
