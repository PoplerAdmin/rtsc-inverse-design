#!/usr/bin/env bash
# Plumbing test of the NEW-scalar-property fine-tune path (no edits to ~/mattergen).
# The label is SYNTHETIC (rtt_test_scalar = number density, atoms/A^3, computed from the CIF);
# it only proves that a custom column flows through csv-to-dataset -> finetune -> generate.
set -ex
FT=${FT:-$HOME/rtsc_targeted_work/T3-mlaudit/ft}
S=$(cd "$(dirname "$0")" && pwd)
MGPY=~/mattergen/.venv/bin/python
export OMP_NUM_THREADS=2 OPENBLAS_NUM_THREADS=2 MKL_NUM_THREADS=2
export MG_EXTRA_PROPS=rtt_test_scalar
cd "$FT"
$MGPY - <<'EOF'
import csv, os
from pymatgen.core import Structure
csv.field_size_limit(10**9)
os.makedirs("newprop_small", exist_ok=True)
for sp in ("train", "val", "test"):
    rows = list(csv.DictReader(open(f"mp20_small/{sp}.csv")))
    keep = ["material_id", "cif"]
    w = csv.DictWriter(open(f"newprop_small/{sp}.csv", "w", newline=""), fieldnames=keep + ["rtt_test_scalar"])
    w.writeheader()
    for r in rows:
        s = Structure.from_str(r["cif"], fmt="cif")
        w.writerow({"material_id": r["material_id"], "cif": r["cif"], "rtt_test_scalar": len(s) / s.volume})
EOF
nice -n 19 $MGPY "$S/mg_newprop.py" csv-to-dataset --csv-folder newprop_small --dataset-name newprop_small --cache-folder cache
export OUTPUT_DIR=$FT/out_newprop_small
cd ~/mattergen
nice -n 19 $MGPY "$S/mg_newprop.py" finetune \
  adapter.pretrained_name=mattergen_base data_module=mp_20 \
  data_module.root_dir="$FT/cache/newprop_small" \
  "hydra.searchpath=[file://$S/mg_conf]" \
  "+lightning_module/diffusion_module/model/property_embeddings@adapter.adapter.property_embeddings_adapt.rtt_test_scalar=rtt_test_scalar" \
  '~trainer.logger' 'data_module.properties=["rtt_test_scalar"]' \
  trainer.max_epochs=2 trainer.accumulate_grad_batches=16 trainer.check_val_every_n_epoch=2
nice -n 19 $MGPY "$S/mg_newprop.py" generate "$FT/../gen_newprop16" --model_path="$OUTPUT_DIR" \
  --batch_size=16 --num_batches=1 --seed=1 --record_trajectories=False \
  "--properties_to_condition_on={'rtt_test_scalar': 0.08}" --diffusion_guidance_factor=2.0
echo NEWPROP_DONE
