"""Run MatterGen CLIs with EXTRA scalar property ids, without editing the ~/mattergen checkout.

MatterGen only accepts property names listed in mattergen.common.utils.globals.PROPERTY_SOURCE_IDS
(README "Fine-tune on your own property data", step 1).  Every module imports that same list
object, so appending to it in-process before the CLI entry point runs is equivalent to editing
globals.py.  The property-embedding YAML (README step 4) is supplied from src/ml/mg_conf/ via a
Hydra search-path override instead of being copied into the checkout.

Run with the MatterGen venv, from the repository root:
  MG_EXTRA_PROPS=rtt_test_scalar ~/mattergen/.venv/bin/python src/ml/mg_newprop.py csv-to-dataset \
      --csv-folder DIR --dataset-name NAME --cache-folder CACHE
  MG_EXTRA_PROPS=rtt_test_scalar ~/mattergen/.venv/bin/python src/ml/mg_newprop.py finetune \
      adapter.pretrained_name=mattergen_base data_module=mp_20 data_module.root_dir=CACHE/NAME \
      "hydra.searchpath=[file://$PWD/src/ml/mg_conf]" \
      '+lightning_module/diffusion_module/model/property_embeddings@adapter.adapter.property_embeddings_adapt.rtt_test_scalar=rtsc_scalar' \
      '~trainer.logger' 'data_module.properties=["rtt_test_scalar"]' ...
  MG_EXTRA_PROPS=rtt_test_scalar ~/mattergen/.venv/bin/python src/ml/mg_newprop.py generate OUT \
      --model_path=RUN_DIR --properties_to_condition_on="{'rtt_test_scalar': 0.08}" ...
"""
import os
import sys

from mattergen.common.utils import globals as g

for p in filter(None, os.environ.get("MG_EXTRA_PROPS", "").split(",")):
    if p not in g.PROPERTY_SOURCE_IDS:
        g.PROPERTY_SOURCE_IDS.append(p)

cmd, sys.argv = sys.argv[1], [sys.argv[1]] + sys.argv[2:]
if cmd == "csv-to-dataset":
    from mattergen.scripts.csv_to_dataset import main
    main()
elif cmd == "finetune":
    from mattergen.scripts.finetune import mattergen_finetune
    mattergen_finetune()
elif cmd == "generate":
    from mattergen.scripts.generate import _main
    _main()
else:
    raise SystemExit(f"unknown command {cmd}")
