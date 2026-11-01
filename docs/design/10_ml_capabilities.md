# 10 — What MatterGen and MatterSim can do on this machine

Agent T3-mlaudit, 2026-10-28. Every number here was **measured** on this machine (RTX 4090 24 GB,
WSL2) unless it is marked *cited* or **UNVERIFIED**. Scripts are in `src/ml/`. Raw outputs are in
`~/rtsc_targeted_work/T3-mlaudit/`, which is scratch and not in the repo. ML results are
**heuristic only**: they rank, filter or defer candidates and never issue a verdict
(contract S2).

## 0. Summary

| Stage | Works? | Measured cost (RTX 4090) |
|---|---|---|
| MatterGen unconditional (`mattergen_base`) | yes | **4.9 s/structure** at batch 64; 5.1 GB GPU |
| MatterGen + `chemical_system` (CFG γ=2) | yes; 58/64 (Li-Mg-H) and 63/64 (Mg-Ir-H) samples are exactly in the requested ternary, 64/64 inside it | **10.0–10.4 s/structure** at batch 64; 10–11 GB. Batch 128 is no faster (1.33 s/it vs 0.64 s/it) and uses 16.4 GB |
| MatterSim-5M batched relax (FIRE + FrechetCellFilter, fmax 0.02) | yes, **only with `TORCH_DISABLE_NATIVE_JIT=1`** (§3.1) | **0.6–0.7 s/structure** (64 × ~9 atoms); 0.21 s/structure for 283 reference cells; 2–8 GB |
| E_hull vs DFT reference (Alex-MP, 846k entries, MP2020) | yes, offline after a one-time 0.87 GB download | seconds per chemical system |
| E_hull vs MatterSim-only hull (self-consistent) | yes | relax of reference entries ~1 min per ternary |
| MatterSim phonons (phonopy, finite displacement) | runs; **not reliable for hydrides** (§3.4) | **2–6 s/structure** (120–243-atom supercells) |
| MatterGen fine-tune on a NEW scalar | **yes, end to end** (plumbing test with a synthetic label) | ~125–170 training samples/s; batch ≤32 (batch 32 = 20.4 GB). Batch ≥64 overflows the 24 GB card and slows down by more than 10× |

Bottom line: with the shipped checkpoints, generation is limited to chemistry, space group,
E_hull, band gap, magnetic density and bulk modulus. A superconductivity-aware generator needs
fine-tuning. The GPU cost of that is small (hours); the binding constraint is labels (§4).
MatterSim works for relaxation and hull pre-screening. Its harmonic phonons gave a **false
dynamical instability** for the best-known ambient-pressure hydride candidate, Mg2IrH6 (§3.4),
so MatterSim phonons must not be used to reject hydrides.

## 1. Installation and checkpoints

- **MatterGen**: source checkout `~/mattergen` at commit `9242366` ("Fix Linux PyPI installation
  dependencies (#252)"), installed editable in its own venv `~/mattergen/.venv` (Python 3.10,
  uv; torch 2.2.1+cu121, pymatgen 2024.10.29, mattersim 1.1.2, lmdb). CLIs:
  `mattergen-generate`, `mattergen-finetune`, `mattergen-train`, `mattergen-evaluate`,
  `csv-to-dataset`.
- **MatterSim**: main venv `~/.rtt/envs/py` (Python 3.12, mattersim 1.2.5, torch 2.14.1+cu130,
  phonopy 4.8.6). Local models are `~/.local/mattersim/pretrained_models/mattersim-v1.0.0-{1M,5M}.pth`.
  Everything here uses 5M.
- **Checkpoints**: `~/mattergen/checkpoints/*/checkpoints/last.ckpt` are Git-LFS pointer files
  (134 bytes; `git lfs` is not installed). Real weights come from Hugging Face
  `microsoft/mattergen` into `~/.cache/huggingface/hub` on first use, about 0.5 GB each, in about
  6 s each. Downloaded so far: `mattergen_base`, `chemical_system`,
  `chemical_system_energy_above_hull`. Available on demand: `mp_20_base`, `space_group`,
  `dft_mag_density`, `dft_band_gap`, `ml_bulk_modulus`, `dft_mag_density_hhi_score`
  (README §"Pretrained models").
- **No pressure conditioning exists.** Every MatterGen model is trained on 0 GPa
  Alex-MP-20 / MP-20 structures with ≤20 atoms per cell and E_hull < 0.1 eV/atom
  (`MODEL_CARD.md`). It excludes noble gases, Tc, Pm and Z > 84.

### Prior experiments

- An earlier MatterGen test (`adapter_runs/G-AMB_B-Mg_sgany_P0_n16_s7_*`, not in the repo): one generation run of
  `chemical_system_energy_above_hull` with `{'chemical_system':'B-Mg','energy_above_hull':0.05}`,
  γ=2. It produced 16 structures in 199 s (12.4 s/structure at batch 16), all 16 in Mg-B. Nothing
  further was done with them.
- The hull cache from the same test: MatterSim "hulls" made only from elemental and MgB2 references, with
  fallback fcc or dimer structures; the H-S lists are empty. They are **not** usable as a
  reference hull.
- An earlier fine-tune export script documented the CSV/YAML fine-tune format but was never
  run. Its format description matches what §4.2 executed. **No fine-tune had been attempted
  before this audit.**

## 2. Conditioning mechanism (MatterGen; Zeni et al., *Nature* 2025, doi:10.1038/s41586-025-08628-5)

Read from the code at `9242366`:

- **Base model**: a diffusion model over atom types, fractional coordinates and lattice, with a
  GemNet-T score network (`mattergen/denoiser.py`), trained unconditionally.
- **Adapters** (`mattergen/adapter.py`, `common/gemnet/gemnet_ctrl.py`): fine-tuning adds one
  `PropertyEmbedding` per new property plus, per message-passing block, an adapt layer and a
  *mixin* linear layer **initialized to zero**. At step 0 the fine-tuned model is therefore exactly
  the base model. The unconditional embedding of an adapter property is forced to zero
  (`ZerosEmbedding`), so the unconditional score is unchanged. The default config has
  `full_finetuning: true` (`conf/adapter/default.yaml`), so all weights train, at lr 5e-6 and for
  200 epochs by default (`conf/finetune.yaml`).
- **Float properties** go through a `StandardScalerTorch` and a `NoiseLevelEncoding` (sinusoidal)
  embedding (`conf/.../property_embeddings/dft_band_gap.yaml`). Labels may be sparse: rows that
  lack the property are dropped by `filter_sparse_properties`.
- **Classifier-free guidance** (`diffusion/sampling/classifier_free_guidance.py`):
  score = γ·s(x|y) + (1−γ)·s(x), i.e. sampling p(x)·p(y|x)^γ. γ=0 is unconditional and γ=1 is
  pure conditional. Each step evaluates both branches, which is why conditional sampling is
  measured at 2.0–2.1× the unconditional cost.
- **Limits** (`MODEL_CARD.md`, cited): fine-tuning needs "at least several thousands of labeled
  structures" with good coverage of the target range. Performance degrades for extreme property
  values with few training examples. That matters here: Tc ≈ 300 K is far outside any DFT label
  set (§4).

## 3. Smoke test results

### 3.1 Commands that work (copy-paste)

```bash
# --- generation (MatterGen venv). mg_generate.sh records wall time + peak GPU memory ---
S=<repo>/src/ml; W=~/rtsc_targeted_work/T3-mlaudit
bash $S/mg_generate.sh $W/gen_uncond64 mattergen_base 64 1 1
bash $S/mg_generate.sh $W/gen_LiMgH64  chemical_system 64 1 1 "{'chemical_system':'Li-Mg-H'}" 2.0
bash $S/mg_generate.sh $W/gen_MgIrH64  chemical_system 64 1 2 "{'chemical_system':'Mg-Ir-H'}" 2.0
#   (equivalent raw call)
~/mattergen/.venv/bin/mattergen-generate OUT --pretrained_name=chemical_system --batch_size=64 \
   --num_batches=1 --seed=1 --record_trajectories=False \
   "--properties_to_condition_on={'chemical_system':'Li-Mg-H'}" --diffusion_guidance_factor=2.0
# all three in one background job: src/ml/run_smoke_generation.sh

# --- MatterSim relax (main venv, GPU) ---
export OMP_NUM_THREADS=2 OPENBLAS_NUM_THREADS=2 MKL_NUM_THREADS=2 TORCH_DISABLE_NATIVE_JIT=1
PY=~/.rtt/envs/py/bin/python
nice -n 19 $PY $S/mattersim_screen.py relax --out $W/relax_LiMgH64 $W/gen_LiMgH64/generated_crystals.extxyz
#   options: --pressure-GPa P (enthalpy relax), --fmax 0.02, --max-steps 1000, --model MatterSim-v1.0.0-5M

# --- E_hull vs DFT reference (MatterGen venv) ---
~/mattergen/.venv/bin/python -c "from mattergen.evaluation.reference.presets import \
  get_reference_dataset_path as g; print(g('reference_MP2020correction.gz'))"     # 0.87 GB, once
zcat <printed path> > $W/ref/reference_MP2020correction.lmdb                          # 3.85 GB, once
~/mattergen/.venv/bin/python $S/ehull_ref.py --lmdb $W/ref/reference_MP2020correction.lmdb \
   --relaxed $W/relax_LiMgH64/relaxed.json --extxyz $W/relax_LiMgH64/relaxed.extxyz --out $W/relax_LiMgH64/ehull.json

# --- E_hull vs MatterSim-only hull ---
~/mattergen/.venv/bin/python $S/ehull_ref.py --lmdb $W/ref/reference_MP2020correction.lmdb \
   --export-ref H-Li-Mg $W/ref/ref_H-Li-Mg.extxyz
nice -n 19 $PY $S/mattersim_screen.py relax --out $W/relax_ref_H-Li-Mg $W/ref/ref_H-Li-Mg.extxyz
$PY $S/ehull_self.py --ref $W/relax_ref_H-Li-Mg/relaxed.json --cand $W/relax_LiMgH64/relaxed.json \
   --out $W/relax_LiMgH64/ehull_self.json
# relax + both hulls for all three sets: src/ml/run_smoke_screen.sh

# --- phonons ---
nice -n 19 $PY $S/mattersim_screen.py phonons --relaxed $W/relax_MgIrH64/relaxed.extxyz \
   --indices 53 49 33 --out $W/ph_MgIrH        # supercell auto (>=9 A, <=250 atoms), 0.01 A, 16^3 mesh
```

**Blocking issue found and fixed.** In `~/.rtt/envs/py`, torch 2.14 routes the `bmm` in
MatterSim's force backward pass to a Triton "native" kernel. That kernel needs a C compiler, and
the WSL image has none, so every CUDA MatterSim force/stress call fails with
`RuntimeError: Failed to find C compiler`. Setting `TORCH_DISABLE_NATIVE_JIT=1` before torch is
imported restores the ATen path. `mattersim_screen.py` sets it itself. **Any other agent using
MatterSim on the GPU in that venv will hit this.** The earlier CPU-only MLFF prescreen did
not.

### 3.2 Generation

| Run | Model | γ | N | Wall | s/structure | Peak GPU |
|---|---|---|---|---|---|---|
| gen_uncond64 | mattergen_base | — | 64 | 314 s | 4.91 | 5.1 GB |
| gen_LiMgH64 | chemical_system Li-Mg-H | 2 | 64 | 643 s | 10.04 | 10.2 GB |
| gen_MgIrH64 | chemical_system Mg-Ir-H | 2 | 64 | 667 s | 10.43 | 11.2 GB |
| gen_MgIrH128 (stopped at 9 %) | same, batch 128 | 2 | — | 1.33 s/it × 1000 steps ≈ 22 min | ≈10.4 | 16.4 GB |
| gen_ft_bandgap16 (fine-tuned, §4.2) | own checkpoint | 2 | 16 | 144 s | 9.0 | — |

The sampler runs 1000 predictor-corrector steps. The GPU is already saturated at batch 64, so
larger batches give no speed-up. Sustained capacity is ≈ 8.6k conditional or 17k unconditional
structures per GPU-day.

The unconditional run gave 64 diverse compositions (1–20 atoms, median 8). The
conditional-run quality checks used pymatgen `StructureMatcher` on the MatterSim-relaxed cells,
with the Alex-MP reference entries of the same system as the known set:

| | in exact ternary | unique | unique and not in Alex-MP |
|---|---|---|---|
| Li-Mg-H | 58/64 | 61 | 58 |
| Mg-Ir-H | 63/64 | 50 | 46 |

Mg-Ir-H output is repetitive: 7× IrMgH2, 5× IrMg2H, 9× Ir2Mg2H6. It rediscovered **Fm-3m
Mg2IrH6** (it matches the Alex-MP entry) and produced two copies of **Mg2IrH5** (I4mm at symprec
0.1). Mg2IrH5 is absent from the Alex-MP reference. That it is the experimentally known compound
is **UNVERIFIED** (I recall an experimental Mg2IrH5 but did not check its structure).

### 3.3 MatterSim relaxation and E_hull

The batched relaxation (`mattersim.applications.batch_relax.BatchRelaxer`) is fast:

| Set | N (atoms) | Relax wall | s/structure | Ionic steps | Converged | Peak alloc |
|---|---|---|---|---|---|---|
| uncond64 | 64 (564) | 37 s | 0.58 | 5025 | 61/64 | 2.1 GB |
| LiMgH64 | 64 (564) | 43 s | 0.68 | 5335 | 62/64 | 4.5 GB |
| MgIrH64 | 64 (574) | 47 s | 0.74 | 5727 | 62/64 | 5.3 GB |
| Alex-MP refs H-Li-Mg | 283 (2126) | 60 s | 0.21 | 10196 | 280/283 | 6.7 GB |
| Alex-MP refs H-Ir-Mg | 118 (1517) | 72 s | 0.61 | 12791 | 114/118 | 8.1 GB |

Sequential single-structure relaxation is ~17 s/structure, measured on 3 structures while the
GPU was shared. Always use the batched path.

**Hull data available offline.** No Materials Project API key is configured, and the earlier
`hull_cache` is not a hull (§1). Microsoft's MatterGen reference set (`reference_MP2020correction.gz`,
845,997 Alexandria + MP PBE(+U) entries with MP2020 corrections, 208,307 chemical systems) is
now downloaded and decompressed to an LMDB, `~/rtsc_targeted_work/T3-mlaudit/ref/`. It is
memory-mapped and queried per chemical system (≈0.5 GB RSS). `mattergen`'s own
`LMDBBackedReferenceDatasetImpl` fails on this lmdb version (it opens the environment twice), so
`ehull_ref.py` has its own reader. The entries are pickled. Unpickling a downloaded file is only
acceptable because it is Microsoft's official data release. The set holds 283 entries for
H-Li-Mg and its subsystems and 118 for H-Ir-Mg, including Mg2IrH6, MgIrH3 and LiMgH3.

Two E_hull definitions were computed (eV/atom, 0 GPa):

| Set | Hull | ≤0 | ≤0.05 | ≤0.1 | median |
|---|---|---|---|---|---|
| uncond64 | DFT ref + MatterSim/MP2020 | 0 | 17 | 43 | 0.068 |
| LiMgH64 | DFT ref + MatterSim/MP2020 | 0 | 2 | 12 | 0.174 |
| LiMgH64 | MatterSim-only | 0 | 10 | 25 | 0.119 |
| MgIrH64 | DFT ref + MatterSim/MP2020 | 5 | 17 | 29 | 0.112 |
| MgIrH64 | MatterSim-only | 4 | 15 | 33 | 0.094 |

(MgIrH64 lowest values: Mg2IrH5 −0.047 on the MatterSim-only hull, which lacks Mg2IrH5;
Mg2IrH6 −0.0003.)

**Calibration of mixing MatterSim with DFT**: each reference entry relaxed with MatterSim, MP2020
applied, compared with its own DFT energy.

- H-Li-Mg: MAE 27 meV/atom (H-containing 43, H-free 23).
- H-Ir-Mg: MAE 64 meV/atom (H-containing 74, H-free 56).
- Systematic: every MP-sourced Mg entry (mp-110, mp-153, …) is +0.10 eV/atom off MatterSim,
  while the Alexandria Mg entries agree to ≤0.01. MgH2 entries scatter from +0.02 to +0.17. The
  reference set is internally inhomogeneous for Mg compounds (probably MP vs Alexandria
  calculation settings, **UNVERIFIED**).

**Recommendation**: rank with the **MatterSim-only hull**, where systematic errors cancel to
first order, and use the DFT-reference hull as a second opinion. Treat |E_hull| < ~0.07 eV/atom
as indistinguishable from 0 at this level. Any finalist needs its own DFT hull. Neither hull
covers P > 0: the DFT reference is a 0 GPa set, and `ehull_ref.py` refuses P > 0. The
MatterSim-only hull works at any pressure if references are relaxed at the same pressure.

### 3.4 MatterSim phonons — throughput and a red flag

| Structure | Supercell (atoms) | Displacements | Wall | min ω (cm⁻¹) | H-weight above 1000 cm⁻¹ |
|---|---|---|---|---|---|
| Fm-3m Mg2IrH6 (gen #53) | 3×3×3 (243) | 54 | 5.1 s | **−363 (L)** | 0.33 |
| Mg2IrH5 I4mm (gen #49) | 3×3×3 (216) | 48 | 2.9 s | ≥0 (stable) | 0.33 |
| MgIrH3 Cm (gen #33) | 3×2×2 (120) | 60 | 2.2 s | −91 | 0.33 |
| LiMgH3 P4₂/mnm (gen #35) | 2×2×2 (160) | 120 | 5.9 s | ≥0 (stable) | 0.45 |

Throughput is ≈ 0.05–0.1 s per supercell force call, so ~2–6 s per structure, or ~15–40k
structures per GPU-day.

**Mg2IrH6 cross-check** (`src/ml/check_mg2irh6.py`; *cited*: Dolui, Conway, Heil, Strobel,
Prasankumar, Pickard, *PRL* 132, 166001 (2024), arXiv:2310.07562):

- **Literature**: PBE a = 6.655 Å, Fm-3m. Imaginary modes were excluded in their screening, and
  SSCHA "do[es] not alter the dynamic stability at ambient pressure". λ ≈ 2.5, Tc 160–175 K.
- **Lattice**:
  - MatterSim relaxing the Alex-MP DFT structure (a_DFT = 6.633 Å) gives a = 6.713 Å.
  - The MatterGen sample relaxed instead to a *second* Fm-3m minimum, a = 6.528 Å
    (Ir–H 1.64 Å vs 1.73 Å), within 0.3 meV/atom of the first.
- **Phonons**: MatterSim finds strong instabilities in **both** geometries:
  - at a = 6.713 Å: Γ −255 cm⁻¹ (triply degenerate) and X −287 cm⁻¹;
  - with the cell fixed at the PBE a = 6.655 Å: Γ −245 and X −279 cm⁻¹.

MatterSim therefore gives a **false negative for dynamical stability** on exactly the compound
class this project targets. Its PES for the [IrH6] complex is also multi-welled where DFT has
one minimum. This fits the literature:

- MatterSim's phonon validation is a max-frequency MAE of 0.87 THz on PhononDB (Yang et al.
  2024, arXiv:2405.04967, *cited*).
- The paper reports no hydride or electron-phonon validation. Its Fig. 3 analysis even excludes
  H-containing materials over anion-correction concerns (*cited*).
- Its pressure range claim is "up to 1000 GPa" (*cited*). No independent hydride benchmark was
  found; Loew et al. arXiv:2412.16551 and Anam et al. arXiv:2509.03401 benchmark generic phonons.

**Rule for the pipeline: MatterSim phonons may promote a hydride or reorder candidates; they may
never reject one.** An imaginary MatterSim mode means "check with DFPT", not "unstable". Mg2IrH6
should serve as a positive control in every ML filter we build.

The "H-weight above 1000 cm⁻¹" column is a crude proxy: it is the |eigenvector|² share of H
coordinates in modes above 1000 cm⁻¹. In octahedral MH6 complexes it is ≈ 1/3 by mode counting
(6 stretches out of 18 H modes), so it carries little information beyond stoichiometry and
should not be used as a label.

## 4. Fine-tuning MatterGen on a superconductivity scalar

### 4.1 What is needed

The pathway is: a scalar label y per structure, then an adapter fine-tune, then sampling with
`{'y': target}` and γ > 1. The MODEL_CARD asks for "several thousands" of labels with
coverage of the target range.

Datasets (*cited* unless marked):

| Dataset | Size | Labels | Fit for this goal |
|---|---|---|---|
| JARVIS-DFT EPC (Choudhary & Garrity, *npj Comput. Mater.* 8, 244 (2022), arXiv:2205.00060) | 1,058 EPC calculations out of 1,736 pre-screened | λ, ω_log, McMillan–Allen–Dynes Tc; 105 with Tc > 5 K | ambient, few hydrides, Tc mostly < 30 K; too small on its own. Exact dataset key in jarvis-tools **UNVERIFIED** |
| Cerqueira et al., *Adv. Mater.* 36, e2307085 (2024), arXiv:2307.10728; Materials Cloud 2023.163 | ~7,000 EPC calculations ("modest convergence") + ML Tc for ~200,000 metals; 541 DFPT-confirmed Tc > 10 K | λ, ω_log, Tc (DFT); Tc (ML) | the most useful: Alexandria-based, so it matches MatterGen's training distribution (≤20 atoms, near hull). Labels nearly all at Tc ≪ 300 K |
| SuperCon (NIMS) | ~30k compositions | experimental Tc, often no structure | composition-only; unsuitable without structure matching |
| In-house | — | A = ∫α²F dω, λ, certified Tc bracket from the `rtt.core` DFPT + certify pipeline | each label costs CPU-hours to CPU-days of QE; tens of labels, not thousands |

### 4.2 Measured: the fine-tune pipeline runs on this GPU

1. **Existing property** (`dft_band_gap`, first 2000/200/200 rows of MP-20): `run_finetune_timing.sh`.

   | Micro-batch (accumulate) | Result |
   |---|---|
   | 128 (×4) or 64 (×8) | **fills 24 GB, then spills into shared host memory (WDDM)**: 100 % utilization at ~85 W, no epoch finished within 13 min. Killed |
   | 32 (×16) | 20.4 GB; 63 it/epoch in 12 s → ≈ 170 samples/s; 10 epochs in 104 s wall |
   | 16 (×32) | 8.6 GB; 125 it/epoch in 16 s → ≈ 125 samples/s; 10 epochs in 115 s |

   After 10 epochs, `loss_val` = 0.40. Generation from the fine-tuned checkpoint with
   `--model_path=<run dir>` works (16 structures in 144 s).
2. **New property id, no edits to `~/mattergen`** (`run_newprop_test.sh`):
   - `src/ml/mg_newprop.py` appends extra ids to `PROPERTY_SOURCE_IDS` in-process (every module
     shares that list object).
   - The embedding YAML lives in `src/ml/mg_conf/` and is found through
     `hydra.searchpath=[file://…/src/ml/mg_conf]`.
   - The full chain `csv-to-dataset → finetune → generate` completed for a **synthetic** label
     (`rtt_test_scalar` = number density). This proves only the plumbing.
   - With 2000 structures and 2 epochs, the conditioned samples do **not** follow the condition;
     no adherence was expected or tested.

   Template command:
   ```bash
   MG_EXTRA_PROPS=rtt_lambda ~/mattergen/.venv/bin/python src/ml/mg_newprop.py finetune \
     adapter.pretrained_name=mattergen_base data_module=mp_20 data_module.root_dir=<cache>/<name> \
     'hydra.searchpath=[file://<repo>/src/ml/mg_conf]' \
     '+lightning_module/diffusion_module/model/property_embeddings@adapter.adapter.property_embeddings_adapt.rtt_lambda=rtt_lambda' \
     '~trainer.logger' 'data_module.properties=["rtt_lambda"]' \
     trainer.accumulate_grad_batches=16      # => micro-batch 32 (20 GB); use 32 for 16 (8.6 GB)
   ```
   This needs `src/ml/mg_conf/.../property_embeddings/rtt_lambda.yaml`, a copy of
   `rtt_test_scalar.yaml` with the name changed. For a real dataset, also set
   `data_module.average_density` to the dataset mean.

### 4.3 GPU-hour estimate

At ~150 samples/s, the cost is N_labels × epochs / 150 s.

- **Cerqueira-size labels**: 7,000 × 200 epochs (the default) = 1.4 M samples ≈ **2.6 GPU-h**.
- **ML-Tc labels on all ~200k Alexandria metals**: 30 epochs ≈ 11 GPU-h.
- **Validation sampling**: ≈ 10 s/structure, so 1,000 conditioned samples ≈ 2.8 GPU-h.

The GPU is **not** the constraint. The constraints are:

1. **Label range.** The target (Tc ≈ 300 K at ~0 GPa, or A ≳ 1300 K per the stretch-plan memo)
   lies far beyond the training labels. CFG then extrapolates, and the model card warns about
   exactly this. Expect samples that move *toward* high λ/ω_log motifs (H-rich, light-element,
   covalent-metallic), not samples that hit the target.
2. **Label quality.** The available DFT EPC data are "modest convergence" isotropic λ with
   McMillan-type Tc. ML-predicted Tc labels inherit the predictor's bias.
3. **Pressure.** All labels and the base model are 0 GPa. That suits the 0–1 GPa goal, but
   high-pressure hydride data cannot be used without a pressure property, and no pressure
   property exists.
4. **Local DFT labels are too few** for fine-tuning: tens, not thousands. They are better used
   as a validation set for whatever label proxy is chosen.

### 4.4 Recommended use in the targeted pipeline

1. **Generation now**: `chemical_system` (γ ≈ 2) over H-rich ternaries chosen by the theory
   documents, plus `chemical_system_energy_above_hull` at E_hull ≤ 0.05. Cost is ~10 s/sample:
   about 1k samples per system in ~3 h.
2. **Screen**: batched MatterSim relax, then the MatterSim-only hull (keep ≤ 0.1 eV/atom), then
   StructureMatcher dedup and novelty checks. Run MatterSim phonons only for *ranking*, never for
   rejection. Mg2IrH6 is the positive control.
3. **Fine-tune later** on Cerqueira et al. λ/ω_log, or on a theory-derived proxy computed for
   Alexandria structures, if the downloadable data have structures and permissive licensing
   (**UNVERIFIED**: not downloaded or inspected in this audit). Condition on the high tail of
   the label distribution, not on an out-of-range 300 K.
4. All verdicts stay with DFPT plus the certified Eliashberg bounds.

## 5. Files

`src/ml/`:

- `mg_generate.sh` — generation with timing and GPU memory logging.
- `run_smoke_generation.sh` — the three runs in §3.2.
- `mattersim_screen.py` — `relax` and `phonons` subcommands.
- `ehull_ref.py` — DFT reference hull via LMDB, plus `--export-ref`.
- `ehull_self.py` — MatterSim-only hull.
- `run_smoke_screen.sh` — relax and hull for the three sets.
- `check_mg2irh6.py` — the §3.4 cross-check.
- `run_batch_scaling.sh` — batch-128 generation test.
- `run_finetune_timing.sh` — `ACC=16 EP=10 bash run_finetune_timing.sh`.
- `mg_newprop.py`, `mg_conf/…/rtt_test_scalar.yaml`, `run_newprop_test.sh` — new-property
  fine-tune path.

## 6. Not verified / open

- Whether the generated I4mm Mg2IrH5 matches the experimental Mg2IrH5 structure.
- Whether MatterSim's Mg2IrH6 instabilities persist with larger supercells or other
  displacement amplitudes. Γ-point instabilities do not depend on supercell size, so they very
  likely persist. A DFPT Γ-point check on Mg2IrH6 would settle MatterSim's reliability for this
  class and is cheap; I recommend it.
- The licence and structure availability of the Cerqueira et al. data (Materials Cloud
  2023.163), and the jarvis-tools key of the Choudhary–Garrity EPC set.
- The MatterSim 1M vs 5M accuracy difference on hydrides; only 5M was used.
- The cause of the +0.10 eV/atom offset of MP-sourced Mg entries inside the Alex-MP reference.
- Downloaded and kept in scratch: HF cache checkpoints (3 × 0.5 GB), the reference set (0.87 GB
  gz + 3.85 GB LMDB) and `mp_20.zip` (12 MB, sha256 matches the LFS pointer). Fine-tune test
  checkpoints were deleted.
