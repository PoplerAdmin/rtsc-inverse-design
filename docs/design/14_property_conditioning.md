# 14 — Conditioning MatterGen on A = ∫α²F dω

Agent D1-labels, 2026-10-29. Every number below is **measured** on this machine unless marked
*cited* or **UNVERIFIED**.

- Scripts: `src/ml/labels/`. Embedding config: `src/ml/mg_conf/.../property_embeddings/rtt_log10A.yaml`.
- Data sources, licences and checksums: `data/LABELS_SOURCES.md`. Label tables and summaries:
  `data/labels/`.
- Raw data, checkpoints and samples stay in scratch, `~/rtsc_targeted_work/D1-labels/`.
- ML output here is heuristic. It never issues a verdict (rule S2).

## 0. Summary

1. **Labels.** 9,307 open DFPT electron–phonon calculations with structures and α²F were
   collected and relabelled. Each one gets A, λ, ω_log and λ⟨ω²⟩, recomputed from its own α²F.
   - Cerqueira–Sanna–Marques: 8,249 entries, CC BY 4.0.
   - JARVIS EPC: 1,058 entries, CC BY 4.0.
   - All entries are at 0 GPa. After the gates (≤ 20 atoms, dynamically stable) **8,675 training
     rows** remain, 634 of them hydrides.
2. **No open structure + α²F label comes anywhere near the target.**
   - Largest A in the whole set: 390 K (JARVIS MgH2, dynamically *unstable*). Largest among
     training rows: 346 K (LiPdH2).
   - Zero entries have A ≥ 500 K. Theorem N1 needs A ≥ 1299 K; the design target is A ≥ 2060 K.
   - The 99th percentile is 150 K. Maximum λ⟨ω²⟩ is 9.4×10⁵ K², against a target of
     ≥ 6×10⁶ K².
   - The design target therefore lies **6–7 standard deviations** above the training mean in
     log10 A (mean 1.64, sd 0.25).
3. **Fine-tune.** `mattergen_base` was fine-tuned on `rtt_log10A` = log10(A/K): 7,835 train /
   409 val / 431 test rows, split by formula, 200 epochs, 1.7 GPU-h. It works end to end.
4. **Adherence** (64 samples per setting, γ = 2, same seed). Inside the label range the
   conditioning moves the output in the expected direction, and the effect is large and
   significant. As the target rises from 30 → 100 → 300 K:
   - light-element atom fraction: 0.07 → 0.33 → 0.63;
   - H-containing samples: 5 % → 22 % → 33 %;
   - composition-kNN label estimate: 33 → 79 → 81 K;
   - calibrated MatterSim phonon proxy: rises;
   - all of these differ from the unconditional set at p < 10⁻³ (Mann–Whitney).

   **Beyond the label range (800 and 1500 K) nothing increases further.** Every proxy plateaus or
   falls back. 1500 K behaves like "the top of the training distribution" (A ≈ 100–300 K), not
   like A ≥ 1300 K. The high targets do produce more A₂MH₆-type complex hydrides (Mg2IrH6
   family) and transition-metal nitrides and borides, plus some metal-free molecular junk (N₂,
   N–O).
5. **Use it as a biased prior, not as a target-hitter.** The adapter shifts generation toward
   light-element, H- and B/C/N-rich metallic chemistries. It cannot certify, rank, or produce
   A ≥ 1300 K. The proxies used to judge adherence are weak (§4.1), so even the in-range
   adherence is only **weak evidence**.

## 1. Data

| Source | Entries | P | α²F | Notes |
|---|---|---|---|---|
| Cerqueira, Sanna, Marques, Materials Cloud 3kbt5-r3n56 (arXiv:2307.10728) | 8,249 | 0 GPa (no `press` in vc-relax) | yes, 10 smearings, QE `a2F.dosK` (Ry) | PBEsol, near-hull Alexandria structures. Biased toward predicted Tc > 5 K (*cited*: ML-screened before DFPT) |
| JARVIS-DFT EPC `supercon_3d` | 1,058 | "0 GPa" | yes, 50 bins (meV) | 432 flagged unstable; 11 with no positive α²F |
| JARVIS `alex_supercon` | 8,253 | — | **no** (λ, ω_log only) | It is the same set as Cerqueira: 8,248 ids match and λ is identical. Used only for a cross-check |

The 50 GB Cerqueira archive was streamed rather than stored (`run_fetch_cerqueira.sh`,
`fetch_cerqueira.py`, ~45 min). Only the label files were kept: 1.3 GB. All 17 stream md5 sums
equal the published ones. Not used, with reasons, in `data/LABELS_SOURCES.md`:

- Bercx et al. 2025: labels are only in a 15 GB AiiDA archive. It is the best next addition.
- 2D JARVIS: outside the bulk domain.
- High-pressure hydride sets: excluded by the pressure rule.

## 2. Labels (`compute_labels.py`)

All four labels are computed from α²F over ω > 0, using the definitions of
`docs/theory/01_eliashberg_targets.md` §0:

- A = ∫α²F dω
- λ = 2∫α²F/ω dω
- ω_log
- λ⟨ω²⟩ = 2∫ω α²F dω

Units are K.

- **Cerqueira smearing choice.** The authors' `McMillan.dat` "integral a2F" equals my A at
  smearing #6 (0.030 Ry) to a median relative error of 6×10⁻⁶ (95th percentile 0.2 %). Smearing #6
  is therefore the primary label, which reproduces the published data.
  - Smearings #4–#8 give a spread of 7 % (median) and 37 % (95th percentile). That spread is the
    label uncertainty, and it is stored as `A_K_smear_min/max`.
- **JARVIS** A is computed on a 50-bin grid. λ agrees with the published λ to a median ratio of
  0.995 (5–95 %: 0.85–1.09).
- **Computed vs derived.** Every training label is *computed* from α²F. The λ, ω_log-only set
  supports only a rigorous **lower bound**, A ≥ λω_log/2. This follows from Jensen,
  ⟨ω⟩ ≥ ω_log, with A = λ⟨ω⟩/2 [DER]. On the computed data the bound is tight to a median
  factor of 1.11 (5–95 %: 1.04–1.55). It is recorded (`A_K_lowerbound`) but not used as a label.
- **Dynamical stability gate.**
  - Cerqueira: QE phonon dispersion minimum ≥ −30 cm⁻¹ (tolerating Γ acoustic-sum-rule noise).
    This drops 191 entries.
  - JARVIS: `stability == stable`. This drops 429.

### Distribution (all entries with finite A)

| Set | n | n with H | A p5 / p50 / p95 / p99 / max (K) | n(A ≥ 300) | n(A ≥ 500, 800, 1000, 1300, 2060) | max λ⟨ω²⟩ (K²) |
|---|---|---|---|---|---|---|
| Cerqueira all | 8,248 | 671 | 16 / 45 / 103 / 147 / 346 | 7 (6 with H) | 0 | 9.4×10⁵ |
| JARVIS all | 1,047 | 19 | 20 / 55 / 156 / 243 / 390 | 3 | 0 | 6.6×10⁵ |
| **Training (gated)** | **8,675** | **634** | 16 / 45 / 104 / 150 / 346 | 7 (4 with H) | **0** | 9.4×10⁵ |

Other properties of the training set:

- **Pressure:** every row is 0 GPa, so nothing was removed by the ≤ 1 GPa rule.
- **Cell size:** every row has ≤ 20 atoms.
- **Hydrides:** 7.3 % of rows.
- **Top A among training rows:** LiPdH₂ 346 K, NaC₆ 331, KCdH₃ 328, NaPtH₃ 313, K₂AuH₄ 307,
  MgB₂ 302, LiMoN₂ 299, NaAuH₃ 297. Listed in `data/labels/summary_log10A.json`.
- **Distance to the target:** the gap between 1300 K and the best known 0 GPa label is a factor
  of **3.8** in A. In standardised log10 A the targets sit at z = +2.1 (300 K), +3.8 (800 K) and
  +6.2 (1500 K).

## 3. Fine-tune

Data: `build_dataset.py` writes the CSVs.

- Split grouped by reduced formula (no formula appears in two splits), 90/5/5, seed 0.
- Label `rtt_log10A`. Average density 0.0694 atoms/Å³.

Command (`run_finetune_A.sh`):

```bash
MG_EXTRA_PROPS=rtt_log10A ~/mattergen/.venv/bin/python src/ml/mg_newprop.py csv-to-dataset \
  --csv-folder $W/ds/csv --dataset-name log10A --cache-folder $W/cache
MG_EXTRA_PROPS=rtt_log10A OUTPUT_DIR=$W/ft_log10A ~/mattergen/.venv/bin/python src/ml/mg_newprop.py finetune \
  adapter.pretrained_name=mattergen_base data_module=mp_20 data_module.root_dir=$W/cache/log10A \
  data_module.average_density=0.06942 \
  'hydra.searchpath=[file://<repo>/src/ml/mg_conf]' \
  '+lightning_module/diffusion_module/model/property_embeddings@adapter.adapter.property_embeddings_adapt.rtt_log10A=rtt_log10A' \
  '~trainer.logger' 'data_module.properties=["rtt_log10A"]' \
  data_module.batch_size.train=16 trainer.accumulate_grad_batches=8 \
  trainer.max_epochs=200 trainer.check_val_every_n_epoch=5
```

- **Settings:** micro-batch 16 × accumulation 8 = effective batch 128; lr 5×10⁻⁶;
  `full_finetuning: true` (default).
- **Wall time:** 6,180 s (1.7 h); peak GPU 21.3 GB, including other agents' use.
- **Why micro-batch 16:** another agent's MatterSim job held 18 GB of the GPU at launch.
  `run_finetune_A.sh` therefore waits for the GPU to free up, and the 32-micro-batch layout was
  dropped.
- **Validation loss:** 0.190 at epoch 4, 0.158 at best epoch 179, 0.168 at epoch 199. The curve is
  noisy (±0.01).
- **Checkpoints:** in scratch, `~/rtsc_targeted_work/D1-labels/ft_log10A/lightning_logs/version_0/checkpoints/`.

| file | sha256 |
|---|---|
| `epoch=179-loss_val=0.16.ckpt` (used for sampling) | 57b7a233db1f4f58c1dd5e32357c17fcf14699eda6f695740d815163a0083069 |
| `last.ckpt` | 5957d557f626cb396a02f5ec033073ffb5feb7dc92083ddcb20133396f89391d |

Generation (`run_generate_A.sh`):

```bash
MG_EXTRA_PROPS=rtt_log10A python src/ml/mg_newprop.py generate OUT --model_path=$W/ft_log10A \
  --checkpoint_epoch=best --batch_size=64 --num_batches=1 --seed=1 --record_trajectories=False \
  "--properties_to_condition_on={'rtt_log10A':3.1761}" --diffusion_guidance_factor=2.0
```

Conditional sampling took 633–1,220 s per 64 structures; the slower runs shared the GPU.

Because every set uses the same seed, all sets have identical atom counts (median 8). The
comparisons are therefore paired in noise and size. One unintended effect: **the fine-tuned
unconditional model is itself different from the base model.** Its samples are lighter (mean
mass 72 vs 103 u) and closer to the training compositions (kNN distance 0.42 vs 0.88). Full
fine-tuning moved the whole distribution toward Alexandria metals. "uncond" below means the
fine-tuned model with no condition.

## 4. Adherence

### 4.1 Proxies, and how weak they are

The true A of a generated structure cannot be computed cheaply: it needs a DFPT run, CPU-days
each. Three proxies were used, each validated on the **held-out test split**, where the DFT A is
known:

| Proxy | Validation on 431 test structures (formula-disjoint) |
|---|---|
| Composition kNN: median log10 A of the 5 nearest training rows in element-fraction space, L1 (`gen_proxies.py`) | Spearman 0.69; MAE 0.129 dex vs 0.178 dex for predicting the mean |
| x_light (H, Li, Be, B, C, N, O, F atom fraction), x_H | Spearman 0.29 and 0.27 |
| MatterSim (5M) relax plus finite-displacement phonons: ω_max, phonon-DOS ω_log, and a 5-fold-CV linear fit log10 A ~ log ω_max + log ω_log + x_light + x_H (`calib_proxy.py`) | ω_max Spearman 0.25; CV-linear Spearman 0.29 and MAE 0.171 vs baseline 0.182 on the 370 MatterSim-stable structures |

**Conclusion:** within this data set, phonon-scale and composition proxies explain little of A.
A ≈ λ⟨ω⟩/2 is dominated by the electronic factor N(E_F)⟨g²⟩, which these proxies do not see.
The kNN proxy is the best of them, but it only interpolates labels of *known* metals. It is
blind to metallicity and novelty, and it cannot return a value above the training maximum.
MatterSim phonons are unreliable for hydrides (doc 10 §3.4). Adherence judged by these proxies
is therefore **weak evidence**.

### 4.2 Results

All sets: 64 samples, γ = 2, seed 1. p-values are two-sided Mann–Whitney tests against
"uncond". "Phonon" columns use MatterSim-stable structures only (|ω_min| < 15 cm⁻¹).

| Target A (K) | x_light mean | frac with H | x_H mean | mean mass (u) | kNN Â median (K) | metal-free samples | MS-stable | ω_max median (K) | phonon-proxy Â median (K) |
|---|---|---|---|---|---|---|---|---|---|
| uncond | 0.157 | 0.08 | 0.026 | 71.6 | 52 | 0 | 48/64 | 408 | 38 |
| 30 | 0.074 (p 0.24) | 0.05 | 0.008 | 76.8 | **33** (p < 10⁻⁴) | 0 | 44/63 | 442 (p 0.49) | 37 (p 0.21) |
| 100 | 0.335 (p 3×10⁻⁴) | 0.22 | 0.089 (p 0.02) | 59.7 | **79** (p < 10⁻⁴) | 2 | 37/64 | 735 (p 0.03) | 44 (p 0.005) |
| 300 | **0.633** (p < 10⁻⁴) | **0.33** | **0.162** (p 3×10⁻⁴) | **44.2** | **81** | 11 | 25/64 | **1138** (p 2×10⁻⁵) | **55** (p 5×10⁻⁶) |
| 800 | 0.565 (p < 10⁻⁴) | 0.20 | 0.088 | 51.4 | 76 | 9 | 26/64 | 968 (p 2×10⁻⁴) | 55 |
| 1500 | 0.509 (p < 10⁻⁴) | 0.28 | 0.111 (p 0.002) | 55.6 | 70 | 3 | 36/64 | 904 (p 7×10⁻⁴) | 53 (p 10⁻⁴) |
| base `mattergen_base` (doc 10) | 0.104 | 0.05 | 0.016 | 102.6 | 41 | 0 | — | — | — |

Reading:

- **Direction is right inside the label range (30 → 100 → 300 K).** Every proxy is monotone:
  light-element fraction, H content, lower mass, kNN Â, MatterSim ω_max and the phonon proxy. The
  A=30 set moves *down* in kNN Â (33 K vs 52 K uncond), so the condition works in both
  directions.
- **Magnitude is far short.** At target 300 K the kNN estimate is 81 K, and the phonon proxy is
  55 K. Both proxies are biased toward the training mean, so they understate. Still, nothing
  suggests the samples reach their target.
- **Out of range (800, 1500 K): saturation.** No proxy goes above its 300 K value. x_light,
  x_H, ω_max and kNN Â all slightly *decrease*. Classifier-free guidance at z = +3.8 and +6.2
  extrapolates the NoiseLevelEncoding embedding into a region it never saw. The result looks like
  a noisy version of "high-A training chemistry", not more extreme chemistry. This matches the
  MODEL_CARD warning (*cited*, doc 10 §2).
- **Chemistry at high targets** (top MatterSim-stable phonon-proxy entries):
  - A=800: H₈Pd₃Ru, Ca₂PdRhH₆, Ca₂CoPtH₆, RuH, PdRhH₂, LaB₂C, Ti₃N₄;
  - A=1500: BHC, TiB₉, NaSr₂RhH₆, K₂PtAuH₆, LiPd₅H₆, MnB₂;
  - A=300: BH₄, Na₂MnBeH₆, BN, ZrH₂, BaIrH₄, Sr₂IrAuH₆, CrB₂.

  The A₂MH₆ / A₂MM′H₆ complex hydrides are the Mg2IrH6 family (Dolui et al., *PRL* 132, 166001
  (2024), *cited* in doc 10). They are a sensible direction, but they appear at 300 K as well, so
  this is not a "1500 K-specific" signal. These samples are not relaxed against a hull and not
  deduplicated against the 0 GPa references here: they are **leads for the B-pipeline screen,
  not candidates**.
- **Junk rate rises with the target:**
  - metal-free samples: 0 → 11/64 at 300 K and 9/64 at 800 K (N₂, N₃O, N₇O, H₃N₅, CaN₇-type
    molecular or polynitrogen cells);
  - MatterSim-stable fraction: 75 % → 39–56 %;
  - relax non-convergence: 0 → 6–8 of 64.

  Every conditioned set must be filtered for metallicity and stability.

## 5. Limits and open problems

1. **The label ceiling is the binding constraint, not the GPU.** No open 0–1 GPa structure+α²F
   datum has A > 400 K. Conditioning at A ≥ 1300 K is extrapolation by a factor ≥ 3.8, and the
   sampler does not extrapolate (§4.2). Ways forward:
   - (a) condition at the top of the range (A ≈ 250–350 K) and rely on downstream screening;
   - (b) add the Bercx et al. EPW set: experimentally known structures, better convergence;
   - (c) add in-house DFPT labels for the high-A tail found by the B-pipeline (tens, not
     thousands), mainly to test the adapter rather than to train it;
   - (d) use high-pressure hydride α²F (A ≈ 1500–2500 K *cited-range*, **UNVERIFIED**) with an
     explicit pressure property. That needs a pressure-conditioned base model, which does not
     exist (doc 10 §1).
2. **Label quality.**
   - "Modest convergence" DFPT (*cited*, Cerqueira); smearing sensitivity 7 % median and 37 % at
     the 95th percentile.
   - PBEsol (Cerqueira) and JARVIS functionals are mixed.
   - The JARVIS 50-bin α²F is coarse.
   - Cerqueira's DFPT set is pre-selected by an ML Tc screen, so it over-represents moderate
     coupling.
   - Dynamical-stability gates differ between sources.
3. **Proxies are weak** (§4.1). The decisive test is DFPT on a sample of conditioned outputs,
   about 5 per target, with the `rtt.core.dfpt` workflow. It is not done. The cost is CPU-days.
4. **Full fine-tuning shifts the unconditional distribution** (§3). An adapter-only run
   (`adapter.full_finetuning=false`) would keep the base prior and is the obvious ablation. Not
   run.
5. **Single seed, single γ = 2, 64 samples per set.** γ and seed sweeps were not run. All
   p-values are uncorrected for multiple comparisons; several in-range effects pass p < 10⁻³ and
   survive any reasonable correction.
6. λ⟨ω²⟩ (the second theory quantity) was computed for every row (`lam_w2_K2`) but **not used**
   for conditioning. It is strongly collinear with A in these data, and its maximum (9.4×10⁵ K²)
   is a factor ≥ 6 below the target, so the same ceiling applies.
7. One A30 structure (index 23) crashed spglib's Niggli reduction in phonopy and was skipped.
   Relaxation did not converge for 2–8 structures per set; those were kept, using their last
   frame.

## 6. Files

- `src/ml/labels/`
  - `fetch_cerqueira.py`, `run_fetch_cerqueira.sh` — stream and extract the archive.
  - `compute_labels.py` — labels.
  - `build_dataset.py` — filters, split and CSVs.
  - `run_finetune_A.sh` — fine-tune.
  - `run_generate_A.sh` — sampling.
  - `calib_proxy.py` — proxy calibration on the test split.
  - `gen_proxies.py` — composition and kNN proxies.
  - `gen_phonon_proxy.py` — MatterSim phonon proxies.
  - `run_adherence.sh` — runs the proxy stage.
- `src/ml/mg_conf/lightning_module/diffusion_module/model/property_embeddings/rtt_log10A.yaml`
- `data/LABELS_SOURCES.md`
- `data/labels/`
  - `labels_all.tsv.gz` — 9,307 rows, every entry.
  - `manifest_log10A.tsv.gz` — training rows with their split.
  - `summary_log10A.json`
  - `adherence_summary.json`
- Scratch (`~/rtsc_targeted_work/D1-labels/`):
  - `raw/` — extracted archive, 1.3 GB, plus the JARVIS files;
  - `labels/*.jsonl` — with structures;
  - `ds/` — CSV splits;
  - `cache/`;
  - `ft_log10A/` — checkpoints and metrics.csv;
  - `gen/<set>/` — samples, relax, phonons;
  - `calib/` — test-split MatterSim proxies;
  - `logs/`.

## Review response (R14, `docs/reviews/R14_conditioning_review.md`)

Accepted corrections; these override the text above where they disagree:

- **Label ceiling depends on smearing.** Labels above use the 0.030 Ry smearing. Recomputed at
  0.005 Ry, the training maximum is 859 K (19 rows ≥ 300 K, 6 ≥ 500 K); no row reaches 1000 K
  under any smearing. The gap to 1299 K is ~1.5×, not 3.8× (small-smearing values may be
  under-converged). Labels should carry the smearing envelope, as the pipeline's own rule does.
- **z-scores corrected:** mean 1.639, sd 0.249 in log10 A give z = +3.38 (300 K), +5.10 (800 K),
  +6.2 (1500 K). MatterGen's sinusoidal encoder repeats every ~2π in z, which explains the
  plateau: no guidance setting can make it extrapolate to 1300 K.
- **"Monotone up to 300 K" is withdrawn.** kNN Â means: 82.7 K (A100) vs 76.5 K (A300),
  p = 0.40. Saturation starts near 100 K; only composition variables keep moving.
- **Doc 12 §3.3 acceptance test not run** (500 samples vs lane G1 at equal GPU cost, judged on
  coarse-DFPT A). Until it passes, the A-adapter is a minor experimental lane, not the backbone.
- **Full fine-tuning memorised training formulas** (18.8 % of unconditioned samples vs 1.6 % for
  base). Next retrain is adapter-only.
- Statistics: ~57 tests (Bonferroni ~9×10⁻⁴), proxies are ~2 independent signals, single seed.

**Pipeline decision.** Generation backbone = shipped chemical-system / energy-above-hull
conditioning plus theory filters (MatterSim hull, DOS(E_F), ω) and coarse DFPT. In-house DFPT
labels feed system reweighting and a ranking surrogate; revisit MatterGen retraining once there
are ~100 in-house labels with A ≥ 300 K. Conditioning on λ⟨ω²⟩ or ω_max is not pursued (no gain;
ω_max stays a filter).
