# 21. Build report: front end S1–S3 (`rtt.source`, `rtt.structure`, `rtt.mlff`)

Agent B2-front, 2026-10-29. Spec: `docs/design/12_architecture.md` rev 2 (§3.1 S1–S3, §3.3,
§3.5, §5.1, §6.1, §7). Codes against B1-core's `rtt.records`, `rtt.config`, `rtt.ledger.store`
(INTERFACES_READY). Labels: **(meas.)** measured on this machine in this build; **UNVERIFIED**
not checked.

## 1. What was implemented

### 1.1 Modules

| Module | Spec item | Contents |
|---|---|---|
| `source/lanes.py` | §3.3 lanes, §5.1 `LaneSpec`, `plan_day` | `CHECKPOINT_PROPERTIES` (the property sets of the released MatterGen checkpoints), `MATTERGEN_ELEMENTS`; `LaneSpec` refuses any conditioning that is not *exactly* the checkpoint's own property set, refuses non-shipped checkpoints, and refuses lane GA unless `allow_adapter=True` + `model_path`; `plan_day(machine, targets, weights)` splits the MatterGen GPU hours (default 15 h, §8.2) over G1 0.55 / G2 0.15 / G3 0.10 / G4 0.15 / G5 0.05, one conditioning value per 64-structure job; G1 E_hull condition drawn U(0, 0.3) per job; G1 chemsys sampled ∝ feedback weights floored at 0.2× mean (§3.3) |
| `source/mattergen.py` | §5.1 `generate(lane, n, seed)` | subprocess `nice -n 19 ~/mattergen/.venv/bin/mattergen-generate` (2 threads), shipped checkpoints only; output dir keyed by (lane spec, n, seed) and reused when finished; parses `generated_crystals.extxyz`; `Origin` gets lane, checkpoint, **checkpoint sha256** (HF blob name), conditioning, γ, batch id, seed, flags `HYPOTHETICAL_STRUCTURE` (+ `OFF_REQUESTED_SYSTEM` / `SUBSYSTEM_OF_REQUESTED`) |
| `source/database.py` | §5.1 `query_alexmp(chemsys, gap_max, ehull_max)` | lane D from the offline Alex-MP LMDB; DFT E_hull ≤ 0.1 on the hull of the same entries; ordered by H_f desc, N asc |
| `source/controls.py` | §3.3 control lane, §5.1 `CONTROLS`, `inject_controls(ledger)` | Mg₂IrH₆ (Alex-MP DFT cell), H3S Im-3m @ 200 GPa (calibration only), MgB₂ (`BLOCKED_PH_CRASH`), Nb, Pb, Al; T\*_s = 80, 150, 20, ½·Tc_exp K; injected through S1 like everything else; `anchor_from_lane_d()` for the future 6–8-atom anchor |
| `source/prefilter.py` | §5.1 `prefilter(rs, tiers) -> Outcome` | S1 removals: `no_pseudo`, `n_atoms_primitive > 20`, `outside_mattergen_domain`; records cost class C1/C2/C3, H_f, elements |
| `source/ingest.py` | (glue) | RawStructure → canonical primitive → sha256 → candidate id → `candidates` row + S1 outcome; `write_outcome` registers the control's surrogate targets text in `inputs` so the R7 trigger accepts it |
| `source/finetune.py` | §5.1 (stub) | `export_labels(ledger, prop)` → MatterGen CSV; `train_adapter` raises (owned by the label agent) |
| `source/external.py`, `source/mgvenv_tool.py` | (bridge) | paths (env/`machine.toml [paths]` overridable), GPU-in-use check, and a standalone helper run with the **MatterGen venv's** python (`lmdb`, `mattergen` live only there): `refquery` and `ehull-dftref` |
| `source/runner.py` | (glue) | `run_lane_job`, `run_database_lane` |
| `structure/canonical.py` | §5.1 | spglib standardized primitive (symprec 1e-3; Niggli fallback), deterministic site order, `structure_sha256` |
| `structure/validate.py` | §5.1 | S2 `PASS | INVALID`: finite/parseable, min distance ≥ 0.5 Å, volume/atom in [1, 100] Å³ |
| `structure/dedupe.py` | §5.1 `find_duplicate(s, ledger)` | grouping by (reduced formula, P), exact sha path, `StructureMatcher(ltol 0.2, stol 0.3, angle 5°)`; `is_novel` vs Alex-MP |
| `structure/stage.py` | S2 runner | validity → duplicate → novelty (metric in the PASS outcome's `value`); controls first, never duplicates themselves |
| `mlff/relax.py` | §5.1 `relax_batch(structs, model, fmax)` | `BatchRelaxer` FIRE + FrechetCellFilter, fmax 0.02, per-structure pressure (grouped), single-point re-evaluation; `TORCH_DISABLE_NATIVE_JIT=1` set before torch |
| `mlff/hull.py` | §5.1 `ehull_self`, `ehull_dftref` | MatterSim-only hull: Alex-MP references of the system relaxed with the same model at the same P, cached per (model, P, chemsys); DFT-reference hull via `mgvenv_tool ehull-dftref` (MP2020 via mattergen `VasprunLike`), 0 GPa only |
| `mlff/features.py` | §5.1 | H_f, vol/H, min H–H, min M–H, motif (H₂ / network / isolated / no H), H neighbours of M, light fraction, contains H/Li/Be |
| `mlff/phonons.py` | §5.1 `ml_phonons(struct) -> PhononSummary` | phonopy finite displacements with MatterSim; min ω, ω_max, ω₂^light (light-element-weighted second moment), harmonic ZPE, ω_log of the phonon DOS, `needs_phonon_check` |
| `mlff/alignn.py` | §5.1 | ALIGNN numbers, recorded only; weight 0 for H/light; inert unless `RTT_ALIGNN_MODELS` is set |
| `mlff/screen.py` | S3 runner + decision | order relax → hulls → features (+ALIGNN) → phonons for hull survivors; `decide_s3` (pure), `apply_rank_cut`, provisional `score_s3`; one `tier_runs` row per candidate written *before* compute (kind `ml`), one JSON artifact registered, one outcome |

### 1.2 Outcomes each stage can emit (and what can trigger them)

| Stage | Outcome | Trigger | Source id |
|---|---|---|---|
| S1 | `PASS` | – | – |
| S1 | `SCREENED_OUT(no_pseudo)` | element without `<El>_ONCV_PBE-x.y.upf` in `tiers [S1] pseudo_dir` | doc12#3.1:S1 |
| S1 | `SCREENED_OUT(n_atoms_primitive)` | primitive cell > `tiers [S1] n_atoms_max` (20) | doc10#1 |
| S1 | `SCREENED_OUT(outside_mattergen_domain)` | noble gas, Tc, Pm, Z > 83 | doc10#1 |
| S2 | `INVALID(finite | min_distance_A | volume_per_atom_A3)` | as named | doc21#S2-defaults (overridable by `tiers [S2]`) |
| S2 | `DUPLICATE(duplicate_of)` | StructureMatcher fit with an earlier S2-PASS candidate at the same P | doc12#3.1:S2 |
| S3 | `SCREENED_OUT(H7-ml:e_hull_ml_self_eV)` | MatterSim-only E_hull > `tiers [S3] ehull_self_max_eV` (0.42) | doc02#5a:H7 + doc10#3.3 |
| S3 | `DEFERRED(rank_cut(rank, cut))` | survivor ranked below `tiers [S3] rank_cut_per_day` (25) by `score_s3`; lane C exempt | doc12#8.3 |
| S3 | `DEFERRED(relax_not_converged)` | E_hull above the window but the relax hit `max_steps` | doc02#5a:H7 |
| S3 | `DEFERRED(ml_failed)` | relax error / non-finite energy | – |
| S3 | `PASS` (+ tags `NO_ML_HULL`, `RELAX_NOT_CONVERGED`, `SCREEN_DISABLED:S3:H7-ml`, `CONTROL_EXEMPT_RANK_CUT`) | otherwise | – |

**Never a removal at S3** (R2; doc 10 §3.4): MatterSim phonons (any min ω), ALIGNN, the
DFT-reference hull, features, novelty. `decide_s3` reads only `e_hull_ml_self_eV` and `error`;
`tests/front/test_s3.py::test_phonons_never_remove` checks this with min ω down to −10⁶ cm⁻¹, and
the end-to-end test runs every control (including Mg₂IrH₆) through S3 with fake MatterSim
phonons of −300 cm⁻¹ and requires PASS. Python `Outcome.__post_init__` and the ledger trigger
`outcomes_allowed` both enforce the §3.1 table.

## 2. Deviations from the spec, and why

1. **Lane D has no band-gap filter.** The Alex-MP reference LMDB stores energy, structure,
   POTCARs and ids only (inspected 2026-10-28). `gap_max` is accepted and recorded
   (`gap_filter_applied: false`), and every lane-D structure is flagged `GAP_UNKNOWN`.
   Metallicity is decided at S4 (H1), where §3.1 puts it anyway. Lane D also includes binary
   subsystems containing H (e.g. MgH₂ for H-Ir-Mg); S2 deduplicates across queries.
2. **Two venvs.** The main venv (Python 3.12, MatterSim) has no `lmdb`; the MatterGen venv
   (Python 3.10) has `lmdb`, `pymatgen`, `mattergen`. Everything that reads the LMDB or uses
   mattergen's `VasprunLike` runs as `mgvenv_tool.py` in the MatterGen venv (subprocess,
   `python -I`, nice 19). It imports nothing from `rtt`.
3. **Dedupe does not use the core `discovery/fingerprint.py`.** `core/` was empty
   when this was written. Grouping by reduced formula + pressure plus StructureMatcher (the same
   matcher doc 10 used) is sufficient at our scale; the fingerprint can be added as a pre-key
   later without changing the interface.
4. **S2 thresholds** (`d_min_A` 0.5, `v_min_A3` 1, `v_max_A3` 100) are not in `tiers.toml`
   (B1 owns it). They are read from `tiers [S2]` if present, else these defaults with source id
   `doc21#S2-defaults`. Proposal for `tiers.toml`:
   ```toml
   [S2]
   source = "doc21#S2-defaults"
   d_min_A = 0.5
   v_min_A3 = 1.0
   v_max_A3 = 100.0
   ```
   `INVALID` is a structural statement (garbage cell), not a physics screen; an H₂ distance
   (0.74 Å) is valid.
5. **`score_s3` is provisional** and lives in `mlff/screen.py` because `rtt.score.s3` did not
   exist. If `rtt.score.s3.score_s3(metrics, targets)` appears it is used automatically. Default:
   0.35·hull + 0.30·H-term + 0.25·ω₂-term + 0.10·light − 0.2·[H₂ molecular]; ALIGNN weight 0.
   It only orders survivors for the rank cut.
6. **Rank cut** is applied per `run_s3` call (one call = one day's batch), not across days.
   Controls are exempt (they must traverse every stage, §3.3).
7. **S3 failures** cannot be `ERROR` (not in the S3 row of §3.1), so a relax failure is
   `DEFERRED(ml_failed)` with the error text; its `tier_runs` row is `failed`.
8. **`ml_bulk_modulus` "high"** for G4 is 300 GPa (doc 02 says only "high"). Provisional;
   every fourth G4 job uses it.
9. **Cost class** for ≤ 5-atom cells that are neither cubic nor hexagonal is C2, and cells > 10
   atoms are C3 (§8.1 does not price either).
10. **Provenance of the source item** (MatterGen file#index, Alex-MP id) is kept as an
    `Origin.flags` entry `source_ref=...`, because `Origin` has no field for it.
11. **Novelty** is stored in the S2 PASS outcome's `value` (`{"novel": bool, "alexmp_match": id}`)
    rather than a separate table.
12. **H3S @ 200 GPa** goes through S3 relaxed at 200 GPa, with no ML hull (item 14); the DFT
    hull is skipped (0 GPa set). Its S3 numbers are calibration data, not a screen of anything.
13. **Controls are never `DUPLICATE`.** They are processed first in S2 and become the
    representatives for later lanes (a generated Fm-3m Mg₂IrH₆ is a duplicate *of the control*).
14. **No ML hull above 1 GPa** (`tiers [S3] ml_hull_max_pressure_GPa`, default 1.0). Only the
    H3S calibration entry is affected; its E_hull is `None` (cannot remove). Reason: in the first
    smoke run the 200 GPa H-S reference hull (72 cells) took 493 s and hit CUDA OOM three times
    (§4); H7 is an ambient-pressure window and the scope is 0–1 GPa.
15. **Unconverged relaxations never remove.** An E_hull above the window from a relax that hit
    `max_steps` gives `DEFERRED(relax_not_converged)`; below the window it is a PASS with tag
    `RELAX_NOT_CONVERGED`. Reason: an unconverged relax overestimates E_hull, and in the smoke
    run the two unconverged structures changed E_hull between identical runs (0.295 → 0.171 and
    0.002 → 0.170 eV/atom).
16. **GPU etiquette.** `relax_batch` and `ml_phonons` call `torch.cuda.empty_cache()` after each
    batch; `smoke_front.py` refuses to start (and to enter S3) when ≥ 4000 MiB are in use.


## 3. Tests (`tests/front/`, CPU unless marked)

`pytest tests/front` (CPU, 2 threads, `CUDA_VISIBLE_DEVICES=`): **37 passed, 3 skipped** in
11 s (meas., final code). Markers registered in `tests/front/conftest.py`:
`gpu` (enable with `RTT_GPU_TESTS=1`) and `external` (`RTT_EXTERNAL_TESTS=1`).

| File | What it checks |
|---|---|
| `test_lanes_mattergen.py` | LaneSpec refuses 8 kinds of unsupported conditioning and GA; every shipped checkpoint accepts exactly its own set; `plan_day` shares, G1 E_hull ∈ [0, 0.3], G2 ∈ S6 list, G4 non-hydride, determinism, weight floor; argv; `generate` with a fake executable (origin fields, flags, reuse without rerun, failure → `MatterGenRunError`) |
| `test_s1_s2.py` | control structures (Mg₂IrH₆ Fm-3m, 9-atom primitive, Ir–H 1.729 Å; H3S 200 GPa override; MgB₂ blocked flag); `inject_controls` idempotent, S1 PASS against surrogate targets; prefilter criteria; supercells judged by the primitive cell; hash invariance (supercell, permutation, origin shift); validity checks; StructureMatcher dedupe and novelty; `run_s2` marks a strained generated fcc Al as `DUPLICATE` of the Al control, an overlapping cell `INVALID`, skips S1-screened candidates; disallowed outcome rejected |
| `test_s3.py` | `decide_s3` hull rule (0.42 boundary kept), missing hull → PASS, error → DEFERRED, unconverged relax never removes, disabled screen; **phonons never remove**; rank cut + control exemption; score monotone; toy MatterSim-only hull; phonon summary; motifs; `compute_s3` and `run_s3` end to end with fake MatterSim hooks (tier_runs, artifacts, outcomes, H3S has no DFT hull, second call is a no-op); far-above-hull → SCREENED_OUT, relax failure → DEFERRED + failed tier_run |
| `test_gpu_external.py` | `gpu`: real MatterSim relax of Al (volume/atom) + phonons (no imaginary modes, ω_max 250–450 cm⁻¹): **passed** (meas., 10 s). `external`: `query_alexmp("H-Ir-Mg")` contains Mg₂IrH₆, all ≤ 0.1, all `GAP_UNKNOWN`; DFT-ref hull bridge on Mg₂IrH₆: **both passed** (meas.) |

Observation from the external test: feeding Mg₂IrH₆'s own *stored* DFT energy back through the
`ehull-dftref` bridge gives E_hull = −0.116 eV/atom, because `VasprunLike` adds the MP2020 H
anion correction (−0.119 eV/atom here) while the stored reference energies already include
their corrections (`correction` field 0). This is the intended behaviour for raw MatterSim
energies (MatterSim predicts uncorrected PBE) and is the same mixing doc 10 §3.3 calibrated
(MAE 64 meV/atom on H-Ir-Mg). It is another reason the DFT-reference hull is a recorded second
opinion only.

## 4. Live smoke run

Script: `tests/front/smoke_front.py` (not collected by pytest). Raw outputs (scratch, not in the
repo): `~/rtsc_targeted_work/B2-front/smoke1/` (cold) and `.../smoke2/` (warm caches, after the
fixes of §2 items 14–16). GPU was free (≈ 130 MiB in use) at both starts.

**Run.** Lane G1, checkpoint `chemical_system_energy_above_hull` (HF blob sha256 recorded),
conditioning `{chemical_system: H-Ir-Mg, energy_above_hull: 0.1}`, γ = 2, batch 16, seed
20261029; plus the six controls; S1 → S2 (with Alex-MP novelty) → S3 (relax, both hulls,
features, ALIGNN recorded as unavailable, phonons).

**Timing (meas.)**

| Step | Cold (smoke1) | Warm (smoke2) |
|---|---|---|
| MatterGen, 16 structures | 254 s (15.9 s/structure at batch 16; doc 10: 12.4 s at batch 16, 10.2 at 64) | reused |
| S1 ingest + controls | < 2 s | 1.9 s |
| S2 (validate, dedupe, novelty; 6 Alex-MP queries) | 18 s | 1.0 s |
| S3 total (22 candidates) | 925 s | **126 s** |
| – candidate relax at 0 GPa (21) | 84 s | 36 s |
| – reference hulls, built once per system and cached: H-Ir-Mg (113 cells), B-Mg (38), Nb, Pb, Al | 106 s, 134 s, 2 s, ~3 s, ~3 s | cached |
| – H-S reference hull at 200 GPa (72 cells) | **493 s, 3 CUDA OOM warnings** | not built (item 14) |
| – phonons (22) | ≈ 75 s | ≈ 75 s |
| peak GPU memory (whole device) | ≥ 13.2 GB, OOM | 7.7 GB |

The cold-run candidate relax (84 s, 4 s/structure) was slower than doc 10's 0.6–0.7 s because
a batch waits for its slowest member: 2 of 16 generated cells (and Nb) ran to `max_steps` = 1000.
In the warm run all relaxations converged and the same batch took 36 s (1.7 s/structure; cells
of 5–20 atoms, batch of 21). Recorded GPU-h for the 22
S3 tier_runs (warm): 0.034.

**Ledger rows (warm run):** 22 candidates (16 G1 + 6 C), 66 outcomes (S1/S2/S3 × 22, all
`PASS`), 22 `tier_runs` (kind `ml`, state `done`, lane G, versions mattersim 1.2.5, torch
2.14.1+cu130, phonopy 4.8.6, ase 3.29.0), 22 artifacts; `export_jsonl` + `manifest.sha256`
written. No `SCREENED_OUT`: every generated structure is inside the 0.42 eV/atom ML window
(max 0.350), and with 16 < 25 the rank cut defers nothing.

**Generated structures (G1, H-Ir-Mg).** 16/16 exactly in the ternary; 16 distinct after S2
(no duplicates); 15 novel vs Alex-MP, 1 known (`H2Ir2Mg2`, i.e. MgIrH, matches an Alex-MP entry).
Cost classes: 15 × C3 (5–20-atom low-symmetry cells, 11 in P1), 1 × C2. MatterSim-only E_hull
(warm run): median 0.093, min −0.0005, max 0.350 eV/atom; 6 of 16 ≤ 0.07 (doc 10's
"indistinguishable from 0" band). DFT-reference minus self hull: −0.01 … +0.08 eV/atom for 15
cells, −0.13 for one. 10 of 16 have imaginary MatterSim modes (down to −417 cm⁻¹); none is removed
for it. Top `score_s3`: H8Ir2Mg2 (Mg2IrH8 stoichiometry, P1, E_hull 0.033, 0.86),
H7Ir2Mg2 (0.80), H12Ir4Mg4 (0.79) — orderings from a placeholder score, not claims.

**Controls (S3 metrics; all PASS S1–S3)**

| Control | SG / N | E_hull self (eV/atom) | E_hull DFT-ref | min ω ML (cm⁻¹) | ω₂^light (K) | Alex-MP match |
|---|---|---|---|---|---|---|
| Mg₂IrH₆ | 225 / 9 | 0.000 | −0.026 | **−363** (reproduces doc 10 §3.4) | 1584 | yes |
| H3S @ 200 GPa | 229 / 4 | – (cold run: −0.005) | – (0 GPa set) | ≥ 0 | 2127 | – |
| MgB₂ | 191 / 3 | 0.000 | 0.035 | ≥ 0 | 833 | yes |
| Nb | 229 / 1 | −0.000 | −0.001 | ≥ 0 | – | yes |
| Pb | 225 / 1 | −0.000 | 0.001 | ≥ 0 | – | **no** (see below) |
| Al | 225 / 1 | 0.003 | 0.017 | ≥ 0 | – | yes |

Mg₂IrH₆ is the critical S3 control: MatterSim gives it a −363 cm⁻¹ mode, exactly the false
instability doc 10 found, and S3 passes it, as R2 requires. Pb fcc at a = 4.95 Å did not match
any Alex-MP Pb entry under StructureMatcher; not investigated (novelty is informational only).
Nb's MatterSim relax did not converge in 1000 FIRE steps in the cold run (it did in the warm
run); E_hull is ~0 either way.

## 5. Open items / not verified

- The 6–8-atom anchor control (§3.3) is not chosen (needs a person's choice from lane D).
- `score_s3` weights are a placeholder until `rtt.score.s3` exists (doc 11 T1a).
- The S2 thresholds and the G4 bulk-modulus value need entries in `tiers.toml` / doc 02.
- Lattice constants of the elemental and MgB₂ controls are textbook values (UNVERIFIED to
  better than ~1 %); every stage relaxes them.
- ALIGNN is not installed in the main venv; `alignn_numbers` reports `available: false`
  unless `RTT_ALIGNN_MODELS` and an `alignn` install are provided.
- No engine wiring: B1's `engine/stages.py` should call `ingest_raw`/`run_lane_job` (S1),
  `run_s2` (S2) and `run_s3` (S3). `run_s3` does not itself check the shared GPU; the scheduler
  (or the caller, as `smoke_front.py` does) must check `gpu_free_for_us()` first.

## 6. Response to RC1 (front-end findings F12, F15, F16)

Review: `docs/reviews/RC1_code_review.md`; tests `tests/review/test_rc1_findings.py`.
After the changes: `tests/front` + `tests/review` → 59 passed, 3 skipped. That includes
F12, F15 and F16, and also the S4 finding F7 once B3's fix landed. B1-core's edits to
`mlff/screen.py` (selection through `rtt.engine.retry.needs_run`) and to `structure/stage.py`
(N6: controls are not representatives) were merged, not overwritten.

### F15 — the S3 ML-hull window is now a measured number

`src/rtt/mlff/calibrate.py` (run on the GPU while it was free; data in
`docs/design/data/21_hull_calibration.json`):

- **Hull error.** Every Alex-MP reference entry of a system is relaxed with MatterSim-5M (the
  S3 self hull). For each converged H-containing entry, err = E_hull(MatterSim-only hull) −
  E_hull(DFT hull of the same reference set).
- **ZPE.** For entries ≤ 20 atoms whose hull decomposition products are also ≤ 20 atoms,
  Δ_ZPE = ZPE(entry) − Σ x_i ZPE(product_i), using harmonic MatterSim ZPE (positive modes
  only).

| System | H entries | MAE (eV/atom) | p95 \|err\| | ZPE pairs |
|---|---|---|---|---|
| H-Ir-Mg | 33 | 0.051 | 0.126 | 3 |
| H-Li-Mg | 40 | 0.027 | 0.093 | 0 |
| H-Ca-Pd | 87 | 0.022 | 0.065 | 8 |
| **pooled** | **160** | 0.030 | **0.092** | 11 (p95 \|Δ\| 0.018) |

Signed error: mean +0.004, range −0.13 … +0.14 eV/atom (meas.).

Window = 0.35 (H7) + √2·0.092 + 0.10 (Mg offset, doc 10 §3.3) + 0.018 (ZPE) = **0.598 →
0.60 eV/atom**. It replaces B1's interim 0.56 in `tiers.toml [S3] ehull_self_max_eV`. It is
marked `provisional = true` and its source string carries the derivation.

`screen.py`'s fallback default is 0.56, the RC1 floor; `tiers.toml` is authoritative.
`tests/front/test_s3.py` now reads the threshold from tiers and asserts ≥ 0.55.

Limits, stated:

- **Coverage.** Only 3 of the 7 planned systems were measured. The run was stopped when
  another agent's MatterGen fine-tune started on the GPU (shared-GPU rule). The remaining four
  were H-Mg-Pt, H-Na-Pd, H-Al-Li and H-La-Pd. Re-run the command when the GPU is free:
  `python -m rtt.mlff.calibrate OUT.json <systems>`.
- **ZPE term is weak.** It rests on 11 pairs, most of them on-hull entries whose
  decomposition is themselves (Δ = 0). It is an under-sampled lower estimate, not a measured
  tail.
- **What the error covers.** The reference entries are hull members on both sides, so the
  √2 factor is the reviewer's allowance for the candidate side. The error is also measured on
  known near-hull compounds, not on generated cells.
- **B1's score still uses 0.42.** `rtt.score.s3` (B1's file) still normalises its hull term
  by 0.42. That affects ordering only, not removal. B1 should read the threshold from tiers.

### F16 — dedupe

- **Tighter matcher.** `StructureMatcher(ltol 0.2, stol 0.1, angle 5°, scale=False)` instead of
  `stol 0.3, scale=True`.
- **Bonding guard.** `structure/dedupe.py::bonding_compatible` never merges two cells when any
  of these holds:
  - their H–H motif class differs (H₂-molecular / short / network / isolated / no H);
  - their shortest H–H distances differ by > 0.2 Å;
  - their shortest (non-H)–H distances differ by > 0.08 Å. This separates doc 10's two Fm-3m
    Mg₂IrH₆ minima, whose Ir–H differ by 0.09 Å.
- **Duplicate of a control: linked, not dropped.** A generated cell that matches a control is
  PASSed at S2 with `value.same_structure_as_control` and tag `SAME_AS_CONTROL`. It is then
  evaluated at the real T\* (B1's N6 change already stopped controls from being
  representatives; this adds the link).
- **Post-relaxation duplicates.** After relaxation, `run_s3` records `relaxed_match_of` in the
  S3 metrics when a relaxed cell matches an earlier relaxed cell. This is informational and
  never an outcome: §3.1 allows DUPLICATE only at S2.

### F12 (front part) — re-entry, rank queue, score

- **Retry.** `run_s3` selects candidates through B1's `rtt.engine.retry.needs_run`.
  `DEFERRED(ml_failed)` (including a whole-batch CUDA OOM) and
  `DEFERRED(relax_not_converged)` re-enter, within `[retry] max_attempts` (3).
  - `Outcome.retry` records the attempt number.
  - A retry doubles the relax `max_steps` per failed attempt (1000 → 2000 → 4000).
  - Test: `test_ml_failure_retried_with_more_steps`.
- **Global rank queue.** `rerank_s3(ledger, targets=, tiers=, rank_cut=)` is the API the
  engine calls daily.
  - `DEFERRED(rank_cut)` candidates are not recomputed. They are ranked from their stored
    `score_s3` together with each new batch.
  - Free slots = the daily cut minus the non-control PASSes still waiting for S4.
  - A promoted candidate gets a new PASS outcome with tag `PROMOTED_FROM_RANK_CUT`.
  - `run_s3` calls it at the end of every batch, and on its own when there is nothing to
    compute. Test: `test_rank_queue_promotes_on_rerank`.
- **No MatterSim phonons in the rank score for H/Li/Be.** `mlff/screen.py::rank_metrics`
  drops the phonon block for any compound containing H, Li or Be, or of unknown composition,
  before the score is computed. Its ω₂ term is then neutral (0.5) in both B1's
  `rtt.score.s3` and the fallback, and the metrics record `score_s3_phonon_term`. Test:
  `test_score_ignores_mattersim_phonons_for_HLiBe`. MatterSim phonons can therefore no longer
  move an H/Li/Be candidate across the rank cut.
- **Not in my files:** the summary line counting DEFERRED by reason and age, and the
  scheduler's back-off. Both belong to B1's engine/report.
