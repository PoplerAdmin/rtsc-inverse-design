# 30. Campaign 1: first real run of `rtt` (targets → MatterGen → named compounds)

Agent E1-campaign, 2026-10-29/30. Run id `campaign1`. Labels: **(meas.)** measured in this run;
**(est.)** estimated; **UNVERIFIED** not checked.

Where things are:

- Run directory: `runs/campaign1/` (privacy-sanitized JSONL export + `manifest.sha256` in `export/`,
  `evidence/<cid>.md` from `rtt evidence`, tables
  used below in `runs/campaign1/tables/`).
- Live ledger and compute scratch (outside the repo): `~/rtsc_targeted_work/E1-campaign/`.
  - `ledger.sqlite`; `work/<candidate>/S4|S5|S6/` (QE); `S1/` (MatterGen output); `S3/s3_artifacts/`.
  - Per-run `PRUNED_MANIFEST.tsv` (sha256 of every pruned file); logs in `logs/`.
- Drivers (new): `scripts/campaign_front.py` (S1 lanes → S2 → S3, day-wide rank cut),
  `scripts/campaign_controls_front.py`, `scripts/campaign_promote.py` (generated-lane quota),
  `scripts/campaign_warm_refs.py`, and the table scripts `scripts/campaign_*_table.py`,
  `campaign_gen_stats.py`, `campaign_relaxed_sg.py`, `campaign_funnel.py`.
- S4–S6 ran through the engine CLI: `rtt --run runs/campaign1 --ledger … --workdir … --ranks N step S4|S5|S6 <cids>`.
  The options `--ledger/--workdir/--ranks` are new (§7).

**Summary.**

- The pipeline ran end to end on real data: 1 040 candidates (528 MatterGen, 506 database, 6 controls).
  - S3: 952 relaxed with MatterSim.
  - S4: 29 candidates on QE (electronic) plus element references.
  - S5: 9 coarse DFPT runs (q 3³ + q 4³).
  - S6: one full-DFPT run launched (NaIrH3), still running when this was written (§5.6).
- **No candidate comes anywhere near the 1 299 K floor.**
  - Best coarse spectral weight of a non-control: A_coarse ≈ 108 K for Pm-3m NaIrH3, i.e.
    8 % of the floor.
  - The Mg₂IrH₆ control reaches 341 K (26 %).
- The evidence level of every named compound is "coarse DFPT, harmonic, isotropic; no full-DFPT
  verdict yet". Most of the 22 evidence items read NOT CHECKED (§6).
- Eleven pipeline bugs were found. Ten were fixed with tests, including four that made results
  wrong or runs hang (§7):
  - an infinite loop in the phonon supercell routine;
  - insulators passing H1 through a Gaussian tail;
  - generated cells reaching QE as P1;
  - a scratch-pruner race that corrupted DFPT results.

---

## 1. Setup

- **Targets.** `targets.toml`: T\* = 300 K, floor a0·T\* = 1 299.44 K (doc 01 §0.4), plus the
  doc 02 §5 hard and soft profile. Controls run at their surrogate T\*_s (doc 12 §3.3).
- **Compute.**
  - QE: ≤ 6 MPI ranks (nice 19), machine total ≤ 12 via the workflow's rank wait.
  - GPU shared with agent D1-labels: start only when < 12 000 MiB in use; MatterGen batch 16.
  - All children: `TORCH_DISABLE_NATIVE_JIT=1`, 2 threads (QE: 1 OpenMP thread per rank).
- **Lane plan.**
  - `plan_day(gpu_h = 2.0, batch = 16)` gave 44 jobs × 16 structures (704). Shares follow
    doc 12 §3.3: G1 0.55, G2 0.15, G3 0.10, G4 0.15, G5 0.05.
  - **G2 and G3 were not run, nor the `ml_bulk_modulus` quarter of G4** (11 jobs, 176
    structures). Their checkpoints (`space_group`, `dft_band_gap`, `dft_mag_density`,
    `ml_bulk_modulus`) are not in the local Hugging Face cache; the repo copies are git-LFS
    pointers.
  - The driver sets `HF_HUB_OFFLINE=1` and skips such jobs with status
    `SKIPPED_CHECKPOINT_NOT_LOCAL` in `lane_plan.json`. Downloading ~4 × 0.5 GB needs the
    user's approval.
  - Executed: **G1 24 jobs, G4 5 jobs, G5 4 jobs = 528 structures.**
- **G1 chemical systems.**
  - `targets.toml [chemsys] list = "targets/chemsys_list.txt"` points to a file that does
    not exist. `load_chemsys_list` therefore falls back silently to a default list that
    also contains donor–donor–H systems (bug B11, not fixed: the file belongs to the theory
    docs).
  - The driver restricted G1 to doc 02 family 1, the 70 systems {Li, Na, Mg, Ca, Sr, Ba, Al,
    Sc, Y, La} × {Rh, Pd, Ir, Pt, Os, Au, Ag} × H. 24 systems were drawn uniformly.
  - E_hull conditioning was drawn U(0, 0.3) per job, γ = 2.
- **A-adapter.** The fine-tuned log10 A checkpoint (R14: experimental minor lane) was **not used**.
- **Lane D.** Alex-MP reference set: the 24 G1 systems + B-C-Li, B-C-Mg, B-Be, B-C-N, with
  E_hull(DFT) ≤ 0.1 eV/atom, giving 506 entries (meas., 106 s).
- **Controls.** Mg₂IrH₆, H3S (200 GPa), MgB₂, Nb, Pb, Al, run through S1–S3. S4 and S5 were run
  for Mg₂IrH₆, Al, Pb and Nb.

## 2. Timings (meas.)

| Step | Work | Wall | Resource |
|---|---|---|---|
| G1 MatterGen (`chemical_system_energy_above_hull`) | 24 jobs, 384 structures | 8 252 s = **2.29 GPU-h** (21.5 s/structure) | GPU shared with D1's generation for the first ~2.5 h |
| G4 MatterGen (`chemical_system`, non-hydride) | 5 jobs, 80 | 1 011 s = **0.28 GPU-h** (12.6 s/structure) | GPU, partly shared |
| G5 MatterGen (`mattergen_base`) | 4 jobs, 64 | 486 s = **0.14 GPU-h** (7.6 s/structure) | GPU alone |
| Lane D query | 26 systems → 506 | 106 s | CPU |
| S2 validity + dedupe + novelty | 992 | 216 s | CPU |
| S3 (relax, two hulls, features, phonons) | 952 | 4.1 GPU-h of tier_runs (D 1.81, G1 1.47, G4 0.41, G5 0.29, C 0.13) | GPU. In the first chunks, ~2.3 h went to building 40+ MatterSim reference hulls (3–8 min per chemical system) |
| S3 hang (bug B1) | – | 1 h 10 min lost | – |
| S4 electronic, 29 candidates | – | 20.8 core-h incl. killed/re-run attempts; 0.02–0.05 core-h for 1–3-atom cells, 0.1–0.6 for 4–9 atoms, 1–3 for 10–14 atoms | 2 ranks per job |
| S4 element references (21 elements, s4 profile) | – | 12.5 core-h | 2 ranks |
| S5 coarse DFPT | 9 runs | see §5 (48 core-h in finished tier_runs) | 2–4 ranks |

MatterGen at batch 16 on a shared GPU took 21.5 s per conditional structure. That is 2× the
10.2 s of doc 10 (batch 64, GPU alone). The unconditioned base model took 7.6 s.

## 3. Funnel (current outcome per stage, meas.; `scripts/campaign_funnel.py`)

| Stage | G1 (384) | G4 (80) | G5 (64) | D (506) | C (6) |
|---|---|---|---|---|---|
| S1 | 384 PASS | 80 PASS | 41 PASS, **23 SCREENED_OUT(no_pseudo)** | 487 PASS, 19 SCREENED_OUT(n_atoms > 20) | 6 PASS |
| S2 | 369 PASS, 15 DUPLICATE | 80 PASS | 41 PASS | 456 PASS, 31 DUPLICATE | 6 PASS |
| S3 | 3 SCREENED_OUT(H7-ml); 366 survivors → **15 PASS** (5 by rank, 10 by quota) | **34 SCREENED_OUT(H7-ml)**; 46 → 0 PASS | 41 → 0 PASS | 1 SCREENED_OUT; 455 → **15 PASS** by rank | 6 PASS |
| S4 | 13 run: 7 PASS, 4 EXCLUDED(H1), 2 ERROR | – | – | 12 run: 8 PASS, 4 EXCLUDED(H1) (3 not run) | 4 PASS |
| S5 | 3 run: BaPtH3 DEFERRED(imaginary modes), NaIrH3 ERROR(admission), LaRhH2 ERROR (run corrupted by bug B4; clean re-run in progress, §5.6) | – | – | 2 run: ScH2 DEFERRED(multiband), LaH2 **SCREENED_OUT(coarse-A)**; Ba₂IrH₆ not run (ranks) | 4 PASS |
| S6 | NaIrH3 running | – | – | – | – |

**Rank cut.** S3's rank cut was applied once, across the whole day: chunks ran with
`rank_cut = 0` and one `rerank_s3(rank_cut = 20)` came at the end. In the first chunk alone the
rank cut would have taken every slot (§7 B9).

- The global cut chose **15 lane-D and 5 G1 candidates**. `score_s3` rewards hull, H fraction,
  light mass and light elements, and has no metallicity term. Lane D's binary hydrides (LaH₃,
  YH₂/₃, AlH₃, SrH₂ …) therefore rank at the top, and 4 of the 15 are insulators that S4 then
  removed.
- So that MatterGen output reached DFT at all, I added a **recorded deviation**:
  `CAMPAIGN_GENERATED_QUOTA`. The top 5 generated candidates by score, plus the 5
  best-scored generated cells with ≤ 6 atoms (cheap enough for S5), were promoted from
  `DEFERRED(rank_cut)` with that tag (`scripts/campaign_promote.py`).
- LiH (G1 subsystem, a known insulator) was skipped.

**S3 outcome mix.** 878 DEFERRED(rank_cut) is the queue, not a removal. Only 38 removals came
from the 0.60 eV/atom ML window: 34 of them are G4 non-hydrides (median E_hull_ML 0.51 eV/atom),
3 G1, 1 D.

## 4. What MatterGen generated (meas.; `scripts/campaign_gen_stats.py`, `campaign_relaxed_sg.py`)

| | G1 (chem.sys + E_hull) | G4 (chem.sys, non-hydride) | G5 (unconditioned) | D (Alex-MP) |
|---|---|---|---|---|
| structures | 384 | 80 | 64 | 506 |
| in the requested system (ternary) | **378/384 exact**, 6 binary subsystems, 0 off-system | 56/80 ternary, 8 subsystems (24 binary) | – | – |
| contain H | 381 | 0 | **2** | 453 |
| H fraction: mean / median | 0.37 / 0.33 | 0 | 0.01 / 0 | 0.41 / 0.50 |
| H_f ≥ 0.5 / ≥ 0.6 | 136 / 62 | 0 | 1 / 1 | 261 / 149 |
| n_atoms median | 10 | 8 | 10 | 10 |
| elements per cell | 3 (378), 2 (6) | 3 (56), 2 (24) | 2–6 (24 four-element, 1 six-element) | 2–3 |
| outside the SG15 set (S1 no_pseudo) | 0 | 0 | **23/64** | 0 |
| space group, raw cell at 0.01 Å: P1 / doc 02 S6 list | 149 / 43 | 30 / 0 | 22 / 7 | 19 / 116 |
| space group, MatterSim-relaxed at 0.1 Å: P1 / S6 list / sg ≥ 75 | **54 / 72 / 121** of 369 | 20 / 5 / 21 of 80 | 8 / 4 / 13 of 41 | 13 / 113 / 211 of 456 |
| novel vs Alex-MP (S2) | **336/369** | 78/80 | 31/41 | 0/456 |
| E_hull,ML median; ≤ 0.1 eV/atom | 0.096; 185/353 | **0.507**; 5/80 | 0.049; 28/35 | 0.055; 367/447 |
| S3 score median / max | 0.58 / 0.81 | 0.38 / 0.69 | 0.46 / 0.66 | 0.67 / 0.83 |

What this says:

- **Conditioning works for composition.** G1 put 98 % of its cells in the requested X-TM-H
  ternary. The unconditioned base model produced H in 2 of 64 cells, and 23 of 64 had elements
  outside the SG15 set.
- **Conditioning does not deliver H-rich cells.**
  - G1 median H_f is 0.33; only 62 of 384 cells reach the A₂MH₆/AMH₃ range (≥ 0.6).
  - Lane D (known compounds) is more H-rich (median 0.50).
  - This is doc 02's warning: H-rich metals are sparse in the training data.
- **Symmetry.**
  - 39 % of raw G1 cells are P1 at 0.01 Å; after MatterSim relaxation, 15 % are P1 and 33 %
    are tetragonal or higher.
  - The raw cells are noisy at the 10⁻³ Å level: P1 at QE's tolerance even when spglib labels
    them Pm-3m at 0.01 Å. This mattered for cost (§7 B3).
- **Non-hydride lane (G4).** Its cells sit far above the hull (median 0.51 eV/atom); 34 of 80
  are outside the ML window. Without the `ml_bulk_modulus` checkpoint the G4 lane is weak.
- **Stability.**
  - G1 cells conditioned on E_hull ∈ U(0, 0.3) came out with median E_hull,ML 0.096 eV/atom;
    50 % are ≤ 0.1.
  - MatterSim phonons find modes below −15 cm⁻¹ in 218 of 365 G1 cells (60 %) and 222 of 454
    D cells (49 %). For H compounds this is not a removal (doc 10 §3.4); Mg₂IrH₆ shows the same
    false instability (−288 cm⁻¹).

## 5. DFT stages

### 5.1 S4 electronic (29 candidates + 4 controls; `scripts/campaign_s4_report.py`)

| cid | lane | formula | sg | N | S4 | gap (eV) | N(E_F) tet. (/eV/atom) | sheets | H_DOS | E_form (eV/atom) |
|---|---|---|---|---|---|---|---|---|---|---|
| c-9a8fa100a45b | C | Mg₂IrH₆ | 225 | 9 | PASS | 0 | 0.386 | 1 | 0.59 | −0.236 |
| c-05c8a97de6cc | C | Nb | 229 | 1 | PASS | 0 | 1.518 | 2 | – | – |
| c-a5b0d5c7da5e | C | Pb | 225 | 1 | PASS | 0 | 0.520 | 2 | – | – |
| c-52c78a594a6c | C | Al | 225 | 1 | PASS | 0 | 0.396 | 2 | – | – |
| c-0ae3f133a308 | G1 | NaIrH₃ | 221 | 5 | PASS | 0 | 0.132 | 2 | 0.41 | −0.146 |
| c-d985155f945a | G1 | BaPtH₃ | 221 | 5 | PASS | 0 | 0.148 | 3 | 0.55 | −0.501 |
| c-f934a43f2a5e | G1 | LaRhH₂ (novel) | 123 | 4 | PASS | 0 | 0.356 | 2 | 0.27 | −0.573 |
| c-336a1451e874 | G1 | BaRhH₆ (novel) | 8 | 8 | PASS | 0 | 0.141 | 1 | 0.57 | −0.244 |
| c-84f99ce466ed | G1 | La₂Rh₂H₇ (novel) | 1 | 11 | PASS | 0 | 0.447 | 3 | 0.39 | −0.536 |
| c-f7c577c86a42 | G1 | LaRh₂H₇ | 65 | 10 | PASS | 0 | 0.277 | 2 | 0.41 | −0.324 |
| c-766e0609329a | G1 | Sr₂Pt₂H₈ | 8 | 12 | PASS | 0 | 0.202 | 4 | 0.53 | −0.455 |
| c-7ce96596b2cf | G1 | SrPtH₆ | 166 | 8 | EXCL(H1) | 2.71 | 0 | 0 | – | −0.443 |
| c-c0432b9f7f83 | G1 | BaPtH₂ | 160 | 4 | EXCL(H1) | 1.68 | 0 | 0 | – | −0.611 |
| c-cad6183c9829 | G1 | Ba₂Pt₂H₈ | 2 | 12 | EXCL(H1) | 2.28 | 0 | 0 | – | −0.544 |
| c-e8a5830ad88c | G1 | NaAuH₂ | 160 | 4 | EXCL(H1) | 2.04 | 0 | 0 | – | −0.193 |
| c-c4f72433527d | G1 | BaIrH₆ | 8 | 8 | ERROR (vc-relax hit nstep) | | | | | |
| c-0c8e5c1acff4 | G1 | La₂Rh₂H₈ | 12 | 12 | ERROR (nscf; re-run stopped for S6 ranks) | | | | | |
| c-37743845859b | D | LaH₂ | 225 | 3 | PASS | 0 | 0.337 | 1 | 0.59 | −0.630 |
| c-6e72932d035d | D | ScH₂ | 225 | 3 | PASS | 0 | 0.242 | 2 | 0.46 | −0.672 |
| c-64571b693470 | D | LaH₃ | 139 | 4 | PASS | 0 | 0.018 | 2 | 0.78 | −0.573 |
| c-4b9f8fcc298f | D | YH₃ | 139 | 4 | PASS | 0 | 0.074 | 2 | 0.56 | −0.542 |
| c-fc90f2183664 | D | Y₂H₆ | 194 | 8 | PASS | 0 | 0.036 | 3 | 0.78 | −0.606 |
| c-539acaa4a6e7 | D | Ba₂IrH₆ | 225 | 9 | PASS | 0 | 0.148 | 1 | 0.60 | −0.576 |
| c-95fcf2dda3de | D | BaRh₂H₇ | 65 | 10 | PASS | 0 | 0.083 | 2 | 0.46 | −0.299 |
| c-1e2338c56800 | D | La₄H₁₀ | 109 | 14 | PASS | 0 | 0.190 | 2 | 0.57 | −0.610 |
| c-35bd5d97a39b | D | AlH₃ | 167 | 8 | EXCL(H1) | 2.36 | 0 | 0 | – | −0.001 |
| c-36e0830f70f8 | D | AlH₃ | 12 | 8 | EXCL(H1) | 3.69 | 0 | 0 | – | +0.030 |
| c-9d7a5acac8e4 | D | AlH₃ | 227 | 16 | EXCL(H1) | 3.19 | 0 | 0 | – | −0.010 |
| c-0e5a5998a1ed | D | BaPtH₆ | 166 | 8 | EXCL(H1) | 2.95 | 0 | 0 | – | −0.473 |

- No candidate was magnetic (H2 never fired; every moment 0).
- H7-heur never fired. The largest E_form is +0.03 eV/atom against the 0.56 window. E_form vs
  elements is only a lower bound on E_hull; the real hull (S7a) was not run.
- **Not run (compute):** BaH₇Rh₂-type c-671064ad9675, Sr₄H₈, La₄H₉, Ba₃Ir₂H₁₂, Ba₄Ir₂H₁₄ (11–20
  atoms). Their S4 was stopped to free ranks for S5/S6; the outcome is missing, not removed.
- The controls Al, Pb, Nb were run before the H7 fix. They are elements, for which there is no
  formation energy anyway.

### 5.2 S5 coarse DFPT (q 3³ + q 4³, doc 13 T2-SCREEN; `scripts/campaign_s5_table.py`)

A = A_mid over admitted σ on q 4³ (A_nom of the rule), compared with the 1 299.44 K floor at
T\* = 300 K. Controls use their surrogate floor.

| cid | formula | lane | S5 outcome | A_coarse,mid (K) | A_coarse,lo | **A/1299** | λ (mid σ, q4³) | ω_log (K) | min ω (cm⁻¹) | core-h q3 + q4 (ranks) |
|---|---|---|---|---|---|---|---|---|---|---|
| c-9a8fa100a45b | Mg₂IrH₆ | C (T\*_s 80 K, floor 346.5 K) | PASS | 341.0 | 303.6 | 0.26 | 0.81 | 586 | −13.5 | 11.0 + 18.7 (4) |
| c-0ae3f133a308 | NaIrH₃ | G1 | ERROR: q 3³ mode-positive admission failed (max negative per-mode α²F fraction 3·10⁻⁶–3·10⁻³ > tol 10⁻⁶) | **107.8** (q 4³ only) | 92.4 | **0.083** | 0.20 | 820 | 0 | 1.17 + 2.13 (2) |
| c-d985155f945a | BaPtH₃ | G1 | DEFERRED(imaginary_modes; also multiband, < 3 σ admitted on q4³) | 85.0 | 85.0 | 0.065 | 0.53 | 95 | **−1 292** | 2.26 + 4.21 (2) |
| c-6e72932d035d | ScH₂ | D | DEFERRED(multiband): rule would remove ((1+e_c)·A_hi = 132 K ≪ 1 299) but n_sheets = 2 | 74.1 | 72.7 | 0.057 | 0.22 | 498 | 0 | 0.50 + 1.17 (2) |
| c-37743845859b | LaH₂ | D | **SCREENED_OUT(coarse-A envelope)**: (1+e_c)·A_hi = 65.2 K < 1 299 K | 32.4 | 31.6 | 0.025 | 0.17 | 224 | 0 | 0.52 + 0.83 (2) |
| c-f934a43f2a5e | LaRhH₂ (novel) | G1 | ledger: ERROR (B4-corrupted run); clean re-run: q 3³ finished, q 4³ running. **q 3³ only, informal (not a ledger outcome):** | 96.1 (q 3³) | 87.9 | **0.074** | 1.23 (soft-mode dominated) | 41 | −54 (matdyn mesh) | q3 ≈ 2.7 h wall (2) |
| c-05c8a97de6cc | Nb | C (floor 20 K) | PASS | 83.5 | 81.4 | – | 1.10 | 101 | −25.6 | 0.25 + 0.59 (2) |
| c-a5b0d5c7da5e | Pb | C (15.6 K) | PASS | 36.3 | 35.8 | – | 1.13 | 62 | 0 | 0.16 + 0.80 (2) |
| c-52c78a594a6c | Al | C (2.6 K) | PASS | 66.3 | 65.8 | – | 0.40 | 318 | 0 | 0.16 + 0.33 (2) |

Candid reading:

- **No S5 candidate has A_coarse above 9 % of the floor.** Even with the screen's full
  allowance (1 + e_c)(1 + δ)(1 + ε_sys)/(1 − ε_ω), the largest A_hi,coarse is 478 K (BaPtH₃,
  inflated by δ_q = 0.59 from its instability).
- NaIrH₃'s 108 K comes from q 4³ alone, because q 3³ failed admission. ScH₂ and BaPtH₃ escape
  removal only through the multiband and imaginary-mode blockers, not because their A is close.
- **Mg₂IrH₆ comes out low.**
  - λ ≈ 0.8 and A ≈ 341 K at the coarse tier. The published values are λ ≈ 2.4–2.8 and
    Tc ≈ 160 K (Dolui et al. 2024; values recalled, **UNVERIFIED** here).
  - Either the s4/coarse relaxation lands in a different minimum (doc 10 found two Fm-3m
    minima with Ir–H differing by 0.09 Å), or coarse q/k under-samples its coupling.
  - Its S5 decision is a PASS against the 346.5 K surrogate floor only through the e_c
    margin: A_nom 351 K vs floor 346.5 K. This is the most important control result of the run
    and needs the S6/A5→6 check.
- **Elemental controls are plausible.** λ(Al) 0.40 and λ(Pb) 1.13 are consistent with the
  known 0.4 / 1.5 (UNVERIFIED). The coarse core verdict (decide mode, evidence only) is
  CERTIFIED_CANDIDATE for Nb and Pb at their surrogate T\*, and UNDECIDED for Al and Mg₂IrH₆.

### 5.3 Cost: measured S5 vs doc 12 §8.1

| Cell | class | S5 measured (q4³ + q3³) | doc 12 estimate |
|---|---|---|---|
| ScH₂, LaH₂ (3 atoms, Fm-3m) | C1 | 1.4–1.7 core-h | 0.53 |
| NaIrH₃, BaPtH₃ (5 atoms, Pm-3m) | C1 | 3.3–6.5 core-h | 0.53 |
| Mg₂IrH₆ (9 atoms, Fm-3m) | C2 | **29.7 core-h** (7.7 h wall on 4 ranks) | 10.5–52.5 |

C1 hydrides with 5d metals cost 6–12× the H3S-based C1 estimate (semicore Ir/Pt, 3 species).
C2 sits inside its band. The scheduler's `COST_CORE_H["S5"]["C1"] = 0.525` should become
≈ 5 core-h (not changed: it is B1's table and only orders the queue).

### 5.4 S6 estimate (doc 12 §8.1 scaling) and decision

Two estimates, from measured q 4³ coarse ph.x cost and the H3S anchor (doc 11 §3.1):

- **Anchored:** full q 6³ single mesh ≈ 45× coarse q 4³ core-h. Two meshes (4³ + 6³) ≈ 1.3× the
  finest. So **S6 ≈ 58 × S5(q 4³)**.
- **Formula:** N_PW × 1.93 (90/60 Ry) · N_k × 3.2 (k 0.25 → 0.17) · (N_irr(6³)/N_irr(4³) + 1)
  ≈ 3 → **≈ 18 × S5(q 4³)**. Lower, because it ignores the k_fit growth of the e-ph
  interpolation.

| Candidate | S5 q 4³ core-h | S6 anchored | S6 formula | wall on 8 cores (0.85 derate), anchored | decision |
|---|---|---|---|---|---|
| NaIrH₃ (best non-control A_coarse) | 2.13 | **124 core-h** | 38 | **18 h** (< 48 h) | **launched** (4 ranks, ≈ 36 h; §5.6) |
| BaPtH₃ | 4.21 | 244 | 76 | 36 h | not launched: −1 292 cm⁻¹ instability; a harmonic S6 cannot leave `undecided` |
| ScH₂ | 1.17 | 68 | 21 | 10 h | not launched (ranks), would be next |
| Mg₂IrH₆ (control) | 18.7 | **1 085 core-h** | 337 | **160 h** (formula: 50 h) | **not launched**: > 48 h on 8 cores. Untested end to end, as the brief allows |

### 5.5 Controls: behaviour and recall

| Control | S1 | S2 | S3 | S4 | S5 | Note |
|---|---|---|---|---|---|---|
| Mg₂IrH₆ | PASS | PASS (not novel) | PASS (E_hull,ML −0.023; ML min ω −288 cm⁻¹, ignored as designed) | PASS (1 sheet, H_DOS 0.59, E_form −0.24) | PASS at T\*_s = 80 K, but only by the e_c margin | coarse λ ≈ 0.8 ≪ published; see §5.2 |
| Nb | PASS | PASS | PASS | PASS (2 sheets) | PASS (NEEDS_PHONON_CHECK: −25.6 cm⁻¹ coarse) | coarse verdict CERTIFIED at T\*_s 4.6 K |
| Pb | PASS | PASS | PASS | PASS | PASS | coarse verdict CERTIFIED at T\*_s 3.6 K |
| Al | PASS | PASS | PASS | PASS | PASS | λ 0.40 |
| H3S (200 GPa) | PASS | PASS | PASS (no ML hull > 1 GPa) | not run | not run | calibration only |
| MgB₂ | PASS | PASS | PASS | not run | not run | blocked (ph.x crash W4) |

Control recall: S1 6/6, S2 6/6, S3 6/6, S4 4/4, S5 4/4 (`report.summary.control_recall`). No
screen false negative at the surrogate floors. With Mg₂IrH₆ passing S5 by under 2 %, though, the
coarse tier is close to removing the one hydride control.

### 5.6 Still running when this was written

- **S6 NaIrH₃** (c-0ae3f133a308): `step S6 --ranks 4`.
  - Log: `~/rtsc_targeted_work/E1-campaign/logs/S6_naih3.log`.
  - Work: `work/c-0ae3f133a308/S6/q4x4x4`, then `q6x6x6`. QE finds 48 symmetry operations.
  - Expected ≈ 30–40 h on 4 ranks. This is the **first real QE run through S6 → core
    certify (bracket) → `decide_t`**, so far tested only on the Al fixture.
  - Expected outcome: `undecided`. Reason: H compound with no anharmonic check (G-VALID 2); a
    full A ≈ 100 K, ≪ 1 299 K, would otherwise give `rejectedProven` (N1-ISO), crosschecked.
- **S5 LaRhH₂** (c-f934a43f2a5e, the smallest novel generated metal), 2 ranks. Log:
  `logs/S5_clean_naih3_larh.log`.
- To rerun or finish:
  `rtt --run runs/campaign1 --ledger ~/rtsc_targeted_work/E1-campaign/ledger.sqlite --workdir ~/rtsc_targeted_work/E1-campaign/work evidence <cid>`.
- When both finish, regenerate the repo export with
  `python scripts/campaign_export.py ~/rtsc_targeted_work/E1-campaign/ledger.sqlite runs/campaign1/export`
  and the tables with the `scripts/campaign_*_table.py` scripts. The export in the repo is a
  snapshot from 2026-10-30 00:48, taken while these two runs were still going.

## 6. Named candidates and their evidence

Order: the §3.6 keys do not apply (no S6 result yet). Ranked instead by the furthest stage
reached, then A_coarse. All are far below target.

| # | Compound | cid | origin | furthest evidence | A_coarse / 1 299 K | Main open items |
|---|---|---|---|---|---|---|
| 1 | **NaIrH₃, Pm-3m** (5 atoms) | c-0ae3f133a308 | G1 (H-Ir-Na, E_hull cond. 0.107); matches an Alexandria hypothetical entry (not novel) | S4 PASS (metal, NM, E_form −0.15); S5 ERROR (q 3³ admission); S6 running | **0.083** | dynamical stability on q 6³; S7 hull, ZPE, anharmonicity; multiband (2 sheets) |
| 2 | **BaPtH₃, Pm-3m** (5) | c-d985155f945a | G1; Alexandria hypothetical | S5 DEFERRED: strongly unstable (−1 292 cm⁻¹) | 0.065 | the cubic cell is not a minimum; a distorted variant would need to be generated |
| 3 | **ScH₂, Fm-3m** (3) | c-6e72932d035d | D (known) | S5 DEFERRED(multiband) | 0.057 | known compound; coarse A is 6 % of the floor |
| 4 | **LaRhH₂, P4/mmm** (4) | c-f934a43f2a5e | G1, **novel** | S4 PASS (metal, 2 sheets, H_DOS 0.27, E_form −0.57); S5 q 3³ done (informal), q 4³ running | ≈ 0.074 (q 3³ only) | a soft mode (ω_log ≈ 41 K, −54 cm⁻¹ on the matdyn mesh) carries λ ≈ 1.2: near-unstable, not a high-A structure |
| 5 | **BaRhH₆, C2/m-like (sg 8)** (8) | c-336a1451e874 | G1, **novel** | S4 PASS (1 sheet, H_DOS 0.57, E_form −0.24) | not run | S5 cost C3 (est. 50–700 core-h) |
| 6 | **Ba₂IrH₆, Fm-3m** (9) | c-539acaa4a6e7 | D (Alexandria) | S4 PASS (1 sheet, H_DOS 0.60) | not run (C2 ≈ 30 core-h like Mg₂IrH₆) | the natural next S5: same family as the Mg₂IrH₆ control |
| – | LaH₂ | c-37743845859b | D | **removed at S5** (coarse-A envelope, 65 K vs 1 299 K, provisional e_c) | 0.025 | – |

Evidence statements: `rtt evidence <cid>` (all 22 items, §3.6 template). For every candidate
above, the following are **NOT CHECKED**:

- full-grid harmonic stability (S6);
- the Γ diagnostic;
- G-VALID 1–5 and the ENVELOPE inputs (no decide inputs below S6);
- the decide outcome;
- the full core verdict and μ\* sensitivity;
- the ZPE-corrected E_hull vs own-DFT references (S7a not run);
- the anharmonic check;
- MD/PIMD.

What is **checked**: structure, lane, novelty, H1, H2, the S4 n_sheets (multiband status, S4
only), the provisional margins used, and control recall. The complete statement for the lead
candidate is reproduced below; the others differ only in the values.

```
# Evidence statement: H3IrNa (c-0ae3f133a308)

Ranking keys (§3.6): decide outcome **None**, certified T_low = 0 K, A_mid = None K (A_lo = None K). T_up and A_hi are ceilings only.

Coarse DFPT (S5, q 3³ + 4³, not a ranking key): A_coarse,mid = 107.8 K, A_coarse,lo = 92.35 K; S5 outcome ERROR; tier_run #1083.

- **structure**: checked (H3IrNa, spacegroup 221, 5 atoms, P = 0 GPa; structure sha256 be514ee64c2ac94d)
- **origin lane**: checked (lane G1)
- **novelty**: checked (S2 PASS, novel=False; S2 outcome PASS, outcome #1112)
- **H1 metallic (S4)**: checked (N_EF_per_eV_atom = 0.1321, gap_eV = 0; S4 outcome PASS, outcome #3064, tier_run #1059)
- **H2 non-magnetic (S4)**: checked (moment_per_atom_muB = 0, abs_moment_muB = 0, moment_muB = 0, magnetic_energy_meV_atom = 0.005742; S4 outcome PASS, outcome #3064, tier_run #1059)
- **harmonic stability, full q grid of the final structure (S6)**: NOT CHECKED (no S6 phonons)
- **Gamma diagnostic**: NOT CHECKED
- **G-VALID 1 non-adiabatic**: NOT CHECKED (no decide inputs)
- **G-VALID 2 stable/anharmonic**: NOT CHECKED (no decide inputs)
- **G-VALID 3 Coulomb sign**: NOT CHECKED (no decide inputs)
- **G-VALID 4 convergence data**: NOT CHECKED (no decide inputs)
- **G-VALID 5 multiband data**: NOT CHECKED (no decide inputs)
- **ENVELOPE inputs**: NOT CHECKED (no decide inputs)
- **decide outcome and firing rule**: NOT CHECKED (no decision)
- **core verdict**: NOT CHECKED (no full/converged verdict)
- **ZPE-corrected E_hull (H7)**: NOT CHECKED (S7a not run)
- **anharmonic check**: NOT CHECKED
- **MD/PIMD persistence (H6)**: NOT CHECKED
- **mu* sensitivity**: NOT CHECKED (no full/converged verdicts)
- **multiband / k-resolved status**: checked (n_sheets=2 (S4 only; no band-/k-resolved check); S4 outcome PASS, outcome #3064, tier_run #1059)
- **provisional margins used upstream**: checked (S4: PASS (PROVISIONAL_THRESHOLD(H2),PROVISIONAL_THRESHOLD(H7)))
- **audit and control recall**: checked (control recall: S1 6/6, S2 6/6, S3 6/6, S4 4/4, S5 4/4)

8 of 22 items checked; 14 NOT CHECKED.
```

**Evidence level, plainly.** The strongest statement this campaign supports is negative.

- Within the harmonic, isotropic, coarse-DFPT approximation, none of the generated or database
  hydrides examined carries more than ~110 K of A = ∫α²F dω. The requirement is ≥ 1 299 K.
- That conclusion rests on the provisional e_c / M_c and on bench-grade calibration (doc 13).
- No candidate has a positive evidence statement toward 300 K.

## 7. Pipeline bugs found

Fixed bugs (each with a test):

| # | Bug | Effect in this run | Fix (file) | Test |
|---|---|---|---|---|
| B1 | `auto_supercell` loops forever when the longest direction already has rep 1 and the atom cap is exceeded | S3 hung 1 h 10 min on a 20-atom cell | shrink only reducible directions (`mlff/phonons.py`) | `tests/front/test_s3.py::test_auto_supercell_terminates_with_long_axis` |
| B2 | H1 took the **max** of the tetrahedron, σ and 2σ Gaussian N(E_F). The 2σ tail (0.27 eV) crosses gaps | 5 insulators (NaAuH₂ gap 2.0 eV, BaPtH₂ 1.7, AlH₃ 3.2, …) PASSed H1; re-decided after the fix | tetrahedron decides when present (`dft/electronic.py::decide_s4`) | `tests/back/test_electronic_belli.py::test_h1_gaussian_tail_does_not_override_tetrahedron` |
| B3 | Generated cells reached QE without symmetry: P1 at QE's tolerance, though labelled e.g. Pm-3m at 0.01 Å | perovskite q 3³ had 14 irreducible q instead of 4, with no k reduction; S5 projected at 40 h → 1–2 h after the fix | `dft_structure()` idealizes to the labelled group via spglib when the cell is only approximately symmetric (`engine/stages.py`, `structure/canonical.py::symmetrized_for_dft`) | `tests/front/test_dft_symmetrize.py` |
| B4 | Stale `dyn*` files of an earlier ph.x run made the per-q ScratchWatcher think q points were finished, and it pruned the **live** wavefunctions/dvscf | BaPtH₃ crash ("pwscf.wfc not found"); NaIrH₃ gave −6 000 cm⁻¹ garbage at X; LaRhH₂ DFPT diverged. All three were re-run clean and gave sane results | move stale dyn files aside before a non-recover ph.x run (`dft/workflow.py::_clear_stale_ph_outputs`) | `tests/back/test_chain_fallback.py::test_stale_dyn_files_moved_aside` |
| B5 | `S4Electronic` never passed `e_form_eV` | H7-heur was NOT CHECKED on every candidate | new `dft/element_refs.py` (elemental references at the s4 profile from Alex-MP ground states; H₂ in a box; cached, file-locked); wired into `S4Electronic.run/decide` | `tests/test_s4_eform.py` |
| B6 | pw.x failures with known robust fixes ended as ERROR (La hydrides: "too many bands are not converged" in nscf; YH₃: BFGS "dE0s is positive") | 5 S4 ERRORs; after the fix LaH₂, LaH₃, YH₃, LaRhH₂, LaRh₂H₇ and La₂Rh₂H₇ passed | one sticky retry with cg diagonalization / damped dynamics, recorded by a `<step>.fallback` marker (`dft/workflow.py::Chain`) | `tests/back/test_chain_fallback.py` |
| B7 | MatterGen `source_ref` stored the absolute home path in ledger rows and exports | privacy: home directory in exports | `~`-relative path (`source/mattergen.py::_tilde`); this run's export is sanitized | `tests/front/test_lanes_mattergen.py::test_source_ref_has_no_home_dir` |
| B8 | Evidence statement read novelty from tier_run metrics, but S2 stores it in the outcome value | novelty always NOT CHECKED | merge the S2 outcome value (`report/evidence.py`); header now also prints A_coarse (S5) | `tests/back/test_report.py::test_novelty_from_s2_outcome_value_and_s5_line` |
| B10 | Smaller fixes | S5 ERROR outcomes said only "error"; S4 tier_runs recorded no wall/core-h; the CLI wrote QE scratch and the ledger into the repo run dir | S5 now names the failing mesh and reason; S4 records `run.wall_s/core_h`; CLI gained `--ledger/--workdir/--ranks` (`RTT_LEDGER/RTT_WORKDIR`) | `tests/test_cli_paths.py` |

Not fixed (noted):

- **B9 — rank cut per call.** `run_s3(rank_cut = N)` called per chunk lets the first chunk take
  all N slots. The campaign worked around it (cut 0 per chunk, one final `rerank_s3`). The
  engine's `S3MLFF.run_batch` should do the same.
- **`score_s3` has no metallicity or composition-class term.** Known insulating binaries
  (AlH₃, SrH₂) dominate the rank cut; 4 of 15 lane-D promotions failed H1. Proposals: an
  S3-level ML gap proxy as a ranking feature, or per-lane rank quotas (what the campaign did by
  hand).
- **B11 — missing `targets/chemsys_list.txt`.** The fallback silently includes donor–donor–H
  systems. The flag `DEFAULT_CHEMSYS` is set but never surfaced.
- **Core mode-positive admission tolerance 10⁻⁶** rejects a smearing whose largest
  negative per-mode α²F fraction is 3·10⁻⁶ (NaIrH₃ q 3³). This is a soundness rule; I did not
  touch it. With the strict rule, a whole S5 becomes ERROR when one mesh fails. A per-mesh
  fallback to the T2-SCREEN one-mesh rule at q ≥ 4³ would keep the stage informative (doc 22
  §10 already allows it).
- **S4 relax `nstep` exhaustion** (BaIrH₆) is not in the fallback patterns.
- **Orphan tier_runs.** Killed jobs leave `state = 'started'` rows: 48 S3 from the hang, 5 S4,
  8 S5. At the end of the run they were closed as `failed` with metrics
  `{"note": "orphaned by a killed process (campaign 1)"}`. The engine should do this on restart.
- **Scheduler cost table.** C1 S5 is underestimated about 10× for 5d hydrides (§5.3).
- **MatterSim reference hulls.** 3–8 min per new chemical system; they dominated S3 wall time.
  G5's unconditioned multi-element cells each trigger a new hull. Proposal: skip the ML hull for
  S1-passing G5 cells with > 3 elements, or share hulls across agents.

Interventions by hand, stated for the record:

- S4/S5 workers were stopped and restarted several times to move ranks to S5/S6. Each restart
  left an orphan tier_run.
- Stale dyn files were moved aside by hand in three running S5 jobs before the B4 fix. Those
  three runs were later discarded (directories `S5_failed_*`) and re-run clean.

## 8. Not done / not verified

- G2, G3 and the G4 bulk-modulus lane (checkpoints not local).
- The A-adapter lane (not used).
- S7 (own-DFT hull, ZPE, anharmonicity, MD); S8.
- The A3→4 / A4→5 audits: no audit sampling was run in this campaign.
- Published Mg₂IrH₆ λ/Tc and the elemental λ values quoted for comparison were recalled
  (**UNVERIFIED**).
- Belli φ / Tc_φ are recorded but not used; Tc_φ for Mg₂IrH₆ came out −25 K, an extrapolation
  artefact of the formula.
- No commits were made: the orchestrator commits.

## 9. Tests after the fixes

`scripts/run_all_tests.sh` plus the new files (`tests/test_cli_paths.py`, `tests/test_s4_eform.py`):
**1 296 passed, 14 skipped, 0 failed** (2026-10-29, meas.). The earlier RC1 F11 failure is gone.
New tests:
- `tests/back/test_chain_fallback.py` (fallback text, sticky retry, no retry on unknown error, stale dyn);
- `tests/front/test_dft_symmetrize.py`;
- one test each in `tests/front/test_s3.py`, `tests/front/test_lanes_mattergen.py`,
  `tests/back/test_electronic_belli.py` and `tests/back/test_report.py`.

## 10. Recommended next steps

1. Let the S6 NaIrH₃ run finish. It is the first real QE run through S6 → certify → decide. Then
   check its ledger chain (verdict, recheck, decision).
2. S5 for Ba₂IrH₆ (Fm-3m, single sheet, H_DOS 0.60). It is the closest family member of the
   Mg₂IrH₆ control and needs ≈ 30 core-h.
3. Resolve the Mg₂IrH₆ coarse λ ≈ 0.8 vs literature discrepancy before trusting any coarse
   screen on hydrides. Options: an S6 on 8 cores (≈ 1 100 core-h est., > 48 h) or an A5→6 audit
   on a smaller relative.
4. Approve the download of the `space_group`, `dft_band_gap`, `dft_mag_density` and
   `ml_bulk_modulus` checkpoints (≈ 0.5 GB each) to run G2, G3 and G4-bulk.
5. Add a metallicity/class term or lane quotas to the S3 rank cut (§7).
