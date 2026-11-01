# 12. Architecture of the targeted pipeline (`rtt`)

Status: **DRAFT v2**, 2026-10-28. Author: agent T5-architecture. v1 was reviewed in
`docs/reviews/R12_architecture_review.md`; the point-by-point response is §9.

This document is meant to be buildable. It fixes the stage order and the outcomes each stage
may emit (§3), the soundness contract (§4), the module list with function signatures (§5), the
candidate record and ledger schemas (§6), the configuration files (§7), and the compute budget
(§8).

Inputs integrated in v2:
- `docs/theory/01_eliashberg_targets.md` (revised). Its §0.4 YAML rule list (G-VALID,
  ENVELOPE, N1-KAPPA, N1-ROWMAX, N1-DOS, N1-ISO, NM-G, RANK-AIM) is the decision layer at the
  DFPT tiers.
- `docs/theory/02_materials_profile.md` (revised): hard conditions H1–H7, soft priors S1–S11,
  the non-hydride channel 5(c), and the generation plan 5(d).
- `docs/theory/03_certificates.md`: Arb certificates for NM-G and N1μ.
- `docs/design/10_ml_capabilities.md`: measured MatterGen and MatterSim costs, the MatterSim
  phonon false instability on Mg₂IrH₆, and the fine-tune path.
- `docs/design/11_scoring_tiers.md`: public ML Tc models get weight 0 for hydrides; coarse
  DFPT at q 3³–4³, judged on A, is the middle tier; coarse imaginary modes never reject.
- `lean/TargetLogic`, described in `lean/README.md`: `decideT` with the outcomes
  `rejectedProven | rejectedAssuming | screenedOut | undecided | pass`.
- `docs/design/13_*` (coarse calibration, **in progress**). Every place that needs its result
  says **⟨doc13⟩**.

Labels: **(meas.)** = read from a run file; **(est.)** = my extrapolation, with the assumption
stated; **UNVERIFIED** = recalled or inferred and not checked.

§1 inventories the initial code base of this project (the certification engine and the
generate-then-verify pipeline around it, written before the targeted design). Its package root
is abbreviated **`base/`** below; the kept parts now live in `src/rtt/core/` (package `rtt.core`).

---

## 0. Summary

**Decisions on the initial subsystems** (details in §1):

| Initial subsystem | Decision |
|---|---|
| `certify/`, `spectrum.py`, `status.py`, `constants.py`, `validity/` | **Keep frozen in `rtt.core`**, sha256-pinned |
| `pipeline.py`: certify post-processing (`certify_spectra`, `apply_input_policy`, `blocking_tags`, `_sigma_exclusion_directions`, `annotate_model`, `ensure_cutoff_assumption`, `cutoff_description`, `trust_base`, `validate_certify_config`, `_certify_params`) | **Keep frozen in `rtt.core`** as `certify_policy.py`, line ranges recorded (R12 #1) |
| `pipeline.py`: orchestration (plan, run, collect, sharding, records) | **Replace** |
| `dfpt/` parsers, units, checks, attached, band_coupling, smearing admission | **Keep frozen in `rtt.core`** |
| `dfpt/` qe_inputs, grids, pseudo, launcher, structure | **Keep frozen in `rtt.core`**, extended only by wrapper code in `rtt.dft` |
| `dfpt/workflow.py` | **Adapt**: profiles, per-q pruning, timing log, re-relax before full DFPT |
| `mlff/`, `generate/mattergen.py` | **Adapt**, using doc 10's measured scripts in `src/ml/` |
| `generate/finetune_export.py` | **Adapt** for the A-adapter (§3.3) |
| `discovery/fingerprint.py`, `sources/{dedupe,fileformats,composition,local}.py` | **Keep frozen in `rtt.core`** |
| `discovery/ledger.py` | **Adapt**: new schema, same append-only triggers |
| `discovery/{loop,policy,goals,executors,report,generators}.py`, `descriptors/`, `report.py`, `cli.py`, `sources/{csp,mp}.py`, cluster templates | **Replace or drop** |
| Lean `TcR1`, `PipelineLogic`, `MultiBand` | Frozen in `lean/`, sha256-pinned |
| Lean `DiscoveryLogic` | Dropped |
| Lean `TargetLogic` (new) | `decideT` is the specification of `rtt.decide` |
| R-LAMW (v1 proposal) | **Deleted** (R12 #2) |

**Stage order and the outcomes each stage may emit** (full table in §3.1):

```
S0  targets            load targets.toml / tiers.toml                        -> (none)
S1  source             MatterGen lanes G1..G5, database lane D, control lane C,
                       composition/size prefilter                             -> PASS | SCREENED_OUT
S2  structure          validity, dedupe, novelty                              -> PASS | INVALID | DUPLICATE
S3  mlff (GPU)         MatterSim relax -> ML hull -> cheap features -> ML phonons (rank only)
                                                                              -> PASS | SCREENED_OUT | DEFERRED(rank)
S4  dft-electronic     pw.x vc-relax, DOS, spin, E_F, n_sheets, Belli phi, heuristic hull
                                                                              -> PASS | EXCLUDED_ASSUMPTION(H1,H2) | SCREENED_OUT | DEFERRED | ERROR
S5  coarse DFPT        q 3^3-4^3, sigma plateau, A_coarse; decideT(tier=coarse)
                                                                              -> PASS | SCREENED_OUT | DEFERRED | ERROR   (never a rejection)
S6  full DFPT          re-relax, two q meshes, >=3 admitted sigma; core certify_spectra
                       + decideT(tier=full)                                   -> verdict evidence + decide outcome
S7  validity closure   own-DFT hull + ZPE (H7), anharmonic check (G-VALID 2 / H6), k- or
                       band-resolved data (N1-MB, 0.75 margin)                -> re-run decideT
S8  converged          EPW protocol, finalists only                           -> decideT(tier=converged)
S9  evidence           statement per finalist; ranking by certified T_low and A_mid
```

**Budget** (§8): ≈ 5k MatterGen structures per GPU-day; about 25 DFT-electronic and 2–6 coarse
DFPT evaluations per day; full DFPT runs at about 7–25 per week for cells of ≤ 5 atoms with
high symmetry (C1), 0.3–1.4 per week for 6–10-atom cubic cells (C2), and fewer than 1 per month
for 6–10-atom cells of moderate symmetry (C3). A hydride cannot reach `rejectedProven` until the anharmonic check (S7)
exists, because of G-VALID item 2. Until then the decisive negative statements for H/Li/Be
compounds are "rejected in the harmonic model" (core evidence) with decide outcome
`undecided`.

---

## 1. Inventory of the initial code base

### 1.1 Python package (`base/`)

Totals: **13 984** LOC of package code and 7 791 LOC of tests, 21 775 LOC together (corrected per
R12 #12). There are 254 `def test_` functions and 442 parametrised cases (435 passed,
13 skipped, on this machine, 2026-10-27). The "Tests" column counts test functions in files that
import the module, so the counts overlap.

| Path (under `base/`) | Purpose | LOC | Tests | Quality notes | Decision + reason |
|---|---|---|---|---|---|
| `certify/classify.py` | the single decision function, a mirror of Lean `classify*`; `recheck` | 525 | ~116 + 756-line correspondence test + oracle (132/132) + forgery suites | the best-tested part | **Keep frozen in `rtt.core`.** Proof-carrying. |
| `certify/{arbcore,cells,rayleigh,multiband,floatcore,search}.py` | Arb primitives, R-CELLS, C-RAYLEIGH, R-MULTI, float search | 201+349+153+323+187+136 | (as above) | float search dominates bracket-mode time | **Keep frozen in `rtt.core`** |
| `spectrum.py`, `status.py`, `constants.py` | atomic α²F, Status/Verdict, a0 = 433148145/10⁸ | 195+58+47 | 140+ | shared types | **Keep frozen in `rtt.core`** |
| `pipeline.py` l. 83–133, 420–650 | `ConfigError`, `validate_certify_config`, `_certify_params`, **`certify_spectra`**, `trust_base`, **`_sigma_exclusion_directions`**, **`blocking_tags`**, **`apply_input_policy`**, `annotate_model`, `cutoff_description`, `ensure_cutoff_assumption` | ≈ 280 | `test_dfpt_smearing_admission.py`, `test_dfpt_soundness.py`, the input-policy parts of `tests/audit*` | the downgrade-only safety layer; its asserts check that status only goes down and the certificate is unchanged | **Keep frozen** as `core/certify_policy.py` (§5.2), including the mode-positive smearing admission fix |
| `pipeline.py` (rest) | plan/run/collect, records, sharding | ≈ 870 | ~69 | monolithic; `--plan-only` bug; logs lost on kill | **Replace** |
| `dfpt/{parsers,units,checks,attached,band_coupling}.py` | QE/EPW parsers, per-mode negative fractions, consistency checks | 549+83+119+122+543 | ~140 | formats checked against QE source | **Keep frozen in `rtt.core`** |
| `dfpt/workflow.py` | relax → scf_fit → scf → ph → q2r → matdyn; resume; `admit_smearings` (l. 511–612) | 757 | 379 + 132 + soundness | works on QE 7.6 (Al, Pb, Nb, H3S) | **Keep frozen** the file. `rtt.dft.workflow` wraps it with profiles; changes go in the wrapper only. |
| `dfpt/{qe_inputs,grids,pseudo,launcher,structure}.py` | input writers, q ⊂ k_fit rule, pseudo table with md5, MPI launcher | 156+97+157+184+204 | (dfpt) | clean | **Keep frozen in `rtt.core`** |
| `dfpt/epw_inputs.py` + `bench_tools/ephproto/` | EPW writer, `filqf`, phself parser, Richardson admission | 169 + ~1 400 | ephproto tests | the EPW protocol is DRAFT v3 | **Adapt later** (S8) |
| `mlff/{mattersim,hull,phonons,_atoms}.py` | relax, ML hull, phonopy phonons | 249+183+182+38 | 7 | 1M by default; crude hull references; single-structure relax (~17 s/structure, doc 10) | **Adapt**: batched 5M relax, MatterSim-only hull (doc 10 §3.3) |
| `generate/mattergen.py` | CLI adapter, `CHECKPOINT_PROPERTIES`, refuses unsupported conditioning | 341 | 18 | verified on real output | **Adapt**: lanes G1–G5 (§3.3) |
| `generate/finetune_export.py` | ledger → MatterGen CSV | 290 | 3 | format matches doc 10 §4.2 | **Adapt** for the A-label adapter |
| `discovery/fingerprint.py`, `sources/{dedupe,fileformats,composition,local}.py` | dedupe and I/O | 219+132+384+143+247 | ~84 | solid | **Keep frozen in `rtt.core`** |
| `discovery/ledger.py` | SQLite with append-only triggers and an origin guard | 459 | soundness suites | good discipline; the guard trusts a self-declared origin string (R12 #4) | **Adapt** (§6.2) |
| `discovery/{loop,policy,goals,executors,interfaces,report,generators}.py` | round loop, Beta bandit | 830+135+296+172+85+271+157 | ~65 | learns from ≤ 1 DFPT outcome a day | **Replace** |
| `stability/*.py` | X-HULL, X-METAL, magnetism, dynamic filters | 462 | 11 | criteria fine | **Adapt** into S4/S7 |
| `descriptors/*.py` | ad-hoc priority weights | 289 | 5 | superseded | **Replace** by doc 11 T1a/T1c in `rtt.score` |
| `validity/__init__.py` | out-of-model flags | 199 | 7 | flags never change status | **Keep frozen in `rtt.core`** |
| `sources/{csp,mp,prototype,synthetic,base}.py`, `report.py`, `cli.py` | | | | | **Drop** csp and mp (no API key; doc 10 uses the offline Alex-MP LMDB); **keep** synthetic for tests; **replace** the rest |

**Frozen-core rule.** The kept files live in `src/rtt/core/` with the relative paths they have
under `base/`; their tests are in `tests/core/`. Every file of both is byte-pinned: its sha256 is
recorded in `src/rtt/core/PINS.json` and checked by `tests/test_core_pins.py`. A change to a core
file requires re-deriving its soundness argument (Lean correspondence, soundness and forgery
suites) and then rewriting the pins with `python tests/test_core_pins.py --update`.

### 1.2 Lean layer

| Library | State in the new repo | What it proves | Role |
|---|---|---|---|
| `TcR1` (529 LOC) | core unchanged | k_B·Tc ≤ A/a0 for finite atomic α²F; the certificate H(a0, η) is a hypothesis | basis of N1-ISO and R-UNIV |
| `PipelineLogic` (1 414) | core unchanged | soundness of the core `classify` (R-UNIV, R-CELLS, C-RAYLEIGH, conjunction, empty list → UNDECIDED) | basis of core `certify_spectra` |
| `MultiBand` (1 175) | core unchanged | Tc ≤ κ/a0; Prop. N is a hypothesis | basis of N1-KAPPA |
| `DiscoveryLogic` (1 406) | not core | ledger state machine of the initial discovery loop | its invariants become SQL triggers (§6.2) |
| `TargetLogic` (new; 49 theorems; axioms only propext, Classical.choice, Quot.sound) | `lean/TargetLogic/` | N1-ROWMAX/N1-DOS algebra, L-MOM1–3 (documentation only), NM-G with `NMHyp` isolated, `decideT` with soundness (`decide_rejectedProven_sound`, `decide_rejectedAssuming_sound`), safety theorems (`not_rejected_of_not_proofGrade`, `decide_of_blocked`, `decide_downgrade`) | **specification of `rtt.decide`** |

Value at the end of the pipeline: the core engine plus `decideT` is the only route to a
removal that is proved rather than heuristic. It costs seconds. Lean proves the logic over
declared fields: `tier` and the envelope inputs. Binding those fields to real artifacts is the
ledger's job (§6.2, R12 #4; `lean/README.md` "Modelling choices" 1).

**Python ↔ Lean correspondence for `decideT`.** `decideT` over ℝ cannot be evaluated in Lean.
`decideCore` depends only on five booleans (GValid, proofGrade, FiresCertified, FiresNM,
BlockedRej), so:

1. **Add to TargetLogic:** `decideB : Bool → Bool → Bool → Bool → Bool → Outcome`, and a theorem
   `decideT_eq_decideB` that states `decideT` equals `decideB` applied to the `decide` of each
   predicate.
2. **Log the truth table:** a `#eval` over all 32 rows goes to `lean/logs/decideB_table.txt`.
3. **Test the Python side:** `tests/test_decide_lean_table.py` checks that `rtt.decide.core_b`
   reproduces the table, and property tests check each predicate (`gvalid`, `fires_*`,
   `blocked_rej`) against the formulas copied verbatim from `Decide.lean`. A drift test
   diff-checks those copied lines, as the core correspondence test in `tests/core/` does.

---

## 2. Known weaknesses of the initial pipeline (unchanged facts; consequences updated)

| # | Weakness | Evidence | v2 consequence |
|---|---|---|---|
| W1 | ph.x is 98–99 % of wall time; no middle tier | Al 3 328 s, Pb 8 268 s (earlier benchmark runs); H3S ph.x 1 h 58 min, end-to-end 122.7 min, user+sys 955 min (benchmark record, read 2026-10-28) | S5 coarse DFPT at q 3³–4³: H3S 271 s / 635 s on 2 ranks (doc 11 §3.1, meas.) |
| W2 | negative α²F at small σ → ERROR | `reports/pb_negative_a2f.md` | core `admit_smearings` and `certify_policy` (A-SIGMA-EXCLUDED blocking) |
| W3 | the smearing spread dominates | H3S: λ 2.40 → 1.51, A 1648 → 1084 K over σ = 0.01–0.10 Ry; 300 K verdict UNDECIDED [143.8, 360.7] K | G-VALID ladder (≥ 3 admitted σ, 2 q meshes); ENVELOPE δ terms; ranking by T_low and A_mid |
| W4 | QE 7.6: ultrasoft/PAW e-ph segfault; MgB2 ph.x SIGSEGV at start (2026-10-28, not diagnosed) | step2 report; benchmark runs | pseudopotential whitelist = SG15 files present; ERROR with a retry policy; MgB2 control blocked until diagnosed |
| W5 | disk, I/O, contention | 14–42 GB scratch per material; WSL crashes; ~15 % per-q penalty | lanes; per-q pruning; disk guard; budgets derated by 15 % (§8) |
| W6 | bracket-mode certify is slow (float search) | Al 1 075 s, Pb 548 s (contended); H3S 16 s | bracket mode only at S6 and later |
| W7 | ML phonon false instabilities | H3S with 1M (initial prescreen); **Mg₂IrH₆ with 5M at Γ and X** (doc 10 §3.4) | MatterSim phonons never remove H/Li/Be candidates |
| W8 | generation was not target-aware | only chemical system + E_hull | lanes G1–G5, database lane, A-adapter (§3.3) |
| W9 | μ* convention; wide brackets | Al T_up ≈ 9× experiment | rejections use μ*-free rules; ranking avoids T_up |
| W10 | multiband | MULTIBAND_RISK; R-MULTI pending audit | G-VALID item 5 and N1-MB rules (doc 01) |
| W11 | harmonic only | ANHARMONIC_H on H3S | G-VALID item 2 + S7 anharmonic check |
| W12 | engineering | monolith; empty logs; V2/V3 not run | small modules; `tier_runs` row written before launch; V2/V3 in the control lane |

---

## 3. Pipeline

### 3.1 Stages, inputs and the outcomes each may emit

Outcome vocabulary (one per candidate × stage; §6.1):
- `PASS`
- `SCREENED_OUT` (heuristic, reversible, with criterion, value, threshold and source id)
- `EXCLUDED_ASSUMPTION` (DFT-level physical filter)
- `DEFERRED` (with a reason: `rank_cut(rank, cut)`, `budget`, `awaiting_dependency`)
- `INVALID`
- `DUPLICATE`
- `ERROR` (with a retry count)
- At S5–S8 only, the decide outcome of `rtt.decide`: `rejectedProven`, `rejectedAssuming`,
  `screenedOut`, `undecided`, `pass`. `screenedOut` and `undecided` are only possible at S5.
- Core certify verdicts (REJECTED_PROVEN, UNDECIDED, CERTIFIED_CANDIDATE) are stored as
  *evidence* in `verdicts` (§6.2). They are never a stage outcome by themselves.

| Stage | Runs | Computes | Removes a candidate when | Emits (and nothing else) | Source of rule/threshold |
|---|---|---|---|---|---|
| **S0** | load `targets.toml`, `tiers.toml`, `machine.toml` | floor = a0·T\* = 1299.44·s1 K; g table; aim box | – | – | doc 01 §0.4; doc 02 §5; doc 11 |
| **S1** source | lanes G1–G5 (MatterGen), D (database), C (controls); composition/size prefilter | cost class (§8.1), H_f, element set | element without an SG15 pseudopotential; N_atoms(primitive) > `n_atoms_max`; element outside MatterGen's domain (noble gases, Tc, Pm, Z > 84) | `PASS`, `SCREENED_OUT` | doc 10 §1; tiers.toml |
| **S2** structure | validity (min distance, density), fingerprint and StructureMatcher dedupe against the ledger, novelty against Alex-MP | `dup_of`, `novel` | invalid; duplicate | `PASS`, `INVALID`, `DUPLICATE` | – |
| **S3** mlff (GPU) | **in this order:** (a) batched MatterSim-5M relax at 0 GPa; (b) MatterSim-only hull, plus the DFT-reference hull as a second opinion; (c) cheap features: H_f (S2), volume per H (S8), H–H and M–H geometry (S4), ALIGNN numbers (recorded; weight 0 for H and light-element candidates, doc 11 §2.3); (d) MatterSim phonons for (c) survivors: min ω, ω_max, H-projected DOS moments, harmonic ZPE (ranking only) | `e_hull_ml_self`, `e_hull_ml_dftref`, features, `score_s3` | `e_hull_ml_self > window_hi + ml_margin` (0.35 + 0.07 eV/atom; H7 and doc 10 §3.3) | `PASS`, `SCREENED_OUT(H7-ml)`, `DEFERRED(rank_cut)` | doc 02 H7/S5; doc 10 |
| **S4** dft-electronic | pw.x vc-relax at 0 GPa (S4 profile); scf; nscf DOS; spin-polarised relax from a ferromagnetic start; band edges → E_F (G-VALID definition); n_sheets; Belli φ and H_DOS (doc 11 T1c); heuristic hull | `gap`, `moment`, `N_EF`, `H_DOS`, `E_F`, `n_sheets`, `phi`, `Tc_phi`, `e_hull_dft_heur` | PBE gap > 0 (**H1**); moment > `moment_tol` (**H2**); heuristic hull above window (no X-HULL here, R12 #8) | `PASS`, `EXCLUDED_ASSUMPTION(H1|H2)`, `SCREENED_OUT(H7-heur)`, `DEFERRED`, `ERROR` | doc 02 H1, H2, H7 |
| **S5** coarse DFPT | profile `coarse` (q 3³ or 4³, k_fit ≥ 24³-equivalent density, ecut ≈ ⅔ converged, σ ladder, mode-positive admission); core `certify_spectra` in decide mode (evidence only); `rtt.decide.decide_t(tier=coarseDFPT)` | A, λ, m, ω_log, ω₂ per σ; σ-plateau values; min coarse ω; coarse ZPE | **only** `A_coarse,hi(σ) < floor / M_c` for every σ on the plateau, with M_c = **⟨doc13⟩** (provisional 2, doc 11 §4) → SCREENED_OUT. Imaginary coarse modes never remove: they set `needs_phonon_check` | `PASS`, `SCREENED_OUT(coarse-A)`, `DEFERRED`, `ERROR`; decide outcome ∈ {`screenedOut`, `undecided`, `pass`} recorded | doc 11 §3.1, §4; ⟨doc13⟩; Lean `decide_coarse` |
| **S6** full DFPT | re-relax at the S6 cutoff and k (or a force/stress check: max force ≤ 1 mRy/bohr, \|stress\| ≤ 0.5 kbar), R12 #14; profile `full` on **two** q meshes (n³ and the next finer), ≥ 3 admitted σ on the finest k; core `certify_spectra` in bracket mode; `rtt.decide.decide_t(tier=fullDFPT)` | ENVELOPE inputs (δ_A, δ_λ, δ_m, ε_ω); per-alternative A, λ, m; certified bracket; full-grid min ω | removed only by decide outcome `rejectedProven` or `rejectedAssuming` (needs G-VALID: harmonic-only H/Li/Be gives `undecided`) | decide outcome ∈ {`rejectedProven`, `rejectedAssuming`, `undecided`, `pass`}; verdict evidence; `ERROR` | doc 01 §0.2, §0.4; Lean `decide_rejectedProven_sound` |
| **S7** validity closure | (a) own-DFT reference hull + harmonic ZPE → H7 (ZPE-corrected E_hull ≤ 0.35; synthesis route stated above 0.1); (b) anharmonic check (G-VALID 2, H6); (c) MD/PIMD kinetic persistence at ≥ 300 K (H6); (d) k- or band-resolved coupling (EPW) when n_sheets > 1 or the N1-ISO margin is < 25 % | updated `Validity` record | (a) ZPE-corrected E_hull > 0.35 → `EXCLUDED_ASSUMPTION(H7)`; otherwise decideT is re-run | `EXCLUDED_ASSUMPTION(H7|H6)`, re-decided outcome | doc 01 G-VALID; doc 02 H6, H7, §4.2 |
| **S8** converged | EPW protocol (EPH_CONVERGENCE_PROTOCOL v3) for finalists | converged alternatives | as S6 | decide outcome with tier = converged | – |
| **S9** evidence | assemble the evidence statement (§3.6) | ranking keys | – | – | – |

Notes on the order (R12 #7, #11):
- ML proxies (ALIGNN, φ geometry) cost seconds and sit in S3. They come before every
  core-hour stage.
- Belli φ needs an SCF, so it is in S4. It shares the S4 scf at no extra cost.
- The Γ-only DFPT tier and "e-ph at Γ" from v1 are **removed**. q = Γ e-ph is not λ, and
  doc 11 §3.1 measured that q 2³ is off by 50–60 %. A Γ dynamical-matrix check is kept as an
  optional diagnostic: it runs for the Mg₂IrH₆ control to settle doc 10's open question about
  MatterSim's reliability. It never removes candidates.
- L-MOM1 is not used as a screen: it is too weak for X–H chemistry (R12 #7). Proxy
  thresholds in S3 use ω₂ from the MatterSim phonon DOS, weighted to light-element modes, as a
  **ranking** feature only: λ_needed = 2·floor/ω₂^proxy. Doc 01 RANK-AIM gives the aim box.

### 3.2 Two decision layers at S5–S8, and how they combine

1. **Core `certify_spectra`.** This is `classify` followed by the input policy. It proves
   statements about Tc_model of the given spectra under envelope e and stores a certificate.
   - Inputs: the admitted alternatives, extended by frequency-scaled copies. For each admitted
     σ, add Ω' = sΩ, p' = p/s for s ∈ {1 − ε_ω, 1 + ε_ω}. This is the Hopfield scaling of
     doc 01 §5.6 (A ∝ ω⁻¹); frequency uncertainty enters as alternatives, as ARCH errata F5
     requires.
   - Envelope: e = (1 + δ_A^q)(1 + ε_sys) − 1. If e ≥ 1, no rejection is possible; this is
     recorded.
   - Per alternative, R-UNIV then tests (1 + e)·A_s', which is ≥ the ENVELOPE A_hi of doc 01,
     because (1 + δ^σ)(1 + δ^q) ≥ 1 + δ^σ + δ^q. The core check is therefore at least as
     conservative as N1-ISO.
   - Mode: decide at S5, bracket at S6–S8.
2. **`rtt.decide.decide_t`.** A Python mirror of Lean `decideT`, over `Inputs`, which are built
   *only* from S4–S7 records (§6.1). It implements G-VALID, ENVELOPE, N1-KAPPA, N1-ROWMAX,
   N1-DOS, N1-ISO and NM-G exactly as in doc 01 §0.4. Comparisons are done in exact rationals
   (`fractions.Fraction` of the IEEE inputs). The g lookup uses the certified table
   `docs/theory/certificates/nm_*.json` (doc 03).

**Combination rule** (contract R1c, §4). The stage outcome is the `decide_t` outcome, with one
extra downgrade. If `decide_t = rejectedProven` and the firing rule is N1-ISO, the core
verdict for the same `tier_run` must be REJECTED_PROVEN after `apply_input_policy`. If it is
not, the outcome becomes `undecided` and the reason `CROSSCHECK_DISAGREE` is recorded. The
multiband rules (N1-KAPPA, N1-ROWMAX, N1-DOS) have no single-band core counterpart and rest
on `decide_t` alone (Lean `n1_rowmax_stable`, `n1_dos_stable`, `stable_of_kappa_bound`). The
`BlockedRej` input of `decide_t` is filled from the core `blocking_tags(spectra)[0]`, so
both layers see the same input-quality tags.

CERTIFIED_CANDIDATE from the core engine is positive evidence. It feeds ranking (T_low)
and the evidence statement. `decideT` has no certification outcome. A certification whose
G-VALID fails is reported as "certified in the harmonic isotropic model; outside the validity
domain (item k)".

### 3.3 Generation plan (only shipped checkpoints, plus one planned adapter)

Measured cost: ≈ 10 s per conditional structure and 5 s per unconditional structure at batch 64;
about 8.6k conditional structures per GPU-day; batch 128 is no faster (doc 10 §3.2). Every lane
uses a shipped checkpoint with **only its own property set**, because the adapter refuses any
other combination (`CHECKPOINT_PROPERTIES`; R12 #9).

| Lane | Share | Checkpoint | Conditioning | Source |
|---|---|---|---|---|
| G1 | ~55 % | `chemical_system_energy_above_hull` | `{chemical_system: X-Y-H from targets.chemsys_list, energy_above_hull: U(0, 0.3)}`, γ = 2 | doc 02 5(d)1, S1, S5, S10 |
| G2 | ~15 % | `space_group` | one group per run from the S6 list {191, 194, 139, 123, 166, 221, 223, 225, 227, 229} | doc 02 5(d)2, S6 |
| G3 | ~10 % | `dft_band_gap` and `dft_mag_density` | {dft_band_gap: 0}; {dft_mag_density: 0} (separate runs) | doc 02 5(d)3, H1, H2 |
| G4 | 10–20 % | `chemical_system`; `ml_bulk_modulus` (high) | non-hydride systems (Li–B–C, Mg–B–C, Be–B, Li–Be–B, B–C–N) | doc 02 5(c) |
| G5 | ~5 % | `mattergen_base` | none (bias control) | doc 02 5(d)5 |
| GA (planned) | replaces part of G1 once validated | fine-tuned adapter on the scalar **A** (property id `rtt_A_K`) | `{rtt_A_K: q_hi}`, where q_hi is the 90–95th percentile of the label distribution, never 1300 K (extrapolation warning, doc 10 §4.3) | labels prepared by another agent (source TBD); doc 10 §4.2 command template; ~2.6 GPU-h for 7k labels |

Rules for the A-adapter:
- **Training.** Labels come from an external DFT set. Our own S5/S6 values are a held-out
  validation set; they are too few to train on (doc 10 §4.3).
- **Acceptance.** The adapter is accepted for production only if, on 500 samples, the
  MatterSim-screened, S4-metallic fraction and the S5 A distribution are not worse than G1 at
  equal GPU cost. This is measured on the first ≤ 10 S5 runs from each lane, by a Mann–Whitney
  test on A_coarse, reported, not a gate.
- **Feedback.** The cheap labels that exist at volume (S4 H_DOS and φ, about 25 a day; S5
  A_coarse, 2–6 a day) are written to `labels` (§6.2) and used to re-weight the G1 chemical
  system list weekly. Each system's sampling weight is proportional to the posterior mean of
  max(A_coarse)/floor over its candidates, with a floor weight of 0.2 so that no system starves.
  This is the replacement for the Beta bandit of the initial discovery loop. It learns from S4/S5, not from S6.

**Database lane D.** Offline Alex-MP reference LMDB, 846k entries (doc 10 §3.3). Query metals
(PBE gap 0) with E_hull ≤ 0.1 eV/atom in the `targets.chemsys_list` systems and in the 5(c)
families, ordered by H_f and N_atoms. These are concrete, known compounds and the baseline the
generator has to beat. They enter at S3 with `origin.lane = "D"`.

**Control lane C** (positive controls and calibration; contracts R5/R6). These run through every
stage at a **surrogate target** T\*_s, whose floor a0·T\*_s replaces 1299.44 K in S5's margin and
in `decide_t`. Each control must *pass* each heuristic stage at T\*_s. Any control that a
screen would have removed is reported as a screen false negative.

| Control | Why | T\*_s | Notes |
|---|---|---|---|
| Mg₂IrH₆ (Fm-3m, 0 GPa) | best predicted ambient hydride (Tc 160–175 K, Dolui et al. 2024); MatterSim phonons fail on it | 80 K | the critical S3 control; also the Γ diagnostic |
| H3S Im-3m (200 GPa) | only H anchor with full DFPT already on disk | 150 K | pressure is set by a control override; outside the 0–1 GPa scope |
| MgB2 (P6/mmm) | multiband; tests G-VALID item 5 and N1-MB | 20 K | blocked until the ph.x crash (W4) is diagnosed |
| Nb, Pb, Al | calibration of coarse vs full; V2/V3 of the pre-registered validation | ½·Tc_exp | spectra exist |
| one 6–8-atom light-element metal from lane D (e.g. a known boride or carbide superconductor; to be chosen from Alex-MP) | the first cost anchor for class C2/C3 (§8) | ½·Tc_exp | needed before any throughput claim (R12 #6, #15) |

### 3.4 Coarse margin and its calibration (placeholder for doc 13)

- S5 removes only by `SCREENED_OUT(coarse-A)` if `A_coarse,hi(σ) < floor / M_c` for every σ on
  the plateau.
  - A_coarse,hi = A_coarse·(1 + ε_sys).
  - M_c = **⟨doc13: margin⟩**. Until doc 13 lands, M_c = 2 with `provisional = true` (doc 11 §4).
    Every outcome made with a provisional value is tagged `PROVISIONAL_MARGIN`, and S9 lists them.
- M_c is used only in the removal direction. It is never passed to the core engine as
  `rel_envelope`. (`Spectrum` requires e < 1, and the envelope covers weights only; R12 #3.)
- The coarse core verdict is stored with `spectrum_tier = 'coarse'` and never drives state.
- Calibration scope, stated as a limit: the references are bench-grade, not converged, and they
  are 1–4-atom cells. M_c relates coarse to bench-grade, and it is re-checked on the first
  6–8-atom control and on every S6 run (`calibration` table, §6.2).

### 3.5 Audits (R12 #5)

Audits use fixed counts sized to downstream capacity, not fractions:

| Audit | Count | Selection | What is measured (continuous) |
|---|---|---|---|
| A3→4 | 10 / week | candidates removed or deferred at S3, sampled with probability ∝ exp(−\|score − cut\|/w); inverse-probability weights stored | E_hull^ML vs DFT heuristic hull; ω₂^proxy vs S5 ω₂; MatterSim phonon instability vs DFPT |
| A4→5 | 2 / week | S4 `SCREENED_OUT(H7-heur)` and `DEFERRED(rank_cut)` near the cut | φ and H_DOS vs A_coarse |
| A5→6 | 1 / month (C1 class only) | S5 `SCREENED_OUT(coarse-A)` with the largest A_coarse/floor | A_coarse vs A_full: feeds ⟨doc13⟩ |
| Controls | all, at T\*_s | lane C | recall of every screen at the surrogate floor |

At T\* = 300 K we expect ≈ 0 positives, so the false-negative *rate* at the real target cannot
be estimated. What the audits deliver are proxy *error distributions* and recall at the
surrogate floors. S9 states exactly this.

### 3.6 Final ranking and evidence statement (goal 3; R12 #10)

**Ranking** is lexicographic over candidates that reached S6, never by T_up:
1. decide outcome (`pass` > `undecided` > `rejectedAssuming` > `rejectedProven`);
2. **certified T_low**: the combined core bracket's lower end (min over all alternatives,
   ARCH F7) at μ\*_hi = 0.13; missing counts as 0;
3. **A_mid**: the median over admitted σ, on the finest q mesh, of A_nom (no envelope);
   A_lo = min over alternatives is reported alongside.

T_up, max σ A and A_hi are reported as ceilings only.

**Evidence statement template** (`rtt.report.evidence`). Every item reads "checked (value,
record id)" or "NOT CHECKED":
- structure, origin lane, novelty;
- H1, H2 (S4);
- harmonic stability on the full q grid of the *final* structure (S6), and the Γ diagnostic;
- G-VALID items 1–5, each with its value;
- ENVELOPE inputs;
- decide outcome and firing rule;
- core verdict, bracket, rule ids, assumptions, flags, and `recheck` result;
- ZPE-corrected E_hull against own-DFT references, plus a synthesis route if E_hull > 0.1 (H7);
- anharmonic check (method, frequency shift, whether α²F was recomputed);
- MD/PIMD persistence (H6);
- μ\* sensitivity: the bracket at μ\*_hi = 0.10, 0.13 and 0.15;
- multiband/k-resolved status;
- provisional margins used upstream;
- audit and control recall summary.

---

## 4. Soundness contract

- **R1** (= contract S1). Only core `classify` produces core REJECTED_PROVEN or
  CERTIFIED_CANDIDATE.
- **R1b** (new, R12 #1). Every core verdict passes through core `certify_spectra`, and
  so through `ensure_cutoff_assumption`, `annotate_model`, `trust_base` and
  `apply_input_policy`. `rtt.certify` calls nothing else. Verdict objects are not built anywhere
  else (enforced by an import-graph test).
- **R1c** (new). Only `rtt.decide.decide_t` produces the decide outcomes `rejectedProven` and
  `rejectedAssuming`, and only at `tier ∈ {fullDFPT, converged}`, as in Lean
  `not_rejected_of_not_proofGrade`. An N1-ISO `rejectedProven` additionally needs agreement
  with the core verdict (§3.2).
- **R2.** Heuristic stages (S1, S3, the heuristic parts of S4, S5's coarse margin, scoring,
  generation feedback) emit only `SCREENED_OUT` or `DEFERRED`. Each carries criterion, value,
  threshold, `source_id` and the provisional flag.
  - **MatterSim phonons and coarse imaginary modes never remove a candidate.**
  - **Public ML Tc/λ models never remove a candidate** and have weight 0 for H and
    light-element chemistries (doc 11).
- **R3.** `EXCLUDED_ASSUMPTION` comes only from S4 (H1, H2) and S7 (H7 with own-DFT references
  and ZPE; H6), each with its criterion. Every filter can be disabled.
- **R4** (R12 #4). A proven core verdict or a removal decide outcome is accepted into the
  ledger only if all of the following hold, enforced by triggers (§6.2):
  - its `tier_run` has `kind = 'dfpt'` and a profile mapped to its `spectrum_tier`;
  - every spectrum's source file is an `artifacts` row whose sha256 matches;
  - `spectrum_digest` is stored per alternative;
  - the core `recheck` returned `ok` before insert.
- **R5** (= contract S4). Every candidate has exactly one current outcome per stage it reached, and
  the summary accounts for all candidates, audits and controls. Counts of REJECTED_PROVEN are
  always split by `spectrum_tier`.
- **R6** (= contract S5, S6). Alternatives combine as a conjunction; flags never upgrade.
- **R7.** Every threshold has a `source_id` (doc section) and an optional `provisional` flag.
  The loader refuses a threshold without one. `targets_sha256` and `tiers_sha256` are stored in
  every outcome.
- **R8** (control lane). A heuristic screen whose control recall at T\*_s drops below 1 is
  disabled automatically (logged), until a person re-enables it.

---

## 5. Package layout and interfaces (`src/rtt/`)

### 5.1 Module list

```
src/rtt/
  __init__.py
  config.py        load_targets(), load_tiers(), load_machine(), sha256_of(); Threshold(value, source_id, provisional)
  records.py       dataclasses of §6.1; to_json()/from_json(); schema "rtt.candidate/2"
  units.py         re-export of core/dfpt/units (K, meV, Ry, cm-1)
  targets/
    spec.py        TargetSpec (floor, s1, s2, g_table, aim_box, H1..H7, S1..S11, chemsys_list)
    gtable.py      g lookup from docs/theory/certificates/nm_*.json (doc 03), conservative 'up' rule
  source/
    mattergen.py   adapted adapter: generate(lane: LaneSpec, n: int, seed: int) -> list[RawStructure]
    lanes.py       LaneSpec, plan_day(machine, targets, weights) -> list[LaneJob]
    database.py    query_alexmp(chemsys, gap_max=0.0, ehull_max=0.1) -> list[RawStructure]
    controls.py    CONTROLS: list[ControlSpec]; inject_controls(ledger) -> None
    prefilter.py   prefilter(rs: RawStructure, tiers) -> Outcome
    finetune.py    export_labels(ledger, prop="rtt_A_K") -> Path; train_adapter(...) (doc 10 §4.2)
  structure/
    canonical.py   canonicalize(rs) -> Structure (spglib primitive, standardized), structure_sha256(s)
    validate.py    validate(s, tiers) -> Outcome
    dedupe.py      wraps core fingerprint + StructureMatcher: find_duplicate(s, ledger) -> str | None
  mlff/
    relax.py       relax_batch(structs, model="MatterSim-v1.0.0-5M", fmax=0.02) -> list[RelaxResult]
    hull.py        ehull_self(struct, energy, chemsys) -> float; ehull_dftref(...) -> float
    features.py    features(struct) -> dict (H_f, vol_per_H, min_HH, motif)
    phonons.py     ml_phonons(struct) -> PhononSummary (min_freq, omega_max, omega2_light, zpe_ml)
    alignn.py      alignn_numbers(struct) -> dict  (recorded only)
  dft/
    profiles.py    PROFILES: {"s4": ..., "coarse": ..., "full": ..., "ref_hull": ...} (ecut, k density, q, sigmas)
    workflow.py    run_profile(struct, profile, workdir, launcher) -> DFPTRun  (wraps core workflow)
    electronic.py  run_electronic(struct, workdir) -> ElectronicResult (gap, moment, N_EF, E_F, n_sheets)
    belli.py       phi_and_hdos(workdir) -> (phi, H_DOS, Tc_phi)   (from src/scoring/phi_elf.py)
    refhull.py     ref_hull(chemsys, profile) -> RefHull (own-DFT, cached by pseudo md5+profile); zpe(...)
    prune.py       prune_q_scratch(workdir, q_done) -> int; disk_guard(path, min_free_gb) -> None
    collect.py     spectra_from_run(run: DFPTRun) -> (list[Spectrum], list[ArtifactRef])  (core parsers + admit_smearings)
  certify/
    run.py         certify(material_id, spectra, ccfg, mode) -> Verdict   (calls core certify_policy.certify_spectra ONLY)
    alternatives.py freq_scaled(spectra, eps_omega) -> list[Spectrum]; envelope_e(delta_A_q, eps_sys) -> float
    recheck.py     recheck_or_raise(verdict, spectra) -> None
  decide/
    inputs.py      build_inputs(record, tier) -> DecideInputs   (Validity, Envelope, MBInputs, tags)
    core.py        gvalid(x), fires_kappa/rowmax/dos/iso/nm(cfg, x), blocked_rej(x), core_b(...), decide_t(cfg, x) -> DecideOutcome
  score/
    s3.py          score_s3(features, phonons, targets) -> float     (doc 11 T1a; ALIGNN weight 0 for H/light)
    s4.py          score_s4(electronic, phi) -> float                (doc 11 T1c)
    rank.py        final_rank(records) -> list[RankRow]               (§3.6)
  engine/
    stage.py       Stage protocol (§5.3)
    stages.py      S1..S9 implementations
    scheduler.py   lanes G/S/F, budgets, rank cuts, audits; next_jobs(ledger, machine) -> list[Job]
    runner.py      run_job(job) -> None   (writes tier_runs 'started' first)
  ledger/
    schema.sql     tables + triggers (§6.2)
    store.py       Ledger(path); insert_*(); current_outcome(cid, stage); export_jsonl(path)
    artifacts.py   register(path) -> ArtifactRef(sha256, bytes); prune events
  report/
    evidence.py    evidence_statement(cid) -> str (Markdown)
    summary.py     daily_summary(ledger) -> str
  cli.py           rtt run|step|status|redecide|evidence|export|controls
  core/    certify/, spectrum.py, status.py, constants.py, validity/, dfpt/, certify_policy.py,
                   discovery_fingerprint.py, sources_{dedupe,fileformats,composition,local}.py, PINS.json
```

The scripts already in `src/ml/`, `src/scoring/`, `src/theory/` and `src/certify/` stay as
research scripts. Their logic is ported, not imported:
- `ehull_self.py`, `ehull_ref.py` → `rtt.mlff.hull`;
- `mattersim_screen.py` → `rtt.mlff.relax` and `rtt.mlff.phonons`;
- `phi_elf.py` → `rtt.dft.belli`;
- `coarse_dfpt.py`, `coarse_calibration.py` → `rtt.dft` and the `calibration` table;
- `src/certify/` (doc 03 certificates) is a data producer for `rtt.targets.gtable`.

### 5.2 Core certify policy (R12 #1)

`core/certify_policy.py` holds `base/pipeline.py` lines 83–133 and 420–650
(`ConfigError`, `validate_certify_config`, `_CERTIFY_KEYMAP`, `_CERTIFY_ORCH_ONLY`,
`_certify_params`, `certify_spectra`, `VERIFIED_PYTHON_FLINT`, `trust_base`, `_TAG_MASS`,
`_EXCL_ITEM`, `_sigma_exclusion_directions`, `blocking_tags`, `apply_input_policy`,
`annotate_model`, `cutoff_description`, `ensure_cutoff_assumption`), byte-identical apart from
imports.
- **Config.** `[certify] allow_unverified_rejection` and `allow_unverified_certification` are
  **forced to false** by `rtt.certify.run`, so that Lean `decideT` (which has no override) and
  the core policy agree.
- **Tests.** Core with it: `test_dfpt_smearing_admission.py`, `test_dfpt_soundness.py`,
  `test_certify*.py`, `test_correspondence.py`, the oracle suite, and the input-policy cases of
  `tests/audit/` and `tests/audit_certify/` (run with `OMP_NUM_THREADS=2`, in parts).

### 5.3 Stage interface

```python
class Stage(Protocol):
    id: str                                   # "S4"
    kind: Literal["heuristic", "physical", "dfpt", "report"]
    lane: Literal["G", "S", "F", "none"]
    requires: tuple[str, ...]                 # stages whose PASS is needed, e.g. ("S3",)
    allowed: frozenset[str]                   # outcomes this stage may emit (table 3.1); enforced by ledger trigger too
    def estimate_cost(self, rec: CandidateRecord, machine: Machine) -> CostEstimate: ...
    def run(self, rec: CandidateRecord, ctx: Context) -> StageResult: ...       # compute; writes tier_runs + artifacts
    def decide(self, res: StageResult, th: Thresholds, tgt: TargetSpec) -> Outcome: ...  # pure
```

`rtt redecide --stage S5` re-evaluates `decide` from stored metrics when `tiers.toml` changes,
without new compute. Each re-decision appends new `outcomes` rows; nothing is overwritten.

---

## 6. Data

### 6.1 Candidate record (`rtt.candidate/2`, JSON; the canonical copy is assembled from the ledger)

```json
{
  "schema": "rtt.candidate/2",
  "candidate_id": "c-3f9a1c2b7d10",
  "structure_sha256": "…",
  "formula": "Mg2IrH6",
  "structure": {"lattice": [[0,3.33,3.33],[3.33,0,3.33],[3.33,3.33,0]], "species": ["Mg","Mg","Ir","H","H","H","H","H","H"],
                "frac_coords": [[0.25,0.25,0.25], "…"], "spacegroup": 225, "n_atoms": 9},
  "pressure_GPa": 0.0,
  "cost_class": "C2",
  "origin": {"lane": "C", "control": "Mg2IrH6", "T_star_surrogate_K": 80.0,
             "generator": null, "checkpoint": null, "checkpoint_sha256": null, "conditioning": null,
             "guidance": null, "batch_id": null, "seed": null, "flags": ["KNOWN_COMPOUND"]},
  "context": {"targets_sha256": "…", "tiers_sha256": "…", "code_git_sha": "…", "core_version": "core-1"},
  "stages": {
    "S3": {"outcome": "PASS", "tier_run_id": 101,
           "metrics": {"e_hull_ml_self_eV": 0.004, "e_hull_ml_dftref_eV": -0.0003, "H_f": 0.667,
                       "vol_per_H_A3": 12.3, "ml_min_freq_cm1": -363.0, "ml_omega2_light_K": 1650.0,
                       "alignn": {"Tc_K": 1.4, "weight": 0.0}},
           "score": 0.71},
    "S4": {"outcome": "PASS", "tier_run_id": 140,
           "metrics": {"gap_eV": 0.0, "moment_muB": 0.0, "N_EF_per_eV_atom": 0.12, "H_DOS": 0.41,
                       "E_F_eV": 0.35, "n_sheets": 2, "phi": 0.42, "Tc_phi_K": 46.0, "e_hull_dft_heur_eV": 0.06},
           "flags": ["DOPING_SENSITIVE"]},
    "S5": {"outcome": "PASS", "tier_run_id": 188, "profile": "coarse", "q": [4,4,4],
           "metrics": {"A_K_by_sigma": [910.0, 880.0, 860.0], "plateau_sigmas": [3,4,5],
                       "lambda_by_sigma": ["…"], "min_freq_cm1": 120.0, "zpe_meV_atom": 145.0},
           "margin": {"M_c": 2.0, "provisional": true, "source_id": "doc11#4"},
           "decide": {"outcome": "pass", "decide_input_id": 77},
           "verdict_ids": [55]},
    "S6": {"outcome": "undecided", "tier_run_id": 230, "profile": "full",
           "decide": {"outcome": "undecided", "reason": "GVALID:anharmonic_check_missing", "decide_input_id": 90},
           "verdict_ids": [61]}
  },
  "flags": ["ANHARMONIC_H", "HYPOTHETICAL_STRUCTURE?"],
  "current": {"stage": "S7", "state": "QUEUED"}
}
```

The numbers above are **illustrative placeholders**, not results (the Mg₂IrH₆ S4–S6 values
have not been computed).

`DecideInputs` (one row in `decide_inputs`; the field names are those of Lean
`TargetLogic.Decide`):

```json
{"tier": "fullDFPT",
 "validity": {"EF": 0.35, "omegaMax": 2500.0, "maxImagModeK": 0.0, "containsHLiBe": true,
              "anharmonicCheck": false, "anharmonicShift": null, "alpha2FRecomputedAnharmonic": false,
              "muStarNegativeSuspected": false, "admittedSmearings": 7, "qMeshes": 2,
              "nSheets": 2, "bandResolvedAvailable": false, "dosFractionsKnown": true, "kResolvedCheckDone": false},
 "envelope": {"Anom": 0.0, "lamNom": 0.0, "mNom": 0.0, "dA": 0.0, "dLam": 0.0, "dM": 0.0, "epsSys": 0.10, "epsOmega": 0.02},
 "mb": {"kappaHi": null, "rowMaxHi": null, "wMin": 0.31},
 "tags": ["A-SIGMA-EXCLUDED(sigma1:lam=…:A=…K)"],
 "cfg": {"a0": "433148145/100000000", "Tstar": 300, "gUp_source": "docs/theory/certificates/SUMMARY.json"},
 "provenance": {"tier_run_ids": [188, 230], "artifact_ids": ["…"], "record_ids": {"S4": 140, "S7a": null}}}
```

Every `validity` field is filled from a named record. A field without a source record is
`null`, and `null` makes G-VALID fail. Missing data therefore always yields `undecided`.

### 6.2 Ledger (SQLite, WAL; `runs/<run_id>/ledger.sqlite`; JSONL export + `manifest.sha256` committed)

```sql
CREATE TABLE meta        (key TEXT PRIMARY KEY, value TEXT NOT NULL);           -- schema, run_id, code sha, core commit+hashes
CREATE TABLE inputs      (id INTEGER PRIMARY KEY, kind TEXT CHECK(kind IN ('targets','tiers','machine')),
                          sha256 TEXT UNIQUE NOT NULL, content TEXT NOT NULL, created TEXT NOT NULL);
CREATE TABLE structures  (sha256 TEXT PRIMARY KEY, json TEXT NOT NULL);
CREATE TABLE candidates  (id TEXT PRIMARY KEY, structure_sha256 TEXT NOT NULL REFERENCES structures,
                          lane TEXT NOT NULL CHECK(lane IN ('G1','G2','G3','G4','G5','GA','D','C')),
                          origin_json TEXT NOT NULL, cost_class TEXT, pressure_GPa REAL NOT NULL, created TEXT NOT NULL);
CREATE TABLE duplicates  (candidate_id TEXT NOT NULL, duplicate_of TEXT NOT NULL, method TEXT NOT NULL);
CREATE TABLE tier_runs   (id INTEGER PRIMARY KEY, candidate_id TEXT NOT NULL REFERENCES candidates,
                          stage TEXT NOT NULL, kind TEXT NOT NULL CHECK(kind IN ('ml','dft','dfpt','epw','md','report')),
                          profile TEXT, attempt INTEGER NOT NULL, state TEXT CHECK(state IN ('started','done','failed')),
                          params_json TEXT NOT NULL, versions_json TEXT NOT NULL, lane TEXT,
                          started TEXT, finished TEXT, wall_s REAL, core_h REAL, gpu_h REAL, metrics_json TEXT);
CREATE TABLE artifacts   (id INTEGER PRIMARY KEY, tier_run_id INTEGER NOT NULL REFERENCES tier_runs,
                          path TEXT NOT NULL, sha256 TEXT NOT NULL, bytes INTEGER NOT NULL, retained INTEGER NOT NULL);
CREATE TABLE outcomes    (id INTEGER PRIMARY KEY, candidate_id TEXT NOT NULL, stage TEXT NOT NULL,
                          outcome TEXT NOT NULL, criterion TEXT, value TEXT, threshold TEXT, source_id TEXT,
                          provisional INTEGER NOT NULL DEFAULT 0, tier_run_id INTEGER REFERENCES tier_runs,
                          targets_sha256 TEXT NOT NULL, tiers_sha256 TEXT NOT NULL, origin TEXT NOT NULL, created TEXT NOT NULL);
CREATE TABLE decide_inputs (id INTEGER PRIMARY KEY, candidate_id TEXT NOT NULL, json TEXT NOT NULL, created TEXT NOT NULL);
CREATE TABLE decisions   (id INTEGER PRIMARY KEY, candidate_id TEXT NOT NULL, stage TEXT NOT NULL,
                          tier TEXT NOT NULL CHECK(tier IN ('ml','coarseDFPT','fullDFPT','converged')),
                          decide_input_id INTEGER NOT NULL REFERENCES decide_inputs, outcome TEXT NOT NULL
                          CHECK(outcome IN ('rejectedProven','rejectedAssuming','screenedOut','undecided','pass')),
                          firing_rule TEXT, crosscheck_verdict_id INTEGER REFERENCES verdicts, origin TEXT NOT NULL, created TEXT NOT NULL);
CREATE TABLE verdicts    (id INTEGER PRIMARY KEY, candidate_id TEXT NOT NULL, tier_run_id INTEGER NOT NULL REFERENCES tier_runs,
                          spectrum_tier TEXT NOT NULL CHECK(spectrum_tier IN ('coarse','full','converged')),
                          status TEXT NOT NULL, T_low_K REAL, T_up_K REAL, T_star_K REAL NOT NULL, mu_hi REAL NOT NULL,
                          rel_envelope REAL NOT NULL, rule_ids TEXT, assumptions TEXT, flags TEXT,
                          certificate_artifact_id INTEGER NOT NULL REFERENCES artifacts, recheck_ok INTEGER NOT NULL,
                          origin TEXT NOT NULL, created TEXT NOT NULL);
CREATE TABLE verdict_spectra (verdict_id INTEGER NOT NULL REFERENCES verdicts, alt_index INTEGER NOT NULL,
                          spectrum_digest TEXT NOT NULL, source_artifact_id INTEGER NOT NULL REFERENCES artifacts);
CREATE TABLE labels      (candidate_id TEXT NOT NULL, name TEXT NOT NULL, value REAL NOT NULL, stage TEXT NOT NULL, tier_run_id INTEGER);
CREATE TABLE calibration (id INTEGER PRIMARY KEY, material TEXT, coarse_run INTEGER, ref_run INTEGER,
                          sigma INTEGER, A_coarse REAL, A_ref REAL, note TEXT);
CREATE TABLE audits      (candidate_id TEXT NOT NULL, audit TEXT NOT NULL, inclusion_prob REAL NOT NULL, created TEXT NOT NULL);
CREATE TABLE queue       (candidate_id TEXT, stage TEXT, priority REAL, state TEXT);   -- derived, mutable, rebuildable
```

Triggers (all insert-only tables: no UPDATE, no DELETE; `queue` is exempt):

- `outcomes_allowed`: `(stage, outcome)` must be in the §3.1 table. For example, S3 cannot
  insert `EXCLUDED_ASSUMPTION`, and S5 cannot insert any rejection.
- `verdicts_bind`:
  - for a proven status: `recheck_ok = 1`;
  - the `tier_runs` row has `kind = 'dfpt'` (or `'epw'` for converged) and `state = 'done'`;
  - the profile ↔ `spectrum_tier` map holds (`coarse` → `'coarse'`, `full` → `'full'`, `epw` →
    `'converged'`);
  - the certificate artifact belongs to the same `tier_run`;
  - `origin = 'rtt.core.certify_policy.certify_spectra'`.
- `verdict_spectra_bind`: the source artifact belongs to the verdict's `tier_run` (or to a
  `tier_run` of the same candidate with kind `dfpt`).
- `decisions_guard`:
  - `rejectedProven` and `rejectedAssuming` require tier ∈ (`fullDFPT`, `converged`), a
    `crosscheck_verdict_id` when `firing_rule = 'N1-ISO'`, and `origin = 'rtt.decide.core.decide_t'`;
  - `screenedOut` requires a non-proof-grade tier;
  - the decide input's `tier` must equal the tier of the tier_runs it references.

The origin strings remain self-declared. The *binding* to hashed DFPT artifacts is what makes a
forged tier label detectable. `rtt export` re-verifies every artifact hash it can still read and
marks pruned ones.

---

## 7. Configuration files

### 7.1 `targets.toml` (owned by the theory docs; values copied from doc 01 §0.4 and doc 02 §5)

```toml
schema = "rtt.targets/1"
[tc]
T_star_K = 300.0                  # s1 = T_star/300
a0 = "433148145/100000000"        # source: core constants.py A0_FRACTION
floor_source = "doc01#0.4:N1-ISO"

[envelope]                        # doc01#0.4:ENVELOPE
eps_sys = { value = 0.10, source = "doc01#0.2b", provisional = false }
eps_omega_min = { value = 0.02, source = "doc01#0.2b" }

[gvalid]                          # doc01#0.4:G-VALID (literal values)
nonadiabatic_max = 0.1
imag_mode_tol_K = 30.0
anharmonic_shift_max = 0.10
min_admitted_sigmas = 3
min_q_meshes = 2
max_delta = 0.5
iso_margin = 0.75                 # FiresISO / FiresNM margin

[nm]
g_table = "docs/theory/certificates/SUMMARY.json"   # doc03 certified; lambda_cap = 100; g_inf = 2.6955e6
source = "doc01#0.3, doc03#2"

[aim]                             # doc01#0.5 RANK-AIM, ranking only
lambda = [2.0, 3.0]
A_K_min = 2060
m_K2_min = 6.0e6

[hard]                            # doc02#5a
H1_metallic = { gap_eV_max = 0.0, source = "doc02#5a:H1" }
H2_nonmagnetic = { moment_muB_tol = 0.05, source = "doc02#5a:H2", provisional = true }
H4_lambda_max = { value = 4.7, flag_above = 3.7, source = "doc02#5a:H4" }
H7_ehull = { window_eV = 0.35, route_required_above_eV = 0.10, zpe_for = ["H"], source = "doc02#5a:H7" }

[soft]                            # doc02#5b, weights used by rtt.score (ranking only)
S2_Hf_floor_H_channel = 0.5
S3_HDOS = [0.4, 0.6]
S6_space_groups = [191, 194, 139, 123, 166, 221, 223, 225, 227, 229]
S10_palette = ["H","Li","Na","Mg","Ca","Sr","Ba","Al","Sc","Y","La","Rh","Pd","Ir","Pt","Os","Au","Ag"]

[chemsys]
list = "targets/chemsys_list.txt"           # enumerated X-Y-H and 5(c) systems; source doc02#5d
```

### 7.2 `tiers.toml` (owned by docs 10, 11 and 13)

```toml
schema = "rtt.tiers/1"
[S1]
n_atoms_max = { value = 20, source = "doc10#1 (MatterGen training <=20 atoms)" }
pseudo_dir = "~/pseudo/sg15"                # whitelist = elements with a file here
[S3]
ml_model = "MatterSim-v1.0.0-5M"
ehull_self_max_eV = { value = 0.42, source = "doc02#5a:H7 + doc10#3.3 (0.35 + 0.07)" }
rank_cut_per_day = { value = 25, source = "this doc #8.3 (S4 capacity)" }
ml_phonons_remove = false                   # doc10#3.4; enforced false for H/Li/Be in code
alignn_weight_H_light = 0.0                 # doc11#2.3
[S4]
profile = "s4"
ehull_heur_max_eV = { value = 0.42, source = "doc02#5a:H7", provisional = true }
[S5]
profile = "coarse"
q = [4, 4, 4]                               # q3^3 allowed for C3 class
plateau_sigma_min_Ry = 0.03                 # doc11#3.1
margin_M_c = { value = 2.0, source = "doc11#4", provisional = true }   # REPLACE with <doc13>
[S6]
profile = "full"
q_meshes = [[4,4,4],[6,6,6]]
force_tol_mRy_bohr = 1.0
stress_tol_kbar = 0.5
certify = { mode = "bracket", mu_hi = 0.13, mu_hi_report = [0.10, 0.13, 0.15], omega_c_factor = 10.0 }
[audit]
A34_per_week = 10
A45_per_week = 2
A56_per_month = 1
[controls]
recall_disable_below = 1.0
```

### 7.3 `machine.toml`

```toml
[lanes.G]  gpu = true;  cpu_threads = 2
[lanes.S]  mpi_ranks = 4; jobs = 2          # two 2-rank jobs (doc 11 coarse runs used 2 ranks)
[lanes.F]  mpi_ranks = 4; jobs = 1
[env]      OMP_NUM_THREADS = 1; TORCH_DISABLE_NATIVE_JIT = 1; nice = 19
[disk]     scratch = "~/rtt_scratch"; min_free_GB = 80; prune_per_q = true
[derate]   contention = 0.85
```

---

## 8. Compute budget

### 8.1 Cost classes and per-stage cost

Cost law, from the benchmark cost notes and doc 11 §3.1: ph.x ∝ N_irr(q) · 3N_at · N_k(q) · N_PW log N_PW ·
N_bands. Anchors (meas.):
- H3S, 4 atoms, O_h: coarse q 3³ 271 s and q 4³ 635 s on 2 ranks; full ph.x 7 100 s on 8 ranks.
- pw.x steps < 1 min for 1–4-atom cells.

Scaling factors (est.):
- **C2** (6–10 atoms, cubic; Mg₂IrH₆ primitive cell has 9 atoms): modes ×2–2.5, N_PW ×3–6
  (volume), bands ×3–4, N_irr(q) and k per q ×1–2 → **×20–100**.
- **C3** (6–10 atoms, ≤ 16 point operations): additionally ×3–6 (q) and ×3–10 (k) → **×150–2000**
  (doc 11; R12 #6).

| Stage | C1 (≤ 5 atoms, cubic/hex) | C2 (6–10 atoms, cubic) | C3 (6–10 atoms, moderate symmetry) | Basis |
|---|---|---|---|---|
| S1 MatterGen | 10 GPU-s | 10 GPU-s | 10 GPU-s | meas. (doc 10) |
| S3 relax + hulls | 0.6–0.7 GPU-s | same | same | meas. |
| S3 phonons | 2–6 GPU-s | 2–6 | 5–15 (est.) | meas. 120–243-atom supercells |
| S4 (relax + scf + DOS + spin + φ) | 0.05–0.2 core-h | 0.5–2 | 1–3 | MgB2/H3S pw.x anchors × N^{2–3} (est.) |
| S5 coarse q 4³ | 0.35 core-h (meas. H3S) | 7–35 core-h | 50–700 core-h (q 3³: ÷2–3) | ×20–100 / ×150–2000 (est.) |
| S6 full, two q meshes (≈ 1.3× the finest) | 20–80 core-h (H3S 16 for one mesh) | 400–2 000 | 3 000–30 000 | same factors (est.) |
| S7a own-DFT reference hull | 1–10 core-h per chemical system, cached | 5–50 | 5–50 | est.; β-B (~105 atoms) excluded, α-B used (R12 #8) |
| S7b anharmonic check | MLFF-MD/TDEP: GPU-hours; SSCHA-DFT: 10²–10⁴ core-h (UNVERIFIED) | | | see §8.4 |
| certify bracket | ≤ 0.3 core-h | ≤ 0.3 | ≤ 0.3 | meas. 16–1 075 s |

### 8.2 Lanes

8 physical cores give 192 core-h/day, or **≈ 163 effective** after the 0.85 contention derate.
- Lane S: 4 cores, ≈ 82 effective core-h/day.
- Lane F: 4 cores, ≈ 82 effective core-h/day.
- Lane G: the GPU, ≈ 24 GPU-h/day, plus 2 CPU threads.

GPU split: about 15 h of MatterGen (≈ 5.4k conditional structures) and about 6 h of MatterSim
(≈ 5k structures with phonons for survivors). Adapter training is about 3 GPU-h when needed.

### 8.3 Daily funnel (survival fractions are assumptions to be measured; `summary.py` reports actuals)

| Stage | Capacity per day | Assumed pass | Flow per day |
|---|---|---|---|
| S1 + S2 | 5.4k (GPU) + lane D (one-off backlog) + controls | 0.8 unique and novel (doc 10: 58/64, 46/64) | 5 400 → ~4 300 |
| S3 | ≥ 5k (GPU) | 0.5 inside the ML window | ~2 150 → **top 25 by score** to S4; the rest are `DEFERRED(rank_cut)` |
| S4 | lane S, 40 core-h: C1 ~200, C2 ~30, C3 ~20 | 0.3 (H1, H2, heuristic hull) | 25 → ~8 |
| S5 | lane S, 42 core-h: C1 ~100, C2 1–6, C3 0.06–0.8 | M_c-dependent ⟨doc13⟩; assume 0.2 | C1/C2 8 → ~2; C3 is the bottleneck |
| S6 | lane F, 82 core-h: C1 1–4/day; C2 0.04–0.2/day (≈ 0.3–1.4/week); C3 ≤ 0.2/week | – | binding |
| S7b | lane F share or GPU | – | binding for hydrides (§8.4) |

What follows from the funnel:
- **C1:** about 7–25 full DFPT runs per week. Few target candidates are this small, apart from
  controls and simple hydrides.
- **C2:** about 1 per week (0.3–1.4).
- **C3:** fewer than 1 per month at full DFPT. Coarse DFPT alone (S5) is already 0.06–0.8 per
  day for C3, so C3 candidates should mostly be bounded by S5 evidence plus decide outcome
  `undecided`. S1 therefore biases toward C1/C2 through the G2 space-group lane, and the
  scheduler prices C3 candidates by expected cost.
- **First weeks:** lane F is reserved for the control lane (Mg₂IrH₆ S5+S6 ≈ 400–2 000 core-h,
  i.e. 5–25 days, est.; the 6–8-atom anchor control) and calibration. No throughput claim is
  made before the Mg₂IrH₆ and anchor costs are measured (R12 #6).
- **Disk:** at most one S6 job in flight (≤ 60 GB with per-q pruning; C2 UNVERIFIED) plus two
  S5 jobs (≤ 10 GB each).

### 8.4 The anharmonic bottleneck (consequence of G-VALID item 2)

For H/Li/Be compounds, `decide_t` returns `undecided` until `anharmonicCheck = true`. Options:

1. **MLFF-MD at 300 K plus TDEP-style effective phonons with MatterSim** (GPU, hours). It
   counts toward G-VALID **only if** MatterSim's harmonic phonons agree with the S6 DFPT phonons
   for *this* compound, within 10 % on H-dominated modes and with no spurious imaginary modes.
   Mg₂IrH₆ shows this can fail (doc 10 §3.4).
2. **SSCHA with DFT forces** (10²–10⁴ core-h, UNVERIFIED), or with a MatterSim model fine-tuned
   on the compound's DFT forces (fine-tuning MatterSim is untested here).

The budget allows option 2 only for the top 1–3 finalists. Option 1 is the default when its
admissibility test passes. Until either has run, the honest result for a hydride is "harmonic
model: core verdict X; decide outcome `undecided` (G-VALID 2)".

---

## 9. Response to review R12

| R12 finding | Response | Where |
|---|---|---|
| **#1 critical**: the input-quality downgrade layer lives in `pipeline.py` | **Accepted.** `certify_spectra` and all helpers core byte-identical as `certify_policy.py` with their tests; contract R1b; the overrides are forced off; `blocking_tags` also feeds `decide_t`'s `BlockedRej` | §1.1, §3.2, §4, §5.2 |
| **#2 critical**: R-LAMW unsound or redundant | **Accepted.** R-LAMW and the classify extension are deleted. L-MOM1–3 are documentation lemmas in `TargetLogic.Moments` with explicit hypotheses (ω_max as a bound, λ > 0) and are not used as screens | §0, §3.1 notes, §1.2 |
| **#3 major**: "@coarse" used as a proof-dressed cut | **Accepted.** S5 removal is `SCREENED_OUT(coarse-A)` with M_c (⟨doc13⟩, provisional 2); coarse verdicts are evidence only (`spectrum_tier='coarse'`); M_c is never passed as `rel_envelope`; the calibration scope is stated; `decide_coarse` in Lean already forbids rejection | §3.1, §3.4, §6.2 |
| **#4 major**: nothing binds the tier to artifacts | **Accepted.** `verdicts` has `spectrum_tier NOT NULL` and a `tier_run` foreign key; `verdict_spectra` stores digests with source artifacts; `recheck_ok` is required; binding triggers; `decisions` are tied to `decide_inputs` provenance; counts are split by tier | §4 R4, §6.2 |
| **#5 major**: audit lane infeasible and mis-specified | **Accepted.** Fixed counts per week; continuous error measures with inverse-probability weights; a surrogate floor and positive-control lane (Mg₂IrH₆, H3S, MgB2, Nb/Pb/Al, a 6–8-atom anchor); "PASS but not selected" = `DEFERRED(rank_cut)` | §3.3, §3.5, R8 |
| **#6 major**: optimistic costs | **Accepted.** Explicit scaling per class C1/C2/C3; full DFPT ≈ 7–25/week (C1), 0.3–1.4/week (C2), < 1/month (C3); 0.85 contention derate; lane F reserved for controls and calibration first; no throughput claim before the anchors are measured | §8 |
| **#7 major**: L-MOM1 screen weak | **Accepted.** The L-MOM1 screen is removed; ω₂-based λ_needed is a ranking feature only; the decisive quantity A is screened at S5 | §3.1 |
| **#8 major**: DFT hull without references or budget | **Accepted.** S4 hull is heuristic (`SCREENED_OUT` only); own-DFT reference hull + ZPE in S7a, cached per chemical system with pseudo md5/profile, costed; X-HULL only there | §3.1, §8.1 |
| **#9 major**: unsupported conditioning; no feedback | **Accepted.** Lanes G1–G5 use only shipped checkpoint/property sets (doc 02 5(d)); planned A-adapter GA with acceptance test; weekly chemsys re-weighting from S4/S5 labels; database lane D | §3.3 |
| **#10 major**: ranking and evidence underspecified | **Accepted.** Ranking by decide outcome, certified T_low, A_mid; T_up as a ceiling only; evidence template covers full-grid stability, ZPE hull, anharmonic check, MD persistence, μ\* sensitivity, synthesis route | §3.6 |
| **#11 major**: ordering; Γ tier weak | **Accepted.** ML proxies moved to S3 before any DFT; S3 internal order relax → hull → features → phonons; Γ e-ph tier removed (Γ phonon only as a diagnostic) | §3.1 |
| **#12 minor**: LOC errors; T3 device | **Accepted.** 13 984 / 7 791; DiscoveryLogic 1 406. T3 costs now come from doc 10's GPU measurements | §1 |
| **#13 minor**: Lean details | Done in TargetLogic by the Lean agent (ω_max as a hypothesis; λ > 0 for L-MOM2/3); the "and nothing more" wording is adopted (λ·ω_log ≥ 2·floor implies only A ≥ floor) | §1.2 |
| **#14 minor**: shared relax | **Accepted.** S6 re-relaxes or checks force and stress at S6 settings | §3.1 |
| **#15 minor**: calibration scope | **Accepted.** MgB2 blocked; a 6–8-atom anchor control added; M_c re-checked on every S6 run | §3.3, §3.4 |

---

## 10. Open problems and what I could not verify

- **⟨doc13⟩:** the coarse margin M_c and its transferability to C2/C3. Doc 11's value of 2
  rests on H3S alone.
- **Anharmonic-check admissibility** of MatterSim MD for hydrides (§8.4). This needs a policy
  decision by the theory owner (doc 01 G-VALID lists "MLFF-MD" without an accuracy condition;
  I added one).
- **C2/C3 cost factors** are scaling estimates (×3 either way). Mg₂IrH₆ S5/S6 is the first
  measurement.
- The **`decideB` theorem and truth-table export** (§1.2) are not yet in `lean/TargetLogic`.
- **ALIGNN numbers** are recorded at weight 0 for H and light elements. Whether they help
  elsewhere is untested.
- **MgB2 ph.x crash** cause; **QE 7.4/7.5 or patched 7.6** for ultrasoft/PAW: both open.
- **A-adapter labels**: source, licence and the A definition used (QE bin rule vs ours) depend
  on the label agent's deliverable. The property id `rtt_A_K` is a proposal.
- **H2 moment tolerance** of 0.05 μ_B per cell is my placeholder (marked provisional).
