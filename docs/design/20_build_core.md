# 20. Build notes: core modules (builder B1-core)

Status: built 2026-10-28/29 against `docs/design/12_architecture.md` rev 2, doc 01 §0.4/§0.6
(post-RL1), doc 03 certificates and `lean/TargetLogic/Decide.lean` **v2** (the RL1 revision,
which landed while this was being built). Python 3.12, `~/.rtt/envs/py`.

## 1. What is implemented

| Module | Content |
|---|---|
| `pyproject.toml` | package `rtt` (src layout), script `rtt = rtt.cli:main`, pytest config (`pythonpath = src`) |
| `rtt/config.py` | `load_targets/tiers/machine`, `sha256_of`, `Threshold(value, source_id, provisional)`, `ConfigError`. **R7:** every numeric leaf in targets/tiers needs a source (inline `{value, source}` or a section-level `source`), else the loader refuses. `Targets.floor = a0·T*` exact; `with_T_star` gives the control-lane surrogate target (own sha256). |
| `rtt/records.py` | §6.1 dataclasses (`StructureData`, `Origin`, `RawStructure`, `Context`, `ArtifactRef`, `CostEstimate`, `Outcome`, `StageResult`, `StageEntry`, `CandidateRecord` = `rtt.candidate/2`), JSON round trip; the §3.1 stage→outcome table (`ALLOWED_OUTCOMES`); the Lean v2 decide types `Prov, IsoData, EpsData, MatData, DosData, Validity, InputTag, DecideInputs, DecideOutcome`. |
| `rtt/units.py` | re-export of core `dfpt/units` |
| `rtt/targets/` | `gtable.py` (certified g lookup from `docs/theory/certificates/nm_*.json` (exact `m0_over_T2`; SUMMARY.json also accepted): `m0·T*²`, "up" rule, monotonicity check on load), `spec.py` (`TargetSpec`) |
| `rtt/ledger/schema.sql` | all §6.2 tables + triggers (list in §3) |
| `rtt/ledger/store.py`, `artifacts.py` | `Ledger` (WAL, append-only), insert helpers, `current_outcome`, `candidate_record` (assembled from the ledger), `export_jsonl` (+ `manifest.sha256`, re-verifies artifact hashes), `summary_counts` (REJECTED_PROVEN split by spectrum tier); artifact hashing, verify/prune *events* (no row is ever changed) |
| `rtt/core/` | `certify/`, `spectrum.py`, `status.py`, `constants.py`, `validity/`, `dfpt/`, `discovery/fingerprint.py`, `sources/{base,composition,dedupe,fileformats,local,synthetic}.py`, and `certify_policy.py`; `PINS.json` (sha256 per file, checked by `tests/test_core_pins.py`) |
| `rtt/certify/` | `run.certify` (the only call into `certify_policy.certify_spectra`; refuses and forces off both `allow_unverified_*`), `alternatives.freq_scaled` (Hopfield copies s ∈ {1−ε_ω, 1+ε_ω}), `envelope_e`, `recheck.recheck_or_raise` (binds run parameters) |
| `rtt/decide/` | `core.py`: exact-rational mirror of Lean `decideT` v2 (`wf`, `gvalid`, `floorR`, `fullSameRun`/`provOK`, `firstRule`, `wouldFire`, `blocksRejection`, `applyPolicy`), `core_b`, verbatim Lean texts for the drift test; `inputs.py`: record → `DecideInputs` mappers (signs → magnitudes, ladder → envelope, band matrix → `MatData` with nominal Perron test vector, spectra tags → `InputTag` with the core `blocking_tags` result); `lean_io.py`: reader for `lean/export/decide_cases.json` |
| `rtt/score/s3.py`, `rank.py` | T1a-style S3 ranking score (ranking only, ALIGNN unused); §3.6 lexicographic final ranking + `rows_from_ledger` |
| `rtt/engine/` | `stage.py` (§5.3 `Stage` protocol, `Context` with surrogate targets, `StageBase`); `stages.py` (S1–S3 delegate to B2's batch runners, S4–S6 per candidate on B3's `rtt.dft`, S7 re-decide, S8 stub, S9 evidence via B3's `report.evidence`; the soundness glue `store_verdict` → `record_dfpt_decision` → `combine`); `scheduler.py` (derived queue, lane budgets from machine.toml × derate, §8.1 cost by class, controls first); `runner.py` |
| `rtt/cli.py` | `rtt run|step|status|redecide|evidence|export|controls|rank|decide-check` |
| `configs/` | `targets.toml`, `tiers.toml`, `machine.toml` (§7, with source ids; see §2) |

## 2. Deviations from the architecture (rev 2) and why

1. **decide layer follows Lean `Decide.lean` v2, not §3.1/§6.1 of rev 2.** After review RL1 the
   Lean spec changed: per-input `Prov(tier, run)` for the groups A, λ, m, ε, matrix, DOS (no
   single `tier`); outcomes `rejectedProven r | screenedOut | undecided | pass` — **no
   `rejectedAssuming`**; NM-G (now [THM-C], doc 03) gives `rejectedProven` recording `[M-ISO, CERT-NM]`; `screenedOut` = "a floor holds but provenance is not proof grade" and can
   occur on full-DFPT records. Consequences here:
   - `DecideInputs` JSON uses the Lean field names (`val, iso, eps, mat, dos, tags`), not
     §6.1's `validity/envelope/mb`. The stage `tier` is kept as ledger metadata only.
   - Lean has no `converged` tier; EPW-derived inputs are mapped to `Prov("fullDFPT")`.
   - Ledger: `decisions_guard_screen` (screenedOut only at non-proof tiers) is **dropped**; a
     `rejectedAssuming` decision is **refused**; `rejectedProven` needs one of the six rules and
     a recorded assumption list. At S6–S8 a decide `screenedOut` becomes stage outcome
     `undecided` (reason `WOULD_FIRE_LOWER_TIER`), since the §3.1 table has no removal there.
   - Python-only addition: a `None` in a field Lean requires (anything but `E_F`) gives
     `undecided` with rule `MISSING` (Lean cannot represent it; §6.1 "missing → undecided").
2. **The 2026-10-29 doc 01 additions are switches** in `targets.toml [decide]`, all **on** (as in
   Lean and the exported cases): `w_lo_sum_check` (Σ w_i(1−δ_w) ≤ 1 in dosWF),
   `moment_consistency` (A² ≤ λm/4·(1 + 1/100) in WELL-FORMED), `nmg_certified` (NM-G records
   `[M-ISO, CERT-NM]`; off = the older NUM-* list). The g table is read exactly from
   `docs/theory/certificates/nm_*.json` (`m0_over_T2`, verdict CERTIFIED, as
   `lean/scripts/gen_nm_table.py`); no λ cap (beyond the last λ_k the λ = ∞ certificate).
3. **tier_runs is not strictly insert-only**: one UPDATE `started → done|failed` filling
   finished/wall/core/gpu/metrics is allowed; every other column is frozen by trigger. (§6.2
   needs both "written before launch" and `state = 'done'` at verdict insert.)
4. **verdicts carry `spectra_json`** (alt index, digest, source artifact id, sha256). The
   `verdicts_bind` trigger checks it with JSON1 *at insert* (every alternative must point to an
   artifact of a done dfpt/epw run of the same candidate with the same sha256, ≥ 1 alternative).
   `verdict_spectra` rows are still written and must match it. Reason: a row-level trigger on
   `verdicts` cannot see `verdict_spectra` rows inserted afterwards.
5. **Removal decisions need a rechecked core verdict for every rule** (not only N1-ISO):
   a full/converged verdict with `recheck_ok = 1` on a tier_run listed in the decide input's
   provenance. Stricter than §6.2; it is what R4 requires for "a removal decide outcome".
6. Extra columns/tables: `outcomes.reason, decision_id, tags`; `stage_outcomes` (the §3.1 table,
   immutable, used by `outcomes_allowed`); `artifact_events` (prune/verify, since artifacts are
   insert-only); `screen_status` (R8 auto-disable log). `ERROR` is allowed at S7 and S8.
7. **The core keeps the package layout** (`core/certify/…`, `sources/dedupe.py`, …)
   instead of the flat names in §5.1, so relative imports work unchanged. The root, `discovery`
   and `sources` `__init__.py` are minimal stubs (they do not import the dropped modules).
   `certify_policy.py` = import header + `base/pipeline.py` **L83–133 and L420–650**.
8. `machine.toml` in §7.3 is not valid TOML (`;` separators); written as valid TOML with a
   `[capacity]` table.
9. `score_s3` accepts both the §5.1 signature `(features, phonons, targets)` and B2's call
   `(s3_metrics, targets)`.

## 3. Ledger triggers (all enforced by SQLite)

Append-only (UPDATE and DELETE abort) on every table except `queue`; `tier_runs` as in §2.3.
`outcomes_allowed` (stage, outcome) ∈ §3.1; `outcomes_inputs_known` (targets/tiers sha256
registered, R7); `outcomes_criterion` (SCREENED_OUT/EXCLUDED_ASSUMPTION need criterion, value,
threshold, source_id; DEFERRED/ERROR need a reason); `outcomes_excluded_criteria` (R3: S4 only
H1/H2, S7 only H6/H7); `outcomes_s5_screen` (S5 removes only by `coarse-A`);
`outcomes_decide_bound` (a decide-vocabulary outcome must reference a decisions row of the same
candidate and stage, equal or downgraded to undecided); `outcomes_rejection_bind` (rejection
needs a done full-DFPT/EPW tier_run of the candidate). `verdicts_bind_*` (origin = core
`certify_spectra`; proven ⇒ `recheck_ok = 1`; done dfpt/epw run of the same candidate with the
profile ↔ spectrum_tier map; certificate artifact of the same run; hashed source artifacts);
`verdict_spectra_bind`. `decisions_guard_*` (origin `rtt.decide.core.decide_t`; decide input of
the same candidate and tier; provenance tier_runs of this candidate and, for proof tiers, a done
run with matching profile; rejection only at fullDFPT/converged; rechecked verdict; N1-ISO needs
a REJECTED_PROVEN crosscheck verdict; rule + assumptions; no rejectedAssuming).

## 4. Test results (2026-10-29, `scripts/pyt.sh`: nice 19, 2 threads, run in parts)

| File | Result | What it checks |
|---|---|---|
| `tests/test_core_pins.py` | 50 passed | sha256 of every core file and of both `certify_policy.py` blocks against `PINS.json` |
| `tests/test_decide_core.py` | 92 passed | every rule (N1-KAPPA, N1-ROWMAX, N1-DOS, N1-ISO, NM-CERT, NM-G) incl. margin clause and k-resolved; WELL-FORMED and G-VALID items; RL1 counterexamples (negative A, signed imaginary modes / shifts, Îµ = 0, mixed provenance, mixed runs) â never rejectedProven; R01 two-band counterexample â undecided / pass, never rejectedProven; low-tier inputs never reject (4 000 random cases); rejection â Lean preconditions (6 000 cases); tag addition downgrade-only (3 000 cases); envelope growth never creates a rejection; the new DOS-sum and moment gates and their switches; g table equals the certificates |
| `tests/test_decide_lean_table.py` | 24 passed, 1 skipped | `core_b` vs an independent transcription of `decideCore`+`applyPolicy`; drift check of 11 verbatim Lean texts; (skipped: `lean/logs/decideB_table.txt` not exported) |
| `tests/test_decide_lean_cases.py` | 565 passed | all **563** cases of `lean/export/decide_cases.json` (outcome, rule, assumptions, inputs-used provenance, gates, per-rule floor/provOK, exact A_hi, Î»_hi, m_hi, floor, NM-CERT floor, g_up, Îº_hi, rowmax, w_lo); export config = rtt config |
| `tests/test_ledger_triggers.py` | 117 passed | append-only on every table; tier_run life cycle; every forbidden (stage, outcome) pair; R2/R3/R7 fields; verdict without DFPT tier_run, unfinished run, profile mismatch, other candidate, foreign certificate, missing recheck, wrong origin, wrong/missing artifact hash; decisions: rejection at coarse/ml tier, wrong origin, tier/provenance mismatch, no rechecked verdict, N1-ISO without REJECTED_PROVEN crosscheck, no rule/assumptions, rejectedAssuming; outcomes not bound to a decision or a full run |
| `tests/test_engine.py` | 14 passed | combination rule; full S6 ledger chain on the Al spectra (verdict â recheck â decision â `rejectedProven`), CROSSCHECK_DISAGREE, tampered hash refused; scheduler; S3 score; Â§3.6 ranking; R7 loader; import-graph (R1b/R1c); `rtt decide-check` |
| `tests/test_certify_e2e.py` | 8 passed | Al and Pb at T* = 300 K: core REJECTED_PROVEN (rechecked) and decide_t rejectedProven by N1-ISO (also with frequency-scaled alternatives); **Al CERTIFIED [0.80282, 10.302] K at T* = 0.5875 K and Pb CERTIFIED [4.5151, 9.4237] K at T* = 3.595 K, Î¼*_hi = 0.13 â equal to the old validation runs**; overrides refused |
| core `test_certify_types, test_dfpt_smearing_admission, test_dfpt_soundness, test_dfpt_parsers` | 45 passed | |
| core `test_certify` | 19 passed | |
| core `test_certify_multiband` | 8 passed | |
| core `test_certify_multiband_review` | 19 passed, 9 skipped | skips need the wMB reference `mb.py`, which is not included |
| core `test_correspondence` | 51 passed | drift-checked against the core `lean/PipelineLogic` |
| core `audit_certify/test_forgery` | 13 passed | |
| core `oracle/test_oracle_vs_certify` | 132 passed | |

Total: **1 157 passed, 10 skipped, 0 failed** (870 new + 287 core).

## 5. Open items / not verified

- **S6 `min_freq_K` is `None`**: the full-grid imaginary-mode check (dyn files of the final
  structure) is not parsed yet, so every real S6 decision is currently `undecided` (MISSING)
  until B3 supplies the value. This is the safe direction, but it blocks rejections.
- S4 → `n_sheets`, `E_F_eV` come from B3's `ElectronicResult.metrics()`; S7 compute
  (own-DFT hull + ZPE, anharmonic check, band-resolved coupling) and S8 (EPW) are not wired;
  S7 only re-decides from records. `rtt redecide --stage S5` is a stub.
- The S6 glue (`S6Full.run` → `certify_and_decide`) is tested only up to the ledger chain on
  the Al fixture (`tests/test_engine.py`); a real two-mesh QE run through it has not been done.
- The run id of all S6 input groups is the finest-mesh tier_run id (one run per Lean
  `fullSameRun`); the coarser mesh enters only through δ^q.
- `decideB_table.txt` is not exported by Lean; `core_b` is checked against an independent
  transcription and, case by case, against the 563 exported Lean cases.
- The scheduler's costs are the §8.1 upper estimates; no spent-budget accounting from
  `tier_runs` yet (callers pass `spent_core_h`).

## 6. Fixes for code review RC1 (2026-10-29)

| Finding | Fix (B1) | Test |
|---|---|---|
| F1 H1/H2 silently off | `targets.toml` `hard.H1_metallic` / `H2_nonmagnetic` carry `value =` (same number as the named key); `S4Electronic.decide` passes `targets.threshold("hard.H1_metallic")` etc. (sourced `Threshold`s) | review F1; `test_engine::test_s4_decide_uses_sourced_thresholds` |
| F12 DEFERRED permanent | `rtt.engine.retry`: DEFERRED is never terminal. Queue states (`rank_cut`, `budget`, `awaiting_dependency`) re-enter without using the budget. Failure deferrals and ERROR re-enter up to `[retry] max_attempts = 3` (tiers.toml, sourced, engineering policy) with back-off `backoff_h·2^(n−1)` (6 h). `run_s3` (one-line change in B2's `mlff/screen.py`) and `S3MLFF.run_batch` select by `needs_run` | review F12; `test_rank_cut_deferral_is_not_a_failure` |
| F13 ERROR terminal | `scheduler.next_stage` re-offers a stage whose current outcome is retryable within budget; `next_jobs` applies the back-off; `runner` writes `Outcome.retry`; `retry.exhausted()` lists used-up candidates for S9 | review F13; `test_retry_error_and_backoff`, `test_terminal_removal_stops` |
| F14 recall overstated | `report.summary.control_recall`: only PASS / pass / undecided count as passed; ERROR and DEFERRED are counted separately (`error`, `deferred` keys, present when non-zero) and are recall failures | review F14 |
| F15 H7-ml margin | `tiers.toml S3.ehull_self_max_eV = 0.56` (0.35 + 2√2·0.074), **provisional**, to be replaced by the A3→4 audit 95th percentile | review F15 |
| N6 controls | Surrogate floor: the runner passes `ctx.targets_for(rec)` to every `decide` (S4, S5) and `record_dfpt_decision` uses it (S6/S7). Dedupe: controls are no longer dedupe representatives (two-line change in B2's `structure/stage.py`), so a generated copy of a control passes S2 and is evaluated at the real T* | `test_control_surrogate_target_reaches_decisions` |

Not B1's: F2–F11 and F17 (B3's S4/S5/collect code) and F16 (B2's dedupe). After these fixes the
review suite has 18 of 19 passing; the remaining F7 (provisional H2 tolerance recorded on PASS)
is in B3's `decide_s4`.

## 7. Wiring of B3's RC1 back-end fixes (2026-10-29)

- S4: `S4Electronic.run` registers `rtt.dft.electronic.electronic_artifacts(res, ledger, tier_run_id)`
  before finishing the tier_run (RC1 N1: S4 removals are re-verifiable from hashed files).
- S5: runs `tiers.toml S5.q_meshes = [[3,3,3],[4,4,4]]` (doc 13 T2-SCREEN), one tier_run per mesh,
  the finest carries the metrics (+ `by_mesh`, + S4 `n_sheets` / `sheet_w_min`). `S5Coarse.decide`
  calls `decide_coarse(finest, metrics_coarser=coarser, e_c=S5.e_c, margin_M_c=S5.margin_M_c,
  n_sheets, w_min)`. Scheduler S5 cost = q4³ cost × 1.5 (q3³ is ÷2–3 of q4³, §8.1).
- S6: `min_freq_K` = min over both meshes of B3's `CoarseMetrics.min_freq_K` (dyn files + matdyn of
  the final structure; None if any mesh lacks it); `max_imag_mode_K = max(0, −min_freq_K)` is
  mapped to `Validity.imagAbsMaxK` (absolute) by `rtt.decide.inputs.build_inputs`. Real S6
  decisions can now leave `undecided` (the earlier open item of §5 is closed in code, untested on QE).
- `scripts/run_all_tests.sh` runs the whole suite in parts. Result (2026-10-29): **1 268 passed,
  14 skipped, 1 failed** — the failure is `tests/review` F11 (B3's `decide_coarse` single-mesh
  fallback with ε_sys = 0.10 returns DEFERRED instead of the PASS the review test expects).
