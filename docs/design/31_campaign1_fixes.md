# 31. Campaign 1 fixes (doc 30 §7 open defects, review R30 F1/F3/F5/F8)

Agent C2-fixes, 2026-10-30. Code changes in `src/rtt/`, each with tests. Labels: **(meas.)** measured,
**(est.)** estimated, **UNVERIFIED** not checked. No QE, GPU or MatterGen runs were made for this
document; all checks are unit tests on fakes or on the campaign-1 numbers quoted in doc 30.

| § | Defect | Fix (files) | Tests |
|---|---|---|---|
| 1 | B9: S3 rank cut applied per call | one deterministic cut per pool, recorded in table `rank_pools` (`mlff/screen.py`, `ledger/schema.sql`, `ledger/store.py`) | `tests/front/test_campaign1_fixes.py` (4 tests) |
| 2 | `score_s3` has no metallicity term; hand quota `CAMPAIGN_GENERATED_QUOTA` | composition metallicity prior incl. valence parity (`score/metallicity.py`, `score/s3.py`, `mlff/features.py`); `[S3.rank_lane_min_share]` in `tiers.toml` | same file (4 tests) |
| 3 | B11: chemsys list missing, silent fallback | list generated from doc 02 (`configs/targets/chemsys_list.txt`, `source/lanes.py`); fallback logs a warning and is recorded in every lane job's source id | same file (2 tests) |
| 4 | one failed admission mesh turns S5 into ERROR | `dft/admission.py`; wired into `S5Coarse.decide` | `tests/back/test_campaign1_fixes_back.py` (4 tests) |
| 5 | vc-relax hitting `nstep` not retried | one recorded continuation in `dft/workflow.py::Chain` | same file (3 tests) |
| 6 | scheduler S5 cost ~10x low for 5d hydrides | refit from doc 30 §5.3 (`engine/scheduler.py`) | `tests/test_campaign1_fixes_engine.py` (2 tests) |
| 7 | R30 F1(a), F1(b), F3, F5 | S3 window domain, S1 generator limits, S2 representative order, `rtt step` gate (`mlff/screen.py`, `source/prefilter.py`, `records.py`, `structure/stage.py`, `cli.py`, `engine/retry.py`) | front + engine files (6 tests) |

## 1. B9: pool-wide rank cut

- `rank_pool(items, cut, …)` ranks the non-control S3 PASSes of a pool by (score desc, candidate id
  asc) and applies the cut once. The result does not depend on input order (test). Controls stay exempt.
- `run_s3(…, chunk_size=N)` computes in chunks. Each chunk's survivors are written as
  `DEFERRED(rank_cut)` queue entries, with reason `rank_pending(pool=…)` and tag `RANK_PENDING`. The pool is
  closed once at the end. `run_s3(…, close_pool=False)` only queues; `close_rank_pool(…)` (= `rerank_s3`)
  then closes the day's or the campaign's pool once. A single call without chunking behaves as before:
  the call is the pool.
- Every close appends a row to the new insert-only table `rank_pools`. The row holds the pool id, the
  cut, the free slots, the S4 waiting count, the lane shares and their source, the score version, the
  full ranking (cid, lane, score, rank, selected, via) and its sha256, and the tiers sha256. Each
  written outcome carries the tag `RANK_POOL:<id>`. An older ledger gets the table automatically
  (`CREATE TABLE IF NOT EXISTS`).
- Queued scores from an older score version are recomputed from the stored metrics at close, so one
  pool is ranked by one score function.
- The engine's `S3MLFF.run_batch` passes `**kw` through, so `chunk_size` / `close_pool` / `pool_id` work
  there unchanged. `scripts/campaign_front.py` (rank_cut 0 per chunk, then one `rerank_s3`) is
  equivalent and was left as the campaign-1 record.

## 2. Metallicity prior and per-lane minimum share

**No offline ML band-gap model exists.**
- The MatterGen `dft_band_gap` checkpoint is not local (doc 30 §1).
- ALIGNN is not installed and has weight 0 for H and light elements (doc 11 §2.3).
- The prior is therefore a composition rule with no fitted parameters (`rtt.score.metallicity`,
  version `metallicity-prior/2`), applied in order:

1. **Valence-electron parity** (R30 F1(c)/F8), on the primitive cell's counts.
   - Valence counts: s+d electrons for groups 1–12; group − 10 for the p block; 3 for the
     lanthanides (4f as core; La checked, others **UNVERIFIED**).
   - An odd count cannot fill spin-degenerate bands, so a non-magnetic band structure is metallic:
     p_metal = 0.95.
   - Parity is used only when the counts are those of a cell. A reduced formula gives no parity.
2. **Formal electron counting** for hydrides of donors and transition metals.
   - H is counted as H⁻. Donors have fixed charges (alkali +1, alkaline earth / Zn / Cd / Hg +2,
     Sc / Y / La / Ln / Al / Ga / In +3).
   - TM oxidation state: ox = (n_H − donor charge) / n_TM; d-count = group − ox.
   - Closed shell → 0.15. This covers a charge-balanced saline hydride, d⁰, d¹⁰, d + 2·(H per TM) = 18,
     and square-planar d⁸ with 4 H.
   - Open shell, fractional ox, or excess electrons → 0.85.
   - Negative ox (intermetallic-like) → 0.80.
   - All-metal compositions → 0.80.
   - Anything else (borides, carbides, mixed anions, H3S) → 0.5 (neutral).
3. The levels are ordinal choices, **not calibrated probabilities**.

**Score.** score = (1 − w_m)·base + w_m·p_metal − 0.2·[H₂ motif], where base is the old weighted sum.
- w_m = `tiers [S3] metallicity_weight` = 0.20 (provisional; a judgment).
- Spread between a closed-shell and an open-shell cell: 0.14 in score.
- `features()` records `composition` and `metallicity_prior` for every S3 candidate.
- The S3 decision does not read the prior. It only reorders the queue (contract R2; tested).

**Check on campaign 1 (in-sample, meas. on the doc 30 §5.1 table).** The rule was written from
textbook counting before this check, but the same table motivated the work, so this is not a
validation.
- All 8 S4 insulators get p < 0.5. These are SrPtH₆, BaPtH₂, Ba₂Pt₂H₈, NaAuH₂, BaPtH₆ and AlH₃.
- 17 or more of the 22 are on the correct side of 0.5 (`test_prior_on_campaign1_s4_outcomes`).
- The misses are LaH₃, YH₃ and Y₂H₆. These are charge-balanced, and PBE puts them at E_F with a
  tiny N(E_F) of 0.02–0.07; experimentally LaH₃ is reported semiconducting (**UNVERIFIED**).
- The other miss is Sr₂Pt₂H₈ (counted square-planar d⁸, S4 metal).

**Per-lane minimum share** replaces `CAMPAIGN_GENERATED_QUOTA`.
- Config: `tiers.toml [S3.rank_lane_min_share]`, `G = 0.40`, sourced, provisional.
- In each pool close, the best floor(0.40 × free slots) candidates of lanes G1–G5 / GA are reserved
  first. The rest go by global rank. Unused share slots fall back to the global order.
- A candidate that got its slot only through the share (global rank > cut) is tagged
  `LANE_MIN_SHARE:G` on its outcome; the `via` column of `rank_pools` records the same. The tagged
  outcome's `source_id` is the share's source.
- Campaign 1 promoted 10 of 25 slots by hand. 0.40 × 25 = 10 reproduces that as configuration.

## 3. B11: chemsys list from the materials profile

`src/rtt/configs/targets/chemsys_list.txt` (the `targets.toml` path resolves against the configs
directory) is generated by `rtt.source.lanes.write_chemsys_list`. It has a source-note header and
two sections:

- **`[family1]`**: 70 systems, donors {Li, Na, Mg, Ca, Sr, Ba, Al, Sc, Y, La} × late 4d/5d metals {Rh,
  Pd, Ir, Pt, Os, Au, Ag} × H. This is doc 02 S10 / priority family 1, the A₂MH₆ / AMH₃ complexes.
  There are no donor–donor–H systems.
- **`[nonhydride]`**: the 5(c) systems doc 02 names: B-C-Li, B-C-Mg, B-Be, B-Be-Li, B-C-N, B-C-Sr
  (SrB₃C₃), B-C-Cl (ClB₂C₈). G4 now reads this section unless `targets.toml chemsys.nonhydride` is set.
- **Not enumerated:** family 2, the light-metal H networks. That family is a structural motif, not a
  composition list. The Li/Be-rich intermetallics are not enumerated either, because doc 02 names no
  system for them.

`load_chemsys_list` returns the flag `CHEMSYS_FILE(sha256=…)`. If the file is missing or empty:
- it logs a WARNING;
- it returns `DEFAULT_CHEMSYS_FALLBACK(<why>)`;
- `plan_day` appends that flag to every G1/G4 job's `source_id`, which reaches each candidate's
  origin meta (`lane_source_id`).

R30 F8 asks for a hard failure instead. I kept a loud, recorded fallback because the in-code
fallback is now the same doc 02 list. If the orchestrator prefers a hard stop, it is a one-line
change in `load_chemsys_list`.

## 4. S5 admission failure on one mesh

The core rule (`MODE_NEGATIVE_TOL_REL = 1e-6`) is unchanged.

`rtt.dft.admission.admission_deferral(res)` returns `DEFERRED(admission_failed)` when both of these
hold:
- every coarse mesh either succeeded or failed **only** the core admission rule (error text
  `smearing admission: …`, flag `SMEARING_ADMISSION_FAIL`);
- at least one mesh failed it.

Crashes, missing or changed files, and mixed errors stay ERROR (the helper returns None).

The outcome records, per mesh:
- the status;
- the core max negative per-mode fraction and the tolerance;
- A_mid and A_lo of the admitted mesh.

It also carries the tags `ADMISSION_FAILED:q3x3x3`, `NO_REMOVAL_WITHOUT_ADMISSION` and
`SINGLE_MESH_ADMITTED`. The source is doc 13 §5: a removal needs admission on both meshes.

**Wiring.** S5-recal has finished, so I wired it myself as the first statement of
`S5Coarse.decide`:

```python
if (o := admission_deferral(res)) is not None: return o
```

`S5Coarse.run` now also puts `collected_by_mesh` in the payload, so the check uses the exact
`Collected.errors` rather than the joined error text.

**Consistency with doc 16's all-σ statistic.** A mesh with failed admission has no admitted spectra,
so no statistic can be formed on it. The deferral never removes, so it cannot conflict with the rule.

**Retry.** `DEFERRED(admission_failed)` remains retryable under the engine's default policy, which
would rerun an identical S5 up to 3 times. This is open (§8).

## 5. vc-relax `nstep` continuation

`Chain.run` makes one continuation when a pw.x relax step fails and its output contains
"maximum number of steps has been reached":
- `pw_continue_text` keeps the input and sets nstep × 2.
- It replaces `CELL_PARAMETERS` (converted to Å from the `alat=` header) and `ATOMIC_POSITIONS`
  with the last blocks QE printed.
- Physics settings are unchanged.

The marker `<step>.nstep_retry` (JSON) records:
- the reason, old and new nstep, and the geometry source;
- the first output's sha256 and a copy of that output (`relax.nstep1.out`);
- the exact new input text.

The marker is sticky, so a finished continuation is adopted on resume. At most one continuation is
made (tested); the cg/damped fallback of B6 can still follow. The S4 artifact list now includes
`relax.nstep_retry`, `relax.nstep1.out` and `relax.fallback`, so they are sha-registered in the
ledger.

**Not covered:** the core EPhWorkflow relax used by S5/S6, which is core code.

## 6. Scheduler S5 cost refit

Measured S5 (q 3³ + q 4³, 60 Ry; doc 30 §5.2–5.3), in core-h:

| cell | atoms | core-h |
|---|---|---|
| ScH₂ | 3 | 1.67 |
| LaH₂ | 3 | 1.35 |
| NaIrH₃ | 5 | 3.30 |
| BaPtH₃ | 5 | 6.47 |
| Mg₂IrH₆ | 9 | 29.7 |
| Al, Nb, Pb | 1 | 0.49–0.96 |

**Model** (`S5_COST_FIT`): core_h = max(1.0, a·N³) × (ecut/60)^1.5.
- a = 0.062, the maximum of core_h/N³ over the five hydrides (0.050, 0.062, 0.026, 0.052, 0.041;
  geometric mean 0.044). The old table also used upper ends.
- N³ fits all five within ×2.4 with no separate 5d factor (tested).
- ecut = max(coarse profile, `tiers [S5] ecut_floor_Ry`). H is at 80 Ry (doc 15/16), giving ×1.54.
- `next_jobs` uses the per-candidate estimate from the cell.
- The class table is the fallback when there is no structure:
  - C1 = 5 atoms at 80 Ry = 11.9 core-h (was 0.525);
  - C2 = 10 atoms = 95 core-h;
  - C3 is unchanged (low symmetry is outside the fit).
- Source: doc30#5.2-5.3; doc31#6.
- Only five hydride points support the fit; it is an upper envelope, not a law (**est.**).

## 7. R30 items

- **F1(b) S1 size limit.**
  - `n_atoms_max = 20` is MatterGen's training limit. For lanes D and C, S1 now gives
    `DEFERRED(cost_class)` with a reason and the tag `GENERATOR_DOMAIN_LIMIT_NOT_APPLIED`. The
    MatterGen element-domain rule is not applied to D/C either (PASS with a tag). Generated lanes
    are unchanged.
  - This needed `('S1','DEFERRED')` in the fixed §3.1 table: `records.ALLOWED_OUTCOMES` and
    `schema.sql` (guard count 36 → 37).
  - Older ledgers are migrated once on open by `Ledger._migrate_stage_outcomes`. It drops the
    insert guard, adds the pair, re-creates the guard from the schema text, and writes a meta row
    `migration:stage_outcomes:S1/DEFERRED` (tested). The 19 campaign-1 removals are not re-decided
    by this code. That needs an explicit re-run of S1 decisions (open).
- **F1(a) S3 window domain.**
  - The 0.60 eV/atom window is calibrated on H-containing entries only (doc 21).
  - For a candidate without H, or with unknown composition, an exceedance is now
    `DEFERRED(uncalibrated_window)`, provisional, with the tag `H7_ML_WINDOW_UNCALIBRATED`.
  - It is non-retryable (`engine/retry.NON_RETRYABLE_DEFERRED`): a rerun is deterministic, and the
    candidate waits for a non-hydride calibration.
  - Hydrides are unchanged. `test_decide_s3_hull_rule` was updated to give its hydride cases `n_H`.
- **F3 S2 representative.**
  - S2 now processes controls first, then lane D, then the generated lanes. Within a group the
    order is fewest atoms, then highest space-group number, then id.
  - A generated copy of a database entry therefore becomes the DUPLICATE, tagged
    `GENERATED_REDISCOVERY` / `REDISCOVERS:<D cid>`; generated stats should count those tags.
  - If D arrives in a later batch than its generated copy, the D row is tagged
    `KNOWN_ENTRY_OF_GENERATED`.
  - Novelty against Alex-MP stays in the PASS value as before.
  - Campaign 1 is not re-cut here.
- **F5 `rtt step` gate.**
  - `cmd_step` refuses a per-candidate stage whose `requires` did not pass (exit 2).
  - `--override-reason "<why>"` runs it and first writes a ledger meta row
    `gate_override:<stage>:<cid>:<n>` holding the unmet requirement, its current outcome, the
    reason and the time.
  - This merged cleanly with S5-recal's `redecide` command; `cli.py` has both.

## 8. Open / not verified

- **Admission retry policy.** `DEFERRED(admission_failed)` is retryable, so an identical S5 can rerun
  up to 3 times. Making it non-retryable is a one-word change in `engine/retry.py`, but it would block
  the case where a settings change (doc 16) should retrigger the run. This is the orchestrator's call.
- The metallicity prior is checked only in-sample. Its levels and w_m = 0.20 are judgments.
- The cost fit has five hydride points and only high-symmetry cells.
- F1(b) and F3 are not applied retroactively to the campaign-1 ledger.
- There is no `ML_BELOW_DFT_HULL` penalty yet (R30 F8, second feature).
- The B3 symmetrization tolerance item (R30 F8) is not addressed.
- **Test results (meas., 2026-10-30):** full `tests/` 1 336 passed, 16 skipped, 0 failed (with
  `OMP_NUM_THREADS=2`, nice 19); 26 new tests in the three `*campaign1_fixes*` files.
