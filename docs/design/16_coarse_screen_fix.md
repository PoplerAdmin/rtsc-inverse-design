# 16 — Coarse screen (T2-SCREEN) fix: diagnosis, new statistic, H cutoff floor, recalibration

Agent S5-recal, 2026-10-30. Answers the soundness problem raised after campaign 1 (doc 30 §5.2)
and review R30 F2: can the S5 coarse screen remove a hydride that the full tier would keep?
Every number below is **computed** from QE outputs that already exist, unless it is marked
otherwise. Scripts: `src/scoring/sigma_ladder_table.py` (per-σ A, λ, ω_log, admission) and
`src/scoring/coarse_screen_recal.py` (old vs new statistic against references). Data:
`docs/design/data/16_mg2irh6_sigma_ladders.csv`, `docs/design/data/16_coarse_screen_recal.csv`.
Moments use the QE bin rule on the signed α²F at E > 0 (doc 13 §2). Admission is the core
`admit_smearings` (mode_positive, tol 1e-6).

## 0. Summary

- **The cause is the σ statistic. The 60 Ry cutoff is not the cause.** At matched σ, the
  campaign-1 coarse run (60 Ry, q 4³, k 8³/16³) and the P1 control (80 Ry, q 3³, k 8³/24³)
  agree within 1 % for σ = 0.02–0.05 Ry: 641/636, 518/520, 439/441, 386/387 K. This confirms R30 F2.
  The 341 K in doc 30 is A_mid, the median over the σ ≥ 0.03 Ry plateau. It is a reporting
  number. The removal statistic never used it.
- **The removal statistic (doc 13 ENVELOPE) was still unsound for this hydride.** It gave
  A_hi,coarse = 1093 K. That is **0.60×** the reference's A_hi lower bound of 1822 K, and
  (1 + e_c)·0.60 = 0.90 < 1. So a hydride with Mg2IrH6's σ-shape and a reference A_hi between
  1299 and 1443 K would have been screened out at T\* = 300 K, although the full tier would not
  reject it. The cause: A_nom sits at the smallest σ admitted on *both* meshes. Admission on
  q 3³ cut σ 0.01–0.05 Ry, so A_nom moved to σ = 0.06 Ry (351 K), while A at the physical
  small-σ end is ≈ 1000 K. Only the d_sig term, which happened to see the q 4³ σ = 0.01 row,
  kept the screen from failing outright. Counterfactual: if q 4³ admission had cut the same
  rows, (1 + e_c)·A_hi falls to **746 K = 0.41×** the reference.
- **The doc 30 "PASS by 2 %" is a misreading.** It compared A_nom (351 K) with the 346.5 K
  surrogate floor. The rule actually compared (1 + e_c)·A_hi = 1639 K with that floor, a 4.7×
  margin.
- **Fix (code + tiers.toml):**
  1. New removal statistic `coarse_envelope_all`: it uses **every** σ row, admitted or not,
     and adds a linear σ → 0 extrapolation, 2A(σ₁) − A(σ₂).
  2. Per-element cutoff floor `S5.ecut_floor_Ry = {H = 80}`, applied to the S5 profile.
     A removal is blocked (DEFERRED `ecut_below_floor`) if any mesh ran below the floor.
  3. `rtt redecide --stage S5` re-applies the current rule to stored metrics and appends a
     new outcome row.
- **Recalibration** (doc 13 method; n = 5 materials, 8 run pairs, 2 hydrides):
  - New A_hi,coarse / A_hi,ref ranges from **1.30** (Mg2IrH6) through 3.36 (H3S) to 11.8 (Nb).
  - Old: 0.60 (Mg2IrH6) … 10.7.
  - e_c stays **0.5, provisional**. The smallest headroom is now 1.5 × 1.30 = 1.95.
  - Hydrides keep the screen (not DEFERRED-only), but only at ≥ 80 Ry. §5 explains why.
- **Controls:**
  - Mg2IrH6 under the new rule: (1 + e_c)·A_hi = 3558 K (60 Ry files). That is 10× its
    346.5 K surrogate floor and 2.7× the 1299 K real target.
  - It also passes at A_target = its own reference A_hi (1822 K), with 1.95× margin.
  - The 80 Ry rtt re-run is running in the background (§6).
- **LaH2 re-opened:** the campaign ledger now has outcome #3081, DEFERRED(ecut_below_floor),
  which supersedes #3080 SCREENED_OUT. Its statistic is still far below target: 90 K against
  1299 K (old rule: 65 K, the ~20× margin R30 F2 cites). DEFERRED is retryable, so the
  scheduler will re-run its S5 at 80 Ry. The ~1.4 core-h run will most likely screen it out
  again, this time legitimately.

## 1. How A_hi,coarse was formed (code as of campaign 1)

`src/rtt/dft/coarse.py`, before this change:

- `coarse_metrics` builds one `SigmaRow` per σ (σ_i = 0.01·i Ry, i = 1..10). Rows that
  mode-positive admission rejects keep A and λ "for reporting" but get no ω_log.
- The plateau is the admitted rows with σ ≥ 0.03 Ry (`plateau_sigma_min_Ry`). It feeds only
  `A_mid_K` (median, reported) and the `no_plateau` blocker. **It never entered the removal
  statistic.**
- `coarse_envelope(fine = q 4³, coarse = q 3³)` (doc 13 §5):
  - s_nom = smallest σ admitted on both meshes;
  - A_nom = A(q 4³, s_nom);
  - d_sig = max |A_i − A_nom| / A_nom over the admitted q 4³ rows;
  - d_q and ε_ω are taken at s_nom;
  - A_hi = A_nom(1 + d_sig + d_q)(1 + ε_sys)/(1 − ε_ω).
- Removal rule: (1 + e_c)·A_hi < A_target, with e_c = 0.5.

Campaign-1 Mg2IrH6 (`~/rtsc_targeted_work/E1-campaign/work/c-9a8fa100a45b/S5`):
- q 3³ admitted σ 6–10 only. Its σ 1–5 rows have max negative per-mode fraction
  1.0·10⁻³ … 5.6·10⁻⁶, above the 10⁻⁶ tolerance.
- q 4³ admitted all ten σ.
- Therefore s_nom = σ₆ (0.06 Ry), A_nom = 351.4 K, d_sig = 1.575 (from q 4³ σ₁ = 904.7 K),
  d_q = 0.067, ε_ω = 0.066, A_hi = 1092.7 K, and (1 + e_c)·A_hi = 1639 K.

## 2. The Mg2IrH6 gap, factor by factor

Per-σ A (K) from the α²F files (`data/16_mg2irh6_sigma_ladders.csv`; * = excluded by
admission):

| σ (Ry) | E1 q3 60 Ry k9/18 | E1 q4 60 Ry k8/16 | P1 q3 80 Ry k8/24 | P1 q3 80 Ry k6/36 |
|---|---|---|---|---|
| 0.005 | – | – | 990.3* | 994.3* |
| 0.01 | 888.3* | 904.7 | 829.5* | 775.0* |
| 0.02 | 670.6* | 640.7 | 635.7* | 595.8* |
| 0.03 | 551.5* | 518.3 | 520.0* | 488.6* |
| 0.04 | 468.6* | 439.3 | 440.8* | 415.1* |
| 0.05 | 411.6* | 385.8 | 387.0* | 365.4* |
| 0.06 / 0.10 | 374.8 / 325.1 | 351.4 / 303.6 | – | – |

At tol 1e-6, neither P1 ladder passes admission on any σ: the max negative fraction is
2·10⁻⁵ … 2·10⁻³. Doc 15 §3 called all ten "passing" under a looser reading. This changes nothing
here, but it means doc 15's moments are A-SIGMA-EXCLUDED-grade.

Decomposition of "coarse A = 341 K" against "A ≈ 990 K":

| Factor | Effect on A | Evidence |
|---|---|---|
| Plateau median (σ ≥ 0.03, A_mid) vs smallest coarse σ (0.01) | **×2.65** (341 → 905) | E1 q4 |
| Ladder floor 0.01 vs 0.005 Ry | ×1.19 (830 → 990) | P1 k8/24 |
| Admission cut on q 3³ → A_nom at σ = 0.06 instead of 0.01 | ×2.57 (905 → 351) in A_nom; partly undone by d_sig | E1 |
| Cutoff 60 vs 80 Ry (with its geometry: a = 6.645 vs 6.652 Å) | **≤ 1 %** at σ 0.02–0.05; +9 % at 0.01 (confounded with q, k) | E1 q4 vs P1 k8/24 |
| k_fit 24³ vs 36³ (80 Ry, q 3³) | −6 … −7 % at σ 0.01–0.05; +0.4 % at 0.005 | P1 two runs |
| q 3³ vs q 4³ (60 Ry) | −2 … +6 % | E1 |

The cutoff is a real convergence issue for SG15 H (doc 15 §1: 60 Ry is 8.7 meV/atom and
−2.7 kbar from 120 Ry; 80 Ry is 1.6 meV/atom and −0.4 kbar). It also moves the relaxed lattice
by 0.1 %. But it does not lower coarse A in this material.

**Cutoff isolated (in progress):** the rtt re-run at 80 Ry uses the *same* k grids as E1
(q 3³: k 9³/18³, q 4³: k 8³/16³). Its relaxed cell is a = 2 × 3.32612 = **6.6522 Å**, compared
with 6.6452 Å at 60 Ry and 6.6521 Å in P1. When it finishes, it gives the pure cutoff effect
on the σ ladder (§6).

## 3. New rule

### 3.1 Cutoff floor (profile + blocker)

`tiers.toml [S5] ecut_floor_Ry = { H = 80.0, source = "doc15#1 …; doc16#3", provisional = true }`.

- `S5Coarse.run` applies `profiles.apply_ecut_floor`. It raises `ecutwfc` to the largest floor
  among the candidate's elements, never lowers it, and records it in `Profile.overrides` and in
  the tier_run params (`ecutwfc_Ry`).
- `decide_coarse(ecut_floors=…)` blocks a removal (DEFERRED `ecut_below_floor`) when any mesh
  ran below the floor.
- It also blocks (DEFERRED `ecut_unknown`) when the run's elements are unknown.
- Only H has a measured floor. Other elements keep the profile's 60 Ry. That is unverified for
  first-row O, N and F, which were not in any campaign so far.

### 3.2 Statistic (`coarse_envelope_all`, removal statistic since doc 16)

```text
rows      = every sigma row with A > 0 on a mesh, admitted or not
a_top(m)  = max( max_i A_i , 2 A(s_1) - A(s_2) )          # s_1 < s_2 the two narrowest sigma
A_nom     = max over meshes of a_top
d_sig     = (A_nom - min_i A(fine, s_i)) / A_nom
s*        = smallest sigma with A and omega_log on both meshes (all rows)
d_q       = |A(fine, s*) - A(coarse, s*)| / A(fine, s*)
eps_w     = |w_log(fine, s*) - w_log(coarse, s*)| / w_log(fine, s*), floor 0.02, cap 0.9
A_hi      = A_nom (1 + d_sig + d_q)(1 + eps_sys)/(1 - eps_w)
SCREENED_OUT(coarse-A envelope (all-sigma))  iff  (1 + e_c) A_hi < A_target
single mesh (q >= 4^3): d_q = 0, eps_w = 0.02;  M_c A_hi < A_target, M_c = 2
```

Why excluded rows are used: mode-positive admission exists to protect *proofs*. A smearing that
is too narrow for the k grid gives per-mode negative α²F. For a heuristic *removal*, the safe
direction is the opposite one. The narrow-σ rows are exactly where a sharp-N(E_F) hydride's A
is largest, so they must raise the statistic, not vanish from it. Noisy-high excluded rows (Al,
§4) only make removal harder.

Why the extrapolation: A(σ) of Mg2IrH6 is convex and decreasing (differences per 0.005 Ry:
160, 114, 80 K). The linear σ → 0 value from the two narrowest coarse rows is 1169 K. That is
above A(0.005) = 994 K and close to what further narrowing would give. For flat ladders
(elements) it equals the narrowest value.

Unchanged: all doc 13 preconditions; never on λ, ω_log, m or imaginary modes; outcomes tagged
PROVISIONAL_MARGIN. The doc 13 value is still printed in the outcome reason
(`doc13 A_hi …`) for comparison.

## 4. Recalibration (doc 13 §2 comparison 2(b))

References:
- Al, Nb, Pb, H3S: doc 13's bench references, A_nom,ref and the A_hi,ref lower bound.
- Mg2IrH6: P1 control, 80 Ry, q 3³, k 6³/36³. A_nom,ref = 994.3 K (σ = 0.005).
  A_hi,ref ≥ 994.3·(1 + 0.633)·1.10/0.98 = **1822 K** (doc 13 formula over all its rows).
- Coarse pairs: q 3³ + q 4³ (doc 13 P1, P2; H3S T4 c2/c3; Mg2IrH6 campaign 1 at 60 Ry).

| Material | Pair | A_hi,old / A_hi,ref | A_hi,new / A_hi,ref | A_hi,new / A_nom,ref | new: A_nom, d_sig, d_q, ε_ω |
|---|---|---|---|---|---|
| Al | P1 | 1.171 | 8.89 | 10.0 | 232, 0.72, 0.52, 0.16 |
| Al | P2 | 1.224 | 9.04 | 10.2 | 244, 0.72, 0.50, 0.14 |
| Nb | P1 | 10.73 | 11.79 | 17.4 | 133, 0.37, 0.08, 0.90 |
| Nb | P2 | 10.45 | 11.50 | 17.0 | 131, 0.37, 0.07, 0.90 |
| Pb | P1 | 1.618 | 4.17 | 4.90 | 77, 0.53, 0.18, 0.18 |
| Pb | P2 | 1.547 | 3.77 | 4.42 | 74, 0.54, 0.18, 0.13 |
| H3S | T4 c2/c3 | 0.899 | 3.36 | 5.06 | 3292, 0.65, 0.45, 0.09 |
| **Mg2IrH6** | E1 60 Ry | **0.600** | **1.302** | 2.39 | 1169, 0.74, 0.02, 0.05 |

The old column reproduces doc 13's `13_coarse_envelope.csv` exactly. That is a check of the
script.

**e_c.**
- No pair now has A_hi,new < A_hi,ref; the worst is 1.30.
- doc 13's rule of thumb was e_c ≈ 3× the worst observed envelope shortfall. Under the new
  statistic there is no observed shortfall, so e_c covers only unobserved materials.
- I keep **e_c = 0.5, provisional**. With it the worst headroom is 1.95×, compared with 1.35×
  claimed by doc 13 and 0.90× actually present.
- This is a judgment, not a quantile: n = 5 materials, 2 of them hydrides.

**Cost in screening power** (campaign-1 S5 runs, new vs old (1 + e_c)·A_hi, all against
1299 K):

| Candidate | Old (K) | New (K) | Margin, new |
|---|---|---|---|
| LaH₂ | 65 | 90–102 | 13–14× |
| ScH₂ | 132 | 249–279 | ≈ 5× |
| NaIrH₃ (q 4³ admitted, q 3³ not) | no envelope | 431 | 3× |
| BaPtH₃ (unstable, DEFERRED anyway) | 477 | 1333–2037 | above target |

- The ranges come from two ways of evaluating the same runs. The higher value re-reads the
  α²F files, which also gives the excluded rows an ω_log. The lower value uses the stored
  pre-doc-16 records, which lack those ω_log values.
- Elements are now 4–12× above their own references instead of 1.2–1.6×. The screen removes
  only candidates far below target, which is what a heuristic removal should do.

## 5. Hydrides: screen or DEFERRED-only?

The brief allows DEFERRED-only for H compounds if no sound e_c exists. I chose to **keep the
screen for hydrides, with the 80 Ry floor as a precondition**, for three reasons:

1. With the new statistic, both hydride pairs overshoot the reference (1.30×, 3.36×) before
   any e_c is applied. The old statistic undershot on both (0.60×, 0.90×).
2. The specific failure mode (admission removing the large-A end) is now structurally
   excluded.
3. DEFERRED-only would leave every hydride retrying S5 (DEFERRED is retryable,
   `rtt.engine.retry`) or force it into S6 at roughly 45× the cost (doc 30 §5.4). That is a
   worse use of compute than a screen whose worst calibrated headroom is 1.95×.

What this does **not** cover, stated plainly:
- n = 2 hydrides. One is at 200 GPa. For the other (Mg2IrH6), the "reference" is itself a
  q 3³ run without δ^q, so A_hi,ref is a lower bound and could be higher.
- A strongly coupled soft mode at a q that 3³/4³ does not sample (a Kohn anomaly) could lower
  coarse A with respect to the converged A. None of the five materials shows this, so it is
  untested.
- Low-symmetry, 6–20-atom hydrides are absent from the calibration.

If the audit lane (doc 12 §3.1: 5 % of SCREENED_OUT pushed to T3) or the first S6 hydride
finds A_hi,coarse / A_hi,full < 1/(1 + e_c) = 0.67, switch H compounds to DEFERRED-only with a
non-retryable reason (`NON_RETRYABLE_DEFERRED` in `rtt.engine.retry`). Do not just raise e_c.

## 6. Mg2IrH6 through the real `rtt` S5 path at 80 Ry (running)

- Command:
  `rtt --run ~/rtsc_targeted_work/S5-recal/run --ledger ~/rtsc_targeted_work/S5-recal/ledger.sqlite --workdir ~/rtsc_targeted_work/S5-recal/work --ranks 4 step S5 c-9a8fa100a45b`.
  The ledger is a SQLite backup copy of the campaign ledger, so the campaign ledger is not
  touched.
- Settings: the profile used is `coarse` with the override `ecutwfc_Ry = 80`. The q 3³ and q 4³
  meshes and k grids match campaign 1.
- Progress: started 2026-10-30 00:57 (local). relax 107 s, scf_fit 101 s and scf 21 s on 4 ranks.
  ph.x on q 3³ is running.
- Estimate: campaign 1 at 60 Ry took 29.7 core-h (7.7 h wall). Scaled by N_pw ∝ ecut^1.5
  (×1.54), this gives ≈ 46 core-h, ≈ 12 h wall. This is an extrapolation, not a measurement.
- The run process loaded the pre-doc-16 `decide_coarse`, so its own S5 outcome uses the doc 13
  statistic. A waiter, `~/rtsc_targeted_work/S5-recal/post_mg.sh`, runs when the job exits.
  It applies `rtt redecide --stage S5 --all c-9a8fa100a45b` (doc 16 rule) to the copy ledger,
  tabulates the 80 Ry σ ladders, and appends the run to the recalibration table. Output goes to
  `~/rtsc_targeted_work/S5-recal/post_mg.out`.
- Expected result, from §2 (≤ 1 % per-σ cutoff effect): (1 + e_c)·A_hi ≈ 3500 K vs the
  346.5 K floor, a ~10× margin. **Not yet verified.**

## 7. LaH₂ (campaign 1) re-opened

- Mechanism: `rtt redecide --stage S5 c-37743845859b` on the campaign ledger (new CLI path,
  `S5Coarse.redecide`). It recomputes `decide_coarse` from the stored per-mesh metrics.
  Records from before doc 16 lack the cutoff, so it reads the cutoff from each mesh's
  `dfpt_run.json` (60 Ry). It then appends an outcome row.
- No SQL row was edited. The append-only triggers stay intact.
- Result: outcome #3081 is DEFERRED, criterion `ecut_below_floor`, tag REDECIDED_DOC16. The
  reason names the superseded row #3080 SCREENED_OUT. #3081 is now current.
- Why: under the new rule a removal needs ≥ 80 Ry for H, and LaH₂ ran at 60 Ry.
- Its statistic is 90 K against 1299 K. A re-run at 80 Ry will almost certainly screen it out
  again, and that removal will then be properly founded.
- No other campaign-1 S5 outcome was a removal.
- A dry run with `--all` on a ledger copy shows:
  - ScH₂ stays DEFERRED, now also blocked by ecut;
  - BaPtH₃ would go from DEFERRED(imaginary_modes) to PASS, because its statistic is now
    above target;
  - every PASS stays PASS.
- `--all` was **not** applied to the campaign ledger. That is an orchestrator decision.

## 8. Code changes and tests

| File | Change |
|---|---|
| `src/rtt/dft/coarse.py` | `a_top`, `coarse_envelope_all` (removal statistic), `ecut_floor_blockers`; `decide_coarse(ecut_floors=…)`; excluded rows now get ω_log; `CoarseMetrics.ecutwfc_Ry`, `.elements`, `.from_metrics`; `to_metrics` stores admission_ok / ecutwfc_Ry / elements |
| `src/rtt/dft/profiles.py` | `ecut_floor_Ry`, `apply_ecut_floor` |
| `src/rtt/engine/stages.py` (S5Coarse only) | floor applied to the profile; ecutwfc in tier_run params; `ecut_floors` passed to decide; `S5Coarse.redecide` |
| `src/rtt/cli.py` | `redecide --stage S5 [cids] [--all] [--dry-run]` (was a stub) |
| `src/rtt/configs/tiers.toml` | `S5.ecut_floor_Ry` (new, provisional); `e_c` and `margin_M_c` source strings now point to doc 16; plateau comment |
| `tests/back/test_coarse_screen_doc16.py` | new, 7 tests: extrapolation, the Mg2IrH6 under-coverage of the doc 13 statistic, the admission-cut counterfactual, cutoff blockers, profile floor + tiers sourcing, metrics round-trip, re-opening of a stored SCREENED_OUT |
| `tests/back/test_collect_coarse.py`, `tests/test_engine.py` | expected values updated to the new statistic; metrics fixture given elements |

`scripts/pyt.sh tests/back tests/review tests/front tests/test_engine.py tests/test_cli_paths.py
tests/test_ledger_triggers.py`: 272 passed, 4 skipped, 1 failed. The failure is
`tests/front/test_s3.py::test_decide_s3_hull_rule`, which belongs to the S3 code (being edited
by another agent) and does not touch S5.

## 9. What remains uncalibrated

- **e_c = 0.5**: judgment on n = 5 materials (2 hydrides). It has no statistical meaning.
- **M_c = 2** (single-mesh fallback): its statistic changed, and it has no single-mesh hydride
  data.
- **The cutoff floor** exists only for H (80 Ry, from total energy, pressure and forces, doc 15).
  The per-σ A effect of 60 → 80 Ry is not yet isolated (§6 pending). The other light elements
  have no measured floor.
- **The Mg2IrH6 reference** is q 3³ at k_fit ≤ 36³. No q 4³/6³ full-grade reference exists
  (doc 15 §4: ≈ 1.4 days on 4 ranks for q 4³ / k_fit 48³).
- **The σ → 0 extrapolation** is linear from two rows. For a σ-ladder much steeper than
  Mg2IrH6's at σ < 0.01 Ry, it could still undershoot. The `d_sig` term adds margin, but by
  construction, not by measurement.
- **Other agents' pending S5 work:** `src/rtt/dft/admission.py` (admission failure →
  DEFERRED) calls into `S5Coarse.decide`, so it should be checked together with this change.
