# 33. Campaign 2 S5 queue (interim)

Agent S5-queue, 2026-10-30. **Interim status:** one result (LaH₂) is in. LiScCoH₆ is running.
Ba₂IrH₆ and LiScIrH₆ have not started. Labels: (meas.) measured, (est.) estimated.

## Setup
- Code: the real `rtt` CLI, `step S5`, doc-16 statistic (`coarse_envelope_all`), H cutoff floor 80 Ry.
  Each job runs on its own SQLite backup copy of the source ledger, under
  `~/rtsc_targeted_work/S5-queue/<job>/ledger.sqlite`. The E1 copy was taken at 07:40 local while
  NaIrH₃ S6 was still writing the E1 ledger. The lane-P copy was taken at the same time. No source
  ledger was written.
- Launcher: `~/rtsc_targeted_work/S5-queue/launch.sh`. It checks the running pw.x/ph.x count
  before each launch (≤ 12 in total) and starts the job with setsid, nice 19 and OMP = 1.
- Scheduler estimate (`estimate_core_h`, est.): 69.6 core-h for each 9-atom cell (Ba₂IrH₆,
  LiScIrH₆, LiScCoH₆, Mg₂IrH₆), 2.58 core-h for LaH₂. All are below 48 h wall on ≥ 2 ranks, so
  none was skipped as compute-limited.

## Results

| cand. | q | σ (Ry) 0.01 … 0.05 A (K) | λ (σ = 0.01) | ω_log (K, σ = 0.01) | (1 + e_c)·A_hi (K) | /1299.44 | rtt decision |
|---|---|---|---|---|---|---|---|
| LaH₂ (c-37743845859b), 80 Ry | 3³ | 36.1, 35.0, 34.7, 34.4, 34.1 | 0.193 | 222 | | | |
| | 4³ | 38.9, 36.8, 35.7, 35.0, 34.3 | 0.212 | 237 | **95.0** | **0.073** | SCREENED_OUT (outcome #3083, copy ledger) |

- LaH₂ statistic: A_nom = 41.1 K (all-σ, σ → 0), d_sig = 0.240, d_q = 0.072, ε_ω = 0.064. The
  doc-13 A_hi is 40.5 K.
- The result supersedes #3081, DEFERRED(ecut_below_floor). The removal now has a proper basis
  under the code rule.
- **Decision: SCREENED_OUT stands** (orchestrator, 2026-10-30).
  - R30 §3.3 allows no campaign-2 removal until precondition P2 holds: the removal rule
    re-derived for hydrides.
  - Doc 16 did that re-derivation: a new all-σ statistic, recalibrated with two hydrides,
    Mg₂IrH₆ (1.30×) and H3S (3.36×). No calibration pair now undershoots its reference.
    Worst headroom with e_c = 0.5: 1.95×.
  - LaH₂'s margin is 13.7× below target, about 7× beyond that worst calibrated headroom.
  - Caveat: P2 as worded in R30 is an S6-lite run of Mg₂IrH₆, and it has not been done.
    The hydride e_c therefore rests on n = 2 hydride calibration points (doc 16 §5, §9).
- Mg₂IrH₆ comparison at matched settings (doc 16, 60 Ry, q 4³): A = 905 K at σ = 0.01 Ry. LaH₂ is
  about 0.04× of the control. The 80 Ry Mg₂IrH₆ data is pending:
  `~/rtsc_targeted_work/S5-recal/post_mg.out`.
- Full per-σ data: `~/rtsc_targeted_work/S5-queue/tables/lah2.csv`.

## Cost (meas.)
- LaH₂ at 80 Ry on 2 ranks: q 3³ took 2707 s (1.50 core-h) and q 4³ took 4810 s (2.66 core-h).
  The total, 4.16 core-h, is **1.6× the scheduler estimate**.
- Campaign 1 at 60 Ry took 1.35 core-h. So 60 → 80 Ry cost ×3.1 here, against the model's ×1.54.
  The machine was loaded with 12 ranks, so wall-time inflation is confounded with the cutoff.

## Not done / NOT CHECKED
- LiScCoH₆ (lane P, c-7c5823b24da2): running on 2 ranks since 09:47 local. Log:
  `~/rtsc_targeted_work/S5-queue/liscoh6/s5.log`.
- Ba₂IrH₆ (c-539acaa4a6e7) and LiScIrH₆ (c-43d572803f9d): not launched. The plan is to start
  Ba₂IrH₆ on 4 ranks when the Mg₂IrH₆ 80 Ry run frees its ranks, using
  `launch.sh ba2irh6 ~/rtsc_targeted_work/S5-queue/ledger_E1copy.sqlite c-539acaa4a6e7 4`.
- The primary campaign-2 question (does any non-control beat Mg₂IrH₆'s A at matched settings) is
  **unanswered**.

## Final results (collected 2026-11-01; the campaign was then closed)

| Candidate | Source | Decision | (1+e_c)·A_hi,coarse | A_nom (all σ, σ→0) | Share of 1299.4 K | Measured cost |
|---|---|---|---|---|---|---|
| Ba₂IrH₆ (c-539acaa4a6e7) | lane D | SCREENED_OUT | 134.7 K | 57.0 K | 0.10 | 88.0 core-h, 23.2 h wall (4 ranks) |
| LiScCoH₆ (c-7c5823b24da2) | lane P | SCREENED_OUT | 395.1 K | 146.6 K | 0.30 | 56.8 core-h, 29.6 h wall (2 ranks) |
| LaH₂ (c-37743845859b) | lane D | SCREENED_OUT | 95.0 K | 41.1 K | 0.07 | 4.2 core-h |
| LiScIrH₆ (c-43d572803f9d) | lane P | not finished | — | — | — | killed on q 3³ when the run machine was restarted |

**Answer to the pre-registered campaign-2 question:** no. No non-control candidate came near
the Mg₂IrH₆ control (A ≈ 994 K at σ → 0, 80 Ry; doc 16). The best, LiScCoH₆, has an upper
envelope of 0.30× the 300 K floor. Ba₂IrH₆, same electron count as Mg₂IrH₆, is about 17× weaker:
the electron count alone does not carry the coupling.

**First real S6 (NaIrH₃, c-0ae3f133a308):** both full meshes (q 4³, q 6³) finished, then the
decision step recorded ERROR. The same mode-positive smearing admission that failed its S5 q 3³
mesh rejected every smearing: the largest negative per-mode α²F fraction was 3×10⁻⁶ at σ = 0.05
Ry against the core tolerance 10⁻⁶ (R30 F5 predicted this). The tolerance is too strict for
hydride spectra; it needs a calibrated value (and a test on the Mg₂IrH₆/H3S references) before
any S6 can complete on a hydride. Open; not changed here because the core file is sha-pinned.

**Cost check:** the 9-atom S5 runs took 57–88 core-h, inside the doc-31 refit's 95 core-h class
estimate, but 1.3–1.9× the per-cell estimate (69.6 core-h) on a fully loaded machine.
