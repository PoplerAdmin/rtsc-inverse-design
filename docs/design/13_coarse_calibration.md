# 13 — Calibration of the coarse-DFPT tier (T2) against bench references: Al, Pb, Nb, H3S

Status: **measured** (all numbers below come from runs whose inputs, outputs and SHA-256 are listed
in `docs/design/data/13_coarse_runs_manifest.csv`). Answers doc 11 §4 ("this margin must be
re-checked as more T3 references … are compared with their own coarse runs") and review R12
finding 3 ("report the calibration honestly").

## 0. Summary

- **A = ∫α²F dω and m = λ⟨ω²⟩ survive coarse DFPT. λ and ω_log do not.** Per smearing σ (at the
  σ values that both the coarse and the reference ladder admit), q 3³–4³ coarse runs give
  A_coarse/A_full = **0.925 … 1.106** and m_coarse/m_full = **0.926 … 1.054** over all 4 materials
  and 3 coarse profiles (14 runs). In the same runs λ_coarse/λ_full = 0.80 … 1.70 and
  ω_log,coarse/ω_log,full = 0.14 … 1.28.
- **The per-σ ratio is not the number that matters for deferral.** Deferral compares the
  coarse run's *own admitted ladder* with a target. The full tier compares the *reference's
  upper envelope* (theory doc 01 §0.4 ENVELOPE) with the same target. The coarse ladder is often
  cut at small σ by mode-positive admission (H3S q 3³ admits only σ ≥ 0.07 Ry; Al admits only σ ≥
  0.02–0.025 Ry), and the reference's A falls steeply with σ (H3S, Nb). So the doc-11 statistic
  max_σ A_coarse falls as low as **0.53 × A_hi,full** (H3S q 3³) and **0.69 × A_hi,full** (Nb
  q 4³). Against the reference's *nominal* A, the worst case is 0.80 (H3S q 3³).
- **The doc-11 factor-2 rule holds on this sample, but with little headroom for a single q 3³
  mesh** (1/0.53 = 1.89 is just below 2). **Recommendation (§5):** screen out on the coarse
  ENVELOPE built from the coarse ladder with two coarse q meshes (3³ and 4³), using
  **(1 + e_c)·A_hi,coarse < A_target with e_c = 0.5** (worst observed A_hi,coarse/A_hi,full =
  0.90, so a margin of 1.5 leaves 1.35× headroom). If only one coarse mesh exists, keep the doc-11
  factor 2 and require q ≥ 4³ (worst 0.69, headroom 1.38). Never screen out on q 2³ alone, or on
  a ladder that admission rejects.
- **Spurious imaginary modes are common at coarse settings.** Nb at q 4³ shows −30 … −32 cm⁻¹ on
  the matdyn interpolation mesh in all 3 profiles. H3S shows −124 cm⁻¹ (q 3³ mesh) and
  −190/−474 cm⁻¹ (q 6³ grid/mesh). Pb (k_fit 16³) shows −8.6 cm⁻¹. None of the four references
  has any mode below −1 cm⁻¹. These runs are exactly the ones where λ inflates (+37 … +122 %) and
  ω_log collapses (×0.14 … 0.48), while A stays within +1 … +32 %.
- Cost on 6 ranks (machine shared with other agents, load average ≈ 12): ph.x 95–894 s at q 3³–4³
  for these 1-atom cells. The q 3³ runs are 13–46× cheaper per core-second than the 8-rank
  references (Nb q 4³ is only 5.5× cheaper).
- The outcome is a heuristic **SCREENED_OUT(T2, "coarse-A envelope", e_c = 0.5, calibration =
  this doc)**. It is never a proof. The calibration relates coarse to **bench-grade** references
  (q 8³ / 6³, k_fit 32³ / 36³), not to converged spectra. It uses four 1–4-atom cubic cells and
  is **not validated for 6–10-atom low-symmetry hydrides** (§6).

## 1. What was run

Chain: `src/scoring/coarse_dfpt.py`, which generalises `coarse_dfpt_probe.sh`: scf_fit (la2F,
k_fit) → scf (k_scf) → ph.x (`electron_phonon='interpolated'`, ldisp, q nq³, same el_ph_sigma
ladder as the reference) → q2r (la2F) → matdyn (la2F, dos, 32³ mesh). The structure and
pseudopotentials are taken from the reference's own `scf_fit.in` (SG15 ONCV PBE 1.2; relaxed bench
cells). The runs used 6 MPI ranks × 1 thread at nice 19, with one job at a time. Each workdir
gets a `coarse_run.json` with settings, per-step wall time, return codes, the lowest frequency
per irreducible q and on the matdyn mesh, and the SHA-256 of every output.

| Material | Reference (earlier benchmark run, 8 ranks) | Reference ph.x wall | σ ladder (Ry) |
|---|---|---|---|
| Al fcc 0 GPa | ecut 60/240, k 16³/32³, q 8³ (29 irr. q) | 54 min | 0.005 … 0.05 |
| Pb fcc 0 GPa (no SOC) | ecut 50/200, k 16³/32³, q 8³ (29) | 2 h 16 min | 0.01 … 0.10 |
| Nb bcc 0 GPa | ecut 120/480, k 16³/32³, q 8³ (29) | 1 h 1 min | 0.01 … 0.10 |
| H3S Im-3m 200 GPa | ecut 90/360, k 18³/36³, q 6³ (16) | 1 h 58 min | 0.01 … 0.10 |

Coarse profiles (Al, Pb, Nb run here; H3S runs reused from `~/rtsc_targeted_work/T4-surrogates/h3s_c1..c4`,
made with `coarse_dfpt_probe.sh` on 2 ranks):

| Profile | ecut | k_scf / k_fit | q meshes |
|---|---|---|---|
| P1 | reference value | 12³ / 24³ | 2³, 3³, 4³ |
| P2 | ≈ 2/3 of reference (Al 40, Pb 34, Nb 80 Ry; ecutrho = 4×) | 12³ / 24³ | 3³, 4³ |
| P3 | reference value | 8³ / 16³ | 2³, 4³ |
| H3S T4c1–c4 | 60 Ry (2/3) | 8³/16³ (c1), 12³/24³ (c2–c4) | 2³, 3³, 4³, 6³ |

ph.x wall time (s, 6 ranks; H3S 2 ranks):

| | q 2³ | q 3³ | q 4³ |
|---|---|---|---|
| Al P1 / P2 / P3 | 39 / – / 14 | 95 / 71 / – | 211 / 166 / 87 |
| Nb P1 / P2 / P3 | 109 / – / 40 | 361 / 200 / – | 894 / 428 / 342 |
| Pb P1 / P2 / P3 | 65 / – / 32 | 203 / 86 / – | 502 / 213 / 161 |
| H3S (2 ranks) | 36 | 271 | 635 (q 6³: 2099) |

pw.x took 3–10 s per step. q2r + matdyn took 20–40 s. Wall times are inflated by
concurrent load from other agents, so treat them as upper values.

## 2. Method

- **Moments** (`src/scoring/coarse_calibration.py`): QE bin rule on the signed a2F_total at
  E > 0, with Δ from the file footer: A = Σ a2F·Δ (K), λ = 2Σ a2F/ω·Δ,
  ω_log = exp[(2/λ)Σ a2F ln ω/ω·Δ], m = λ⟨ω²⟩ = 2Σ ω·a2F·Δ (K²). These are the same quantities the
  core pipeline's `A-SIGMA-EXCLUDED` tag uses. λ agrees with QE's footer to < 0.5 %, except in
  runs with imaginary mesh modes. There the footer, which also sums the E ≤ 0 bins, differs by 8 % (Nb P2 q4) up to a factor 3.3 (Nb P1 q4); H3S c4 differs by 16 %.
- **Smearing admission:** `admit_smearings` and `mode_negative_fractions` are taken
  unchanged in logic from the core (mode-positive admission)
  (`rtt/core/dfpt/workflow.py`, `rtt/core/dfpt/parsers.py`; tolerance 1e-6). They are applied to
  coarse and reference ladders separately. Admitted reference ladders: Al σ₂…σ₁₀, Pb σ₃…σ₁₀,
  Nb and H3S σ₁…σ₁₀.
- **Three comparisons**, each answering a different question:
  1. *Per-σ ratio* X_coarse(σ_i)/X_full(σ_i) over σ admitted in both. Answers "is the coarse
     spectrum close to the full one?"
  2. *Doc-11 statistic* max over the coarse run's own admitted σ of A_coarse, divided by
     (a) the reference's nominal A (smallest admitted σ) and (b) a lower bound on the reference's
     ENVELOPE value A_hi,full ≥ A_nom(1+δ_A^σ)(1+ε_sys)/(1−0.02). This is a lower bound because
     the bench has one q mesh, so δ^q = 0 and ε_ω is at its floor. (b) is the comparison that
     decides whether a coarse screen-out could remove something the full-tier N1-ISO rule would
     not reject.
  3. *Coarse ENVELOPE*: theory doc 01 §0.2(b)/§0.4 applied to the coarse ladder, with q 4³ as
     "finest" and q 3³ as the second mesh (`src/scoring/coarse_envelope.py`):
     A_hi,c = A_nom,c(1+δ^σ+δ^q)(1+0.10)/(1−ε_ω), with ε_ω = max(q-spread of ω_log, 0.02) capped
     at 0.9, and X_nom = value at the smallest σ admitted on both coarse meshes. The theory doc
     does not define X_nom explicitly; this is my reading.

## 3. Results

### 3.1 Per-σ ratios coarse/full (min … max over σ admitted in both)

| Run | adm. σ (coarse) | A | λ | ω_log | m | lowest freq. (cm⁻¹): q grid excl. Γ / matdyn mesh |
|---|---|---|---|---|---|---|
| Al P1 q2 | 2–10 | 0.95 … 1.67 | 0.99 … 1.76 | 0.94 … 0.97 | 0.91 … 1.56 | 144 / 0 |
| Al P1 q3 | 4–10 | 0.93 … 1.04 | 0.94 … 1.03 | 0.98 … 1.00 | 0.93 … 1.04 | 135 / 0 |
| Al P1 q4 | 5–10 | 0.97 … 1.03 | 0.94 … 1.00 | 1.05 … 1.06 | 0.98 … 1.05 | 115 / 0 |
| Al P2 q3 | 4–10 | 0.97 … 1.08 | 1.02 … 1.12 | 0.95 … 0.97 | 0.93 … 1.04 | 132 / 0 |
| Al P2 q4 | 5–10 | 1.02 … 1.08 | 1.02 … 1.08 | 1.02 … 1.02 | 0.98 … 1.05 | 113 / 0 |
| Al P3 q2 | 2–10 | 0.65 … 0.94 | 0.71 … 0.97 | 0.93 … 0.98 | 0.59 … 0.91 | 144 / 0 |
| Al P3 q4 | **none** (admission fails) | – (all σ: 0.51 … 1.09) | – | – | – | 117 / 0 |
| Nb P1 q2 | 1–10 | 0.85 … 1.40 | 0.67 … 1.12 | 1.35 … 1.37 | 0.98 … 1.60 | 133 / 0 |
| Nb P1 q3 | 1–10 | 0.94 … 0.99 | 0.81 … 0.86 | 1.24 … 1.27 | 0.99 … 1.04 | 111 / 0 |
| Nb P1 q4 | 1–10 | 1.01 … 1.02 | **1.37 … 1.43** | **0.39 … 0.45** | 0.98 … 0.99 | 72 / **−32.4** |
| Nb P2 q3 | 1–10 | 0.94 … 0.98 | 0.80 … 0.85 | 1.25 … 1.28 | 0.99 … 1.04 | 112 / 0 |
| Nb P2 q4 | 1–10 | 0.99 … 1.01 | **1.47 … 1.58** | **0.14 … 0.19** | 0.98 … 0.99 | 74 / **−30.7** |
| Nb P3 q2 | 1–10 | 0.95 … 1.37 | 0.86 … 1.27 | 1.15 … 1.17 | 0.98 … 1.37 | 109 / 0 |
| Nb P3 q4 | 1–10 | 1.03 … 1.11 | **1.51 … 1.70** | **0.36 … 0.41** | 0.93 … 0.99 | 65 / **−29.5** |
| Pb P1 q2 | 1–10 | 1.01 … 1.09 | 1.07 … 1.17 | 0.95 … 0.96 | 0.93 … 0.99 | 34 / 0 |
| Pb P1 q3 | 2–10 | 0.97 … 1.02 | 1.00 … 1.03 | 0.95 … 0.98 | 0.97 … 1.04 | 28 / 0 |
| Pb P1 q4 | 1–10 | 0.98 … 1.00 | 0.98 … 1.00 | 1.00 … 1.02 | 0.97 … 0.99 | 25 / 0 |
| Pb P2 q3 | 2–10 | 0.92 … 0.97 | 0.88 … 0.91 | 1.05 … 1.07 | 0.97 … 1.04 | 32 / 0 |
| Pb P2 q4 | 1–10 | 0.94 … 0.96 | 0.87 … 0.90 | 1.08 … 1.11 | 0.97 … 0.99 | 28 / 0 |
| Pb P3 q2 | 1–10 | 0.99 … 1.05 | 0.99 … 1.06 | 1.00 … 1.01 | 0.96 … 1.02 | 34 / 0 |
| Pb P3 q4 | 1–10 | 0.97 … 1.02 | 1.23 … 1.34 | **0.43 … 0.48** | 0.97 … 1.03 | 21 / **−8.6** |
| H3S c1 q2 | **none** (non-contiguous) | – (all σ: 0.45 … 0.71) | – | – | – | 321 / −0.9 |
| H3S c2 q3 | **7–10 only** | 1.09 … 1.10 | 1.37 … 1.38 | 0.74 | 0.97 | 285 / **−124** |
| H3S c3 q4 | 1–10 | 0.99 … 1.08 | 1.10 … 1.20 | 0.81 … 0.85 | 0.93 … 1.04 | 243 / −0.9 |
| H3S c4 q6 | 1–10 | 1.13 … 1.32 | **1.90 … 2.22** | **0.29 … 0.31** | 0.95 … 1.13 | **−190** / **−474** |

The reference meshes' lowest frequencies are 0 (Al, Pb, Nb) and −0.9 cm⁻¹ (H3S).
Full per-σ values are in `data/13_coarse_calibration_per_sigma.csv`.

**Aggregates, q 3³–4³ (14 runs with an admitted ladder):** A 0.925 … 1.106; m 0.926 … 1.054;
λ 0.80 … 1.70; ω_log 0.14 … 1.28. **q 2³:** A 0.65 … 1.67, and admission fails outright for H3S.

### 3.2 Deferral-relevant statistics

| Material | max_σ A_coarse / A_nom,full (q3; q4) | max_σ A_coarse / A_hi,full lower bound (q3; q4; max(q3,q4)) | coarse ENVELOPE A_hi,c (q3+q4) / A_hi,full lower bound |
|---|---|---|---|
| Al (P1 / P2) | 1.03; 1.03 / 1.08; 1.07 | 0.92; 0.91; 0.92 / 0.96; 0.95; 0.96 | 1.17 / 1.22 |
| Nb (P1 / P2) | 0.94; 1.02 / 0.94; 1.01 | 0.64; 0.69; 0.69 / 0.63; 0.68; 0.68 | 10.7 / 10.5 (ε_ω capped: soft-mode ω_log spread) |
| Pb (P1 / P2) | 1.02; 1.61 / 0.98; 1.54 | 0.87; 1.37; 1.37 / 0.83; 1.31; 1.31 | 1.62 / 1.54 |
| H3S (c2, c3) | 0.80; 1.08 | **0.53**; 0.72; 0.72 | **0.90** |

(A_hi,full lower bounds: Al 76.3 K, Nb 180.2 K, Pb 42.5 K, H3S 2482 K; reference A_nom: 67.5,
121.8, 36.1, 1648 K.) Source: `data/13_coarse_envelope.csv`.

Reading:
- Against the reference *nominal* A, max_σ A_coarse is never below 0.80, so a 20 % shortfall is
  the worst seen. This is the number doc 11 implicitly assumed.
- Against what the full tier actually certifies with (A_hi,full, which includes the reference's
  own σ spread and ε_sys), the shortfall reaches 47 % (H3S q 3³). The reason is not that coarse A
  is wrong at a given σ (per-σ it is +9 … +10 %). Admission removes the small-σ end of the coarse
  ladder, where A is largest. Nb's 31–37 % shortfall has the same origin: the reference A drops
  32 % from σ₁ to σ₉, and the coarse max sits at σ₁ but is not raised by the reference's spread.
- The coarse ENVELOPE (two coarse meshes) carries its own σ- and q-spreads, so it tracks
  A_hi,full much better (worst 0.90, H3S). When the coarse phonons are soft (Nb q 4³, H3S q 6³),
  the ω_log q-spread drives ε_ω to its cap. A_hi,c then becomes huge, and nothing is screened
  out. That is the safe direction.

### 3.3 False imaginary modes

| Run | q grid (non-Γ) | matdyn mesh | Reference mesh | Effect on A / λ / ω_log (per σ) |
|---|---|---|---|---|
| Nb q 4³ (P1, P2, P3) | +65 … +74 | −29.5 … −32.4 cm⁻¹ (−42 … −47 K) | 0 | +1 … +11 % / +37 … +70 % / −55 … −86 % |
| Pb P3 q 4³ (k_fit 16³) | +21 | −8.6 cm⁻¹ (−12 K) | 0 | −3 … +2 % / +23 … +34 % / −52 … −57 % |
| H3S c2 q 3³ | +285 | −124 cm⁻¹ | −0.9 | +9 … +10 % / +37 % / −26 % |
| H3S c4 q 6³ | −190 | −474 cm⁻¹ | −0.9 | +13 … +32 % / +90 … +122 % / −69 … −71 % |

All of them are spurious: the bench references at denser q and k are stable. Except for H3S c4,
they appear only after Fourier interpolation (the dyn files on the q grid are positive). The Nb
and H3S-c4 modes exceed the 30 K noise threshold of theory doc 01 G-VALID item 2. Two
consequences:
1. A coarse imaginary mode is no grounds for rejection or deferral (as doc 11 already says). It
   triggers a T3 phonon check.
2. In this sample, coarse soft modes **raised** A, never lowered it, which is the safe direction
   for deferral. They wreck λ and ω_log. So any coarse rule on λ, ω_log or the NM-G moment rule
   (which needs λ_hi) is unsafe at coarse settings. Coarse m is fine (0.93–0.99 in those runs),
   but NM-G also needs λ_hi.

### 3.4 Other observations

- **k_fit 16³ is too coarse for Al.** At q 4³, no smearing up to 0.05 Ry passes per-mode
  positivity, so the run yields no admitted ladder and cannot screen anything. Al at q 2³ / k_fit
  16³ admits, but underestimates A by up to 35 % per σ. Outside q 2³, no admitted
  per-σ A ratio in this sample is below 0.925.
- **ecut 2/3 (P2) costs little accuracy in A**: relative to P1 at the same q, Pb −4 … −5 %,
  Al +4 … +5 %, Nb −1 … −2 %. It also cuts ph.x time by 22–58 %. λ and ω_log move more.
- **Refining q is not monotone** for λ / ω_log (Nb q 3³ → 4³: λ ratio 0.83 → 1.40). This is as
  doc 11 found for H3S. A is monotone-ish within ±10 %.

## 4. Empirical envelope

Coarse/full ratios (q 3³–4³, k_fit 24³, any ecut ≥ 2/3; n = 14 runs, 4 materials):

| Quantity | per-σ ratio range | Usable as a T2 screen? |
|---|---|---|
| A | 0.925 … 1.106 | yes, via the coarse ENVELOPE (§3.2) |
| m = λ⟨ω²⟩ | 0.926 … 1.054 | robust per σ, but its only consumer (NM-G) needs λ_hi → no |
| λ | 0.80 … 1.70 | no |
| ω_log | 0.14 … 1.28 | no |
| doc-11 statistic max_σ A_coarse / A_hi,full | 0.53 (q3) … ; ≥ 0.69 at q 4³ | only with factor ≥ 2 and q ≥ 4³ |
| coarse ENVELOPE A_hi,c / A_hi,full | ≥ 0.90 | yes, with e_c = 0.5 |

## 5. Recommendation

**Rule T2-SCREEN (replaces doc 11 §4 item 3's "Allowed DEFER rule" and doc 12 T6's
"REJECTED_PROVEN@coarse" outcome):**

```yaml
- name: T2-SCREEN
  kind: heuristic screen (never a proof; outcome SCREENED_OUT, not REJECTED_PROVEN)
  inputs: coarse alpha^2F ladders at q 3^3 AND q 4^3 (k_fit >= 24^3-equivalent density,
          ecut >= 2/3 of the T3 value), mode-positive admission applied to each
  preconditions:
    - admission succeeded on both meshes, >= 3 admitted sigma on q 4^3 (else: no screen-out)
    - no q 2^3 data used
  statistic: A_hi_c = ENVELOPE(A) on the coarse ladder (theory doc 01 §0.4, eps_sys = 0.10,
             eps_omega from the coarse q-spread of omega_log, capped at 0.9)
  test: SCREENED_OUT(T2, criterion="coarse-A envelope", e_c=0.5, calib="doc 13")
        if (1 + e_c) * A_hi_c < A_target   # A_target = 1299.44 K * Tc_target/300 K
  fallback (only one coarse mesh, q >= 4^3): SCREENED_OUT if 2 * max_sigma_adm A_coarse < A_target
  never: screen on lambda, omega_log, m, or on coarse imaginary modes
  audit lane: as doc 12 §3.1 (default 5 % of SCREENED_OUT pushed to T3)
```

Why e_c = 0.5: the worst observed A_hi,c / A_hi,full is 0.90, so (1+e_c)·0.90 = 1.35 ≥ 1 leaves
35 % headroom. The worst per-σ undershoot at q 3³–4³ is −7.5 %. e_c = 0.5 is the value doc 12
already used as "provisional". It is now backed by n = 4 materials, not one. It does **not** come
from a statistical bound: with four materials no defensible quantile exists. It is a judgment of
"about 3× the worst observed shortfall of the envelope (10 %)", and should be re-derived when the
first T3 finalists give hydride references (§6).

Why the doc-11 single-statistic rule needs q ≥ 4³: at q 3³ it reached 0.53 × A_hi,full (H3S),
which a factor 2 only barely covers. At q 4³ its worst is 0.69 (Nb), so factor 2 leaves 38 %
headroom.

**How the coarse envelope feeds the ENVELOPE rule:**
1. The coarse A_hi,c is computed with the *same* ENVELOPE formula as the full tier, from the
   coarse candidate's own ladder. This follows R12 finding 3(b): no special e at the coarse tier
   inside `classify`.
2. e_c is applied **only** outside `classify`, in the screen-out direction, as the factor
   (1+e_c) on A_hi,c (R12 3(c)). It is never passed as `rel_envelope`. That keeps the
   `Spectrum.rel_envelope < 1` constraint irrelevant (R12 noted e_c ≥ 1 cannot be expressed
   there).
3. A coarse spectrum may still be passed to `classify` with the full-tier e for *evidence*. The
   result is stored with `spectrum_tier = coarse` and never drives state (R12 3(a)). The state
   transition is `SCREENED_OUT(T2, "coarse-A envelope", e_c = 0.5, source = doc 13 +
   coarse_run.json SHA-256)`.
4. For candidates that pass, the q 4³ coarse ladder can serve as the "second q mesh" for the T3
   ENVELOPE (δ^q = |X(q_T3) − X(q 4³)|/X_nom). This is conservative: it overstates δ^q relative to
   the second-finest converged mesh. Measured at the reference's nominal σ, δ_A^q(P1 q 4³ vs
   bench q 8³/6³) is 0.5 % (Pb), 2.4 % (Nb), 7.9 % (H3S) and 24 % (Al, at a σ the coarse Al
   ladder does not admit). δ_λ^q at the same σ is 0.4 % (Pb), 18 % (H3S), 20 % (Al) and 37 % (Nb,
   driven by the spurious soft mode). Values near or above 0.5 make G-VALID return UNDECIDED.
   That is the correct outcome, and the T3 run should then add an intermediate mesh.

## 6. Limitations and open problems

- **References are bench-grade, not converged** (R12 W3). The bench meshes (q 8³ / 6³, one mesh
  each) give no δ^q for the reference, so A_hi,full is a lower bound. A larger true A_hi,full would
  make the coarse screen *less* safe than measured.
- **Four materials, all 1–4-atom cubic cells.** Three are 1-atom elements with nearly
  free-electron or d-band Fermi surfaces. Low-symmetry 6–10-atom hydrides have many more
  branches, flatter H modes and more Kohn-anomaly risk. The numbers here do not transfer without
  checking the first T3 finalists against their own coarse runs.
- **MgB2 is missing** (bench ph.x crashed). Two-gap, strongly anisotropic coupling on few q points
  is the case most likely to break a coarse-A rule, and it is not covered here.
- **Pb without SOC** (SG15 scalar-relativistic). The coarse/full ratio is internally consistent,
  but the reference itself is not physical Pb.
- **H3S coarse runs differ from the rest** (ecut 2/3, matdyn 16³, 2 ranks). They come from an
  earlier probe and were not re-run here.
- Wall times were measured under shared load (load average ≈ 12 on 16 logical CPUs) and are
  upper values.
- **X_nom is not defined in theory doc 01.** I took it as the value at the smallest σ admitted on
  both meshes. A different choice (for example the σ-plateau) changes δ^σ and therefore A_hi,c.
- Not verified: whether the Nb q 4³ instability is a real Kohn-type soft spot that the q 8³
  reference misses. It appears only after interpolation, and the dense reference is stable, so I
  treat it as spurious.

## 7. Files

| Path | Content |
|---|---|
| `src/scoring/coarse_dfpt.py` | coarse chain runner (pw template or ASE structure; timings, SHA-256, min frequencies) |
| `src/scoring/coarse_calibration.py` | moments, `admit_smearings` (same logic as the core), per-σ and per-run comparison |
| `src/scoring/coarse_calibration_chain.sh` | the 18-run chain used here (6 ranks, one job at a time) |
| `src/scoring/coarse_envelope.py` | ENVELOPE on coarse ladders vs reference nominal / A_hi lower bound |
| `docs/design/data/13_coarse_calibration_per_sigma.csv` | per run, per σ: admission, A, λ, ω_log, m, ratios |
| `docs/design/data/13_coarse_calibration_per_run.csv` | per run: settings, ph.x wall, min frequencies, ratio ranges |
| `docs/design/data/13_coarse_envelope.csv` | §3.2 statistics |
| `docs/design/data/13_coarse_runs_manifest.csv` | SHA-256 of every a2F.dos, lambda, ph.in, scf_fit.in (coarse and reference) + step wall times |

Reproduce (WSL; coarse workdirs under `~/rtsc_targeted_work/K1-calibrate/runs/`, references
copied to `…/ref/`): `src/scoring/coarse_calibration_chain.sh` (run from that directory) calls `coarse_dfpt.py`
for every profile; then `coarse_calibration.py` and `coarse_envelope.py` as in their docstrings.
