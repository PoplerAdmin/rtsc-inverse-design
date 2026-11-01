# R12. Adversarial review of `docs/design/12_architecture.md` (DRAFT v1)

Reviewer: agent R12-review-architecture, 2026-10-28. Read-only. No compute was run beyond
`wc`/`grep` on source files. Module paths below name the certified core package, which the design
places under `src/rtt/core/` (package `rtt.core`). The package version reviewed includes the
mode-positive smearing admission fix for Pb.

Severity: **critical** = can produce a wrong proven statement or silently breaks a soundness
property; **major** = invalidates a design decision, the budget, or one of the three goals;
**minor** = wrong fact or missing detail with limited impact.

---

## A. Math of L-MOM1..3 (check 1): the inequalities are correct; one proposed use of them is not

Setting (as in the design and in `TcR1.Spectrum`): finite atomic measure, weights p_j ≥ 0 at
Ω_j ∈ (0, ω_max], λ = Σ 2p_j/Ω_j > 0, A = Σ p_j, probability weights w_j = (2/λ)·p_j/Ω_j.

| Claim | Verdict | Proof / remark |
|---|---|---|
| ⟨ω⟩ := Σ w_jΩ_j = 2A/λ | true | Σ (2/λ)(p_j/Ω_j)Ω_j = 2A/λ. Matches `Spectrum.lambda_`, `omega2`, `omega_log` in `rtt/core/spectrum.py` (same moment convention). |
| **L-MOM1** A ≤ λ·ω_max/2 | **true** | A = Σ (p_j/Ω_j)·Ω_j ≤ ω_max·Σ p_j/Ω_j. Equality iff all weight sits at ω_max (Einstein). Needs p_j ≥ 0. |
| **L-MOM2** A ≤ λ·ω₂/2 | **true**, and sharper than L-MOM1 | ⟨ω⟩ ≤ (Σ w_jΩ_j²)^{1/2} = ω₂ (Jensen / power mean), ω₂ ≤ ω_max. Cauchy–Schwarz form A² = (Σ √(p_j/Ω_j)·√(p_jΩ_j))² ≤ (λ/2)(λω₂²/2) is also correct. Equality iff Einstein. |
| **L-MOM3** A ≥ λ·ω_log/2 | **true** (Jensen direction correct) | ln concave ⇒ Σ w_j ln Ω_j ≤ ln Σ w_jΩ_j ⇒ ω_log ≤ ⟨ω⟩. Needs λ > 0 (ω_log undefined otherwise; `Spectrum.omega_log` returns NaN). Equality iff Einstein. |
| 2·a0·T\* = 2598.9 K at T\* = 300 K | **true** | a0 = 433148145/10⁸ exactly (`rtt/core/constants.py`, `A0_FRACTION`); 2·4.33148145·300 = 2598.88887 K = 223.96 meV. A_min = a0·T\* = 1299.444435 K, identical to the constants file's "Corollary (screening filter)". |
| Core bound statement | consistent | `constants.py`: "k_B Tc ≤ A/a0 (strictly: every T ≥ A/a0 is stable)", A = Σ weight_K in kelvin, k_B = ħ = 1; Lean `TcR1.stable_of_ge`/`main` with H(a0, 1e-9) as hypothesis. R-UNIV in `certify/classify.py` tests A_hi = (1+e)·A < a0·T\* exactly in rationals. The design's use of A and a0 matches. |

So L-MOM1 gives a correct necessary condition **for the model's α²F**: Tc_model ≥ T\* ⇒ λ·ω_max ≥
2598.9 K, with λ and ω_max of the *same non-negative spectrum*. The problems are in how the design
proposes to use it (findings 2 and 7).

---

## B. Findings

### 1. [critical] The input-quality downgrade layer lives in `pipeline.py`, which the design replaces

Evidence. The design keeps `certify/`, `spectrum.py`, `status.py`, `constants.py`, `validity/`
in the core (§0, §3.5 `core/`) and **replaces** `pipeline.py` (§0, §1.1). But the post-classify
policy that makes a proven status safe for real QE input is the `pipeline.py` certify
post-processing (now `rtt/core/certify_policy.py`):
`certify_spectra` (≈ l. 443–465) calls `classify`, then `ensure_cutoff_assumption`, `annotate_model`
(adds M-ISO-SINGLEBAND, DISCRETISATION, the F9 Coulomb-block check), `trust_base` (python-flint
version pin), and `apply_input_policy` → `blocking_tags` / `_sigma_exclusion_directions`
(l. 484–585). These **downgrade** REJECTED_PROVEN to UNDECIDED when an input carries
`A-FILE-INCOMPLETE`, `A-FILE-COMPLETE(unverified)`, `A-DROP-NONPOS(m>0)` (missing mass lowers A,
anti-conservative for an upper bound), or an `A-SIGMA-EXCLUDED` tag whose excluded λ or A lies above
the admitted envelope. The Pb smearing-admission fix states this explicitly ("pipeline: the tag
blocks a verdict only when an excluded smearing's λ or A falls outside the envelope-widened
range"). So the Pb fix is half in `dfpt/workflow.py` (core) and half in `pipeline.py` (dropped).
§1.1's statement that the fix is carried by keeping `dfpt/` + `certify/` in the core is wrong.

Consequence. A new `rtt.certify` wrapper written from the module list in §3.5 issues
REJECTED_PROVEN for a Pb-like ladder where an excluded smearing has larger A than every admitted
one, or for a truncated a2F file. That is exactly the failure the core pipeline was hardened against.
R1 ("only core classify produces proven statuses") does not catch it, because the status
*is* produced by classify; the guarantee depends on the downgrade that follows.

Required change. Keep `certify_spectra` and its helpers (`ensure_cutoff_assumption`,
`annotate_model`, `trust_base`, `apply_input_policy`, `blocking_tags`,
`_sigma_exclusion_directions`, `validate_certify_config`, `_certify_params`) as a separate
hash-pinned module (e.g. `core/certify_policy.py`, source lines recorded in
`CORE.txt`), together with the tests that exercise them (`test_dfpt_smearing_admission.py`,
the input-policy parts of the soundness suites). Make `rtt.certify` a pure call to that function.
Add an explicit contract R1b: "every proven status passes through the core input policy".

### 2. [critical] R-LAMW is proposed for exactly the inputs on which L-MOM1 is false

Evidence. §1.3: "A natural source is a coarse DFPT run (T6) when its α²F **fails admission** but its λ
and ω_max are well defined". Admission (`admit_smearings`, mode_positive) fails precisely when
per-mode α²F has negative weight. L-MOM1 needs p_j ≥ 0 (and the Lean `TcR1.Spectrum` has
`p_nonneg`). Counterexample with signed weights: p₁ = −0.01 K at Ω₁ = 10 K, p₂ = 1 K at
Ω₂ = ω_max = 1000 K. Then λ = 2(−0.001 + 0.001) = 0, A = 0.99 K, and A > λ·ω_max/2 = 0. Small
negative mass at low ω (W2: Pb/Al negative mass sits at ω ≤ 15–31 K) depresses λ strongly while
barely changing A, so L-MOM1 fails in the anti-conservative direction for rejection. Second
problem: a coarse-run λ is not a *bound* on the full or converged λ ("rigorous pair (λ_hi,
ω_max,hi)" has no source anywhere in the pipeline). EPW per-mode λ_qν are also discretised
estimates, not bounds.

Consequence. R-LAMW is either unsound (with failed-admission or coarse inputs) or redundant
(with an admitted, non-negative spectrum, A is available directly and R-UNIV is strictly sharper
than R-LAMW, since A ≤ λω_max/2 always). There is no input for which it adds a sound rejection.

Required change. Delete R-LAMW and the planned `classify` extension. Keep L-MOM1–3 as Lean lemmas
documenting the target derivation only (they are useful for docs/theory and for proxy thresholds).
State the hypotheses in Lean explicitly: `∀ j ∈ S, Ω j ≤ ωmax` (the TcR1 structure has no ω_max
field) and `0 < λ` for L-MOM3. If anything uses ω_max with a QE spectrum, it must be the maximum
of the spectrum actually passed to classify (grid points, including the broadening tail of
`a2F.dos`), not the phonon ω_max.

### 3. [major] The "@coarse" verdict is used as a funnel cut, which makes it a heuristic dressed as a proof

Evidence. §3.1 table, T6: "Outcome on failure: REJECTED_PROVEN@coarse", with the margin rule
(1+e_c)·max_σ A_coarse < A_min, e_c = 0.5 "provisional"; §4.3 assumes a T6 pass fraction of
0.1–0.3 on "A_coarse within margin of A_min", i.e. T6 *does* remove candidates. A REJECTED_PROVEN
on the coarse spectrum proves nothing about the full or converged spectrum; the removal is valid
only to the extent the empirical e_c is. That is the definition of a heuristic screen (R2). If
e_c is passed into classify as `rel_envelope`, the verdict silently becomes "proved, assuming an
empirical envelope" (`E-ENV(e)` is recorded, but it then looks the same as a DFT-error envelope).
Also: `Spectrum` requires rel_envelope ∈ [0, 1) (`spectrum.py` `__post_init__`), so a calibrated
e_c ≥ 1 cannot be expressed at all; and the envelope is weights-only (ERRATA F5), so it is
meaningful for R-UNIV (depends on A only) but not for the certification direction.

Required change. (a) The T6 *tier outcome* that removes a candidate is `SCREENED_OUT(T6,
criterion = "coarse-A margin", e_c, source id)`; the coarse certify verdict is stored as evidence in
`verdicts` with spectrum_tier = coarse and never drives state. (b) Call classify on coarse spectra
with the same e as at the full tier, never with e_c. (c) e_c is only ever used in the rejection
direction. (d) Report the calibration honestly: the planned calibration set (Al, Pb, Nb, H3S;
MgB2 crashes) is four 1–4-atom cells, and the "full" reference values are themselves
unconverged (W3), so e_c relates coarse to bench-grade, not to converged, and is not
transferable to 6–10-atom hydrides without validation on the first T7 finalists.

### 4. [major] Nothing binds `spectrum_tier` or the spectrum's provenance to a proven status; the origin guard is a self-declared string

Evidence. The earlier ledger guard (`discovery/ledger.py` l. 123–126, 312) checks
`origin = 'rtt.core.certify.classify'`, a string the caller writes. The core `Verdict` has no
`spectrum_tier` field; the design adds it in the wrapper outside the core (§3.5 "spectrum_tier
tagging"). Nothing prevents (i) a caller tagging a coarse run "full", (ii) a Spectrum built from
ML or proxy data (e.g. a synthetic Einstein spectrum from λ_ML, ω_log^ML) being passed to classify
and stored as REJECTED_PROVEN with a valid origin string. The `classify` result is a proof about
whatever list of Spectrum objects it receives.

Required change. In `verdicts`: `spectrum_tier NOT NULL CHECK IN ('coarse','full','converged')`,
plus a NOT NULL foreign key to the `tier_runs` row (kind = DFPT, profile) and to `artifacts` rows
whose sha256 equals the files the spectra were parsed from; a trigger refuses a proven status
unless that join exists and the tier_run profile maps to the stated spectrum_tier. Store
`spectrum_digest` per alternative (classify already computes it) and run the core `recheck`
before insert. Report queries (T9) select by spectrum_tier ∈ {full, converged} and by the join,
not by status alone; summary counts of REJECTED_PROVEN must be split by tier.

### 5. [major] The audit lane cannot run at the stated rate and would not measure what it claims

Evidence (arithmetic from §3.1 and §4.3). T3 screens out ≈ 2 250–2 375 of 2 500 per day; 5 % is
≈ 115/day pushed to T4, whose total capacity is ≈ 40/day. At T6 → T7, 5 % of ≈ 3 rejects is 0.15
candidates/day × 100–500 core-h = 15–75 core-h/day, i.e. 15–80 % of lane F. "Budget permitting"
therefore means the T3 audit almost never runs, or the audits displace the main funnel. The
125–250 T3 PASS candidates/day that miss the "top 40" cut are neither screened out nor audited
(their outcome is undefined in §3.3's vocabulary; DEFERRED is the only fit).
Statistics: (i) pushing a reject *one tier further* measures disagreement between two
heuristics, not the false-negative rate with respect to the target (A ≥ A_min on a full spectrum);
the composed FN rate needs audits that reach T7. (ii) The design itself expects ≈ 0 true positives
at 300 K (§4.3 "Sanity anchor"). With zero positives in the sample, a binary FN rate is 0/n and
uninformative; even with a non-zero positive rate p, bounding the FN rate needs on the order of
3/p audited items per tier (rule of three) — far beyond capacity.

Required change. (a) Audit a fixed *number* per tier per week sized to downstream capacity, not a
fraction. (b) Measure the proxy's *continuous* error (e.g. λ_ML·ω vs DFT λ·⟨ω⟩; ML ω_max vs DFT)
on audited items, oversampling near the threshold, with inverse-probability weights recorded so
the estimate is unbiased. (c) Use a surrogate target with real positives for recall measurement
(e.g. A ≥ a0·T for T = 40–100 K, where known 0-GPa superconductors exist) and a positive-control
set of known superconductors (MgB2, Nb, NbN/Nb3Sn-class, known hydride/boride/carbide metals) run
through T1–T6 to show the screens do not drop them. (d) Define the outcome of "PASS but not
selected" (DEFERRED with rank and cut value).

### 6. [major] Per-tier costs are extrapolated optimistically; the funnel output (2–6 full DFPT/week) is not supported

Evidence.
- T7: the only H-containing anchor is H3S Im-3m (4 atoms, 48 point operations, 16 irreducible q,
  ≈ 16 core-h). For a 6–10-atom cell of "moderate symmetry" the design's own model
  (cost ∝ N^{3–4} × Σ_q N_k(q)) gives ×8–16 from N (4 → 8) and roughly ×5–10 from symmetry
  (more irreducible q, each with a larger IBZ of k): ≈ 600–2 500 core-h (est., same ×3 caveat as
  §5), not 100–500. At 96 core-h/day that is 0.3–1 T7 per **week**, before T8 (also lane F, 1–3×
  T7) and before the audit share (finding 5). The read-out "about 2–6 per week" holds only for
  ≤ 4-atom high-symmetry cells, which §4.3 explicitly does not assume.
- T5/T6 use per-q anchors of 10–46 s on 2 ranks for **1–4-atom elemental/binary cells at
  high-symmetry q** (step2 report §4). Scaling N^{3–4} from MgB2 (3 atoms) to 12 atoms gives ×64–256:
  0.026 core-h × 256 ≈ 6.6 core-h for Γ alone, above the stated 0.2–2; T6 with 4–30 irreducible q
  plus lower-symmetry k sets reaches tens to > 100 core-h, versus the stated 1–10. Not binding in
  the funnel as drawn (arrivals are only 3–4/day), but T6 cannot be the "cheap decisive" tier the
  design wants at that cost.
- T3: 2 500 arrivals/day at the stated 7–60 s each on 2 CPU threads is 4.9–42 CPU-hours of a
  single 2-thread worker, i.e. 2 880/day at 30 s and 1 440/day at 60 s (6–10-atom cells with
  ≤ 160-atom supercells are at the upper end). The "≥ 3 000" capacity needs the GPU or more threads.
- Lanes: S (4 cores) + F (4 cores) use all 8 physical cores; lane G adds 1–2 more CPU threads plus
  MatterGen's host work, and CONTEXT allows earlier validation re-certifications to run concurrently.
  W5 itself measures ~15 % per-q slowdown under this kind of oversubscription; 192 core-h/day is an
  upper bound, not a planning number.

Required change. Re-derive §4.2/§4.3 with the N and symmetry factors applied explicitly (show the
formula per tier), give the funnel for two cell classes (≤ 4 atoms cubic/hexagonal; 6–10 atoms
moderate symmetry), budget T8, audits and calibration out of lane F, and run T3 on the GPU. Make
the first-week measurement of one 6–8-atom hydride at the coarse and full tiers a prerequisite for
any throughput claim. Expect the honest number for 6–10-atom cells to be below one full DFPT per
week.

### 7. [major] The L-MOM1 screen at T3/T5 is barely decisive for the intended chemistry

Evidence. λ_needed = 2598.9 K / ω_max. ω_max of any compound with covalent X–H bonds is set by the
stretch modes (C–H, B–H, N–H ≈ 2 900–3 400 cm⁻¹ ≈ 4 200–4 900 K; UNVERIFIED typical values),
giving λ_needed ≈ 0.55–0.65, which almost every hydride "passes" while those modes contribute little
to λ. Doc 02 (`docs/theory/02_materials_profile.md` §5 row 9) asks for ω_max(H) ≥ 1 500 K, which
by L-MOM1 corresponds to λ ≥ 1.73 *at most*. ⟨ω⟩ is typically a fraction of ω_max, so L-MOM1 is
loose by that factor. The assumed T3 pass fraction 0.05–0.1 therefore cannot come from the
"theory-anchored" test; it would come from E_hull and the ranking cut.

Required change. State in §3.1 that L-MOM1 at T3/T5 is weak and why. Use the sharper quantities
in the proxy screens: ⟨ω⟩ = 2A/λ or ω₂ (L-MOM2) estimated from the ML phonon DOS weighted toward
the light-element modes, i.e. a screen on λ_needed = 2·A_min/ω₂^proxy, still heuristic. Ask doc
01/10 for a calibrated λ·⟨ω⟩ (= 2A) predictor, since A is the quantity the decisive test uses.

### 8. [major] DFT E_hull at T4 has no reference set and no budget

Evidence. T4 computes "DFT E_hull" (X-HULL, an EXCLUDED_ASSUMPTION). A hull is only meaningful
against competing phases computed with the same functional, pseudopotentials (SG15 NC, forced by
W4) and cutoffs. MP energies (VASP PAW, with corrections) cannot be mixed in. For a ternary such
as Li–B–H the reference set is every MP-stable phase of the system (elements, binaries, ternaries),
some with large cells (α-B 12 atoms; β-B ≈ 105 atoms; LiBH4 24 atoms in the conventional cell),
at 80–120 Ry NC cutoffs. None of this appears in §4.2/§4.3. `mlff/hull.py` has the same issue for
ML (the design addresses that one).

Required change. Add a reference-phase sub-tier per chemical system (cached, versioned by
pseudopotential md5 + cutoff + k density), with its cost; restrict T1 chemical systems to those
whose references are affordable; record the reference set sha in each T4 record. Until it exists,
X-HULL is a heuristic (SCREENED_OUT), not EXCLUDED_ASSUMPTION.

### 9. [major] The T1 conditioning in the design is not supported by any released MatterGen checkpoint, and nothing feeds e–ph targets back into generation (goal 2)

Evidence. §3.3 example: `"conditioning": {"chemical_system": ..., "dft_band_gap": 0.0,
"space_group": 225}`. The earlier adapter's own table (`generate/mattergen.py`
l. 64–75, `CHECKPOINT_PROPERTIES`) lists single-property checkpoints plus two pairs
(`dft_mag_density_hhi_score`, `chemical_system_energy_above_hull`); no checkpoint conditions on
chemical_system + band gap + space group, and the adapter raises `UnsupportedConditioningError` for
such a request (l. 154–190). No released property relates to phonons or e–ph coupling. The
fine-tuning path (`finetune_export.py`) is deferred "until hundreds of DFT labels" exist; with
2–6 (realistically < 1, finding 6) α²F per week, e–ph labels will never reach that. So the
pipeline's link from "theory targets" to "meaningful MatterGen output" is only: chemical system,
one other property, and post-hoc filtering. The bandit was removed (W8) but no replacement
learning signal is specified.

Required change. (a) Correct the example and §3.6 contract to the real checkpoint/property
combinations (or plan a multi-property fine-tune and budget its GPU time). (b) Specify the feedback
mechanism explicitly: a cheap label available at scale (e.g. T3 ML light-element phonon weight /
ω₂^proxy, T4 DFT H-projected DOS at E_F at ~40/day) used to fine-tune or to steer chemical-system
sampling, with the number of labels per week and when the first fine-tune happens. (c) Include a
database lane (known MP/OQMD/ICSD metals in the target chemistry) as a baseline: these are already
"concrete compounds", cheaper to justify than MatterGen hypotheticals, and they tell whether
generation adds anything.

### 10. [major] Final ranking and evidence statement are underspecified for goal 3

Evidence. §4.3: "best-supported candidates will be the candidates with the largest certified
T_up / largest A". T_up is an upper bound; ranking by it (and by max_σ A) rewards wide brackets and
large envelopes. The evidence statement (T9) does not include: anharmonic/quantum-nuclear
stability of H sublattices at 0 GPa (doc 02 row 9 explicitly asks for SSCHA for H; the design only
carries the ANHARMONIC_H flag), DFT phonon stability on the full q grid of the *final* structure,
zero-point energy in the hull, sensitivity to μ\*, and synthesis route/metastability window.

Required change. Rank finalists by a conservative quantity (min over admitted alternatives of A,
or certified T_low, or an Eliashberg Tc at fixed μ\* = 0.10–0.13 with the spread across
smearings), report T_up only as the ceiling. Add an optional T8b "H-anharmonic check" (SSCHA or at
least frozen-phonon anharmonicity of the dominant H modes, plus ZPE-corrected hull) for any
finalist whose λ is carried by H modes, and list these items in the T9 template as checked /
not checked.

### 11. [major] Stage ordering: a seconds-cost heuristic sits behind a core-hour tier, and the Γ tier is weaker than stated

Evidence. The ML λ/ω_log predictor (seconds, heuristic) is in T5, after T4 vc-relax (0.3–3 core-h,
the stated binding tier at ~40/day). It does not need DFT inputs. Also T5's "e-ph at q = Γ": the
q = Γ e–ph contribution is not an estimate of λ (acoustic γ vanish, intraband double-delta at q = 0
is ill-defined at finite smearing); only the Γ dynamical matrix is informative. Γ-only stability
misses zone-boundary/CDW instabilities, which are the typical failure of strongly coupled metals;
"catches many unstable structures" is not supported by any measurement (the §4.4 item 2 comparison
is the right test, but the ordering is decided before it).

Required change. Move the ML λ (and any ML A/⟨ω⟩) proxy into T3, ahead of T4, so it reduces T4
arrivals. Drop e–ph at Γ from T5 (keep the Γ phonon check, or fold it into T6 if §4.4 item 2 shows
Γ catches few instabilities). Order T3 internally: relax → ML hull → ML proxies → ML phonons (the
most expensive part) only for survivors.

### 12. [minor] Wrong or unverifiable inventory facts

- "7 280 LOC of package code": `wc -l` over the package's `**/*.py` gives **13 984**; tests
  give 7 791; 13 984 + 7 791 = 21 775, so the "21 775 including tests" total is right and the
  package figure is wrong.
- Lean LOC: TcR1 529, PipelineLogic 1 414, MultiBand 1 175 confirmed; DiscoveryLogic is 1 406
  (not 1 402).
- Per-file LOC in the §1.1 table spot-checked (classify 525, workflow 757, parsers 549, mattergen
  341, pipeline 1 150, ledger 459, loop 830, hull 183, phonons 182, mattersim 249): all correct.
  254 `def test_` confirmed. `admit_smearings`/`mode_positive` present in `dfpt/workflow.py`
  l. 511–612: confirmed.
- T3 anchor "2 CPU threads": `reports/step3_mattersim_prescreen.md` §2 gives "~1.75 cores, 2.1 GB
  RSS" but does not state the device for the Mg-B run (the report's fix concerns CUDA execution).
  Mark the device as UNVERIFIED or re-measure.
- T7 anchors Al 7.4 core-h and Pb ≈ 16 core-h agree with PERFORMANCE_DRAFT (3 328 s × 7.89,
  8 268 s × 6.84).

### 13. [minor] L-MOM Lean effort and statement details

L-MOM1 is a few lines over `Finset.sum_le_sum`. L-MOM3 needs weighted Jensen for `Real.log`
(concavity is in Mathlib, but the weighted finite form with the λ normalisation, exp/log
round-trip and λ > 0 side condition is more like 50–150 lines, est.). The TcR1 `Spectrum` has no
λ, ω_max, ω_log: define them in `TargetLogic`, with ω_max as a hypothesis bound rather than a
`Finset.sup'` (avoids the non-empty side condition). The design's sentence "λ·ω_log ≥ 2·a0·T\* is a
sufficient condition for passing the A-necessary condition" is correct but should add "and
nothing more": it implies A ≥ A_min, not Tc ≥ T\*.

### 14. [minor] Shared relax between tiers

§3.1/§1.1 plan "one relaxed structure and one scf" shared by T4–T7. A vc-relax at T4 settings (lower
k, possibly lower cutoff) leaves residual stress/forces at T7 settings, which shows up as spurious
soft or imaginary modes in DFPT. Require a re-relax (or a force/stress check with a threshold) at
the T7 cutoff and k grid before T7.

### 15. [minor] Calibration item cost and scope (§4.4)

"T6 on Al, Pb, Nb, H3S, MgB2 … 1–2 core-h total": plausible for these cells, but MgB2 ph.x crashes
(W4), leaving four materials, all ≤ 4 atoms. Add at least one 6–8-atom metal with light-element
modes before trusting e_c or the Γ-vs-grid agreement for the target class (see finding 3, 6).

---

## C. Summary of required changes (priority order)

1. Keep the `pipeline.py` certify post-processing in the core (input policy, annotations) with its tests;
   add contract R1b (finding 1).
2. Delete R-LAMW; keep L-MOM1–3 as documentation lemmas with explicit hypotheses (2, 13).
3. T6 removal = SCREENED_OUT with e_c; coarse verdicts are evidence only; never pass e_c as
   rel_envelope (3).
4. Bind proven statuses to DFPT tier_run + artifact hashes + recheck; spectrum_tier NOT NULL (4).
5. Redesign audits (fixed counts, continuous error, surrogate threshold, positive controls) (5).
6. Redo the budget with explicit N and symmetry scaling; expect < 1 full DFPT/week for 6–10 atoms (6).
7. Replace the L-MOM1 proxy with ⟨ω⟩/ω₂-based proxies; move ML λ to T3 (7, 11).
8. Add DFT hull reference sub-tier, real MatterGen conditioning combos, a feedback/label plan, a
   database lane, conservative ranking, and an H-anharmonicity check (8, 9, 10).

## D. What I could not verify

- No compute was run; all scaling estimates in findings 6 and 11 are mine (est.), from the
  design's own cost law and the step2/PERFORMANCE anchors.
- Whether newer MatterGen releases ship multi-property checkpoints beyond the adapter's table
  (the table is pinned to MatterGen commit 92423660) was not checked on-line.
- Typical X–H stretch frequencies in finding 7 are recalled (UNVERIFIED).
- The bench-2026-10-28 numbers for Nb/H3S quoted by the design were not re-read from the bench
  directory (only the reports under `reports/`).
