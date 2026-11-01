# R01 — Adversarial review of `docs/theory/01_eliashberg_targets.md`

Reviewer: agent R1-review-eliashberg. Scope: the document, `src/theory/*.py` and `docs/theory/data/`.
This is a review. I did not edit the target document or its code.

Verification scripts are in `docs/reviews/r01_scripts/`:

| script | what it does |
|---|---|
| `indep_solver.py` | Independent linearized ME solver. It uses the full Matsubara axis n ∈ [−N, N−1], keeps the λ(0) terms in Z and in the kernel (no hand cancellation), uses a non-symmetric eigenproblem, and lets you pick the μ* cutoff convention. It shares no code with `src/theory`. |
| `spot_checks.py` | Re-derives table values with the independent solver. |
| `adversarial_opt.py` | Differential-evolution search over 2–4-atom spectra for counterexamples to N1μ, N2 and N3. |
| `sanity_n3.py` | Checks that the N3 objective used by the optimizer is non-trivial. |
| `reentrance_check.py` | Tests the step "Tc ≥ 300 K ⇒ A ≥ 300·a_E(300 K)" across the discontinuities of the block cutoff N_c(T). |

Everything ran in WSL with nice 19 and 1–2 threads; the commands are in each file's docstring. I re-ran
`validate_solver.py`, and its output matches `data/validation.txt`.

Tags: **critical** means a stated hard reject rule can discard a real Tc ≥ 300 K material, or a claim
is wrong in a way that changes the pipeline. **major** means the rule's direction or scope is
under-specified, or a proof needs repair. **minor** covers citations, typos and presentation.

---

## Summary

The solver is correct. The headline numbers reproduce independently: a* = 4.33148, the Table 2, 4
and 6 entries, and C∞. The extreme-point reductions behind N1μ, N2 and N3 are sound. I found no
counterexample to N1μ, N2 or N3 by continuous optimisation (§F10).

The problems are in how the results are **handed to the pipeline as hard reject rules**:

- **F1 (critical).** The floors are theorems for a *single isotropic* α²F. A DFPT run gives the
  Fermi-surface-averaged α²F. For anisotropic or multiband materials, N2 and N3 evaluated on that
  average are not necessary conditions; a simple two-band case breaks them by a factor of 2. N1 on
  the averaged A is also invalid, and the doc's κ fix needs data the pipeline will usually not have.
- **F2/F3 (major).** The "1.2 × nominal" uncertainty recipe is ill-defined for ω₂, and the
  "λ credibly ≤ 3" gate is stated in the unsafe direction. A unified, direction-safe form exists:
  λ⟨ω²⟩_hi < g(λ_hi).
- **F4 (major).** The reject floor is not gated on the validity of the model: non-adiabatic,
  low-E_F or anharmonic candidates.
- **F5 (minor).** The Gao et al. article number is wrong (it is 8253, not 8256).
- **Further minor findings.** One misattributed Esterlis claim, a typo in the §0 summary row, and the
  KAY journal references.

---

## Findings

### F1 [critical] N2 and N3 (and N1 with A) are not necessary conditions for the band-averaged α²F that DFPT produces

**Evidence.** §5.1 of the doc already concedes that the band-averaged A is not a valid bound. It says
Tc·a0/A_avg is unbounded and gives a +9 % MgB2-like model. It also says that N2 and N3 "have no
multiband proof yet". Even so, §0.1 lists "λ⟨ω²⟩_hi < 2.7×10⁶ K² ⇒ reject" and "ω₂,hi < 1075 K ⇒
reject" with no multiband qualifier. The κ prescription for A applies only "with multiband data".

A standard QE/EPW isotropic run gives exactly the averaged α²F = Σ_ij (N_i/N) α²F_ij.

**Explicit counterexample.** This is the decoupled limit, which continuity extends to weak
interband coupling.

1. Take two bands with equal DOS. Band 1 has an Einstein α²F_11 with λ = 3 and Ω = 1075.3 K. Band 2
   is uncoupled.
2. Band 1 alone gives Tc = 300.0 K at μ* = 0. This is the doc's own Table 4 entry, which I
   reproduced with the independent solver (`spot_checks.py` item 3).
3. The averaged α²F has λ = 1.5, ω₂ = 1075 K, λ⟨ω²⟩ = 1.73×10⁶ K² and A = 806 K.

All three floors therefore reject a 300 K superconductor:

- N2: 1.73×10⁶ < 2.70×10⁶.
- N3: λ = 1.5 ≤ 3, and ω₂ = 1075 K is below Ω_E(1.5) = 1774 K.
- N1 with A: 806 < 1299.

With a DOS fraction w for the coupled band, λ⟨ω²⟩ and A are under-estimated by the factor w, without
limit. Only κ (here κ = A_11 = 1613 K) gives the right answer.

The realistic size of the effect is smaller than this toy case: the doc's MgB2-like model gives
+9 % in A. It is not zero for multigap hydrides or borocarbides, and the floors have no margin at all.

**Required change.**

1. Restrict N2 and N3 explicitly to "single isotropic α²F (one-band, or proven weakly anisotropic)".
   Until a multiband version exists, either drop them from the hard-reject list or apply them only
   with an explicit anisotropy margin.
2. For N1 without band-resolved data, use the rule "reject if (1+δ_an)·A_hi < 1299.4 K", where δ_an
   is stated and justified, for example δ_an ≥ 0.15–0.25 when the Fermi surface has several sheets
   with different character, or when the gap anisotropy from an anisotropic run is not small.
3. Alternatively, require κ, which needs the band-resolved α²F_ij, before any hard rejection within
   ~25 % of the floor.
4. Add the decoupled two-band case above to §5.1 as the reason.

### F2 [major] The "1.2 × nominal" uncertainty envelope is ill-defined for ω₂ and too coarse for A and λ⟨ω²⟩

**Evidence.** §0.1 says to apply the floors to "1.2 × the nominal value for a DFPT α²F with ≈20 %
λ error". This has three problems:

1. **ω₂ = sqrt(λ⟨ω²⟩/λ) does not change when α²F is multiplied by a constant.** "ω₂,hi = 1.2 ×
   nominal" therefore corresponds to no error model. ω₂'s uncertainty comes from frequency errors,
   such as anharmonic shifts of H modes or the q-mesh. It does not come from the λ amplitude.
2. **The three moments scale differently with the two physically distinct errors.** For an
   amplitude error (1+ε_a) and a frequency stretch (1+ε_ω) at fixed electronic matrix elements
   (λ ∝ 1/ω² per mode):
   - A ∝ (1+ε_a)/(1+ε_ω);
   - λ⟨ω²⟩ ∝ (1+ε_a), with no frequency dependence at fixed g².

   So frequency errors do not enter λ⟨ω²⟩ at all, and enter A with the opposite sign from a naive
   "stretch".
3. **The 20 % figure is asserted.** In hydrides, λ from DFPT routinely changes by 20–50 % with the
   smearing/k-mesh and with anharmonicity. Errea et al. (2015, 2020), which the doc cites, report
   tens of percent.

**Required change.**

1. State the error model as separate amplitude and frequency envelopes, taken from the candidate's
   own convergence data (smearing sweep, q-mesh) and an anharmonic flag.
2. Propagate them to A_hi, λ⟨ω²⟩_hi and λ_hi, ω₂,hi separately.
3. Make the default ε_a at least the spread seen in the smearing sweep, not a fixed 20 %.

### F3 [major] The N3 gate "λ credibly ≤ 3" is in the unsafe direction; use the unified form λ⟨ω²⟩ < g(λ_hi)

**Evidence.** N3 needs the *true* λ ≤ 3. If DFPT under-estimates λ (nominal 2.8, true 3.6), the
gate fires but N3 does not apply.

The direction-safe statement follows from the doc's own results. Define
g(λ) := λ·Ω_E(λ, μ*=0)², the minimum λ⟨ω²⟩ at fixed λ (§3.3). From `two_atom_fixed_lambda.csv`:

| λ | g(λ) (K²) |
|---|---|
| 1 | 6.85×10⁶ |
| 1.5 | 4.72×10⁶ |
| 2 | 4.01×10⁶ |
| 2.5 | 3.67×10⁶ |
| 3 | 3.47×10⁶ |
| 4 | 3.24×10⁶ |
| ∞ | 2.70×10⁶ (= N2) |

g is decreasing in λ. Then:

- **Theorem-level statement (modulo the monotonicity of g):** Tc ≥ 300 K ⇒ λ⟨ω²⟩ ≥ g(λ) ≥ g(λ_hi)
  for any λ_hi ≥ λ_true.
- **The rule:** reject if λ⟨ω²⟩_hi < g(λ_hi).

This rule:

- uses the upper envelope for both quantities, so it is safe in both directions;
- reduces to N2 when λ_hi is unknown (λ_hi = ∞);
- is never weaker than "ω₂ < 1075 K with λ ≤ 3" when λ is actually ≤ 3;
- is robust to a spurious soft mode that inflates λ. The soft mode lowers ω₂, which under the
  fixed-1075 K form could cause a wrong rejection. Here it only lowers g, which loosens the test.

One more point on the monotonicity steps:

- KAY III proves that f(λ) = Tc/Ω is strictly increasing for λ ≥ λ* ≈ 0.767 [CIT, arXiv:2409.02121
  v2]. That makes "Ω_E(λ) decreasing", the step N3 uses, a published theorem for μ* = 0 and λ ≥ 0.77.
  The doc lists it only as numerics.
- g decreasing is equivalent to Tc/(Ω√λ) increasing. That is slightly *stronger* than KAY
  Conjecture 1 and is not proven there. KAY III does not state it; it proves only the monotonicity of
  Tc/Ω. Numerically it holds (validation §2).

**Required change.**

1. Replace the N3 row and the "λ credibly ≤ 3" bullet with the tabulated g(λ) and the rule
   "λ⟨ω²⟩_hi < g(λ_hi) ⇒ reject". Interpolate in λ only downward, i.e. use the g value of the next
   *larger* tabulated λ.
2. Cite KAY III for the monotonicity of Tc/Ω.
3. State that the monotonicity of g is a separate numerical assumption.

### F4 [major] The hard-reject floor is not gated on the validity of the model

**Evidence.**

- §5.3 asks the materials layer to flag λΩ/E_F > 0.1, and §5.4 asks for a harmonic/anharmonic flag.
  §0.1 does not say what happens to a *flagged* candidate.
- N1–N3 are theorems of isotropic, adiabatic ME with μ* ≥ 0. Vertex corrections in the
  non-adiabatic regime (Pietronero–Strässler–Grimaldi 1995, cited) are argued in that literature to
  *raise* Tc. Spin-fluctuation or plasmon contributions are also outside the model.
- A "hard reject" applied outside the model's domain is not a theorem.

**Required change.** §0.1 should say that hard rejection by N1–N3 is valid only when all of these
hold:

- λΩ_max/E_F ≤ 0.1, or a bandwidth criterion;
- the α²F is from a dynamically stable structure (harmonic, or anharmonic when flagged);
- the material is single-band or κ is used (F1).

Otherwise the candidate should be "UNDECIDED/out-of-model", not rejected.

### F5 [minor] Gao et al. article number, and two paraphrase issues

**Evidence.** arXiv:2502.18281 gives the journal reference **Nat Commun 16, 8253 (2025)**, DOI
10.1038/s41467-025-63702-w. The doc (§5.5) has 8256, which is wrong. The sibling
`02_materials_profile.md` has the correct 8253.

The paraphrases:

1. The paper says the hypothetical compounds with the highest λ are "≈3.3" (ClB2C8 and Al2OsH7), but
   it also says Li₂AgH₆ has λ ≈ 4. The source is internally inconsistent, and the doc repeats both
   without noting it.
2. "Practical ceiling about 100–120 K" is stronger than the source. The source says one can design
   hypothetical compounds reaching 100–120 K, and that higher Tc is not forbidden but realisation is
   "extremely unlikely".

The other numbers check out against the arXiv HTML:

- ω_log rarely above 1800 K;
- Li₂AgH₆: Tc 108.8 K, E_hull 0.319 eV/atom;
- Li₂AuH₆: 116.1 K, 0.172 eV/atom.

**Required change.** Change 8256 to 8253, note the λ ≈ 3.3 vs 4 inconsistency, and soften
"ceiling".

### F6 [minor] Esterlis citations: one claim misattributed, one coupling mislabelled

**Evidence.**

1. **PRB 99, 174516 (2019)** ("Pseudogap crossover in the electron-phonon system", arXiv:1811.10622).
   Its abstract does not state "breakdown at λ ∼ 1, with Tc peaking near that value". It describes a
   metal → pseudogap crossover, SC at small λ₀ and CDW at large λ₀, with no numerical threshold. The
   "Tc peaks near λ ∼ 1" statement probably comes from Esterlis, Kivelson & Scalapino, npj Quantum
   Materials 3, 59 (2018), "A bound on the superconducting transition temperature". I did not check
   that paper; it is UNVERIFIED by me.
2. **PRB 97, 140501(R) (2018)** is correctly cited. The coupling ≈ 0.4 is the *bare* dimensionless
   λ₀ of the 2D Holstein model, not the renormalised λ used everywhere else in the doc. The doc
   writes "λ ≈ 0.4" without saying so.

**Required change.** Re-attribute or mark claim 1 UNVERIFIED, and write λ₀ (bare) for claim 2.

### F7 [minor] KAY citations are accurate, but the journal references are missing or incomplete

**Evidence.** Checked against arXiv:2409.02121 v2 and arXiv:2409.00532:

- C∞ = (1/2π)·𝔤(2)^{1/2} = 0.1827262477….
- The rigorous upper bound is ≈ 2.034 C∞.
- KAY III Conjecture 1 is "Tc(λ, Ω) ≤ (Ω/2π)·sqrt(𝔤(2)λ) for all λ, Ω". This matches the doc's use.
- Corollary 1, Tc ≥ (Ω/2π)·sqrt(λ−1) for λ > 1, is attributed to AD as well. This matches the doc.
- The doc's c₂ = 1/(2π sqrt(B_c)) implies 𝔤(2) = 1/B_c = 1.31814, which is consistent.

Both papers are J. Stat. Phys. (2025):

- Part III: DOI 10.1007/s10955-025-03469-y, art. 93;
- Part II: DOI 10.1007/s10955-025-03468-z, art. 94.

The volume came back as "92" from the abstract page, which is probably 192; check it. The doc cites
Part III only as an arXiv number.

The KAY results are for μ* = 0. The doc correctly applies them only there.

**Required change.** Add the journal references, and add the proven monotonicity of Tc/Ω (F3).

### F8 [minor] Typo in the §0 summary row for N1μ

**Evidence.** The §0 N1μ row says μ* = 0.10: "A ≥ 1678–1706 K", with "range = ω_c ∈ [5000, 50000]".
Table 2 and `einstein_abs.csv` give **1671.5** (ω_c = 5000 K), 1706.0, 1703.3 and 1677.5. The range is
therefore **1672–1706 K**. The μ* = 0.13 and 0.20 rows are correct.

**Required change.** Fix the range.

### F9 [minor] μ* cutoff convention: symmetric block form, which is the permissive choice; state it

**Evidence.** `build_M` subtracts 2μ* on the N_c × N_c block, i.e. θ(ω_c−|ω_n|)θ(ω_c−|ω_m|). The
textbook (Allen–Mitrović) form puts θ only on the summed index m. With the independent solver
(`spot_checks.py` item 6):

- the "asym" form gives Tc lower by 0.13–0.25 % at ω_c = 10 Ω (λ = 0.43–3, μ* = 0.10–0.13);
- it gives Tc lower by **2.6 %** at ω_c = 5000 K, μ* = 0.10 (Table 2 Einstein optimum, N_c = 3:
  300.0 K → 292.3 K).

The symmetric form is therefore the more permissive one, which is the correct direction for
necessary conditions. N1, and the μ* = 0 versions of N2 and N3, do not depend on it.

On the "convention gap" with AD (×1.3–1.9 at λ ≈ 0.4): most of it is the Morel–Anderson
renormalisation between ω_c = 10 ω_max and the phonon scale implicit in AD. With ln 10 ≈ 2.3,
μ* = 0.13 at 10 ω_max acts like about 0.10 at ω_max. The gap is not a defect of the solver, but §4.3
should say this explicitly. Using `mu_rescale` to map literature μ* values (usually quoted for AD or
for ω_c ≈ 5–10 ω_max) onto the doc's cutoff would make the aim box less ambiguous.

**Required change.** In §1, state that the block-symmetric μ* form is a choice that is permissive
relative to the textbook form, and give the 0.2–2.6 % magnitude. In §4.3, add the Morel–Anderson
explanation of the AD gap.

### F10 [info — passes] Adversarial optimisation found no counterexample to N1μ, N2 or N3

`adversarial_opt.py` runs scipy differential evolution with polishing over K = 2, 3, 4 atoms. Shapes
and frequency ratios are free (up to e⁸ ≈ 3000 for N3), and the amplitude is solved for the first
instability at 300 K. It minimises:

- ω₂ at fixed λ (N3);
- λ⟨ω²⟩ (N2);
- A (N1μ).

Each best point was rechecked with truncation multiple 40. `sanity_n3.py` confirms that the
objective separates shapes, e.g. ω₂ = 1169 K for r = 3, f = 0.5 vs 1075 K for Einstein, at λ = 3.

Results (ratio = best found / claimed bound; ≥ 1 means the bound holds):

| problem | K = 2 | K = 3 | K = 4 | minimiser found |
|---|---|---|---|---|
| N3 ω₂, λ = 3, μ* = 0 (bound 1075.33 K) | 1.00000 | 1.00000 | 1.00000 | atoms collapse to one frequency (Einstein) |
| N3 ω₂, λ = 2, μ* = 0 (1416.82 K) | 1.00000 | 1.00000 | 1.00000 | Einstein |
| N3 ω₂, λ = 3, μ* = 0.10, ω_c = 20000 K (1273.31 K) | 1.00000 | 1.00000 | 1.00000 | Einstein |
| N2 λ⟨ω²⟩, μ* = 0 (2.6955×10⁶) | 1.00010 | 1.00010 | — | pushed to the Ω lower box edge (20 K): the soft-mode limit, approached from above |
| N2 λ⟨ω²⟩, μ* = 0.13, ω_c = 20000 K (3.634×10⁶) | 1.00015 | 1.00015 | — | same, soft-mode limit |
| N1μ A, μ* = 0.13, ω_c = 20000 K (1811.33 K) | 1.00000 | 1.00000 | — | Einstein at Ω ≈ 2134 K |
| N1μ A, μ* = 0.13, ω_c = 5000 K (1781.61 K) | 1.00000 | 1.00000 | — | Einstein at Ω ≈ 2264 K |

Raw output: `~/rtsc_targeted_work/R1-review-eliashberg/adv2.out`. For N3 the optimizer had r up to
e⁸ ≈ 3000 and could place a soft atom carrying up to 99 % of λ. It still returned Einstein. This is
consistent with the argument in F3: a soft atom adds λ but no λ⟨ω²⟩, so at fixed total λ it costs
ω₂.

These results support the "minimiser is Einstein" claims beyond the r ≤ 20 grid. They are float
evidence, not a certificate (the doc says this in §6).

### F11 [info — passes] Solver correctness

1. **Gap equation.** `build_M` implements the even sector
   M_nm = λ(|n−m|)[n≠m] + λ(n+m+1) − 2μ*θθ − δ_nm[(2n+1) + 2Σ_{k≤n}λ(k)], with M_nn = λ(2n+1) −
   (2n+1) − 2S_n. This is exactly what follows from Z_n with
   |ω̃_n|/πT = (2n+1) + λ(0) + 2Σ_{k=1}^{n} λ(k): the λ(0) self-term cancels on the diagonal. The code
   achieves this with `fill_diagonal(M, L[2n+1])`, which overwrites the L[0] + L[2n+1] diagonal.
   This is correct. My independent solver keeps λ(0) on both sides and works on the full ± axis. It
   reproduces the reduced form (next item), which also confirms that the even sector is the relevant
   one for these spectra.
2. **Reproductions with the independent solver:**
   - a* = 4.3314817 (doc 4.3314821; certified [4.3314814, 4.3314815]);
   - Table 2 μ* = 0.13, ω_c = 20000 K Einstein optimum → Tc = 300.000 K;
   - Table 4 (λ = 3, μ* = 0) → 300.000 K;
   - Table 4 (λ = 2, μ* = 0.10) → 300.000 K;
   - Table 6 r = 128: Tc/ω_log = 1.2970, Tc/A = 0.2275, Tc/sqrt(λ⟨ω²⟩) = 0.1146 (all match);
   - Einstein Tc/(Ω√λ) = 0.176385 / 0.182098 / 0.182663 at λ = 10 / 100 / 1000 (match), approaching
     C∞ = 0.1827262.
3. **Truncation.** The doc's 30→60 test (≤ 1.4×10⁻⁷) is confirmed in the re-run. The soft-kernel
   B_c is identical at N = 200 and N = 400.
4. **Root finding.** `tc()` scans down from the certified upper bound A/4.33148 and takes the first
   sign change. It would miss only two roots inside one 0.8× bracket, which is harmless for the
   upper-bound use.
5. **Possible weakness, not an error.** `scale_for_300` and `lam_c` use brentq on a bracket that
   assumes monotonicity in the amplitude. If ρ were non-monotone in the amplitude, they could return
   a non-minimal root. My optimiser scans the amplitude upward from below, and the two agree.

### F12 [info — passes, but the proof should state the step] N1μ: Tc as a supremum, with an absolute cutoff

**Evidence.**

- The extreme-point proof in §3.1 is correct: M is affine in α²F, ρ is convex, and ρ < 0 at both ends
  in ω.
- It proves "unstable at T ⇒ A ≥ T·a_E(T)". The table is computed only at T = 300 K. Tc is defined as
  sup{T unstable}, so Tc ≥ 300 K needs T·a_E(T) ≥ 300·a_E(300) for every T ≥ 300 K.
- With an absolute ω_c, the block size N_c(T) drops discontinuously, at T = 318.3 K for ω_c = 5000 K,
  and at 303.2, 335.1 and 374.5 K for ω_c = 20000 K. a_E(T) jumps down at each of those points.

`reentrance_check.py` scans T ∈ [300, 420] K, including both sides of every step, for μ* = 0.10 and
0.13 and ω_c = 5000, 10000 and 20000 K. In every case the minimum of T·a_E(T) is at T = 300 K. The
jumps are much smaller than the linear factor T/300.

**Required change.** Add this step to the proof, with the scan as [NUM] evidence. The same remark
applies to the μ* > 0 rows of N2, whose c₂ depends on ω_c/T.

### F13 [info — passes] The ω_log non-bound example is right

Reproduced: Tc/ω_log = 1.2970 at λ = 2, r = 128, μ* = 0. It stays large at μ* > 0: 0.905 at
μ* = 0.10 and 0.819 at μ* = 0.13, with ω_c = 10 ω_max. The mechanism is simple. The soft atom adds λ
and pulls ω_log down by √r, but contributes λ(k≥1) ∝ (Ω_L/2πkT)² ≈ 0. So
Tc/ω_log ≈ 0.1146·√r → ∞. The claim "no shape-independent ω_log bound at fixed λ" is correct.

### F14 [minor] Smaller points

1. **Rounding at the A floor.** "A_hi < 1300 K ⇒ reject" rounds the certified 1299.44 K *up*, so it
   rejects the interval [1299.44, 1300). The effect is negligible, but a "theorem-backed" rule should
   use 1299.4 K, or round down.
2. **Aim box A.** "A ≥ 2050 K" is described as 1.1 × (1860–1970). From Table 1, μ* = 0.13, the
   hydride rows are 1875 (H3S-like) and 1969 (LaH10-like), so 1.1× gives 2063–2166 K. Either the
   source range or the result is slightly off.
3. **Mixed cutoff conventions in §3.1.** "Every other shape tested needs more A (Table 1)" compares
   Table 1 (ω_c = 10 ω_max) with Table 2 (absolute ω_c). The statement is true by the proof itself,
   but the cross-table comparison is not like-for-like.
4. **Citations checked:**
   - Bergmann & Rainer, Z. Phys. 263, 59 (1973). The abstract says "positive contribution [from] all
     frequency regions", most effective "slightly above 2πTc". The doc is correct.
   - Carbotte, RMP 62, 1027 (1990). The abstract calls the asymptotic limit "of theoretical interest
     only". The doc is correct.
   - Allen & Dynes, PRB 12, 905 (1975). The abstract gives Tc ≈ 0.15 (λ⟨ω²⟩)^{1/2} *at μ* = 0.1*,
     says Tc is not limited by phonon frequencies, and calls the "λ = 2 limit" spurious. The doc
     quotes the scaling correctly. It could add that AD's 0.15 at μ* = 0.1 is 7 % below the doc's
     c₂(0.10) = 0.1617, which is the convention gap again (F9).
   - Yuzbashyan & Altshuler, PRB 106, 054518 (2022), λ_c ≈ 3.0–3.7, negative specific heat near
     3.7. The doc is correct.
   - The doc does not cite a Semenok λ_c paper, so there was nothing to check.

---

## What I could not verify

- Leavens' 1975 constant. The doc already marks it UNVERIFIED.
- The full text of Carbotte, and of Esterlis–Kivelson–Scalapino npj 2018.
- The exact J. Stat. Phys. volume numbers for KAY II and III.
- A certified (interval) version of any N1μ, N2 or N3 number. All of my checks are float checks.
