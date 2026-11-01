# R03: Adversarial review of `docs/theory/03_certificates.md` and `src/certify/`

Reviewer: agent R3-review-certificates. Date: 2026-10-28.
Scripts and logs: `docs/reviews/r03_scripts/`. Independent solver: `docs/reviews/r01_scripts/indep_solver.py`.
It uses the full Matsubara axis, keeps λ(0) in both Z and the kernel, and solves a non-symmetric eigenproblem.

## Verdict

**I found no error that would produce a false "proof".**

I checked the following line by line, and I re-derived each analytic bound myself:

- the Lagrangian reduction (§2.2);
- the end lemmas;
- the shift and Coulomb additions to the cell lemma;
- the Schur and tail bookkeeping in `cellcert.py`;
- the N1μ segment logic;
- the lookup and T-supremum arguments.

All of them are correct.

`verify.py --quick` and the full `verify.py` both pass on my run (results below). The adversarial float search finds no counterexample. The worst case is the Einstein point itself, as the theory predicts.

The findings below are bookkeeping errors, statements of scope that need tightening, and robustness gaps in the checker. None is critical. Two are worth fixing before the rule is used to issue rejections: findings 2 and 3.

## Re-checks run

| Check | Result | Log |
|---|---|---|
| `verify.py --quick` (every 10th cell, all ends, tiling, Rayleigh) | `ALL CERTIFICATES VERIFIED`, exit 0 | `r03_scripts/verify_quick.log` |
| `verify.py` (full, every stored cell) | 40/40 files PASS, `ALL CERTIFICATES VERIFIED`, exit 0 | `r03_scripts/verify_full.log` |
| Adversarial NM-G spectra, m = 0.9999·g_lo(λ_k), λ ≤ λ_k, T ∈ {300, 330, 450} K, λ_k ∈ {0.8, 1, 1.5, 2, 3, 5, 10, 100} | No violation. Worst top eigenvalue 0.999999 (λ_k = 100). Every maximiser collapses onto the Einstein atom. | `r03_scripts/adversarial_nm.log` |
| Sharpness: Einstein(λ_k) at m = 1.001·g_lo | Unstable for every λ_k tested (top eigenvalue 1.00001–1.00025) | same |
| Large-λ truncation convergence (N = 80/300/700) | Converged to < 1e-10. Einstein λ = 1e4 is stable at 0.99999·m₀(∞) and unstable at 1.0001·m₀(∞). | `r03_scripts/large_lambda_conv.log` |
| N1μ: A = 0.9999·A_floor, Einstein and 2-atom spectra, T grid ≥ 300 K including both sides of every N_c jump | No violation. Max top eigenvalue 0.99985. At A = 1.002·A_float the solver finds instability (1.0009–1.0011), so the floors are tight. | `r03_scripts/adversarial_nm.log` |
| Corollary Tc ≤ 0.18272634·√m | Holds. Einstein Tc/(C√m) = 0.9653 (λ = 10), 0.99656 (λ = 100), 0.99966 (λ = 1000): it approaches the constant from below. | same |
| μ* in the *asymmetric* cutoff convention, 150 random cases | Never raised the top eigenvalue (float only) | `r03_scripts/asym_mu_check.log` |

The tested shapes were:

- Einstein atoms at λ_k and at 0.95, 0.8, 0.5 and 0.2·λ_k;
- two-atom spectra with frequency ratios 1.5–1000 and weight splits 0.02–0.98;
- soft + hard pairs with λ < λ_k;
- box spectra;
- Nelder–Mead maximisation over two- and three-atom spectra on the constraint surface.

## Proof check (item 1)

1. **Lagrangian reduction (§2.2). Correct.**
   - The even-sector form is
     Q⁺(x) = Σ_{n≠m}[l(|n−m|)+l(n+m+1)]x_nx_m + Σ_n[l(2n+1) − 2S_n − (2n+1)]x_n².
     I re-derived it from Z_n and Δ_n with y_n = Δ_n/(2n+1). The λ(0) terms cancel, and the matrix agrees with `arbcore.even_block` and with the independent solver.
   - Q⁺ is affine in α²F. The constant part −D is integrated against a probability measure ν, so it is unchanged.
   - The atom E_ω = (c/ω)δ_ω gives l[E_ω](k) = m/(ω²+(2πk)²), λ_ω = m/ω² and ∫λ_ω dν = λ. I checked all three.
   - Affinity in m′ at fixed w turns (C) at m₀ plus (C₀) at m′ = 0 into a bound for all m′ ∈ [0, m₀].
   - The sign is right. The step needs β ≥ 0, and then −β(λ_k − λ) ≤ 0 when λ ≤ λ_k. All stored β are ≥ 0 (see finding 4).
   - Integrability is fine: each matrix entry is bounded by m/(4π²k²), and x has finite support.
   - "Every μ* ≥ 0, any cutoff" follows from (F1) in the symmetric-block convention (finding 3).
2. **Scaling (§2.1). Correct.** α̂ is the push-forward divided by T. Then λ is invariant and m scales as m/T².
3. **Ends (§2.3). Correct.**
   - The unit-weight Schur row sum of the full form is S_n + L + (L − S_n) − 2S_n − (2n+1) = 2L − 2S_n − (2n+1), which is ≤ 2L − 1.
   - L ≤ m₀/24 for all w.
   - L = (m₀/4w)(coth(w/2) − 2/w) ≤ m₀/(4w), because coth u ≤ 1 + 1/u.
   - The small-w bound is increasing in w, so its value at w_s bounds all smaller w.
   - The large-w bound drops −βm₀/w² ≤ 0, which is legitimate.
4. **Cell additions. Correct.**
   - s″ = −6βm₀/w⁴ ≤ 0, so the tangent line at the cell centre is an upper bound. The Taylor model stays linear in t, so checking t = ±h suffices.
   - The tail uses sup_cell s = s(w₂), and s is increasing. For s < 0 this only enlarges γ_m, which is correct.
   - The Coulomb block is constant, enters A0 only, and is supported on n < N_c ≤ N = 120.
   - I re-derived the tail lemma: row ≤ T(m+1) + T(m+N+1) − (S_m − S_{m−N}) − (2m+1) ≤ 2T(m+1) − (2m+1).
   - I re-derived the far-coupling term: the row sums are T(N+K−n) + T(N+K+n+1). The column sums are ≤ T(K+1) ≤ cmax, so the code's choice is conservative.
   - The derivative formulas for both families are correct: l″ = m₀(6w²−2b²)/(w²+b²)³ and l″ = −4aw(3b²−w²)/(w²+b²)³.
   - The closed forms for L and the integral tail bounds are correct.
5. **T-supremum. Correct.** λ does not depend on T, m/T² is decreasing, and Theorem NM holds separately at each T.
6. **Lookup (§2.5). Sound, apart from a wording gap (finding 5).**
7. **N1μ (§4). Correct.**
   - N_c(T) = j exactly on (ω_c/((2j+1)π), ω_c/((2j−1)π)]. `Nc_of` = ⌊(ω_c/(πT)+1)/2⌋ is the right count.
   - I recomputed the segment lists by hand for all 8 (μ*, ω_c) pairs. They match the JSONs, and every segment with T_lo(j) < T_u = A_floor/a0 carries a certificate.
   - T ≥ T_u gives A/T ≤ a0, so (F2) applies as invoked.
   - The monotonicity-in-a argument and the reduction of a general α²F to Einstein kernels at fixed A are correct.
   - I re-derived the large-w weighted Schur bounds:
     - rows n < n₀: ε − 1 + 2(n₀−1)ε + a/H. The doc's (a + (n₀−1)ε)/H is larger, so it is conservative.
     - rows n ≥ n₀: 2L − (2n₀+1) + 2(H−1)n₀ε, using 2L ≤ a.
8. **Code against proof. Matches.**
   - `MFamily.shift` = β(λ_k − m₀/w²) and `dshift` = 2βm₀/w³.
   - The `ends_check` formulas match §2.3.
   - `rigorous_negdef` is a valid PD certificate: diagonal scaling by exact powers of two, then ‖E‖_∞ < s, where E = S − sI − LLᵀ is an Arb ball.
   - `tiled` demands exact float equality of adjacent cell endpoints, and the w_range endpoints are the same values the end lemmas use.
   - Rounding goes the right way:
     - m₀ and A_floor are decimal floors;
     - a_j is a decimal ceiling, with a_j·T_lo ≥ A_floor re-checked in Arb;
     - C_cert = 0.18272634 is ≥ 1/√29.95005976 = 0.1827263391.
   - The downstream rule (`src/rtt/targets/gtable.py`) reads the exact decimal `m0_over_T2` as a Fraction, not the float `g_lo_K2`.

## Findings

### 1. Cell and segment counts in the document are wrong. Severity: low (bookkeeping).

- §0 and §3 say "1809 certified cells". The JSONs sum to **1793**, which is also the sum of the doc's own cells column.
- §4.3 says "All **38** segment certificates closed". There are **32** segment certificates (2+2+4+7+2+2+4+9).

**Required change:** correct both numbers. Have `summarize.py` print the totals so the doc can quote them.

### 2. The checker is the generator: there is no independent implementation of the Arb inequalities. Severity: medium (acknowledged in §5.3, but it is the main residual risk).

`verify.py` re-runs `cellcert.certify_cell`, `ends_check` and `ends_A`, which are the same functions that produced the certificates. A shared bug would pass silently. Examples:

- a sign error in a tail bound;
- a wrong index in `d`;
- a wrong `shift_sup` direction.

I found no such bug by reading, and the float adversarial tests are consistent with the claims. Still, a "computer-assisted theorem" used to issue proved rejections should have at least one of the following:

- (a) a second, independently written cell checker (for example mpmath intervals, or a different blocking N/K), run on all cells;
- (b) a Lean statement of the Lagrangian reduction, Lemma E and the shift/Coulomb additions. These are short; §5.3 lists them as paper-only.

**Required change:** add (a) before relying on proved rejections, or label the rejections "THM-C (single implementation)".

### 3. "Every μ* ≥ 0 and every cutoff" is proved only in the symmetric-block convention. Severity: medium (scope statement).

(F1) is a quadratic-form statement: μ* adds −2μ*(Σ_{n<N_c}x_n)². That form exists only for the symmetric Coulomb block.

In the textbook asymmetric convention (cutoff on the summed index only), the eigenproblem is not symmetric, and "μ* only lowers Q⁺" is not what (F1) proves. My float check of 150 random cases never saw μ* raise the top eigenvalue, but that is not a proof.

**Required change:** in §0, §2.2 and the Corollary, write "every μ* ≥ 0 and every cutoff ω_c, symmetric-block convention (01 §1)". Alternatively, add a proof for the asymmetric convention. At μ* = 0 the issue disappears, and NM-G as a rule only uses μ* = 0 as a lower bound.

### 4. `verify.py` does not check β ≥ 0. Severity: low (checker robustness).

The proof needs β ≥ 0 (§2.2, step "using (2) and β ≥ 0"). `verify.py` re-checks βλ_k < 1 − η through the m = 0 endpoint, but never checks the sign.

A JSON with β < 0 would verify, but it would certify the statement "λ ≥ λ_k" instead of "λ ≤ λ_k". All current β are ≥ 0, because `propose_beta` searches [0, 0.95/λ].

**Required change:** in `verify_nm`, assert `arb(beta) >= 0`. Also assert that `upper["lambda"] == lambda` for finite λ_k, and that `g_lo_K2 <= 90000*m0`.

### 5. Lookup proof (§2.5) equates "Tc ≥ Tc_target" with "some T ≥ Tc_target is unstable". Severity: low (wording).

Tc = sup{T unstable} need not be attained. The rule is still sound, because it uses a strict inequality. If m_hi < m₀T*², then by continuity m_hi < m₀T² for all T in some interval [T* − ε, ∞). So every unstable T must lie below T* − ε, which gives Tc ≤ T* − ε < T*.

**Required change:** replace the sentence with this two-line argument. §5.1 already notes the reading, but the proof of §2.5 should not depend on it.

### 6. The N1μ checker does not consider the segment j = 0 (T > ω_c/π, N_c = 0). Severity: low (latent; not triggered).

`verify_n1mu` loops j = j300 … 1. For all 8 stored pairs this is harmless:

- the j = 1 segment is already covered by the universal bound;
- ω_c/π ≥ 1591 K > T_u ≈ 386–412 K.

For a small ω_c (ω_c < π·T_u ≈ 1300 K) the gap would be real.

**Required change:** add an explicit check `wc/pi >= T_u`, or a j = 0 certificate.

### 7. KAY III corollary: consistent with the literature, but over-labelled. Severity: low (citation and statement).

What I verified from the full text of arXiv:2409.02121v2 (Kiessling, Altshuler, Yuzbashyan, J. Stat. Phys. 192, 93 (2025)):

- C∞ = (1/2π)√𝔤(2) = 0.1827262477.
- Their kernel is λΩ²/(Ω² + (ω_n − ω_m)²), which is the same as l(k) here, at μ* = 0.
- Their **Conjecture 1** is the *strict* Einstein bound Tc(λ, Ω) < C∞Ω√λ for all λ > 0.
- Their rigorous Einstein bound has constant ≈ 2.034·C∞ asymptotically.
- They prove no general-spectrum bound of this form.
- They define Tc as the crossing λ = L_E(ϖ), unique for λ ≥ λ* < 0.767. The doc's sup-definition is at least as large, so the doc's bound implies one in their definition.

So the doc's corollary has three properties:

- (i) It is consistent: C_cert/C∞ − 1 = 5.0×10⁻⁷, and float Einstein ratios approach 1 from below.
- (ii) It proves a version of Conjecture 1 relaxed by the factor (1 + 5.0×10⁻⁷), not the strict conjecture.
- (iii) It extends to all spectra through the elementary affinity at fixed m. That is new relative to KAY III, but it should be stated as an easy consequence of the reduction, not as a deep result.

The large-λ asymptote 0.183√(λ⟨ω²⟩) for arbitrary shapes is classical and should be cited: Allen & Dynes, Phys. Rev. B 12, 905 (1975).

**Required change:**

- Replace "[UNVERIFIED]: only the abstract was checked" by the facts above: same kernel, μ* = 0, and the Tc-definition difference.
- Say "proves Conjecture 1 up to a factor 1 + 5.1×10⁻⁷ (non-strict)".
- Cite Allen–Dynes 1975 for the asymptotic constant.

### 8. "Why a β exists" (§2.2) is heuristic. Severity: informational.

The equivalence "a supporting line exists iff (C)" is correct, because the shift is a scalar multiple of I. But the existence of the supporting line is only observed, not proved. This does not affect soundness: each certificate exhibits its own β.

Note also the small margins for some λ_k. The minimum scaled cell margin is about 6×10⁻⁸ for λ = 1.1 and about 1×10⁻⁷ for λ = 0.8 and 0.9. This is fine for 128-bit Arb, but it explains why δ = 10⁻⁶ failed at λ = 0.8. A δ of 10⁻⁵ everywhere would cost nothing in practice and add robustness.

## What I could not verify

- The cited core results: (F1) odd-sector and μ*-monotonicity, and (F2) the a0 = 4.33148145 universal certificate. I took them as given. My float checks with A ≤ A_floor at T up to 520 K are consistent with (F2) but do not test it.
- KAY II (the general dispersive paper). I could not confirm whether it already contains the reduction to Einstein atoms at fixed m.
- Floating-point Cholesky and Arb behaviour in python-flint itself. These are trusted, as in the certified core.
