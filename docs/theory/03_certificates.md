# 03 — Interval-arithmetic certificates for NM-G and N1μ

Status: theory deliverable C1, revised after review R03 (§7).

This file does two things:

- It replaces the numerical steps (a) "Einstein minimiser" and (b) "g decreasing" of
  `01_eliashberg_targets.md` §0.3/§3.3 with one certified theorem per grid value of λ.
- It certifies the N1μ floors (01 §3.1) at μ* ∈ {0.10, 0.13} and ω_c ∈ {5000, 10000, 20000, 50000} K,
  including the T-supremum step.

**Scope of every statement:** the isotropic linearized ME model of 01 §1 with the **symmetric
Coulomb-block convention**. For μ* = 0 the convention does not matter.

Tags as in 01:

- **[THM-C]**: theorem, computer-assisted. Arb ball arithmetic, checked by two independent checkers,
  `src/certify/verify.py` and `src/certify/verify_indep.py`.
- **[THM]**: analytic proof given here.
- **[CIT-CORE]**: proved earlier for the certified core `rtt.core` (sources given).

**Code** is in `src/certify/`:

- Generator: `arbcore.py`, `families.py`, `cellcert.py`, `cover.py`, `rayleigh.py`, `run_nm.py`,
  `run_n1mu.py`.
- Checkers: `verify.py` (re-runs the generator's routines) and `verify_indep.py` (independent, §6).
- Tools: `summarize.py`, `controls.py`.

`arbcore.py`, `cellcert.py`, `cover.py` and `rayleigh.py` generalise the routines of the wE1
certificate of the core (`core`, `cell`, `run_cells`, `sharp`).

**Certificates** are in `docs/theory/certificates/`:

- `nm_<λ>.json` (32 files) and `n1mu_mu<μ*>_wc<ω_c>.json` (8 files);
- `SUMMARY.json`;
- the checker logs `verify_log.txt`, `verify_indep_log.txt` and `controls_log.txt`;
- per-cell generation logs in `logs/`.

**Re-check in one command each** (WSL, repo root; python-flint 0.9.0, numpy, scipy; one thread):

```
OMP_NUM_THREADS=1 nice -n 19 python src/certify/verify.py          # generator-code re-check (~12 min)
OMP_NUM_THREADS=1 nice -n 19 python src/certify/verify_indep.py    # independent re-check (~15 min)
```

**`verify.py`** does the following:

- re-reads the exact parameters;
- checks the preconditions:
  - β ≥ 0 (β = 0 at λ = ∞);
  - the Rayleigh certificate is for the same λ;
  - the reported g_lo ≤ 90000·m₀ (exact rational comparison);
  - the N1μ segment j = 0 is covered;
- re-evaluates the end lemmas;
- re-certifies every stored cell, without splitting;
- checks that the cells tile the stated range exactly;
- re-evaluates the Rayleigh certificates.

**`verify_indep.py`** shares no code with the generator (§6).

Both exit with code 0 iff everything passes.

---

## 0. Results

**NM-G is a computer-assisted theorem on the whole grid** [THM-C]:

- Every tabulated λ_k (0.8 … 100) and λ = ∞ is certified.
- 32 certificates, **1661 cells**, **0 open boxes**.
- Each certificate was re-checked by both checkers.
- Step (a) (Einstein minimiser) and step (b) (g decreasing) are no longer needed.
- Every certified g_lo is 1.0×10⁻⁵ relative below the float g of 01 §0.3 (slack δ = 10⁻⁵).

For each row the certified statement is Theorem NM(λ_k, m₀) of §2.2:

> Every α²F with λ ≤ λ_k and λ⟨ω²⟩ ≤ g_lo(λ_k)·(T/300 K)² is stable at T, for every μ* ≥ 0 and
> every cutoff ω_c (symmetric-block convention).

The Rayleigh certificate of §2.4 gives the upper side: Einstein(λ_k, Ω_hi) is unstable at 300 K.
Hence Ω_E(λ_k) ∈ [Ω_lo, Ω_hi], and the sharp g(λ_k) ≤ g_hi.

| λ_k | g_lo (K², certified; use in NM-G) | g_hi (K²) | Ω_lo (K) | Ω_hi (K) | float g (01 §0.3) | β (Lagrange) | cells |
|---|---|---|---|---|---|---|---|
| 0.8 | 9.452681e6 | 9.454666e6 | 3437.419 | 3437.780 | 9.452776e6 | 0.97561 | 68 |
| 0.9 | 7.871489e6 | 7.873143e6 | 2957.381 | 2957.691 | 7.871568e6 | 0.79547 | 64 |
| 1.0 | 6.850483e6 | 6.851923e6 | 2617.343 | 2617.618 | 6.850552e6 | 0.66044 | 59 |
| 1.1 | 6.147092e6 | 6.148384e6 | 2363.951 | 2364.200 | 6.147154e6 | 0.55661 | 55 |
| 1.2 | 5.637762e6 | 5.638946e6 | 2167.518 | 2167.746 | 5.637819e6 | 0.47506 | 53 |
| 1.3 | 5.254201e6 | 5.255305e6 | 2010.396 | 2010.607 | 5.254254e6 | 0.40994 | 51 |
| 1.4 | 4.956119e6 | 4.957161e6 | 1881.512 | 1881.709 | 4.956169e6 | 0.35717 | 50 |
| 1.5 | 4.718445e6 | 4.719436e6 | 1773.592 | 1773.779 | 4.718493e6 | 0.31385 | 49 |
| 1.6 | 4.524867e6 | 4.525818e6 | 1681.678 | 1681.855 | 4.524913e6 | 0.27785 | 48 |
| 1.7 | 4.364371e6 | 4.365288e6 | 1602.273 | 1602.441 | 4.364415e6 | 0.24767 | 49 |
| 1.8 | 4.229276e6 | 4.230164e6 | 1532.840 | 1533.001 | 4.229318e6 | 0.22212 | 45 |
| 1.9 | 4.114078e6 | 4.114942e6 | 1471.497 | 1471.652 | 4.114120e6 | 0.20028 | 49 |
| 2.0 | 4.014738e6 | 4.015582e6 | 1416.817 | 1416.965 | 4.014779e6 | 0.18151 | 47 |
| 2.25 | 3.817620e6 | 3.818422e6 | 1302.582 | 1302.719 | 3.817658e6 | 0.14462 | 48 |
| 2.5 | 3.671310e6 | 3.672082e6 | 1211.827 | 1211.954 | 3.671347e6 | 0.11789 | 46 |
| 2.75 | 3.558519e6 | 3.559267e6 | 1137.544 | 1137.664 | 3.558555e6 | 0.097924 | 44 |
| 3.0 | 3.468968e6 | 3.469697e6 | 1075.324 | 1075.437 | 3.469003e6 | 0.082609 | 46 |
| 3.25 | 3.396173e6 | 3.396887e6 | 1022.241 | 1022.348 | 3.396208e6 | 0.070624 | 44 |
| 3.5 | 3.335852e6 | 3.336553e6 | 976.269 | 976.371 | 3.335886e6 | 0.061071 | 46 |
| 3.75 | 3.285060e6 | 3.285750e6 | 935.957 | 936.056 | 3.285093e6 | 0.053326 | 46 |
| 4.0 | 3.241711e6 | 3.242392e6 | 900.238 | 900.332 | 3.241744e6 | 0.046964 | 48 |
| 4.5 | 3.171646e6 | 3.172313e6 | 839.530 | 839.618 | 3.171678e6 | 0.037231 | 45 |
| 5.0 | 3.117483e6 | 3.118138e6 | 789.618 | 789.701 | 3.117515e6 | 0.030232 | 49 |
| 6.0 | 3.039235e6 | 3.039874e6 | 711.716 | 711.790 | 3.039266e6 | 0.021076 | 49 |
| 7.0 | 2.985444e6 | 2.986072e6 | 653.064 | 653.132 | 2.985474e6 | 0.015525 | 50 |
| 8.0 | 2.946200e6 | 2.946819e6 | 606.857 | 606.920 | 2.946230e6 | 0.011909 | 51 |
| 10.0 | 2.892781e6 | 2.893389e6 | 537.846 | 537.902 | 2.892810e6 | 0.0076422 | 52 |
| 15.0 | 2.824200e6 | 2.824793e6 | 433.912 | 433.958 | 2.824229e6 | 0.0034082 | 57 |
| 20.0 | 2.790995e6 | 2.791581e6 | 373.564 | 373.603 | 2.791023e6 | 0.0019203 | 59 |
| 50.0 | 2.732967e6 | 2.733542e6 | 233.793 | 233.818 | 2.732995e6 | 0.00030815 | 64 |
| 100.0 | 2.714107e6 | 2.714677e6 | 164.745 | 164.763 | 2.714135e6 | 7.7111e-05 | 66 |
| ∞ | **2.695481e6** | 2.696047e6 | – | – | 2.695508e6 | 0 | 64 |

Notes on the table:

- **Exact value:** the certified value is g_lo = 90000·m₀ exactly, with the decimal m₀ stored in
  each JSON as `m0_over_T2`. Downstream code should read m₀, not the float.
- **Rounding:** g_lo is rounded *down* here (floor, 7 significant digits). Ω_lo is rounded down and
  Ω_hi up. The JSON floats use directed rounding.
- **g_hi:** this is the Einstein value at Ω_hi, which is Rayleigh-certified unstable. So g_hi is an
  upper bound on the sharp g. The ratio g_hi/g_lo − 1 ≈ 2×10⁻⁴ only reflects the chosen
  10⁻⁴ step in Ω for the upper certificate.
- **Corollary (λ = ∞) [THM-C]:**
  - Within the model, every α²F satisfies Tc ≤ C_cert·sqrt(λ⟨ω²⟩), for every μ* ≥ 0 and every
    cutoff (symmetric block).
  - C_cert = 1/sqrt(29.94979021) = 0.18272717 = C∞·(1 + 5.0×10⁻⁶).
  - See §2.5 for the relation to KAY III and Allen–Dynes.

**Rule changes, applied in 01 §0.3–0.4:**

- Use the g_lo column with the same lookup: the smallest tabulated λ_k ≥ λ_hi, or ∞.
- Mark NM-G as [THM-C], with recorded assumptions M-ISO and CERT-NM.
- Drop NUM-EINSTEIN-MIN and NUM-G-MONOTONE.
- Any sub-grid of the certified grid is also sound.

**N1μ is certified** at μ* ∈ {0.10, 0.13} × ω_c ∈ {5000, 10000, 20000, 50000} K, including the
T-supremum: 32 segment certificates, 2864 cells, 0 open boxes. The table is in §4.3.

---

## 1. Model, notation, and the facts used from earlier work

The model is that of 01 §1:

- isotropic linearized ME;
- the even sector Q⁺_T;
- the symmetric Coulomb block −2μ*(Σ_{n<N_c}x_n)²;
- ρ(T) = sup_x Q⁺_T(x)/|x|²;
- T is unstable iff ρ(T) ≥ 0;
- Tc := sup{T unstable}.

Units are k_B = ħ = 1. For a spectrum α²F we write λ = 2∫α²F/ω, m := λ⟨ω²⟩ = 2∫ω α²F and A = ∫α²F.

The even-sector form is

  Q⁺(x) = Σ_{n≠m}[l(|n−m|) + l(n+m+1)] x_n x_m + Σ_n [l(2n+1) − 2S_n − (2n+1)] x_n²,

with S_n = Σ_{k=1}^n l(k). R03 re-derived it from Z_n and Δ_n.

Facts used without re-proof:

- **(F1)** [CIT-CORE] The odd sector is never unstable, and μ* ≥ 0 only lowers Q⁺. Sources: the
  multiband proof wMB Thm 3, and `src/rtt/core/certify/PROOFS.md` R-MU/R-ODD. This is a
  statement about quadratic forms, so it holds in the **symmetric-block convention**: there "stable
  for μ* = 0 in the even sector" implies stable for every μ* ≥ 0 and every cutoff. In the
  asymmetric textbook convention the eigenproblem is not symmetric and (F1) does not apply as
  stated. R03's float check of 150 random cases saw no counterexample, but this is **not proved**.
- **(F2)** [CIT-CORE] The universal certificate. For a ≤ a0 = 4.33148145, every Einstein kernel
  l(k) = 2aw/(w² + (2πk)²) satisfies Q⁺ ≤ −η|x|² for all w > 0 (the wE1 certificate).
  Hence every α²F with A/T ≤ a0 is stable at T, for every μ* ≥ 0 and every N_c.
- **(F3)** [THM; wE1 certificate Lemmas 1–6, reproduced in `cellcert.py`] The cell method. Take a
  one-parameter kernel family l_w(k) that is positive and decreasing in k. A cell [w₁, w₂] is
  certified by four ingredients:
  1. a first-order Taylor expansion of the N×N leading block, with a rigorous |l''| remainder;
  2. the tail-block bound z^T C z ≤ −Σ(2m+1−2T(m+1)) z_m² (Lemma 4);
  3. a Schur-complement coupling with K explicit columns and an analytic remainder (Lemma 5);
  4. verified positive definiteness of −P(±h) (Lemma 6).

**Additions to the cell method made here** (each has a one-line proof; R03 re-derived all of them).

1. **Diagonal shift.** For the NM family a term s(w)·|x|² is added, with s(w) = β(λ − m₀/w²).
   - It enters the N×N block on the diagonal.
   - It enters the tail diagonal through γ_m ≤ 2m+1 − 2T̄(m+1) − sup_cell s − η, which is Lemma 4
     applied to C + sI.
   - Since s'' = −6βm₀/w⁴ ≤ 0, we have s(w_c + t) ≤ s(w_c) + t s'(w_c). So the Taylor model of the
     shift needs no remainder.
2. **Coulomb block.** For N1μ the constant term −2μ*(Σ_{n<N_c} x_n)² is added to the leading
   block (N_c ≤ N). It does not touch the tail or the coupling blocks.
3. **The kernel l(k) = m₀/(w² + (2πk)²).**
   - Closed form for its sum: L = Σ_{k≥1} l(k) = (m₀/2)(coth(w/2)/(2w) − 1/w²).
   - Tail bound: T(j) ≤ l(j) + min(m₀ arctan(w/2πj)/(2πw), m₀/(4π²j)). The second form is valid on
     balls containing w = 0.
   - Derivatives: l' = −2m₀w/(w²+b²)² and l'' = m₀(6w²−2b²)/(w²+b²)³.

**Controls** (`src/certify/controls.py`; output in `certificates/controls_log.txt`; all OK):

- The Arb block equals the float matrix of `src/theory/eliashberg.build_M` to 1.4×10⁻¹⁴.
- The cell certifier *fails*, as it must, in three cases:
  - m₀ = 1.0001·g at w ≈ w_E;
  - λ = ∞ with m₀ = 1.0001·g(∞) on [0, 2⁻⁶];
  - finite λ without the shift (β = 0) on [0, 2⁻⁶].
- It succeeds for 0.9999·g on the same cells.

R03 also ran adversarial float searches over 2–3-atom, box and soft+hard spectra. It found no
violation, and every maximiser collapsed onto the Einstein atom.

---

## 2. Rule NM-G: the certified theorem

### 2.1 Scaling (μ* = 0) [THM]

Let T > 0, and let α̂ be the image of α²F under ω ↦ ω/T, divided by T. Then

- l_T[α²F](k) = l_1[α̂](k) for all k;
- λ[α̂] = λ[α²F];
- m[α̂] = m[α²F]/T².

To see this, substitute ω = Tu in l_T(k) = ∫2ω dα²F/(ω² + (2πkT)²). So α²F is stable at T iff α̂
is stable at T = 1, and only the pair (λ, m/T²) matters.

### 2.2 Theorem NM(λ_k, m₀) [THM-C]

> For every finite positive measure α²F on (0, ∞) with λ[α²F] ≤ λ_k and m[α²F] ≤ m₀·T², the
> even-sector form at temperature T satisfies Q⁺_T(x) ≤ −η|x|², with η = 10⁻⁷. Hence α²F is stable
> at T for every μ* ≥ 0 and every cutoff ω_c, in the symmetric-block convention (F1).

For λ_k = ∞ the hypothesis is only m ≤ m₀T².

**Proof.**

1. By §2.1 we may take T = 1. If m = 0 then α²F = 0 and Q⁺ = −Σ(2n+1)x_n², so there is nothing to
   prove.
2. Otherwise put c = m/2 = ∫ω dα²F and define the probability measure ν = ω α²F/c. Let
   E_ω := (c/ω)δ_ω. This is the Einstein atom at ω with moment m and coupling λ_ω = m/ω².
3. Then l[α²F](k) = ∫ l[E_ω](k) dν(ω) for every k. So for every finitely supported x

     Q⁺[α²F](x) = ∫ Q⁺[E_ω](x) dν(ω).                                  (1)

   The bare term −Σ(2n+1)x_n² integrates to itself because ν(ℝ₊) = 1. The interchange of sum and
   integral is fine: each entry is bounded by m/(4π²k²), and x has finite support. Also

     ∫ λ_ω dν(ω) = ∫ (m/ω²)(ω/c) dα²F = 2∫ dα²F/ω = λ[α²F] ≤ λ_k.        (2)

4. The matrix of Q⁺[E_ω] is m·Φ(ω) − D, where Φ is built from 1/(ω² + (2πk)²) and D = diag(2n+1).
   The certificate (§2.3) proves, for all w > 0 and a fixed β with 0 ≤ β < 1/λ_k:

     (C)  m₀Φ(w) − D + β(λ_k − m₀/w²) I ≤ −η I,
     (C₀) −D + βλ_k I ≤ −η I.

5. For fixed w and x, the quantity x^T[m'Φ(w) − D + β(λ_k − m'/w²)I]x is affine in m'. So (C) and
   (C₀) give it ≤ −η|x|² for all m' ∈ [0, m₀], in particular for m' = m. Hence

     Q⁺[E_ω](x) ≤ −η|x|² − β(λ_k − λ_ω)|x|²   for every ω.

6. Integrate this with (1), and use (2) and β ≥ 0:

     Q⁺[α²F](x) ≤ −η|x|² − β(λ_k − λ[α²F])|x|² ≤ −η|x|².  ∎

Remarks on the proof:

- β is a Lagrange multiplier for the constraint on λ. For λ_k = ∞ take β = 0; then (2) is not
  needed.
- The proof needs no reduction to two-atom spectra, no statement that Einstein spectra are the
  minimisers, and no monotonicity of g.
- The condition β ≥ 0 is essential. With β < 0 the same computation would certify the statement
  for λ ≥ λ_k instead. Both checkers enforce β ≥ 0.

**Existence of β (heuristic, R03 finding 8; this does not affect soundness).** Let ψ(λ') be ρ of
the Einstein atom with moment m₀ and coupling λ'. Because the shift is a scalar multiple of I,
(C) says exactly that the line λ' ↦ −η + β(λ' − λ_k) lies above ψ. So a suitable β exists iff ψ
has a supporting line at λ_k.

We have not proved that such a line exists. We only observe it: `run_nm.propose_beta` finds one
for every λ_k. Each certificate exhibits its own β, so the theorem does not depend on this
observation.

That the relaxation loses nothing at the Einstein point is why g_lo can be taken within δ of the
Einstein value. Without the shift (β = 0) the certificate must fail for finite λ_k: the soft-mode
limit ω → 0 is unstable at m₀ = g(λ_k) > g(∞). The controls confirm this.

### 2.3 How (C) is certified [THM-C]

The w-axis is split into (0, w_s] ∪ [w_s, W] ∪ [W, ∞).

- **Cells on [w_s, W].**
  - Dyadic cells, adaptively bisected.
  - Each cell is certified with F3 plus the shift (N = 120, K = 600, 128-bit Arb).
  - For λ_k = ∞ the cells start at w_s = 0.
- **Small w** (β > 0, w ≤ w_s).
  - By Lemma 3 and the unit-weight Schur test, x^T(m₀Φ − D)x ≤ (2L − 1)|x|².
  - L ≤ m₀Σ1/(2πk)² = m₀/24.
  - So the left side of (C) is ≤ m₀/12 − 1 + βλ_k − βm₀/w². This is increasing in w, so it is
    bounded by its value at w_s, which is checked < −η in Arb.
- **Large w** (w ≥ W).
  - L = (m₀/4w)(coth(w/2) − 2/w) ≤ m₀/(4w), because coth u ≤ 1 + 1/u.
  - So the left side of (C) is ≤ m₀/(2W) − 1 + βλ_k, which is checked < −η.
- **(C₀).** βλ_k − 1 < −η is checked in Arb.

The parameters for each λ_k are stored in its JSON:

- m₀, as an exact decimal;
- β, as an exact dyadic;
- w_s and W;
- every cell.

Floating point only *proposes* β, w_s, W and the cell splits.

### 2.4 Upper side: Einstein brackets [THM-C]

For each finite λ_k we have a dyadic w_hi ≈ w_E(1 + 10⁻⁴) and a dyadic vector x with support
n < 200–300. In Arb these give Q⁺(x) > 0 for the Einstein spectrum (λ_k, Ω = w_hi·T). So
Einstein(λ_k, 300·w_hi K) is unstable at 300 K.

Define Ω_E(λ) := inf{Ω : Einstein(λ, Ω) is unstable at some T ≥ 300 K}. Then:

- Theorem NM gives Ω_E ≥ Ω_lo := 300·sqrt(m₀/λ_k).
- The Rayleigh certificate gives Ω_E ≤ Ω_hi.

For λ = ∞ the Rayleigh certificate is for an Einstein atom at w = 2⁻⁶, with m slightly above
g(∞).

### 2.5 Soundness of the lookup rule [THM]

The rule (01 §0.4) is: REJECT if m_hi < g_lo(λ_up)·s². Here λ_up is the smallest tabulated
λ_k ≥ λ_hi (λ_up = ∞ if λ_hi > 100 or unknown), s = T* / 300 K, and T* = Tc_target.

**Claim.** If the rule fires and λ_true ≤ λ_hi and m_true ≤ m_hi, then Tc < T*.

**Proof.**

1. Write m₀ := m₀(λ_up), so that g_lo(λ_up)s² = m₀T*².
2. Since m_true ≤ m_hi < m₀T*², there is ε > 0 with m_true < m₀T² for all T ≥ T* − ε.
3. Also λ_true ≤ λ_up. So Theorem NM makes every T ≥ T* − ε stable.
4. Hence every unstable T is < T* − ε, and Tc = sup{T unstable} ≤ T* − ε < T*. ∎

This argument does not need Tc to be attained. It also holds when the unstable set is empty, in
which case Tc = sup ∅ is taken as 0.

Only the certified grid values are used. **No monotonicity of g, Ω_E or Tc is needed**, because
the theorem already covers every λ ≤ λ_k. The "≤" comes from β ≥ 0 in (2). The T-supremum step is
exact, because m/T² decreases in T.

**Remarks (proved, not needed for the rule).**

- **G is non-increasing [THM].** Let G(λ) := inf{m/T² : some α²F with coupling λ is unstable at T}.
  To see that G is non-increasing, take α²F unstable with coupling λ' < λ.
  - Scaling α²F by s > 1 makes ρ > 0 strictly, because Q_s(x) = s·q(x) − Σ(2n+1)x_n² and
    q(x) > 0 at near-maximizers.
  - Then add a soft atom of coupling λ − sλ' at frequency ε. This changes the form by at most
    4L_ε|x|² → 0 (row-sum bound) and adds (λ − sλ')ε² → 0 to m.
- **Einstein Tc is non-decreasing in λ at fixed Ω [THM].** The same scaling argument applies,
  because l(k) = λ·φ_k(Ω, T) with φ independent of λ. This was the statement proposed in the brief;
  with Theorem NM it is not needed.
- **Soft-mode bound for all spectra [THM-C].** The λ = ∞ certificate says that every α²F satisfies
  Tc ≤ C_cert·sqrt(λ⟨ω²⟩), with C_cert = 0.18272717 = C∞(1 + 5.0×10⁻⁶).
  - **Relation to KAY III.** This is Kiessling, Altshuler & Yuzbashyan, J. Stat. Phys. 192, 93
    (2025), arXiv:2409.02121. R03 checked the full text: their kernel λΩ²/(Ω² + (ω_n − ω_m)²) is our
    l(k) at μ* = 0, and C∞ = (1/2π)√𝔤(2) = 0.1827262477. They define Tc as the crossing
    λ = L_E(ϖ), which is unique for λ ≥ λ* < 0.767. Our sup-definition gives a Tc at least as large,
    so our bound implies the bound in their definition.
  - **What this proves.** The certificate proves their Conjecture 1 (Tc < C∞Ω√λ for Einstein
    spectra) up to a factor 1 + 5.0×10⁻⁶, in non-strict form. It is not the strict conjecture.
    (An earlier δ = 10⁻⁶ certificate gave the factor 1 + 5×10⁻⁷.)
  - **Extension to all spectra.** The extension to every α²F is an elementary consequence of the
    affinity at fixed m (eq. (1)). It is not a deep result.
  - **Prior work.** The large-λ asymptote Tc ≈ 0.183·sqrt(λ⟨ω²⟩) for arbitrary shapes is
    classical: Allen & Dynes, Phys. Rev. B 12, 905 (1975), already non-rigorously.

---

## 3. Coverage (NM-G)

- **Coverage.** Every λ_k is closed on all of w ∈ (0, ∞). This uses cells (44–68 per λ_k, 1661 in
  total), the analytic ends on both sides, and the m = 0 endpoint. **No open boxes.**
- **Slack.** The relative slack in m₀ is δ = 10⁻⁵ for every λ_k. This followed R03 finding 8,
  which asked for robustness.
  - The minimum scaled cell margins are about 10⁻⁶ for most λ_k.
  - Two are smaller: 1.1×10⁻⁷ at λ = 0.8, and 7×10⁻⁹ at λ = 2.75 (a less well centred β).
  - These are far above the 128-bit working precision.
- **Cell budget.** N = 120, K = 600, η = 10⁻⁷. Each λ_k took 10–18 s on one core.
- **No two-atom search.** No two-atom parameter box is searched. The Lagrangian argument bounds
  every spectrum directly, so the two-atom reduction of 01 §3.3 is no longer used.

---

## 4. N1μ: μ*-dependent A floors [THM-C]

### 4.1 Statement

> For μ* and ω_c as tabulated, with the symmetric Coulomb block (01 §1), every α²F with
> A ≤ A_floor is stable at **every** T ≥ 300 K.

Hence Tc ≥ 300 K ⇒ A > A_floor, provided μ* is at least the stated value. The rule depends on a
*lower* bound on μ*, so 01 keeps it as a ranking rule, not a reject rule.

### 4.2 Proof

N_c(T) = #{n ≥ 0 : (2n+1)πT ≤ ω_c} equals j on the interval T ∈ (ω_c/((2j+1)π), ω_c/((2j−1)π)].
For j = 0 this interval is T > ω_c/π. Split T ≥ 300 K into pieces:

- **T ≥ T_u := A_floor/a0** (a0 = 4.33148145). Here A/T ≤ a0, which is stable by (F2). In
  particular the segment j = 0 is covered: ω_c/π ≥ T_u is checked in Arb by both checkers, and it
  holds for every listed ω_c.
- **Each j with T_lo(j) := max(300, ω_c/((2j+1)π)) < T_u.** On that segment A/T < a_j, where
  a_j := A_floor/T_lo(j) rounded up to a decimal. Arb checks a_j·T_lo(j) ≥ A_floor. The cells and
  ends certify that the Einstein kernel l(k) = 2a_j w/(w²+(2πk)²), with the Coulomb block −2μ*
  on the j×j block, satisfies Q⁺ ≤ −η|x|² for all w > 0.
- **Monotone in a.** Q_a = a·q_w(x) − Σ(2n+1)x_n² − 2μ*(Σ_{n<j}x_n)². If q_w(x) ≥ 0, then
  Q_a ≤ Q_{a_j} for a ≤ a_j. Otherwise Q_a ≤ −|x|². So every a ≤ a_j is covered.
- **General α²F.** Affinity at fixed A (01 §3.1) extends the Einstein statement to every α²F. The
  Coulomb and bare terms do not depend on the spectrum.

Ends:

- **Small w:** Q ≤ (a·Lang(w/2) − 1)|x|², which is increasing in w (wE1 Lemma 7).
- **Large w ≥ W = 2048:** a weighted Schur test with h_n = H for n < n₀ and h_n = 1 for n ≥ n₀.
  It uses l(k) ≤ ε := 2a/W, L ≤ a/2 and Σ_{m≠n}K_nm ≤ 2L + S_n. This gives:
  - rows n < n₀: R ≤ ε − 1 + 2(n₀−1)ε + (a + (n₀−1)ε)/H;
  - rows n ≥ n₀: R ≤ a − (2n₀+1) + 2(H−1)n₀ε.

  We take n₀ = ⌈a/2⌉ + 1 and H = ⌈2a⌉, and check both bounds < −η in Arb. This generalizes wE1
  Lemma 8, which needed a < 5.

This also certifies the T-supremum step of 01 §3.1 for these (μ*, ω_c); it was [NUM] there.

### 4.3 Certified floors (δ = 10⁻⁴ below the float value; Tc_target = 300 K)

| μ* | ω_c (K) | **A_floor (K), certified** | float (01 Table 2) | N_c(300 K) | segments j certified by cells |
|---|---|---|---|---|---|
| 0.10 | 5000 | **1671.3055** | 1671.47 | 3 | 3, 2 |
| 0.10 | 10000 | **1705.7914** | 1705.96 | 5 | 5, 4 |
| 0.10 | 20000 | **1703.1664** | 1703.34 | 11 | 11–8 |
| 0.10 | 50000 | **1677.3115** | 1677.48 | 27 | 27–21 |
| 0.13 | 5000 | **1781.4284** | 1781.61 | 3 | 3, 2 |
| 0.13 | 10000 | **1823.4421** | 1823.62 | 5 | 5, 4 |
| 0.13 | 20000 | **1811.1467** | 1811.33 | 11 | 11–8 |
| 0.13 | 50000 | **1768.8811** | 1769.06 | 27 | 27–19 |

- **Coverage.** There are 32 segment certificates (2+2+4+7+2+2+4+9), with 2864 cells (86–97 each),
  and **no open boxes**. All other T, including j = 0, lie in T ≥ T_u and are covered by (F2).
- **Only these four ω_c.** A_floor is **not monotone** in ω_c, so an intermediate cutoff needs its
  own run, which takes about 1 min.
- **Small cutoffs.** For ω_c < π·T_u ≈ 1300 K the j = 0 segment would need its own certificate.
  Both checkers detect this case.
- **μ* must be a lower bound.** The μ* used must be a lower bound on the true μ*, at the stated
  cutoff and in the symmetric-block convention.

---

## 5. Remaining gaps

1. **Model level.** Everything holds within 01 §1 only:
   - isotropic, linearized, adiabatic ME;
   - **symmetric Coulomb block** (for μ* > 0; finding 3);
   - Tc := sup of unstable T;
   - a single α²F, so M-ISO and the G-VALID gates stay.
2. **Steps taken from the certified core [CIT-CORE].** These are (F1), the odd sector and μ*
   monotonicity, and (F2), the universal a0 certificate used in N1μ.
3. **No Lean formalization of the new lemmas.** The following are proved on paper only. R03
   re-derived each of them, but they are not machine-checked:
   - Lemma E (the ends);
   - the Lagrangian reduction (§2.2);
   - the shift and Coulomb additions to the cell lemma;
   - the weighted large-w lemma (§4.2).
4. **Numerical trust.** The numerics are now checked by two independent implementations (§6). Both
   rely on python-flint 0.9.0/Arb, so a bug in Arb itself or in python-flint is not excluded.
5. **N1μ coverage.** N1μ is certified only at the listed (μ*, ω_c), and it still needs a lower bound
   on μ*.
6. **Grid.** NM-G is certified only at the tabulated λ_k. The lookup is sound (§2.5) but
   conservative by up to one grid step. Each extra grid point costs about 15 s.
7. **Not covered.** Multiband and anisotropic NM (01 §6.3) are not addressed here.

---

## 6. The independent checker `verify_indep.py` (R03 finding 2)

The script was written from the statements in §1–§4. It imports none of `arbcore`, `cellcert`,
`families`, `cover`, `run_nm`, `run_n1mu` or `verify`. It uses only python-flint (`arb`,
`arb_mat`), numpy (for float midpoints, never in a verdict) and the JSON files.

Differences from the generator, chosen on purpose:

- **Positive-definiteness proof.** It uses an exact-interval **LDLᵀ factorisation in Arb**: every
  pivot must have a positive lower bound. The generator instead uses a float Cholesky factor plus
  an Arb residual check.
- **Truncation.** N = 110 and K = 700 (generator: 120 and 600), and η = 10⁻⁹.
- **Tail sums.** These are bounded only by integral comparison; it does not use the closed form
  for L.
- **Far coupling** (m ≥ N+K). It uses Cauchy–Schwarz with column sums ≤ 2N·l(K+1).
- **End lemmas.** These are re-derived independently. For the A family at small w it uses
  Lang(u) ≤ u/3 instead of evaluating coth.
- **Segment bookkeeping.** N_c(300 K) is recomputed in Arb, the j = 0 segment is checked, and
  a_j·T_lo ≥ A_floor is verified.
- **Bisection.** A stored cell that fails may be bisected, up to 8 levels. A refinement certifies
  the same statement, and the number of bisections is reported.
- **Not re-checked.** It does not re-check the Rayleigh (upper-side) certificates; no rule uses them.

Result: see `docs/theory/certificates/verify_indep_log.txt` (summary in §7).

---

## 7. Response to review R03

| finding | change |
|---|---|
| 1 counts | Corrected. After regeneration at δ = 10⁻⁵: NM 32 certificates / 1661 cells; N1μ 32 segment certificates / 2864 cells. `summarize.py` now prints the totals. |
| 2 single implementation | Added `verify_indep.py` (§6). All 40 certificate files re-checked with it: see the result line below. |
| 3 Coulomb convention | Every statement with "every μ* ≥ 0 / every cutoff" is now scoped to the symmetric-block convention (header, §0, §1 F1, §2.2, §4, §5). |
| 4 checker preconditions | `verify.py` and `verify_indep.py` check β ≥ 0 (= 0 at λ = ∞), the Rayleigh λ = certificate λ, and g_lo ≤ 90000·m₀ by exact rational comparison. This check showed that the earlier JSON floats were rounded to nearest and could exceed the exact value by 1 ulp, so the generator now rounds g_lo and Ω_lo down and g_hi up. |
| 5 lookup wording | §2.5 now uses the strict-inequality and ε argument and no longer equates Tc ≥ T* with the existence of an unstable T ≥ T*. |
| 6 j = 0 segment | Both checkers now require ω_c/π ≥ A_floor/a0. It holds for all 8 files. |
| 7 KAY wording | §2.5 states the verified facts: same kernel at μ* = 0, the Tc-definition difference, and "proves Conjecture 1 up to a factor 1 + 5.0×10⁻⁶ (non-strict)". The extension to all spectra is described as elementary, and Allen–Dynes 1975 is cited. |
| 8 β existence / margins | The heuristic status of the existence of β is stated in §2.2. All NM certificates were regenerated at δ = 10⁻⁵. |

Independent-checker result (`certificates/verify_indep_log.txt`): **40/40 files PASS, "INDEPENDENT CHECK (40 files): ALL VERIFIED"**. All 1661 NM cells and 2864 N1μ cells closed with **0 bisections**, in about 15 min on one core.

Generator-code re-check (`certificates/verify_log.txt`): **40/40 PASS, "ALL CERTIFICATES VERIFIED"**, including the new precondition checks and the j = 0 check.

---

## 8. Reproduction

```
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
# generate (~10 min NM per process, ~10 min N1mu, one core each)
nice -n 19 python src/certify/run_nm.py --deltas 1e-5,1e-4 --max_time 600 --lams inf,100.0,50.0,20.0,15.0,10.0,8.0,7.0,6.0,5.0,4.5,4.0,3.75,3.5,3.25,3.0,2.75,2.5,2.25,2.0,1.9,1.8,1.7,1.6,1.5,1.4,1.3,1.2,1.1,1.0,0.9,0.8 --out docs/theory/certificates
for mu in 0.10 0.13; do for wc in 5000 10000 20000 50000; do
  nice -n 19 python src/certify/run_n1mu.py --mu $mu --wc $wc --out docs/theory/certificates; done; done
nice -n 19 python src/certify/summarize.py
# verify (two independent checkers) and controls
nice -n 19 python src/certify/verify.py
nice -n 19 python src/certify/verify_indep.py
nice -n 19 python src/certify/controls.py
```
