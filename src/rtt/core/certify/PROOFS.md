# rtt.core.certify: proofs of the rigorous checks

Team CERT. This file proves every analytic inequality used by the Arb checks in
`arbcore.py`, `cells.py`, `rayleigh.py` and `classify.py`. For each rule it also names the Lean
theorems (library `PipelineLogic` in `lean/`; the correspondence is tested by `tests/core/test_correspondence.py`) that
turn the checked inequalities into the verdict.

Units: k_B = ħ = 1, so all energies are in kelvin. Notation:

- Nominal spectrum: atoms Ω_j > 0 with weights p_j ≥ 0. Envelope e ∈ [0, 1): every spectrum p'
  with (1−e)p_j ≤ p'_j ≤ (1+e)p_j (Lean `InEnv`, ERRATA F5: weights only).
- g_T(k, ω) = 2ω/(ω² + (2πkT)²) and λ(k; T) = Σ_j p_j g_T(k, Ω_j) for k ≥ 1.
- ω_n = 2n+1 and S_L(n) = Σ_{k=1}^{n} L(k). λ(0) never appears (ERRATA F2).
- f_p(m) = (2m+1)^{−p}, p ∈ {2, 3}. The test-vector tail is v_m = C f_p(m) for m ≥ N.

**Trust base.** We trust:

- the correctness of Arb/FLINT through python-flint 0.9.0: `arb` arithmetic, `arb.pi`,
  `arb_poly` multiplication, which returns balls containing the exact coefficients, and certain
  comparisons;
- exact conversion of IEEE doubles and Python ints to `arb`/`fmpq`;
- exact `fmpq` arithmetic;
- CPython's exact comparison of doubles and ints, used for the structural and tiling checks.

Floats (numpy/scipy) only *propose* data. A check's verdict depends only on its certificate
and on the spectrum. **Do not run with `python -O`**: no `assert` is used for decisions, but keep
the interpreter standard.

---

## 0. Rules, what is checked, and the Lean theorems used

| Rule | Python check | Lean theorems that turn the checks into the conclusion |
|---|---|---|
| R-UNIV | `A_hi_q < a0_q · T*_q` in exact rationals, where `A_hi = (1+e)Σp_j` is exact | `R_UNIV`, `univ_stable`, `stable_of_A_le`, `TcR1.stable_of_ge` (hypothesis `H(a0, η)` = PROOF.md Prop. 6) |
| R-CELLS (default `mhat`) | `check_cell` per cell (Thm C), `Tiles` with exact double equality, `A_hi ≤ a0·t_end` exactly (ERRATA D1) | `CellOK` (def), `cell_stable`, `cells_stable`, `R_CELLS`, `envelope_form_bound`, `schur_row_test_inf`, `rows_of_eventually`, `rowSum_mono`, `lam_box`, `cover` |
| R-MU | none (Lemma 3) | `R_MU`, `stable_of_even_bound` |
| R-ODD | none | `TcR1.L6` via `stable_of_even_bound` |
| C-RAYLEIGH | `check_rayleigh` (Thm R) | `rayleigh_cert`, `boxMin_le`, `boxMin_attained`, `sel_mul_eq_min`, `Qp_anti_mu`, `rayleigh_unstable` |
| classify / S5 | `classify_evidence` (literal transcription) | `classifyAlt`, `classify`, `classify_rejected_sound`, `classify_rejected_sound'`, `classify_certified_sound`, `classify_eq_rejected_iff`, `classify_eq_certified_iff`, `classify_cases`, `evidence_consistent` |
| R-CELLS variant `cform` (opt-in) | `check_cell` with `variant="cform"` (Thm C2) | **not formalised.** Needs a new Lean lemma, "c-form envelope" (statement in §4). Verdicts that use it record the assumption `R-CELLS[cform]` |

`classify` derives `AltEvidence(univFires, cellsBelowTarget, rayleighAtOrAboveTarget)` for each
alternative spectrum from the stored certificate. It uses only the checks above, so each `True` is
backed by `UnivCert` / `CellsCert` / `RayleighCert` exactly as defined in CORRESPONDENCE_LEAN §P7.

- `cellsBelowTarget` requires a non-empty tiling, because Lean `Tiles _ _ [] = False`, and
  `T_up < T*` strictly (ERRATA F4).
- `rayleighAtOrAboveTarget` requires the certificate's `mu_hi` and `omega_c` to equal the run's
  values, so the same cutoff rule is used (ERRATA F3).
- `classify([]) = UNDECIDED` (ERRATA F1).
- The combined bracket is `[min T_low, max T_up]` (ERRATA F7).

When the cells do not certify, an alternative's `T_up` is `t_end`. This is valid by
`univ_stable`, given the D1 check.

---

## 1. Coupling function: enclosures

**Lemma T0 (decay).** For T ≥ T1 > 0, k ≥ 1 and p' in the envelope:

    λ'(k; T) ≤ (1+e) Σ_j p_j Ω_j / (2π² T1² k²) =: b / k².

The same bound holds for `hi(k) = (1+e)λ(k; T1)`.

*Proof.* We have g_T(k, ω) ≤ 2ω/(2πkT)² ≤ 2ω/(2πkT1)², and p'_j ≤ (1+e)p_j. ∎

**Lemma M (monotonicity, Lean `g_antitone`, `lam_box`).** For k ≥ 1, g_T(k, ω) is decreasing in
T > 0 and decreasing in k. Hence, for T ∈ [T1, T2] and p' in the envelope,

    lo(k) := (1−e)λ(k; T2) ≤ λ'(k; T) ≤ hi(k) := (1+e)λ(k; T1),

and hi and lo are both decreasing in k. ∎

**Lemma S1 (moment series for λ, `ArbSpectrum.lam`).** Let x = 2πkT, Ω_max = max_j Ω_j,
r² = Ω_max²/x² < 1, t_j = Ω_j/Ω_max ≤ 1 and μ_s = Σ_j p_j t_j^{2s+1}. Then for every R ≥ 1,

    λ(k; T) ∈ (2Ω_max/x²) · [ Σ_{s<R} (−1)^s μ_s r^{2s} ± μ_R r^{2R} ].

*Proof.* For each atom, g = (2Ω_j/x²) Σ_{s≥0} (−q_j)^s with q_j = Ω_j²/x² = t_j² r² < 1. This is
an alternating series with decreasing terms, so the remainder after R terms has absolute value
≤ (2Ω_j/x²) q_j^R. We have 2p_j Ω_j^{2s+1}/x^{2s+2} = (2Ω_max/x²) p_j t_j^{2s+1} r^{2s}. Summing
over atoms with weights p_j ≥ 0 gives the claim. ∎

The code checks `r² < 1` in Arb. It uses the series only when 2πkT ≥ 2Ω_max (float proposal; the
lemma needs only r² < 1). It chooses R so that r^{2R} ≤ 2^{−110}, and adds the remainder as a
ball. For smaller k it uses the direct sum Σ_j 2p_jΩ_j/(Ω_j² + x²) in Arb.

---

## 2. The row identity (c-form)

**Lemma Rw.** Let L(k) ≥ 0 with L(k) ≤ b/k² for k ≥ 1, and let v be a bounded sequence. Then for
every n ≥ 0 the following series converge absolutely and

    Σ_{m≥0, m≠n} L(|n−m|) v_m + Σ_{m≥0} L(n+m+1) v_m − 2 S_L(n) v_n  =  Σ_{k≥1} c_{n,k}(v) L(k),

    c_{n,k}(v) = v_{n−k} + v_{n+k} − 2 v_n   (1 ≤ k ≤ n),
    c_{n,k}(v) = v_{n+k} + v_{k−n−1}         (k > n).

*Proof.* Collect the terms by k. In the Toeplitz sum, m = n−k (for k ≤ n) and m = n+k each carry
L(k). In the Hankel sum, m = k−n−1 ≥ 0 iff k ≥ n+1, and the term m = n is the diagonal entry
L(2n+1)v_n. The term −2S_L(n)v_n = −Σ_{k≤n} 2L(k)v_n. Absolute convergence follows from
|v| ≤ V and Σ L(k) ≤ b ζ(2), which justifies the rearrangement. (This is wA Lemma 5 and wE2 F5
with a general L.) ∎

With `Mq L` from Lean, the full row sum is Σ_m (Mq L)_{nm} v_m = Σ_k c_{n,k}L(k) − ω_n v_n. The
Lean `Mhat(hi, lo)` row is

    Σ_m Mhat_{nm} v_m = Σ_k c_{n,k}(v) hi(k) − ω_n v_n + 2 v_n (S_hi(n) − S_lo(n)).          (Rw')

---

## 3. Theorem C: a certified `mhat` cell gives Lean `CellOK`

Cell certificate: T1, T2 (doubles, 0 < T1 ≤ T2 < ∞), N ≥ 1, head u_0..u_{N−1} > 0 (doubles),
C > 0, η > 0, p ∈ {2, 3}, N_e ≥ max(2N, 2), M_v ≥ N_e + 2, K := N_e + M_v. The test vector is
v_m = u_m for m < N and v_m = C f_p(m) for m ≥ N. So v_n > 0 for every n, which Lean's `CellOK`
requires.

**Theorem C.** If `check_cell(cert)["ok"]`, then for every n ≥ 0

    Σ_{m≥0} Mhat(hi, lo)_{nm} v_m ≤ −η v_n,     hi = (1+e)λ(·;T1),  lo = (1−e)λ(·;T2).

Since Mhat's off-diagonal entries are ≥ 0 and v > 0, partial row sums increase (Lean
`rowSum_mono`). So every truncated row sum is ≤ −η v_n (`rows_of_eventually`), which is exactly
`CellOK κ sp e T1 T2` with κ = 2π. Lean `cell_stable` then gives Stable(T) for every T ∈ [T1, T2],
every envelope spectrum, every μ* ≥ 0 and every cutoff.

*Proof, head rows n < N_e.* The check evaluates, in Arb,

    row_n = TH_n + R_n − 2 v_n S_lo(n) − ω_n v_n,   then tests  row_n + η v_n < 0.

1. TH_n = P1[n] + P2[M_v−1−n] + P2[M_v+n]. Here P1 = Λ·V and P2 = Λ·V^rev, with
   Λ(z) = Σ_{k=0}^{K} hi(k) z^k, hi(0) := 0, V(z) = Σ_{m<M_v} v_m z^m and V^rev its reversal.
   - P1[n] = Σ_{m≤n} hi(n−m) v_m.
   - P2[M_v−1−n] = Σ_{k≥1} hi(k) v_{n+k} over n+k < M_v.
   - P2[M_v+n] = Σ_{m<M_v} hi(n+m+1) v_m.

   All indices k involved are ≤ N_e + M_v − 1 ≤ K. So TH_n is the exact sum of the
   Mhat off-diagonal terms and the Hankel diagonal term hi(2n+1)v_n over the columns m < M_v.
2. **Lemma R1 (column remainder).** Σ_{m≥M_v} [hi(m−n) + hi(n+m+1)] C f_p(m) ≤ R_n, where
   R_n := 2C f_p(M_v) · Tsuf(M_v − n) and Tsuf(J) := Σ_{k=J}^{K} hi(k) + b/K.
   *Proof.* f_p decreases, so f_p(m) ≤ f_p(M_v). Since hi decreases (Lemma M) and
   n+m+1 ≥ m−n ≥ M_v−n ≥ 2, both sums are ≤ C f_p(M_v) Σ_{k≥M_v−n} hi(k). Finally
   Σ_{k>K} hi(k) ≤ Σ_{k>K} b/k² ≤ b/K (Lemma T0). ∎
   **Lemma R2 (second column bound; the check uses min(R1, R2)).** With J = M_v − n ≥ 2,
   Σ_{m≥M_v} [hi(m−n) + hi(n+m+1)] C f_p(m) ≤ 2bC·2^{−p} [J^{−(p+2)} + J^{−(p+1)}/(p+1)].
   *Proof.* By Lemma T0, hi(m−n), hi(n+m+1) ≤ b/(m−n)². Also 2m+1 ≥ 2(m−n), so
   f_p(m) ≤ 2^{−p}(m−n)^{−p}. The summand is therefore ≤ 2bC 2^{−p} j^{−(p+2)} with j = m−n ≥ J.
   Since j^{−(p+2)} is decreasing, Σ_{j≥J} j^{−(p+2)} ≤ J^{−(p+2)} + ∫_J^∞ x^{−(p+2)} dx =
   J^{−(p+2)} + J^{−(p+1)}/(p+1). ∎
3. The remaining diagonal part, −2S_lo(n) − ω_n, is computed exactly. The sum runs over
   k = 1..n (ERRATA F2).

So the full row Σ_m Mhat_{nm} v_m ≤ row_n, and the check gives ≤ −η v_n.

*Proof, tail rows n ≥ N_e.* **Lemma T1(p).** Let p ∈ {2, 3}, N ≥ 1, n ≥ N_e ≥ 2N,
U0 = Σ_{m<N} u_m, r = U0/C and Z := (2N)^{1−p}/(2(p−1)) ≥ Z_N := Σ_{j≥N} f_p(j). For any
L with 0 ≤ L(k) ≤ b/k²:

    Σ_k c_{n,k}(v) L(k) ≤ b C [ 2p(p+1)/(n+1)^{p+1} + 4Z/n² ] + b U0/(n−N+1)²
                        + b C f_p(2n+1)/n + b (U0 + C Z)/(n+1)².

*Proof.* Split the k-sum into five pieces.

- (a1) 1 ≤ k ≤ ⌊n/2⌋. Here n−k ≥ n/2 ≥ N, so every index is in the tail, and
  c_{n,k} = C Δ²_k f_p(n) = C ∫_{−k}^{k} (k−|t|) f_p''(n+t) dt ≤ C k² f_p''(n−k). Since
  f_p'' = 4p(p+1)(2x+1)^{−p−2} is decreasing, f_p''(n−k) ≤ f_p''(n/2) = 4p(p+1)(n+1)^{−p−2}.
  With k² L(k) ≤ b and ⌊n/2⌋ ≤ n/2, the piece is ≤ 2p(p+1) b C n/(n+1)^{p+2}
  ≤ 2p(p+1) b C/(n+1)^{p+1}.
- (a2) n/2 < k ≤ n−N. Here c_{n,k} = C[f(n−k) + f(n+k) − 2f(n)] ≤ C f(n−k), because
  f(n+k) ≤ f(n). Also L(k) ≤ b/k² < 4b/n², and Σ_{n−k ≥ N} f(n−k) ≤ Z_N. So the piece is
  ≤ 4bCZ/n².
- (b) n−N < k ≤ n. Put m = n−k < N. Then c_{n,k} = u_m + C f(n+k) − 2C f(n) ≤ u_m, so
  L(k) c_{n,k} ≤ L(k) u_m (L ≥ 0, u_m > 0). With L(n−m) ≤ b/(n−N+1)², the piece is
  ≤ bU0/(n−N+1)².
- (c1) k > n, term v_{n+k} = C f(n+k) ≤ C f(2n+1). Using Σ_{k>n} b/k² ≤ b/n, the piece is
  ≤ bC f(2n+1)/n.
- (c2) k > n, term v_{k−n−1} = v_m with m ≥ 0. Here L(n+1+m) ≤ b/(n+1)², so the piece is
  ≤ b(U0 + CZ)/(n+1)².

The integral bound uses Z_N ≤ ∫_{N−1/2}^{∞} (2x+1)^{−p} dx = (2N)^{1−p}/(2(p−1)), because
(2x+1)^{−p} is convex, so each term is at most its unit-interval integral. ∎

Apply (Rw') with L = hi. Lemma T0 gives 0 ≤ hi ≤ b/k². Also

    S_hi(n) − S_lo(n) ≤ dL := Σ_{k=1}^{N_e} (hi(k) − lo(k)) + b/N_e,

since lo ≥ 0 for e < 1 and Σ_{k>N_e} hi(k) ≤ b/N_e. Write D(n) := 2n+1 − η − 2dL. The row
inequality Σ_m Mhat_{nm} v_m ≤ −η v_n then follows from

    G(n) := ((2n+1)/D(n)) · b · [ 2p(p+1)(2n+1)^{p−1}/(n+1)^{p+1} + 4Z(2n+1)^{p−1}/n²
            + r(2n+1)^{p−1}/(n−N+1)² + (2n+1)^{p−1}/((4n+3)^p n) + last(n) ] ≤ 1,

where last(n) = 4(r+Z) for p = 3, which bounds the increasing term (r+Z)(2n+1)²/(n+1)², and
last(n) = (r+Z)(2n+1)/(n+1)² for p = 2. To see this, divide Lemma T1 plus 2v_n dL by
C f_p(n) D(n) = C(2n+1)^{−p} D(n), with D(n) > 0.

**Monotonicity.** For n ≥ 2N ≥ 2, each bracketed term is nonincreasing in n:

- 2p(p+1)(2n+1)^{p−1}/(n+1)^{p+1}: the log-derivative is (p − 4n − 3)/((2n+1)(n+1)) < 0.
- (2n+1)^{p−1}/n²: obvious for p ≤ 3.
- ((2n+1)/(n−N+1))^{·}: the derivative of (2n+1)/(n−N+1) is (1−2N)/(n−N+1)² < 0.
- (2n+1)^{p−1}/((4n+3)^p n): the log-derivative is 2(p−1)/(2n+1) − 4p/(4n+3) − 1/n < 0.
- last(n): constant for p = 3; (2n+1)/(n+1)² decreases for p = 2.

Also (2n+1)/D(n) = (2n+1)/(2n+1−c) with c ≥ 0 is nonincreasing. Hence G(n) ≤ G(N_e) for every
n ≥ N_e. The check verifies D(N_e) > 0 and G(N_e) < 1 in Arb. Since D(n) increases, D(n) > 0
for every n ≥ N_e. ∎

This settles every row, so Theorem C holds.

**Cover (Lean `Tiles`, `cover`, ERRATA D1).** `evidence_for_alternative` checks with exact double
equality:

- `cells[0].T1 == T_up`;
- `cells[i].T2 == cells[i+1].T1`;
- `cells[-1].T2 == t_end`, with at least one cell;
- `A_hi ≤ a0 · t_end` in exact rationals.

Lean `cells_stable` then gives Stable(T) for every T ≥ T_up, every envelope spectrum, every
μ* ≥ 0 and every cutoff. So Tc_model ≤ T_up.

---

## 4. Theorem C2: the opt-in `cform` cell. Hand proof, not in Lean

**Lemma C2a.** For any L with lo ≤ L ≤ hi and 0 ≤ L ≤ b/k²,

    Σ_k c_{n,k} L(k) ≤ Σ_k c_{n,k} hi(k) + Σ_{k≤n} c⁻_{n,k} (hi(k) − lo(k)),   c⁻ = max(−c, 0).

*Proof.* cL = c·hi − c(hi − L), and −c(hi − L) ≤ c⁻(hi − L) ≤ c⁻(hi − lo). For k > n,
c_{n,k} ≥ 0 because v > 0. ∎

**Lemma C2b.** For 1 ≤ k ≤ n,

    c_{n,k} = Σ_{|j|<k} (k−|j|) δ²v(n+j),   where δ²v(m) = v_{m−1} + v_{m+1} − 2v_m.

Hence c⁻_{n,k} ≤ Σ_{|j|<k} (k−|j|) ν(n+j), where ν = (δ²v)⁻. Since the tail f_p is convex,
ν(m) = 0 for m ≥ N+1. Summing,

    Σ_{k≤n} c⁻_{n,k} w(k) ≤ Σ_i ν(i) W_n(|i−n|),   W_n(d) = Σ_{k=d+1}^{n} (k−d) w(k),  w = hi − lo.

*Proof.* The first identity holds by induction on k (telescoping). For the second, exchange the
finite sums. ∎

The check uses ν_i = |δ²v(i)| as a ball whenever δ²v(i) ≥ 0 is not certain. This is a valid upper
bound. It computes W_n(d) exactly from suffix sums for d < n, and W_n(d) = 0 for d ≥ n.

**Theorem C2.** If `check_cell(cert)` with `variant = "cform"` returns ok, then for every L in the
cell box and every n, Σ_k c_{n,k}(v) L(k) − ω_n v_n ≤ −η v_n. Head rows use Lemmas C2a, C2b and
R1. Tail rows use Lemma T1, which bounds Σ c⁺ L ≥ Σ c L for every 0 ≤ L ≤ b/k². By Lemma Rw this
is the full row of the *true* matrix Mq(λ'_T). The Schur test (Lean `schur_row_test_inf`, via
`rowSum_mono` on Mq, whose off-diagonal entries are ≥ 0) then gives Stable(T).

**Missing Lean lemma (for the Lean team).** State and prove:

    CellOK_cform κ sp e T1 T2 :≡ 0 ≤ T1 ∧ ∃ η > 0, ∃ v > 0,
      ∀ L, (lo ≤ L ≤ hi) → ∀ n, Σ' m, Mq L n m * v m ≤ −η v n.

This is an infinite-row statement. It needs Lemma Rw (tsum rearrangement) and Lemmas C2a/C2b.
Then `cell_stable` follows for `CellOK_cform` exactly as for `CellOK`.

Status: implemented and tested. It is not used by default: with the default `row_check = "mhat"`
every verdict is Lean-backed. In tests it did not beat `mhat`. The reason is structural: any
smooth eigenvector-like test vector is asymptotically critical for Lemma T1, because the head
mass of the true eigenvector satisfies 8bU_tot/C → 1. So the tail must be slightly heavier
than the head's continuation, which creates a concave kink at m = N. On wide cells the kink's
correction W_n is large.

---

## 5. Theorem R: C-RAYLEIGH

Certificate: T (double, > 0), μ_hi ≥ 0, ω_c > 0, N_c and x ∈ ℝ^{N_x} (doubles, not all 0).

**Cutoff.** N_c(T) = #{n ≥ 0 : (2n+1)πT ≤ ω_c}. The left side is increasing in n. The check
verifies in Arb that (2(N_c−1)+1)πT ≤ ω_c (if N_c ≥ 1) and (2N_c+1)πT > ω_c, and that N_c equals
the claimed value. Equality is impossible for T > 0 rational because π is irrational, so 128-bit
Arb decides in practice. If it does not, the check fails.

**Form.** For finitely supported x (support < N_x), Q⁺_T[p'](x) = C(x) − 2μ(Σ_{n<min(N_c,N_x)} x_n)²
+ Σ_j p'_j F_j, where

    F_j = Σ_{k=1}^{2N_x−1} g_T(k, Ω_j) γ_k,
    γ_k = 2Σ_n x_n x_{n+k} [k ≤ N_x−1] + Σ_{n+m=k−1} x_n x_m − 2Σ_{n≥k} x_n².

This is Lemma Rw summed against x, i.e. Lean `Qp_env_eq` and `L1`. Here γ_k collects the
Toeplitz pairs (both orders), the Hankel pairs, and −2Σ_n x_n² S_n = −2Σ_k L(k) Σ_{n≥k} x_n². It
is computed exactly with `arb_poly` products (autocorrelation and square).

**Lemma S2 (per-atom series).** For k with 2πkT ≥ 2Ω_max (checked: q_k = Ω_max²/x_k² ≤ 1/4),

    g_T(k, Ω_j) = (2/Ω_max) Σ_{s<R} (−1)^s t_j^{2s+1} q_k^{s+1} ± (2/Ω_max) t_j^{2R+1} q_k^{R+1},

with the same alternating-series proof as Lemma S1. Hence

    Σ_{k≥k_s} g γ_k ∈ (2/Ω_max) [ Σ_s (−1)^s t_j^{2s+1} Γ_s ± t_j^{2R+1} E ],
    Γ_s = Σ_k γ_k q_k^{s+1},   E = Σ_k |γ_k| q_k^{R+1}.

R = 40, so the remainder is ≤ 4^{−40}·E relative.

**Box minimum (Lean `boxMinQ`, `sel_mul_eq_min`, ERRATA F6).** The check computes

    boxMinQ = −Σ ω_n x_n² − 2μ_hi(Σ_{n<min(N_c,N_x)} x_n)² + Σ_j min((1−e)p_j F_j, (1+e)p_j F_j)

with Arb `min`. The ball of the min of two balls contains the min of every pair of points, so this
is valid even when F_j's ball contains 0. **Theorem R:** if boxMinQ ≥ 0 is certain, then Lean
`rayleigh_cert` (with `boxMin_le` and `Qp_anti_mu`) gives Unstable(Q⁺_{T,μ,N_c(T)}[p']) for every
p' in the envelope and every μ ≤ μ_hi. So T is not stable and Tc_model ≥ T, for the block-cutoff
rule with this ω_c (ERRATA F3, recorded in `Verdict.assumptions`). ∎

---

## 5b. R-MULTI (ARCH §4c; wMB Prop. N **audited SOUND**; opt-in via `rejection_model="multiband"`, not the default)

> **PITFALL (non-conservative). Read before supplying band data.** 𝒜_ij must be in the wMB
> Theorem 1 convention, where band i's mass enhancement is the **row sum**,
> λ_i = Σ_j λ_ij, and 𝒜_ij = ∫α²F_ij. Some codes report band-resolved α²F_ij as
> "contributions" normalised to sum to the isotropic α²F, i.e. W/N_tot. With such data κ comes
> out **too small**, and a superconductor can be **wrongly REJECTED**.
> `multiband.normalisation_refusal` therefore refuses any data whose provenance has a
> normalisation/convention key mentioning "contribution", "isotropic" or "N_tot", or an
> assumption mentioning "contribution". R-MULTI then does not fire. The refusal is covered by
> the digest, so it cannot be stripped from a stored certificate.
>
> **Envelope semantics.** `rel_envelope = e` means true 𝒜_ij ∈ [(1−e), (1+e)] × nominal 𝒜_ij
> (only the upper side is used). If an error estimate instead says nominal = (1−δ)·true, the
> correct value is e = δ/(1−δ), not δ. The envelope covers 𝒜, not N.

Code: `multiband.py`. Input: `rtt.core.dfpt.band_coupling.BandCoupling` (team DFT). It is read by
duck typing, and the check uses only `kind`, `A_K`, `N`, `row_sum_upper_K` and `rel_envelope`.
The float helpers are never used. Lean (library `MultiBand`, hypotheses H(a0, η) and PropN):
`R_UNIV_MB`, `R_UNIV_MB_rowsum`, `sub_of_K`, `sub_env`, `thm1_even`, `thm3_odd`.

**Route "eig" (`R_UNIV_MB` with `sub_of_K` and `sub_env`).**

1. W_ij = N_i A_ij is computed as exact rationals.
2. If W ≠ Wᵀ exactly, W̄ = (W + Wᵀ)/2 is used, and the verdict carries the assumption
   "W symmetrized". The verdict is then about the model with W̄. For `entrywise_upper_bound`
   data this is still an upper bound on the true symmetric W:
   W_true = W_trueᵀ ≤ min(W, Wᵀ) ≤ W̄.
3. K_ij = W̄_ij/(r_i r_j), with r_i = √N_i as Arb balls, so K is enclosed entrywise.
4. Given doubles u_i > 0, the check certifies (1+e)·(Ku)_i ≤ κ_hi u_i for every i in Arb.

Then K u ≤ κ u with κ = κ_hi/(1+e). Lean `sub_of_K` gives 𝒜v ≤ κv with v = u/r > 0, and
`sub_env` gives 𝒜'v ≤ κ_hi v for every 𝒜' ≤ (1+e)𝒜 entrywise. (The proposal is
u = (κ_c − K)^{−1}·1 with κ_c slightly above the float λ_max(K). It is positive even for
reducible K. Untrusted.)

**Route "cw" (direct Collatz–Wielandt vector on 𝒜; used only when no DOS N is supplied).**
Given doubles v_i > 0, the check certifies (1+e)·(𝒜v)_i ≤ κ_hi v_i in exact rationals.
`R_UNIV_MB` takes any v > 0 with 𝒜v ≤ κv, and `sub_env` handles the envelope. The recorded
assumption is "detailed balance with some N > 0", because the theorem's model needs one.

**Route "rowsum" (`R_UNIV_MB_rowsum`, v = 1).** The check verifies
(1+e)·max_i Σ_j A_ij ≤ κ_hi in exact rationals, or uses the stated `row_sum_upper_K` for kind
`row_sum_upper_bound`. Upper-bound kinds are safe for rejection, since κ and the row sums are
monotone in 𝒜.

**Decision.** `fires` requires κ_hi < a0·T*, strict and in exact rationals. The certified
multiband upper bound is T_up = κ_hi/a0, rounded up to a double. This holds because, by
`R_UNIV_MB` applied with any T' > κ_hi/a0 in place of T*, every T > κ_hi/a0 is stable.

**Data kinds (referee R3).** For `kind = "row_sum_upper_bound"` only the declared R_i
(`row_sum_upper_K`) are used, in both the proposal and the check. Any eig, cw or A_K route is
rejected for that kind, because A_K has no defined meaning there.

**Robustness (referee R4).** A non-finite or non-positive T* returns "not fired". It never
raises an exception.

**Recorded assumptions.**

- "M-MULTIBAND (isotropic per band; wMB Thm 1; Prop. N audited SOUND)";
- "A_ij in the lambda_i = sum_j lambda_ij convention (not W/N_tot contributions)";
- "Coulomb N_i mu*_ij PSD";
- "E-ENV-MB(e)";
- "W symmetrized" when it applies;
- the data kind;
- every assumption carried in `BandCoupling.provenance["assumptions"]`, for example
  "M-BAND-ISO". It is propagated, never dropped;
- (referee R1) every `provenance["caveats"]` entry as "DATA CAVEAT: …", the stated derived bound
  as "DATA BOUND: …", and for EPW data "EPW: omega_max bounds all fine-grid phonon frequencies;
  modes with wf <= eps_acustic excluded". Upper-bound kinds are called upper bounds only under
  these recorded conditions;
- "detailed balance with the same N_i for the true model; the envelope covers A, not N" for the
  eig route, and "detailed balance with some N > 0" for the cw and rowsum routes.

**Binding of data to the verdict (referee R2).** The digest covers kind, 𝒜, N, R, e, the
refusal and every data assumption. At recheck the payload is rebuilt from the *supplied*
`BandCoupling`, and the evidence is re-derived from it. Every assumption it implies must appear
in the stored `Verdict.assumptions`. So stripping, for example, "M-BAND-ISO" makes recheck fail.

**Combination (`combine_multiband`).**

- `rejection_model="isotropic"` (default): the status is the Lean-mirrored single-band status.
  Every REJECTED verdict carries the assumption "M-ISO-SINGLEBAND" and the flag
  "MULTIBAND_RISK" (ARCH §4c).
- `rejection_model="multiband"`:
  - REJECTED iff the list of band couplings is non-empty, R-MULTI fires on every element (the S5
    conjunction), and the single-band status is not CERTIFIED;
  - CERTIFIED iff the single-band status is CERTIFIED and R-MULTI does not reject. Such a verdict
    carries "M-ISO-SINGLEBAND", because the lower bound is single-band;
  - otherwise UNDECIDED. In particular, without band data a rejection is impossible.

There is no Lean `classify` for the multiband combination yet. The logic above is a direct
transcription of ARCH §4c.

---

## 6. R-UNIV and the a0 constant

A_hi = (1+e)·Σ_j p_j is computed as an exact rational (`fmpq`) from the doubles. a0 = 433148145/10⁸
is exact. The check is `A_hi < a0·T*` in exact rationals. Lean `R_UNIV` needs `(1+e)A ≤ A_hi`,
which holds with equality, and `A_hi < a0 T*`. The hypothesis H(a0 = 4.33148145, η = 10⁻⁹) is
PROOF.md Proposition 6 (wE1 certificate, audited SOUND).

---

## 7. Proposal heuristics (not part of any proof)

These are documented only for reproducibility. None is trusted.

- **Cell head.** u solves (M_NN + η'I)u = −t, with M = Mhat(hi, lo) in floats truncated at N,
  η' = |ρ_N|/2, and t the float tail drive. The margin is η = |ρ_N|/4. N starts at
  max(48, 3·Ω_max/(2πT1)) and doubles on failure up to `N_cap`.
- **Sizes.** N_e is the smallest n ≥ 2N with float G(n) ≤ 0.9. M_v makes the float remainder
  ≤ 2% of η·v_n.
- **Cover.** Top-down from t_end.
  - Phase 1 tries the single decision cell [T*(1−2⁻³⁰), t_end] before any float Tc is
    computed. This is cheap when Tc ≪ T*.
  - Phase 2 then removes a fraction s of the distance to the float Tc of the (1+e)-scaled
    spectrum at each step. s grows ×1.6 after a success and shrinks ÷3 after a failure. Phase 2
    stops at Tc_float·(1 + rel_tol), after the time budget, or when s < 10⁻³.
  - `mode="decide"` skips phase 2.
- **Float Tc.** A Cholesky test of −M is used, first at truncation 3·k0 (scan and bisection to
  10⁻³), then at 10·k0. The truncated ρ is nondecreasing in N, so a coarse-unstable T is
  fine-unstable.
- **Rayleigh.** x is the float top eigenvector at the (1−e)-scaled spectrum with μ_hi and N_c(T).
  T = Tc_float_lo·(1 − d) with d = rel_tol, 4·rel_tol, …
