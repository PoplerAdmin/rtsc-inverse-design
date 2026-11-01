/-
P1 (weighted Schur / Collatz–Wielandt row test) and P2 (interval envelope).
-/
import PipelineLogic.Common

open Finset

namespace PipelineLogic

open TcR1

noncomputable section

/-- AM–GM with weights: `M x_n x_m ≤ (M/2)(x_n² v_m/v_n) + (M/2)(x_m² v_n/v_m)` when either
`n = m`, or `M ≥ 0` (any `n, m`), for `v_n, v_m > 0`. -/
theorem term_le (Mnm xn xm vn vm : ℝ) (hvn : 0 < vn) (hvm : 0 < vm) (hM : 0 ≤ Mnm) :
    Mnm * xn * xm ≤ Mnm / 2 * (xn ^ 2 * vm / vn) + Mnm / 2 * (xm ^ 2 * vn / vm) := by
  have key : xn ^ 2 * vm / vn + xm ^ 2 * vn / vm - 2 * (xn * xm)
      = (xn * vm - xm * vn) ^ 2 / (vn * vm) := by
    field_simp; ring
  have h0 : 0 ≤ (xn * vm - xm * vn) ^ 2 / (vn * vm) := by positivity
  have h1 : 2 * (xn * xm) ≤ xn ^ 2 * vm / vn + xm ^ 2 * vn / vm := by linarith
  have h2 := mul_le_mul_of_nonneg_left h1 hM
  nlinarith [h2]

/-- **P1 (finite).** Weighted Schur / Collatz–Wielandt row test on `range N`.
`M` symmetric with nonnegative off-diagonal entries, `v_n > 0`, and
`∑_{m<N} M_nm v_m ≤ -η v_n` for all `n < N` imply
`∑_{n,m<N} M_nm x_n x_m ≤ -η ∑_{n<N} x_n²` for every `x`. -/
theorem schur_row_test (M : ℕ → ℕ → ℝ) (v : ℕ → ℝ) (η : ℝ) (N : ℕ)
    (hsym : ∀ n ∈ range N, ∀ m ∈ range N, M n m = M m n)
    (hoff : ∀ n ∈ range N, ∀ m ∈ range N, n ≠ m → 0 ≤ M n m)
    (hv : ∀ n ∈ range N, 0 < v n)
    (hrow : ∀ n ∈ range N, rowSum M v N n ≤ -η * v n) (x : ℕ → ℝ) :
    ∑ n ∈ range N, ∑ m ∈ range N, M n m * x n * x m ≤ -η * nsq N x := by
  set A : ℕ → ℕ → ℝ := fun n m => M n m / 2 * (x n ^ 2 * v m / v n) with hAdef
  set B : ℕ → ℕ → ℝ := fun n m => M n m / 2 * (x m ^ 2 * v n / v m) with hBdef
  -- termwise bound
  have hterm : ∀ n ∈ range N, ∀ m ∈ range N, M n m * x n * x m ≤ A n m + B n m := by
    intro n hn m hm
    by_cases h : n = m
    · subst h
      have hvn := hv n hn
      simp only [hAdef, hBdef]
      field_simp
      ring_nf
      exact le_refl _
    · exact term_le _ _ _ _ _ (hv n hn) (hv m hm) (hoff n hn m hm h)
  have h1 : ∑ n ∈ range N, ∑ m ∈ range N, M n m * x n * x m
      ≤ ∑ n ∈ range N, ∑ m ∈ range N, (A n m + B n m) :=
    sum_le_sum (fun n hn => sum_le_sum (fun m hm => hterm n hn m hm))
  -- the B-part equals the A-part by symmetry
  have hBA : ∑ n ∈ range N, ∑ m ∈ range N, B n m = ∑ n ∈ range N, ∑ m ∈ range N, A n m := by
    rw [sum_comm]
    apply sum_congr rfl; intro n hn
    apply sum_congr rfl; intro m hm
    simp only [hAdef, hBdef]
    rw [hsym m hm n hn]
  -- the A-part is a weighted sum of row sums
  have hArow : ∀ n ∈ range N, ∑ m ∈ range N, A n m = (x n ^ 2 / v n) / 2 * rowSum M v N n := by
    intro n hn
    unfold rowSum
    rw [mul_sum]
    apply sum_congr rfl; intro m _
    simp only [hAdef]
    field_simp
  have h2 : ∑ n ∈ range N, ∑ m ∈ range N, (A n m + B n m)
      = ∑ n ∈ range N, (x n ^ 2 / v n) * rowSum M v N n := by
    simp only [sum_add_distrib]
    rw [hBA, ← two_mul, sum_congr rfl hArow, mul_sum]
    apply sum_congr rfl; intro n _; ring
  have h3 : ∑ n ∈ range N, (x n ^ 2 / v n) * rowSum M v N n
      ≤ ∑ n ∈ range N, (x n ^ 2 / v n) * (-η * v n) := by
    apply sum_le_sum; intro n hn
    have : 0 ≤ x n ^ 2 / v n := div_nonneg (sq_nonneg _) (hv n hn).le
    exact mul_le_mul_of_nonneg_left (hrow n hn) this
  have h4 : ∑ n ∈ range N, (x n ^ 2 / v n) * (-η * v n) = -η * nsq N x := by
    unfold nsq
    rw [mul_sum]
    apply sum_congr rfl; intro n hn
    have := (hv n hn).ne'
    field_simp
  linarith

/-- Truncated row sums increase with the truncation once the diagonal is included
(off-diagonals `≥ 0`, `v ≥ 0`). Hence a bound on the full infinite row sum (the limit of these
partial sums) bounds every partial sum `N > n`. -/
theorem rowSum_mono (M : ℕ → ℕ → ℝ) (v : ℕ → ℝ) (hoff : ∀ n m, n ≠ m → 0 ≤ M n m)
    (hv : ∀ n, 0 ≤ v n) {n N N' : ℕ} (hn : n < N) (hN : N ≤ N') :
    rowSum M v N n ≤ rowSum M v N' n := by
  unfold rowSum
  apply sum_le_sum_of_subset_of_nonneg (range_subset_range.mpr hN)
  intro m _ hm
  have : n ≠ m := by
    intro h; subst h; exact hm (mem_range.mpr hn)
  exact mul_nonneg (hoff n m this) (hv m)

/-- A row bound valid for all sufficiently large truncations (i.e. for the infinite row sum)
holds for every truncation `N > n`. -/
theorem rows_of_eventually (M : ℕ → ℕ → ℝ) (v : ℕ → ℝ) (η : ℝ)
    (hoff : ∀ n m, n ≠ m → 0 ≤ M n m) (hv : ∀ n, 0 ≤ v n)
    (hev : ∀ n, ∃ N0, ∀ N, N0 ≤ N → rowSum M v N n ≤ -η * v n) :
    ∀ N, ∀ n < N, rowSum M v N n ≤ -η * v n := by
  intro N n hn
  obtain ⟨N0, hN0⟩ := hev n
  calc rowSum M v N n ≤ rowSum M v (max N N0) n :=
        rowSum_mono M v hoff hv hn (le_max_left _ _)
    _ ≤ -η * v n := hN0 _ (le_max_right _ _)

/-- **P1 (infinite matrix, truncation form).** If the infinite symmetric matrix `M` has
nonnegative off-diagonal entries, `v_n > 0` for all `n`, and every truncated row satisfies
`∑_{m<N} M_nm v_m ≤ -η v_n` (`n < N`), then for every `N` and `x`,
`∑_{n,m<N} M_nm x_n x_m ≤ -η |x|²`. -/
theorem schur_row_test_inf (M : ℕ → ℕ → ℝ) (v : ℕ → ℝ) (η : ℝ)
    (hsym : ∀ n m, M n m = M m n) (hoff : ∀ n m, n ≠ m → 0 ≤ M n m) (hv : ∀ n, 0 < v n)
    (hrow : ∀ N, ∀ n < N, rowSum M v N n ≤ -η * v n) (N : ℕ) (x : ℕ → ℝ) :
    ∑ n ∈ range N, ∑ m ∈ range N, M n m * x n * x m ≤ -η * nsq N x :=
  schur_row_test M v η N (fun n _ m _ => hsym n m) (fun n _ m _ h => hoff n m h)
    (fun n _ => hv n) (fun n hn => hrow N n (mem_range.mp hn)) x

/-! ## P2: interval envelope -/

/-- **P2 (row comparison).** If `lo ≤ L ≤ hi` pointwise and `v ≥ 0`, the true truncated row
sums of `Mq L` are bounded by those of the envelope `Mhat hi lo`. -/
theorem rowSum_le_envelope (L lo hi : ℕ → ℝ) (hlo : ∀ k, lo k ≤ L k) (hhi : ∀ k, L k ≤ hi k)
    (v : ℕ → ℝ) (hv : ∀ m, 0 ≤ v m) (N n : ℕ) :
    rowSum (Mq L) v N n ≤ rowSum (Mhat hi lo) v N n := by
  unfold rowSum
  exact sum_le_sum (fun m _ => mul_le_mul_of_nonneg_right (Mq_le_Mhat L lo hi hlo hhi n m) (hv m))

/-- **P2 + P1.** For every coupling sequence `L` in the box `[lo, hi]` (with `L(k) ≥ 0` for
`k ≥ 1`, which holds for every spectrum), if the *envelope* rows satisfy
`∑_{m<N} Mhat_nm v_m ≤ -η v_n` for all `N`, `n < N`, with `v > 0`, then the μ* = 0 even form
satisfies `C(x) + F⁺_L(x) ≤ -η |x|²` for all `N`, `x`. -/
theorem envelope_form_bound (L lo hi : ℕ → ℝ) (hL : ∀ k, 1 ≤ k → 0 ≤ L k)
    (hlo : ∀ k, lo k ≤ L k) (hhi : ∀ k, L k ≤ hi k) (v : ℕ → ℝ) (hv : ∀ n, 0 < v n) (η : ℝ)
    (hrow : ∀ N, ∀ n < N, rowSum (Mhat hi lo) v N n ≤ -η * v n) (N : ℕ) (x : ℕ → ℝ) :
    C N x + Fp L N x ≤ -η * nsq N x := by
  rw [form_eq_Mq]
  apply schur_row_test_inf (Mq L) v η (Mq_symm L) (Mq_offdiag_nonneg L hL) hv _ N x
  intro N n hn
  exact le_trans (rowSum_le_envelope L lo hi hlo hhi v (fun m => (hv m).le) N n) (hrow N n hn)

/-- From an even-sector μ* = 0 bound `Q⁺_T(x) ≤ -η|x|²` (η > 0) for all `x`, both sectors are
stable for every `μ* ≥ 0` and every cutoff (Lemma 3 = R-MU, and TcR1.L6 for the odd sector). -/
theorem stable_of_even_bound {ι : Type*} (κ : ℝ) (sp : Spectrum ι) (T η : ℝ) (hη : 0 < η)
    (h : ∀ (N : ℕ) (x : ℕ → ℝ), C N x + Fp (lam κ sp T) N x ≤ -η * nsq N x)
    (μ : ℝ) (hμ : 0 ≤ μ) (Nc : WithTop ℕ) : Stable κ sp μ Nc T := by
  constructor
  · apply not_unstable_of_bound _ η hη
    intro N x
    have := L3 μ hμ Nc N x
    have := h N x
    unfold Qp; linarith
  · apply not_unstable_of_bound _ η hη
    intro N x
    have h6 := L6 κ sp ⊤ T N x
    have h0 := h N (fun n => |x n|)
    rw [nsq_abs] at h0
    unfold Qp at h6
    rw [muTerm_zero] at h6
    linarith

/-- **R-MU.** Even-sector stability at μ* = 0 implies it at every μ* ≥ 0 (same cutoff). -/
theorem R_MU {ι : Type*} (κ : ℝ) (sp : Spectrum ι) (Nc : WithTop ℕ) (T μ : ℝ) (hμ : 0 ≤ μ)
    (h0 : ¬ Unstable (Qp κ sp 0 Nc T)) : ¬ Unstable (Qp κ sp μ Nc T) := by
  intro hU
  apply h0
  intro ε hε
  obtain ⟨N, x, hx, hlt⟩ := hU ε hε
  refine ⟨N, x, hx, ?_⟩
  have := L3 μ hμ Nc N x
  unfold Qp at hlt ⊢
  rw [muTerm_zero]
  linarith

end

end PipelineLogic
