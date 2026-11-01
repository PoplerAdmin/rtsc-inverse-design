/-
MultiBand.Odd: Theorem 3 of `wMB/PROOF.md` (the odd multiband sector is never unstable).

The proof here is elementary (a Schur/Gershgorin row-sum test), not the symbol-positivity
argument of wMB §4.5: the odd kernel `b_nm = l(|n-m|) - l(n+m+1)` is entrywise `≥ 0`
(`l` decreasing) and every truncated row sum satisfies `∑_{m<N} b_nm ≤ D_n`.
-/
import MultiBand.Identities

open Finset

namespace MultiBand

open TcR1

noncomputable section

variable {B : ℕ} {α : Type*}

/-- Prefix sums `F(K) = ∑_{k<K} l(k)`. -/
def pre (l : ℕ → ℝ) (K : ℕ) : ℝ := ∑ k ∈ range K, l k

theorem pre_mono (l : ℕ → ℝ) (hl : ∀ k, 0 ≤ l k) {a b : ℕ} (h : a ≤ b) : pre l a ≤ pre l b := by
  unfold pre
  exact sum_le_sum_of_subset_of_nonneg (range_subset_range.mpr h) (fun k _ _ => hl k)

theorem pre_nonneg (l : ℕ → ℝ) (hl : ∀ k, 0 ≤ l k) (a : ℕ) : 0 ≤ pre l a :=
  sum_nonneg (fun k _ => hl k)

theorem pre_succ (l : ℕ → ℝ) (K : ℕ) : pre l (K + 1) = pre l K + l K := by
  unfold pre; rw [sum_range_succ]

theorem sum_shift (l : ℕ → ℝ) (n N : ℕ) :
    ∑ m ∈ range N, l (n + m + 1) = pre l (n + N + 1) - pre l (n + 1) := by
  induction N with
  | zero => simp
  | succ N ih =>
    rw [sum_range_succ, ih]
    have : n + (N + 1) + 1 = (n + N + 1) + 1 := by omega
    rw [this, pre_succ l (n + N + 1)]; ring

theorem sum_adist_le (l : ℕ → ℝ) (n : ℕ) :
    ∀ N, N ≤ n + 1 → ∑ m ∈ range N, l (adist n m) = pre l (n + 1) - pre l (n + 1 - N) := by
  intro N
  induction N with
  | zero => intro _; simp
  | succ N ih =>
    intro h
    rw [sum_range_succ, ih (by omega)]
    have h1 : adist n N = n - N := by unfold adist; omega
    have h2 : n + 1 - N = (n + 1 - (N + 1)) + 1 := by omega
    have h3 : n + 1 - (N + 1) = n - N := by omega
    rw [h1, h2, pre_succ l (n + 1 - (N + 1)), h3]; ring

theorem sum_adist_ge (l : ℕ → ℝ) (n : ℕ) :
    ∀ N, n + 1 ≤ N → ∑ m ∈ range N, l (adist n m) = pre l (n + 1) + pre l (N - n) - l 0 := by
  intro N hN
  induction N, hN using Nat.le_induction with
  | base =>
    rw [sum_adist_le l n (n + 1) le_rfl]
    have : n + 1 - n = 0 + 1 := by omega
    rw [Nat.sub_self, this, pre_succ]
    simp [pre]
  | succ N hN ih =>
    rw [sum_range_succ, ih]
    have h1 : adist n N = N - n := by unfold adist; omega
    have h2 : N + 1 - n = (N - n) + 1 := by omega
    rw [h1, h2, pre_succ l (N - n)]; ring

theorem Ssum_eq_pre (l : ℕ → ℝ) (n : ℕ) : Ssum l n = pre l (n + 1) - l 0 := by
  unfold Ssum pre
  rw [sum_range_succ']; ring

/-- Truncated row sums of the odd kernel are bounded by the Z-diagonal:
`∑_{m<N} (l(|n-m|) - l(n+m+1)) ≤ l(0) + 2 ∑_{k=1}^n l(k)` whenever `l ≥ 0`. -/
theorem odd_row_le (l : ℕ → ℝ) (hl : ∀ k, 0 ≤ l k) (n N : ℕ) :
    ∑ m ∈ range N, (l (adist n m) - l (n + m + 1)) ≤ l 0 + 2 * Ssum l n := by
  rw [sum_sub_distrib, sum_shift, Ssum_eq_pre]
  have hp1 : l 0 ≤ pre l (n + N + 1) := by
    have := pre_mono l hl (show 1 ≤ n + N + 1 by omega)
    simpa [pre] using this
  rcases le_total N (n + 1) with h | h
  · rw [sum_adist_le l n N h]
    linarith [pre_nonneg l hl (n + 1 - N)]
  · rw [sum_adist_ge l n N h]
    linarith [pre_mono l hl (show N - n ≤ n + N + 1 by omega)]

/-- The odd kernel is entrywise nonnegative for `w > 0`. -/
theorem bm_nonneg (tp w : ℝ) (hw : 0 < w) (n m : ℕ) :
    0 ≤ lw tp w (adist n m) - lw tp w (n + m + 1) := by
  have : adist n m ≤ n + m + 1 := by unfold adist; omega
  linarith [lw_anti tp w hw this]

theorem bm_symm (tp w : ℝ) (n m : ℕ) :
    lw tp w (adist n m) - lw tp w (n + m + 1) = lw tp w (adist m n) - lw tp w (m + n + 1) := by
  rw [adist_comm' n m, Nat.add_comm n m]

/-- Per-atom odd λ-part is `≤ 0`: for a symmetric nonnegative `W` and `w > 0`,
`∑_ij W_ij [φ_iᵀ B⁻_w φ_j - φ_iᵀ D_w φ_i] ≤ 0`. -/
theorem odd_atom_le (tp w : ℝ) (hw : 0 < w) (W : Fin B → Fin B → ℝ)
    (hW0 : ∀ i j, 0 ≤ W i j) (hW : ∀ i j, W i j = W j i) (N : ℕ) (φ : Fin B → ℕ → ℝ) :
    ∑ i, ∑ j, W i j * (Bm tp w N (φ i) (φ j) - Dq tp w N (φ i) (φ i)) ≤ 0 := by
  set b : ℕ → ℕ → ℝ := fun n m => lw tp w (adist n m) - lw tp w (n + m + 1) with hb
  have hb0 : ∀ n m, 0 ≤ b n m := fun n m => bm_nonneg tp w hw n m
  have hbs : ∀ n m, b n m = b m n := fun n m => bm_symm tp w n m
  -- R(y) = ∑_n (∑_m b_nm) y_n²
  set R : (ℕ → ℝ) → ℝ := fun y => ∑ n ∈ range N, ∑ m ∈ range N, b n m * y n ^ 2 with hR
  -- AM-GM, termwise
  have hAM : ∀ i j, Bm tp w N (φ i) (φ j) ≤ (R (φ i) + R (φ j)) / 2 := by
    intro i j
    have e2 : R (φ j) = ∑ n ∈ range N, ∑ m ∈ range N, b n m * φ j m ^ 2 := by
      simp only [hR]
      rw [sum_comm]
      apply sum_congr rfl; intro n _; apply sum_congr rfl; intro m _
      rw [hbs]
    rw [e2]
    simp only [hR]
    unfold Bm
    rw [← sum_add_distrib, div_eq_mul_inv, sum_mul]
    apply sum_le_sum; intro n _
    rw [← sum_add_distrib, sum_mul]
    apply sum_le_sum; intro m _
    have := mul_le_mul_of_nonneg_left
      (show φ i n * φ j m ≤ (φ i n ^ 2 + φ j m ^ 2) / 2 by nlinarith [sq_nonneg (φ i n - φ j m)])
      (hb0 n m)
    simp only [hb] at this ⊢
    nlinarith [this]
  have hRD : ∀ i, R (φ i) ≤ Dq tp w N (φ i) (φ i) := by
    intro i
    simp only [hR]
    unfold Dq
    apply sum_le_sum; intro n _
    rw [← sum_mul]
    have h1 := odd_row_le (lw tp w) (fun k => lw_nonneg tp w hw.le k) n N
    have h2 : Dn tp w n * φ i n * φ i n = Dn tp w n * φ i n ^ 2 := by ring
    rw [h2]
    unfold Dn
    exact mul_le_mul_of_nonneg_right h1 (sq_nonneg _)
  have hsw := swap_sum W hW (fun i _ => R (φ i))
  have step : ∑ i, ∑ j, W i j * (Bm tp w N (φ i) (φ j) - Dq tp w N (φ i) (φ i))
      ≤ ∑ i, ∑ j, W i j * ((R (φ i) + R (φ j)) / 2 - R (φ i)) := by
    apply sum_le_sum; intro i _; apply sum_le_sum; intro j _
    apply mul_le_mul_of_nonneg_left _ (hW0 i j)
    linarith [hAM i j, hRD i]
  have hzero : ∑ i, ∑ j, W i j * ((R (φ i) + R (φ j)) / 2 - R (φ i)) = 0 := by
    have e : ∑ i, ∑ j, W i j * ((R (φ i) + R (φ j)) / 2 - R (φ i))
        = (1 / 2) * ∑ i, ∑ j, W i j * R (φ j) - (1 / 2) * ∑ i, ∑ j, W i j * R (φ i) := by
      rw [mul_sum, mul_sum, ← sum_sub_distrib]
      apply sum_congr rfl; intro i _
      rw [mul_sum, mul_sum, ← sum_sub_distrib]
      apply sum_congr rfl; intro j _
      ring
    rw [e, hsw]; ring
  linarith

/-- **Theorem 3 (odd sector).** For `T > 0` and `N_i ≥ 0`:
`Q⁻(φ) ≤ -∑_i N_i C(φ_i) ≤ -|x|²`. No hypothesis on the spectrum beyond symmetry and
nonnegativity, and no use of `H` or Proposition N. -/
theorem thm3_odd (tp : ℝ) (sp : MSpec B α) (dos : Fin B → ℝ) (hdos : ∀ i, 0 ≤ dos i)
    (T : ℝ) (hT : 0 < T) (N : ℕ) (φ : Fin B → ℕ → ℝ) :
    Qminus tp sp dos T N φ ≤ -nrm2 dos N φ := by
  unfold Qminus nrm2
  have h1 : ∑ s ∈ sp.S, ∑ i, ∑ j, sp.W s i j *
      (Bm tp (sp.Ω s / T) N (φ i) (φ j) - Dq tp (sp.Ω s / T) N (φ i) (φ i)) ≤ 0 := by
    apply sum_nonpos; intro s hs
    exact odd_atom_le tp _ (div_pos (sp.Ω_pos s hs) hT) (sp.W s) (sp.W_nonneg s hs)
      (sp.W_symm s hs) N φ
  have h2 : (1 / T) * ∑ s ∈ sp.S, ∑ i, ∑ j, sp.W s i j *
      (Bm tp (sp.Ω s / T) N (φ i) (φ j) - Dq tp (sp.Ω s / T) N (φ i) (φ i)) ≤ 0 :=
    mul_nonpos_of_nonneg_of_nonpos (by positivity) h1
  have h3 : ∑ i, dos i * nsq N (φ i) ≤ ∑ i, dos i * Cw N (φ i) :=
    sum_le_sum (fun i _ => mul_le_mul_of_nonneg_left (Cw_ge_nsq N _) (hdos i))
  linarith

end

end MultiBand
