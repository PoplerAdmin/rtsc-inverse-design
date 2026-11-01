/-
MultiBand.Algebra: bilinear-form algebra for the single-band kernels
(link to TcR1's `Fp`/`Qe`, gauge scaling, polarization, parallelogram law).
-/
import MultiBand.Defs

open Finset

namespace MultiBand

open TcR1

noncomputable section

/-- Generic bilinear form `yᵀ K z` on `n, m < N`. -/
def bil (K : ℕ → ℕ → ℝ) (N : ℕ) (y z : ℕ → ℝ) : ℝ :=
  ∑ n ∈ range N, ∑ m ∈ range N, K n m * y n * z m

/-- Kernel of `G⁺_w`: `l(|n-m|) + l(n+m+1) - δ_nm D_{w,n}`. -/
def KG (tp w : ℝ) (n m : ℕ) : ℝ :=
  (lw tp w (adist n m) + lw tp w (n + m + 1)) - if n = m then Dn tp w n else 0

theorem adist_comm' (n m : ℕ) : adist n m = adist m n := by unfold adist; omega

theorem KG_symm (tp w : ℝ) (n m : ℕ) : KG tp w n m = KG tp w m n := by
  unfold KG
  by_cases h : n = m
  · subst h; rfl
  · have h' : m ≠ n := fun h'' => h h''.symm
    simp only [h, h', ite_false, adist_comm' n m, Nat.add_comm n m]

theorem Dq_eq_sum (tp w : ℝ) (N : ℕ) (y z : ℕ → ℝ) :
    Dq tp w N y z = ∑ n ∈ range N, ∑ m ∈ range N,
      (if n = m then Dn tp w n else 0) * y n * z m := by
  unfold Dq
  apply sum_congr rfl; intro n hn
  have : ∀ m ∈ range N, (if n = m then Dn tp w n else 0) * y n * z m
      = if n = m then Dn tp w n * y n * z n else 0 := by
    intro m _; split_ifs with h
    · subst h; rfl
    · ring
  rw [sum_congr rfl this, sum_ite_eq (range N) n, ite_eq_left hn]

theorem Gp_eq_bil (tp w : ℝ) (N : ℕ) (y z : ℕ → ℝ) :
    Gp tp w N y z = bil (KG tp w) N y z := by
  unfold Gp Bp bil KG
  rw [Dq_eq_sum, ← sum_sub_distrib]
  apply sum_congr rfl; intro n _
  rw [← sum_sub_distrib]
  apply sum_congr rfl; intro m _
  ring

theorem bil_smul (K : ℕ → ℕ → ℝ) (N : ℕ) (c d : ℝ) (y z : ℕ → ℝ) :
    bil K N (fun n => c * y n) (fun n => d * z n) = c * d * bil K N y z := by
  unfold bil
  rw [mul_sum]; apply sum_congr rfl; intro n _
  rw [mul_sum]; apply sum_congr rfl; intro m _
  ring

/-- Polarization for a symmetric kernel: `K(ψ,ψ) - K(χ,χ) = K(a,b)` with
`ψ = (a+b)/2`, `χ = (a-b)/2`. -/
theorem bil_polar (K : ℕ → ℕ → ℝ) (hK : ∀ n m, K n m = K m n) (N : ℕ) (a b : ℕ → ℝ) :
    bil K N (fun n => (a n + b n) / 2) (fun n => (a n + b n) / 2)
      - bil K N (fun n => (a n - b n) / 2) (fun n => (a n - b n) / 2) = bil K N a b := by
  unfold bil
  have h1 : ∑ n ∈ range N, ∑ m ∈ range N, K n m * ((a n + b n) / 2) * ((a m + b m) / 2)
      - ∑ n ∈ range N, ∑ m ∈ range N, K n m * ((a n - b n) / 2) * ((a m - b m) / 2)
      = (1 / 2) * ∑ n ∈ range N, ∑ m ∈ range N, K n m * a n * b m
        + (1 / 2) * ∑ n ∈ range N, ∑ m ∈ range N, K n m * b n * a m := by
    rw [← sum_sub_distrib, mul_sum, mul_sum, ← sum_add_distrib]
    apply sum_congr rfl; intro n _
    rw [← sum_sub_distrib, mul_sum, mul_sum, ← sum_add_distrib]
    apply sum_congr rfl; intro m _
    ring
  have h2 : ∑ n ∈ range N, ∑ m ∈ range N, K n m * b n * a m
      = ∑ n ∈ range N, ∑ m ∈ range N, K n m * a n * b m := by
    rw [sum_comm]
    apply sum_congr rfl; intro n _
    apply sum_congr rfl; intro m _
    rw [hK m n]; ring
  rw [h1, h2]; ring

/-- `Gp` in TcR1 form: `G⁺_w(y,y) = Fp (l_w) y` (the `l_w(0)` terms cancel). -/
theorem Gp_self (tp w : ℝ) (N : ℕ) (y : ℕ → ℝ) :
    Gp tp w N y y = Fp (lw tp w) N y := by
  rw [Fp_spec_form (lw tp w) (lw tp w 0)]
  unfold Gp Bp Dq Dn
  congr 1
  · apply sum_congr rfl; intro n _
    apply sum_congr rfl; intro m _
    split_ifs with h
    · rw [h]
    · rfl
  · apply sum_congr rfl; intro n _; ring

/-- The single-band Einstein form of wMB, `Q^E_{a,w}(y) = -C(y) + a G⁺_w(y)`, is TcR1's `Qe`
(for every real `a`, in particular `a < 0`). -/
theorem Qe_eq (tp a w : ℝ) (N : ℕ) (y : ℕ → ℝ) :
    Qe tp a w N y = -Cw N y + a * Gp tp w N y y := by
  rw [Qe_linear, Gp_self]
  unfold C Cw; ring

theorem Cw_eq (N : ℕ) (y : ℕ → ℝ) : Cw N y = -C N y := by unfold C Cw; ring

theorem Cw_ge_nsq (N : ℕ) (y : ℕ → ℝ) : nsq N y ≤ Cw N y := by
  rw [Cw_eq]; linarith [C_le N y]

theorem Cw_nonneg (N : ℕ) (y : ℕ → ℝ) : 0 ≤ Cw N y :=
  le_trans (nsq_nonneg N y) (Cw_ge_nsq N y)

theorem Cw_smul (N : ℕ) (c : ℝ) (y : ℕ → ℝ) : Cw N (fun n => c * y n) = c ^ 2 * Cw N y := by
  unfold Cw; rw [mul_sum]; apply sum_congr rfl; intro n _; ring

/-- Parallelogram law: `C(ψ) + C(χ) = ½(C(a) + C(b))`. -/
theorem Cw_parallelogram (N : ℕ) (a b : ℕ → ℝ) :
    Cw N (fun n => (a n + b n) / 2) + Cw N (fun n => (a n - b n) / 2)
      = (Cw N a + Cw N b) / 2 := by
  unfold Cw
  rw [← sum_add_distrib, ← sum_add_distrib, div_eq_mul_inv, sum_mul]
  apply sum_congr rfl; intro n _; ring

/-- Gauge scaling: `G(φ_i, φ_j) = v_i v_j G(θ_i, θ_j)`. -/
theorem Gp_gauge (tp w : ℝ) (N : ℕ) (vi vj : ℝ) (hi : vi ≠ 0) (hj : vj ≠ 0) (y z : ℕ → ℝ) :
    Gp tp w N y z = vi * vj * Gp tp w N (fun n => y n / vi) (fun n => z n / vj) := by
  rw [Gp_eq_bil, Gp_eq_bil, ← bil_smul]
  congr 1
  · funext n; field_simp
  · funext n; field_simp

/-- Polarization for `G⁺`. -/
theorem Gp_polar (tp w : ℝ) (N : ℕ) (a b : ℕ → ℝ) :
    Gp tp w N (fun n => (a n + b n) / 2) (fun n => (a n + b n) / 2)
      - Gp tp w N (fun n => (a n - b n) / 2) (fun n => (a n - b n) / 2) = Gp tp w N a b := by
  simp only [Gp_eq_bil]
  exact bil_polar _ (KG_symm tp w) N a b

/-- Nonnegativity of `l_w(k)` for `w ≥ 0`. -/
theorem lw_nonneg (tp w : ℝ) (hw : 0 ≤ w) (k : ℕ) : 0 ≤ lw tp w k := by
  unfold lw ge; positivity

theorem Dn_nonneg (tp w : ℝ) (hw : 0 ≤ w) (n : ℕ) : 0 ≤ Dn tp w n := by
  unfold Dn Ssum
  have : 0 ≤ ∑ k ∈ range n, lw tp w (k + 1) := sum_nonneg (fun k _ => lw_nonneg tp w hw _)
  linarith [lw_nonneg tp w hw 0]

/-- `l_w` is antitone in `k` for `w > 0`. -/
theorem lw_anti (tp w : ℝ) (hw : 0 < w) {k k' : ℕ} (h : k ≤ k') : lw tp w k' ≤ lw tp w k := by
  unfold lw ge
  have hk : (k : ℝ) ≤ k' := by exact_mod_cast h
  have h0 : (0 : ℝ) ≤ k := Nat.cast_nonneg k
  have hsq : (tp * k) ^ 2 ≤ (tp * k') ^ 2 := by
    have : (k : ℝ) ^ 2 ≤ (k' : ℝ) ^ 2 := by nlinarith
    have e1 : (tp * k) ^ 2 = tp ^ 2 * (k : ℝ) ^ 2 := by ring
    have e2 : (tp * k') ^ 2 = tp ^ 2 * (k' : ℝ) ^ 2 := by ring
    rw [e1, e2]; exact mul_le_mul_of_nonneg_left this (sq_nonneg tp)
  apply div_le_div_of_nonneg_left (by linarith) (by positivity) (by linarith)

end

end MultiBand
