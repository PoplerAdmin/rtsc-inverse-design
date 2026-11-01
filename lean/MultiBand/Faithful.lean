/-
MultiBand.Faithful: the multiband even form written with the physical couplings
`λ_ij(k) = ∫ 2ω dα²F_ij(ω) / (ω² + (2πkT)²)` (wMB §1), showing that `Qplus` is
`∑_i N_i φ_iᵀ [pairing with the full matrix λ_ij − Z-diagonal with row sums ∑_j λ_ij] φ`,
and the identification `|x|² = ∑_i N_i |φ_i|²` with `x_i = √N_i φ_i`.
-/
import MultiBand.Algebra

open Finset

namespace MultiBand

open TcR1

noncomputable section

variable {B : ℕ} {α : Type*}

/-- `λ_ij(k) = ∑_s α²F_{ij,s} g_T(k, Ω_s)`, `α²F_ij = W_ij / N_i`. -/
def lamMB (tp : ℝ) (sp : MSpec B α) (dos : Fin B → ℝ) (T : ℝ) (i j : Fin B) (k : ℕ) : ℝ :=
  ∑ s ∈ sp.S, (sp.W s i j / dos i) * g tp T k (sp.Ω s)

/-- Pairing minus Z-diagonal for a coupling sequence `L` (uncancelled `L(0)`):
`∑_{nm} (L(|n-m|) + L(n+m+1)) y_n z_m - ∑_n (L(0) + 2∑_{k=1}^n L(k)) y_n²`. -/
def Fbil (L : ℕ → ℝ) (N : ℕ) (y z : ℕ → ℝ) : ℝ :=
  ∑ n ∈ range N, ∑ m ∈ range N, (L (adist n m) + L (n + m + 1)) * y n * z m
    - ∑ n ∈ range N, (L 0 + 2 * Ssum L n) * y n ^ 2

theorem Fbil_add (L1 L2 : ℕ → ℝ) (N : ℕ) (y z : ℕ → ℝ) :
    Fbil (fun k => L1 k + L2 k) N y z = Fbil L1 N y z + Fbil L2 N y z := by
  unfold Fbil Ssum
  have h1 : ∑ n ∈ range N, ∑ m ∈ range N,
      (L1 (adist n m) + L2 (adist n m) + (L1 (n + m + 1) + L2 (n + m + 1))) * y n * z m
      = ∑ n ∈ range N, ∑ m ∈ range N, (L1 (adist n m) + L1 (n + m + 1)) * y n * z m
        + ∑ n ∈ range N, ∑ m ∈ range N, (L2 (adist n m) + L2 (n + m + 1)) * y n * z m := by
    rw [← sum_add_distrib]; apply sum_congr rfl; intro n _
    rw [← sum_add_distrib]; apply sum_congr rfl; intro m _; ring
  have h2 : ∑ n ∈ range N,
      (L1 0 + L2 0 + 2 * ∑ k ∈ range n, (L1 (k + 1) + L2 (k + 1))) * y n ^ 2
      = ∑ n ∈ range N, (L1 0 + 2 * ∑ k ∈ range n, L1 (k + 1)) * y n ^ 2
        + ∑ n ∈ range N, (L2 0 + 2 * ∑ k ∈ range n, L2 (k + 1)) * y n ^ 2 := by
    rw [← sum_add_distrib]; apply sum_congr rfl; intro n _
    rw [sum_add_distrib]; ring
  rw [h1, h2]; ring

theorem Fbil_zero (N : ℕ) (y z : ℕ → ℝ) : Fbil (fun _ => 0) N y z = 0 := by
  simp [Fbil, Ssum]

theorem Fbil_sum {ι : Type*} (S : Finset ι) (f : ι → ℕ → ℝ) (N : ℕ) (y z : ℕ → ℝ) :
    Fbil (fun k => ∑ s ∈ S, f s k) N y z = ∑ s ∈ S, Fbil (f s) N y z := by
  classical
  induction S using Finset.induction_on with
  | empty => simp only [sum_empty]; exact Fbil_zero N y z
  | insert a S ha ih =>
    simp only [sum_insert ha]
    rw [Fbil_add (f a) (fun k => ∑ s ∈ S, f s k), ih]

theorem Fbil_smul (c : ℝ) (L : ℕ → ℝ) (N : ℕ) (y z : ℕ → ℝ) :
    Fbil (fun k => c * L k) N y z = c * Fbil L N y z := by
  unfold Fbil Ssum
  rw [mul_sub, mul_sum, mul_sum]
  congr 1
  · apply sum_congr rfl; intro n _; rw [mul_sum]; apply sum_congr rfl; intro m _; ring
  · apply sum_congr rfl; intro n _; rw [← mul_sum]; ring

theorem Bp_Dq_eq_Fbil (tp w : ℝ) (N : ℕ) (y z : ℕ → ℝ) :
    Bp tp w N y z - Dq tp w N y y = Fbil (lw tp w) N y z := by
  unfold Bp Dq Fbil Dn
  congr 1
  apply sum_congr rfl; intro n _; ring

/-- `N_i λ_ij(k) = (1/T) ∑_s W s i j l_{Ω_s/T}(k)` (from TcR1's `L4_pointwise`). -/
theorem dos_lam (tp : ℝ) (sp : MSpec B α) (dos : Fin B → ℝ) (hdos : ∀ i, dos i ≠ 0) (T : ℝ)
    (hT : 0 < T) (i j : Fin B) (k : ℕ) :
    dos i * lamMB tp sp dos T i j k
      = ∑ s ∈ sp.S, (1 / T) * sp.W s i j * lw tp (sp.Ω s / T) k := by
  unfold lamMB
  rw [mul_sum]
  apply sum_congr rfl; intro s hs
  have hd := hdos i
  have h4 := L4_pointwise tp (sp.W s i j) T (sp.Ω s) k hT (sp.Ω_pos s hs)
  have e : dos i * (sp.W s i j / dos i * g tp T k (sp.Ω s)) = sp.W s i j * g tp T k (sp.Ω s) := by
    field_simp
  rw [e, h4]
  unfold lw; ring

/-- **Faithfulness.** For `T > 0`, `N_i ≠ 0`:
`Q⁺(φ) = -∑_i N_i C(φ_i)
   + ∑_i N_i ∑_j [ ∑_{n,m} (λ_ij(|n-m|) + λ_ij(n+m+1)) φ_{n,i} φ_{m,j}
                    - ∑_n (λ_ij(0) + 2∑_{k=1}^n λ_ij(k)) φ_{n,i}² ] - 2∑_ij N_i μ*_ij S_i S_j`.
Summing the last bracket over `j` gives the Z-renormalization with the row sum
`∑_j λ_ij`, while the pairing uses the full matrix `λ_ij` (wMB §1). -/
theorem Qplus_lambda_form (tp : ℝ) (sp : MSpec B α) (dos : Fin B → ℝ) (hdos : ∀ i, dos i ≠ 0)
    (μ : Fin B → Fin B → ℝ) (Nc : WithTop ℕ) (T : ℝ) (hT : 0 < T) (N : ℕ)
    (φ : Fin B → ℕ → ℝ) :
    Qplus tp sp dos μ Nc T N φ
      = -∑ i, dos i * Cw N (φ i)
        + ∑ i, dos i * ∑ j,
            (∑ n ∈ range N, ∑ m ∈ range N,
                (lamMB tp sp dos T i j (adist n m) + lamMB tp sp dos T i j (n + m + 1))
                  * φ i n * φ j m
              - ∑ n ∈ range N,
                (lamMB tp sp dos T i j 0 + 2 * Ssum (lamMB tp sp dos T i j) n) * φ i n ^ 2)
        - coul dos μ Nc N φ := by
  unfold Qplus
  congr 2
  have h : (1 / T) * ∑ s ∈ sp.S, ∑ i, ∑ j, sp.W s i j *
      (Bp tp (sp.Ω s / T) N (φ i) (φ j) - Dq tp (sp.Ω s / T) N (φ i) (φ i))
      = ∑ i, ∑ j, ∑ s ∈ sp.S, (1 / T) * (sp.W s i j *
      (Bp tp (sp.Ω s / T) N (φ i) (φ j) - Dq tp (sp.Ω s / T) N (φ i) (φ i))) := by
    simp only [mul_sum]
    rw [sum_comm]
    apply sum_congr rfl; intro i _
    rw [sum_comm]
  rw [h]
  apply sum_congr rfl; intro i _
  rw [mul_sum]
  apply sum_congr rfl; intro j _
  have e1 : ∀ s ∈ sp.S, 1 / T * (sp.W s i j *
      (Bp tp (sp.Ω s / T) N (φ i) (φ j) - Dq tp (sp.Ω s / T) N (φ i) (φ i)))
      = Fbil (fun k => (1 / T) * sp.W s i j * lw tp (sp.Ω s / T) k) N (φ i) (φ j) := by
    intro s _
    rw [Bp_Dq_eq_Fbil, Fbil_smul]; ring
  rw [sum_congr rfl e1, ← Fbil_sum]
  have e2 : (fun k => ∑ s ∈ sp.S, (1 / T) * sp.W s i j * lw tp (sp.Ω s / T) k)
      = fun k => dos i * lamMB tp sp dos T i j k := by
    funext k; rw [dos_lam tp sp dos hdos T hT i j k]
  rw [e2, Fbil_smul]
  rfl

/-- `|x|² = ∑_i N_i |φ_i|²` with `x_i = r_i φ_i`, `r_i² = N_i` (i.e. `r_i = √N_i`). -/
theorem nrm2_eq_x (dos r : Fin B → ℝ) (hr : ∀ i, r i ^ 2 = dos i) (N : ℕ)
    (φ : Fin B → ℕ → ℝ) :
    nrm2 dos N φ = ∑ i, nsq N (fun n => r i * φ i n) := by
  unfold nrm2 nsq
  apply sum_congr rfl; intro i _
  rw [mul_sum]; apply sum_congr rfl; intro n _
  rw [← hr i]; ring

end

end MultiBand
