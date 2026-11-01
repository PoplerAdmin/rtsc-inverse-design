/-
MultiBand.Identities: Lemma 1.1 (Z-split, pair-breaking term `P ≥ 0`) and
Lemma 4.1 (Perron-gauge identity) of `wMB/PROOF.md`.
-/
import MultiBand.Algebra

open Finset

namespace MultiBand

open TcR1

noncomputable section

variable {B : ℕ} {α : Type*}

/-- Swapping the band indices in a double sum against a symmetric weight. -/
theorem swap_sum (W : Fin B → Fin B → ℝ) (hW : ∀ i j, W i j = W j i)
    (f : Fin B → Fin B → ℝ) :
    ∑ i, ∑ j, W i j * f j i = ∑ i, ∑ j, W i j * f i j := by
  rw [sum_comm]
  apply sum_congr rfl; intro i _
  apply sum_congr rfl; intro j _
  rw [hW j i]

theorem Dq_symm (tp w : ℝ) (N : ℕ) (y z : ℕ → ℝ) : Dq tp w N y z = Dq tp w N z y := by
  unfold Dq; apply sum_congr rfl; intro n _; ring

theorem Dq_sq (tp w : ℝ) (N : ℕ) (a b : ℕ → ℝ) :
    ∑ n ∈ range N, Dn tp w n * (a n - b n) ^ 2
      = Dq tp w N a a - 2 * Dq tp w N a b + Dq tp w N b b := by
  unfold Dq
  rw [mul_sum, ← sum_sub_distrib, ← sum_add_distrib]
  apply sum_congr rfl; intro n _; ring

/-- Per-atom Z-split: for a symmetric weight matrix `W`,
`∑_ij W_ij [φ_iᵀBφ_j - φ_iᵀDφ_i] = ∑_ij W_ij [φ_iᵀBφ_j - φ_iᵀDφ_j]
   - ½ ∑_ij W_ij ∑_n D_n (φ_{n,i} - φ_{n,j})²`. -/
theorem Z_split (W : Fin B → Fin B → ℝ) (hW : ∀ i j, W i j = W j i) (tp w : ℝ) (N : ℕ)
    (Bf : Fin B → Fin B → ℝ) (φ : Fin B → ℕ → ℝ) :
    ∑ i, ∑ j, W i j * (Bf i j - Dq tp w N (φ i) (φ i))
      = ∑ i, ∑ j, W i j * (Bf i j - Dq tp w N (φ i) (φ j))
        - (1 / 2) * ∑ i, ∑ j, W i j * ∑ n ∈ range N, Dn tp w n * (φ i n - φ j n) ^ 2 := by
  have hswap := swap_sum W hW (fun i j => Dq tp w N (φ i) (φ i))
  have e1 : ∑ i, ∑ j, W i j * ∑ n ∈ range N, Dn tp w n * (φ i n - φ j n) ^ 2
      = ∑ i, ∑ j, W i j * Dq tp w N (φ i) (φ i)
        - 2 * ∑ i, ∑ j, W i j * Dq tp w N (φ i) (φ j)
        + ∑ i, ∑ j, W i j * Dq tp w N (φ j) (φ j) := by
    simp only [Dq_sq, mul_add, mul_sub, sum_add_distrib, sum_sub_distrib, mul_sum]
    congr 1; congr 1
    apply sum_congr rfl; intro i _; apply sum_congr rfl; intro j _; ring
  rw [e1, hswap]
  simp only [mul_sub, sum_sub_distrib]
  ring

/-- **Lemma 1.1 (even).** `Q⁺(φ) = -∑_i N_i C(φ_i) + (1/T) ∑_ij ∫dW_ij G⁺_w(φ_i, φ_j) - P(φ)
- Coulomb` (any `μ*`; with `μ* = 0` this is the display of wMB Lemma 1.1). -/
theorem lemma11_even (tp : ℝ) (sp : MSpec B α) (dos : Fin B → ℝ) (μ : Fin B → Fin B → ℝ)
    (Nc : WithTop ℕ) (T : ℝ) (N : ℕ) (φ : Fin B → ℕ → ℝ) :
    Qplus tp sp dos μ Nc T N φ
      = -∑ i, dos i * Cw N (φ i)
        + (1 / T) * ∑ s ∈ sp.S, ∑ i, ∑ j, sp.W s i j * Gp tp (sp.Ω s / T) N (φ i) (φ j)
        - Pterm tp sp T N φ - coul dos μ Nc N φ := by
  unfold Qplus Pterm Gp
  have h : ∀ s ∈ sp.S, ∑ i, ∑ j, sp.W s i j *
      (Bp tp (sp.Ω s / T) N (φ i) (φ j) - Dq tp (sp.Ω s / T) N (φ i) (φ i))
      = ∑ i, ∑ j, sp.W s i j *
          (Bp tp (sp.Ω s / T) N (φ i) (φ j) - Dq tp (sp.Ω s / T) N (φ i) (φ j))
        - (1 / 2) * ∑ i, ∑ j, sp.W s i j *
          ∑ n ∈ range N, Dn tp (sp.Ω s / T) n * (φ i n - φ j n) ^ 2 :=
    fun s hs => Z_split (sp.W s) (sp.W_symm s hs) tp _ N _ φ
  rw [sum_congr rfl h, sum_sub_distrib, ← mul_sum]
  ring

/-- **Lemma 1.1 (odd).** Same with `B⁻`, `G⁻`, no Coulomb term. -/
theorem lemma11_odd (tp : ℝ) (sp : MSpec B α) (dos : Fin B → ℝ) (T : ℝ) (N : ℕ)
    (φ : Fin B → ℕ → ℝ) :
    Qminus tp sp dos T N φ
      = -∑ i, dos i * Cw N (φ i)
        + (1 / T) * ∑ s ∈ sp.S, ∑ i, ∑ j, sp.W s i j * Gm tp (sp.Ω s / T) N (φ i) (φ j)
        - Pterm tp sp T N φ := by
  unfold Qminus Pterm Gm
  have h : ∀ s ∈ sp.S, ∑ i, ∑ j, sp.W s i j *
      (Bm tp (sp.Ω s / T) N (φ i) (φ j) - Dq tp (sp.Ω s / T) N (φ i) (φ i))
      = ∑ i, ∑ j, sp.W s i j *
          (Bm tp (sp.Ω s / T) N (φ i) (φ j) - Dq tp (sp.Ω s / T) N (φ i) (φ j))
        - (1 / 2) * ∑ i, ∑ j, sp.W s i j *
          ∑ n ∈ range N, Dn tp (sp.Ω s / T) n * (φ i n - φ j n) ^ 2 :=
    fun s hs => Z_split (sp.W s) (sp.W_symm s hs) tp _ N _ φ
  rw [sum_congr rfl h, sum_sub_distrib, ← mul_sum]
  ring

/-- **Lemma 1.1, `P ≥ 0`** (for `T > 0`). -/
theorem Pterm_nonneg (tp : ℝ) (sp : MSpec B α) (T : ℝ) (hT : 0 < T) (N : ℕ)
    (φ : Fin B → ℕ → ℝ) : 0 ≤ Pterm tp sp T N φ := by
  unfold Pterm
  apply mul_nonneg (by positivity)
  apply sum_nonneg; intro s hs
  apply sum_nonneg; intro i _
  apply sum_nonneg; intro j _
  apply mul_nonneg (sp.W_nonneg s hs i j)
  apply sum_nonneg; intro n _
  exact mul_nonneg (Dn_nonneg tp _ (div_pos (sp.Ω_pos s hs) hT).le n) (sq_nonneg _)

/-- The Coulomb term is `≥ 0` when `U_ij = N_i μ*_ij` is PSD. -/
theorem coul_nonneg (dos : Fin B → ℝ) (μ : Fin B → Fin B → ℝ) (hU : CoulPSD dos μ)
    (Nc : WithTop ℕ) (N : ℕ) (φ : Fin B → ℕ → ℝ) : 0 ≤ coul dos μ Nc N φ := by
  unfold coul
  have := hU (fun i => Sc Nc N (φ i))
  linarith

/-! ## Perron gauge -/

theorem mw_nonneg (sp : MSpec B α) (v : Fin B → ℝ) (hv : ∀ i, 0 < v i) (s : α) (hs : s ∈ sp.S)
    (i j : Fin B) : 0 ≤ mw sp v s i j := by
  unfold mw
  exact mul_nonneg (mul_nonneg (hv i).le (hv j).le) (sp.W_nonneg s hs i j)

/-- Steps 1–2 of Lemma 4.1: `W_ij G(φ_i,φ_j) = m_ij [G(ψ_ij) - G(χ_ij)]`. -/
theorem gauge_G (tp w : ℝ) (sp : MSpec B α) (v : Fin B → ℝ) (hv : ∀ i, v i ≠ 0) (N : ℕ)
    (φ : Fin B → ℕ → ℝ) (s : α) (i j : Fin B) :
    sp.W s i j * Gp tp w N (φ i) (φ j)
      = mw sp v s i j * (Gp tp w N (psi v φ i j) (psi v φ i j)
          - Gp tp w N (chi v φ i j) (chi v φ i j)) := by
  have hp := Gp_polar tp w N (theta v φ i) (theta v φ j)
  unfold psi chi
  rw [hp, Gp_gauge tp w N (v i) (v j) (hv i) (hv j) (φ i) (φ j)]
  unfold mw theta; ring

/-- Step 3 of Lemma 4.1 (left part):
`∑_ij m_ij [C(ψ_ij) + C(χ_ij)] = ∑_i C(θ_i) v_i ∑_j W^tot_ij v_j`. -/
theorem gauge_C_sum (sp : MSpec B α) (v : Fin B → ℝ) (N : ℕ) (φ : Fin B → ℕ → ℝ) :
    ∑ s ∈ sp.S, ∑ i, ∑ j, mw sp v s i j * (Cw N (psi v φ i j) + Cw N (chi v φ i j))
      = ∑ i, Cw N (theta v φ i) * v i * ∑ j, (∑ s ∈ sp.S, sp.W s i j) * v j := by
  have h1 : ∀ s ∈ sp.S, ∑ i, ∑ j, mw sp v s i j * (Cw N (psi v φ i j) + Cw N (chi v φ i j))
      = ∑ i, ∑ j, mw sp v s i j * Cw N (theta v φ i) := by
    intro s hs
    have hpar : ∀ i j, Cw N (psi v φ i j) + Cw N (chi v φ i j)
        = (Cw N (theta v φ i) + Cw N (theta v φ j)) / 2 := by
      intro i j; unfold psi chi; exact Cw_parallelogram N _ _
    simp only [hpar]
    have hsw := swap_sum (mw sp v s) (fun i j => by unfold mw; rw [sp.W_symm s hs i j]; ring)
      (fun i _ => Cw N (theta v φ i))
    have e : ∑ i, ∑ j, mw sp v s i j * ((Cw N (theta v φ i) + Cw N (theta v φ j)) / 2)
        = (1 / 2) * ∑ i, ∑ j, mw sp v s i j * Cw N (theta v φ i)
          + (1 / 2) * ∑ i, ∑ j, mw sp v s i j * Cw N (theta v φ j) := by
      simp only [mul_sum, ← sum_add_distrib]
      apply sum_congr rfl; intro i _; apply sum_congr rfl; intro j _; ring
    rw [e, hsw]; ring
  rw [sum_congr rfl h1, sum_comm]
  apply sum_congr rfl; intro i _
  rw [sum_comm, mul_sum]
  apply sum_congr rfl; intro j _
  rw [sum_mul, mul_sum]
  apply sum_congr rfl; intro s _
  unfold mw; ring

/-- `∑_j W^tot_ij v_j = N_i (𝒜 v)_i`. -/
theorem row_eq (sp : MSpec B α) (dos : Fin B → ℝ) (hdos : ∀ i, dos i ≠ 0) (v : Fin B → ℝ)
    (i : Fin B) :
    ∑ j, (∑ s ∈ sp.S, sp.W s i j) * v j = dos i * ∑ j, Amat sp dos i j * v j := by
  rw [mul_sum]; apply sum_congr rfl; intro j _
  have := hdos i
  unfold Amat; field_simp

/-- `C(θ_i) v_i² = C(φ_i)`. -/
theorem Cw_theta (v : Fin B → ℝ) (hv : ∀ i, v i ≠ 0) (N : ℕ) (φ : Fin B → ℕ → ℝ) (i : Fin B) :
    Cw N (theta v φ i) * v i ^ 2 = Cw N (φ i) := by
  have hvi := hv i
  have : φ i = fun n => v i * theta v φ i n := by
    funext n; unfold theta; field_simp
  rw [this, Cw_smul]
  ring

/-- Step 3 of Lemma 4.1 (exact, eigenvector case):
`∑_ij m_ij [C(ψ_ij) + C(χ_ij)] = κ ∑_i N_i C(φ_i)`. -/
theorem gauge_C_eig (sp : MSpec B α) (dos : Fin B → ℝ) (hdos : ∀ i, dos i ≠ 0)
    (v : Fin B → ℝ) (hv : ∀ i, v i ≠ 0) (κ : ℝ)
    (heig : ∀ i, ∑ j, Amat sp dos i j * v j = κ * v i) (N : ℕ) (φ : Fin B → ℕ → ℝ) :
    ∑ s ∈ sp.S, ∑ i, ∑ j, mw sp v s i j * (Cw N (psi v φ i j) + Cw N (chi v φ i j))
      = κ * ∑ i, dos i * Cw N (φ i) := by
  rw [gauge_C_sum, mul_sum]
  apply sum_congr rfl; intro i _
  rw [row_eq sp dos hdos v i, heig i, ← Cw_theta v hv N φ i]
  ring

/-- Step 3, Collatz–Wielandt version: if `𝒜 v ≤ κ v` (`v > 0`, `N_i > 0`) then
`∑_ij m_ij [C(ψ_ij) + C(χ_ij)] ≤ κ ∑_i N_i C(φ_i)`. -/
theorem gauge_C_sub (sp : MSpec B α) (dos : Fin B → ℝ) (hdos : ∀ i, 0 < dos i)
    (v : Fin B → ℝ) (hv : ∀ i, 0 < v i) (κ : ℝ)
    (hsub : ∀ i, ∑ j, Amat sp dos i j * v j ≤ κ * v i) (N : ℕ) (φ : Fin B → ℕ → ℝ) :
    ∑ s ∈ sp.S, ∑ i, ∑ j, mw sp v s i j * (Cw N (psi v φ i j) + Cw N (chi v φ i j))
      ≤ κ * ∑ i, dos i * Cw N (φ i) := by
  rw [gauge_C_sum, mul_sum]
  apply sum_le_sum; intro i _
  rw [row_eq sp dos (fun i => (hdos i).ne') v i,
    ← Cw_theta v (fun i => (hv i).ne') N φ i]
  have h0 : 0 ≤ Cw N (theta v φ i) * v i * dos i :=
    mul_nonneg (mul_nonneg (Cw_nonneg N _) (hv i).le) (hdos i).le
  have := mul_le_mul_of_nonneg_left (hsub i) h0
  nlinarith

/-- **Lemma 4.1 (Perron-gauge identity).** For `𝒜 v = κ v` (`κ ≠ 0`, `v_i ≠ 0`, `N_i ≠ 0`,
`T ≠ 0`), with `a = κ/T`, `w = Ω/T`:
`Q⁺(φ) = (1/κ) ∑_ij ∫ dm_ij [Q^E_{+a,w}(ψ_ij) + Q^E_{-a,w}(χ_ij)] - P(φ) - Coulomb`.
(With `μ* = 0` the Coulomb term vanishes and this is wMB's display.) -/
theorem lemma41 (tp : ℝ) (sp : MSpec B α) (dos : Fin B → ℝ) (μ : Fin B → Fin B → ℝ)
    (Nc : WithTop ℕ) (T : ℝ) (N : ℕ) (φ : Fin B → ℕ → ℝ) (v : Fin B → ℝ) (κ : ℝ)
    (hκ : κ ≠ 0) (hT : T ≠ 0) (hv : ∀ i, v i ≠ 0) (hdos : ∀ i, dos i ≠ 0)
    (heig : ∀ i, ∑ j, Amat sp dos i j * v j = κ * v i) :
    Qplus tp sp dos μ Nc T N φ
      = (1 / κ) * ∑ s ∈ sp.S, ∑ i, ∑ j, mw sp v s i j *
          (Qe tp (κ / T) (sp.Ω s / T) N (psi v φ i j)
            + Qe tp (-(κ / T)) (sp.Ω s / T) N (chi v φ i j))
        - Pterm tp sp T N φ - coul dos μ Nc N φ := by
  rw [lemma11_even]
  have hQ : ∑ s ∈ sp.S, ∑ i, ∑ j, mw sp v s i j *
          (Qe tp (κ / T) (sp.Ω s / T) N (psi v φ i j)
            + Qe tp (-(κ / T)) (sp.Ω s / T) N (chi v φ i j))
      = -∑ s ∈ sp.S, ∑ i, ∑ j, mw sp v s i j * (Cw N (psi v φ i j) + Cw N (chi v φ i j))
        + (κ / T) * ∑ s ∈ sp.S, ∑ i, ∑ j, mw sp v s i j *
            (Gp tp (sp.Ω s / T) N (psi v φ i j) (psi v φ i j)
              - Gp tp (sp.Ω s / T) N (chi v φ i j) (chi v φ i j)) := by
    simp only [mul_sum, ← sum_neg_distrib, ← sum_add_distrib]
    apply sum_congr rfl; intro s _; apply sum_congr rfl; intro i _
    apply sum_congr rfl; intro j _
    rw [Qe_eq, Qe_eq]; ring
  have hG : ∑ s ∈ sp.S, ∑ i, ∑ j, sp.W s i j * Gp tp (sp.Ω s / T) N (φ i) (φ j)
      = ∑ s ∈ sp.S, ∑ i, ∑ j, mw sp v s i j *
            (Gp tp (sp.Ω s / T) N (psi v φ i j) (psi v φ i j)
              - Gp tp (sp.Ω s / T) N (chi v φ i j) (chi v φ i j)) := by
    apply sum_congr rfl; intro s _; apply sum_congr rfl; intro i _
    apply sum_congr rfl; intro j _
    exact gauge_G tp _ sp v hv N φ s i j
  rw [hQ, hG, gauge_C_eig sp dos hdos v hv κ heig N φ]
  field_simp

end

end MultiBand
