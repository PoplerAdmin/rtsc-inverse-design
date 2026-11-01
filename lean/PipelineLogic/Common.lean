/-
PipelineLogic.Common: shared notions for the pipeline decision logic.

* `Stable κ sp μ Nc T`: both sectors have `ρ < 0` (negation of TcR1's `Unstable`).
* `InEnv e sp sp'`: `sp'` lies in the relative weight envelope of the nominal spectrum `sp`
  (same atoms, same frequencies, masses in `[(1-e) p_j, (1+e) p_j]`).
* `Mq L`: the infinite symmetric matrix of `C + Fp L` (the μ* = 0 even form).
* `Mhat hi lo`: its interval upper envelope (ARCH §4, R-CELLS).
-/
import TcR1.Basic

open Finset

namespace PipelineLogic

open TcR1

noncomputable section

/-- `Stable(T)`: `ρ⁺(T) < 0` and `ρ⁻(T) < 0`, i.e. neither sector is `TcR1.Unstable`. -/
def Stable {ι : Type*} (κ : ℝ) (sp : Spectrum ι) (μ : ℝ) (Nc : WithTop ℕ) (T : ℝ) : Prop :=
  ¬ Unstable (Qp κ sp μ Nc T) ∧ ¬ Unstable (Qm κ sp T)

/-- `sp'` is in the relative envelope `e` of `sp`: same atoms and frequencies, and
`(1-e) p_j ≤ p'_j ≤ (1+e) p_j` for every atom. -/
def InEnv {ι : Type*} (e : ℝ) (sp sp' : Spectrum ι) : Prop :=
  sp'.S = sp.S ∧ (∀ j ∈ sp.S, sp'.Ω j = sp.Ω j) ∧
    ∀ j ∈ sp.S, (1 - e) * sp.p j ≤ sp'.p j ∧ sp'.p j ≤ (1 + e) * sp.p j

/-- Matrix of the μ* = 0 even form `C(x) + Fp L x`:
`M_nm = L(|n-m|) + L(n+m+1)` (n ≠ m), `M_nn = L(2n+1) - 2∑_{k=1}^n L(k) - (2n+1)`. -/
def Mq (L : ℕ → ℝ) (n m : ℕ) : ℝ :=
  Kp L n m - if n = m then (2 * (n : ℝ) + 1) else 0

/-- Interval upper envelope of `Mq` (ARCH R-CELLS): off-diagonals `hi(|n-m|) + hi(n+m+1)`,
diagonal `hi(2n+1) - 2∑_{k=1}^n lo(k) - (2n+1)`. -/
def Mhat (hi lo : ℕ → ℝ) (n m : ℕ) : ℝ :=
  if n = m then hi (2 * n + 1) - 2 * Ssum lo n - (2 * (n : ℝ) + 1)
  else hi (adist n m) + hi (n + m + 1)

/-- Truncated weighted row sum `∑_{m<N} M_nm v_m`. -/
def rowSum (M : ℕ → ℕ → ℝ) (v : ℕ → ℝ) (N n : ℕ) : ℝ := ∑ m ∈ range N, M n m * v m

theorem nsq_pos_of_ne {N : ℕ} {x : ℕ → ℝ} (hx : ∃ n ∈ range N, x n ≠ 0) : 0 < nsq N x := by
  obtain ⟨n, hn, hxn⟩ := hx
  unfold nsq
  have h1 : 0 < x n ^ 2 := by positivity
  have h2 : x n ^ 2 ≤ ∑ m ∈ range N, x m ^ 2 :=
    single_le_sum (f := fun m => x m ^ 2) (fun m _ => sq_nonneg (x m)) hn
  linarith

/-- The even form at μ* = 0 is the quadratic form of `Mq`. -/
theorem form_eq_Mq (L : ℕ → ℝ) (N : ℕ) (x : ℕ → ℝ) :
    C N x + Fp L N x = ∑ n ∈ range N, ∑ m ∈ range N, Mq L n m * x n * x m := by
  unfold C Fp Mq
  have h : ∀ n ∈ range N, ∑ m ∈ range N, (Kp L n m - if n = m then (2 * (n : ℝ) + 1) else 0)
      * x n * x m = ∑ m ∈ range N, Kp L n m * x n * x m - (2 * (n : ℝ) + 1) * x n ^ 2 := by
    intro n hn
    have : ∀ m ∈ range N, (Kp L n m - if n = m then (2 * (n : ℝ) + 1) else 0) * x n * x m
        = Kp L n m * x n * x m - (if n = m then (2 * (n : ℝ) + 1) * x n * x m else 0) := by
      intro m _; split_ifs <;> ring
    rw [sum_congr rfl this, sum_sub_distrib, sum_ite_eq (range N) n, ite_eq_left hn]
    ring
  rw [sum_congr rfl h, sum_sub_distrib]
  ring

theorem adist_comm (n m : ℕ) : adist n m = adist m n := by unfold adist; omega

theorem Mq_symm (L : ℕ → ℝ) (n m : ℕ) : Mq L n m = Mq L m n := by
  unfold Mq Kp
  by_cases h : n = m
  · subst h; rfl
  · have h' : m ≠ n := fun h'' => h h''.symm
    simp only [h, h', ite_false, adist_comm n m, Nat.add_comm n m]

theorem Mq_offdiag_nonneg (L : ℕ → ℝ) (hL : ∀ k, 1 ≤ k → 0 ≤ L k) (n m : ℕ) (h : n ≠ m) :
    0 ≤ Mq L n m := by
  unfold Mq Kp
  simp only [h, ite_false, sub_zero]
  have h1 : 1 ≤ adist n m := by unfold adist; omega
  exact add_nonneg (hL _ h1) (hL _ (by omega))

/-- Entrywise: the interval envelope dominates the true matrix. -/
theorem Mq_le_Mhat (L lo hi : ℕ → ℝ) (hlo : ∀ k, lo k ≤ L k) (hhi : ∀ k, L k ≤ hi k)
    (n m : ℕ) : Mq L n m ≤ Mhat hi lo n m := by
  unfold Mq Mhat Kp
  by_cases h : n = m
  · simp only [h, ite_true]
    have hs : Ssum lo m ≤ Ssum L m := by
      unfold Ssum; exact sum_le_sum (fun k _ => hlo (k + 1))
    linarith [hhi (2 * m + 1)]
  · simp only [h, ite_false, sub_zero]
    linarith [hhi (adist n m), hhi (n + m + 1)]

end

end PipelineLogic
