/-
TargetLogic.Moments: moment inequalities L-MOM1..3 for a finite atomic NONNEGATIVE spectrum
(`TcR1.Spectrum`: weights `p_j ≥ 0` at frequencies `Ω_j > 0`). Documentation lemmas only
(review R12 finding 2: they must not be used as a rejection rule; with signed weights, e.g.
a non-admitted QE smearing, L-MOM1 is false).

Definitions (`01_eliashberg_targets.md` §0, finite atomic case `α²F = ∑_j p_j δ(ω - Ω_j)`):
* `A = ∑_j p_j` (`TcR1.Spectrum.A`),
* `λ = 2 ∑_j p_j / Ω_j` (`lamS`),
* `m = λ⟨ω²⟩ = 2 ∑_j p_j Ω_j` (`momS`),
* `ω₂ = √(m / λ)` (`omega2`),
* `ω_log = exp((2/λ) ∑_j (p_j/Ω_j) ln Ω_j)` (`omegaLog`).
ω_max is not a field of `TcR1.Spectrum`; it enters as the hypothesis `∀ j ∈ S, Ω_j ≤ ω_max`.
-/
import TcR1.Basic
import Mathlib.Algebra.BigOperators.Field
import Mathlib.Analysis.SpecialFunctions.Log.Basic
import Mathlib.Analysis.SpecialFunctions.Sqrt

open Finset

namespace TargetLogic

open TcR1

noncomputable section

variable {ι : Type*}

/-- `λ = 2 ∑_j p_j / Ω_j`. -/
def lamS (sp : Spectrum ι) : ℝ := 2 * ∑ j ∈ sp.S, sp.p j / sp.Ω j

/-- `m = λ⟨ω²⟩ = 2 ∑_j p_j Ω_j`. -/
def momS (sp : Spectrum ι) : ℝ := 2 * ∑ j ∈ sp.S, sp.p j * sp.Ω j

/-- `ω₂ = √⟨ω²⟩ = √(m/λ)`. -/
def omega2 (sp : Spectrum ι) : ℝ := Real.sqrt (momS sp / lamS sp)

/-- `ω_log = exp((2/λ) ∑_j (p_j/Ω_j) ln Ω_j)`. -/
def omegaLog (sp : Spectrum ι) : ℝ :=
  Real.exp ((2 / lamS sp) * ∑ j ∈ sp.S, sp.p j / sp.Ω j * Real.log (sp.Ω j))

theorem A_nonneg (sp : Spectrum ι) : 0 ≤ sp.A := sum_nonneg sp.p_nonneg

theorem lamS_nonneg (sp : Spectrum ι) : 0 ≤ lamS sp :=
  mul_nonneg (by norm_num)
    (sum_nonneg (fun j hj => div_nonneg (sp.p_nonneg j hj) (sp.Ω_pos j hj).le))

theorem momS_nonneg (sp : Spectrum ι) : 0 ≤ momS sp :=
  mul_nonneg (by norm_num)
    (sum_nonneg (fun j hj => mul_nonneg (sp.p_nonneg j hj) (sp.Ω_pos j hj).le))

/-- **L-MOM1.** `A ≤ λ ω_max / 2` when every frequency is `≤ ω_max`. -/
theorem lmom1 (sp : Spectrum ι) (ωmax : ℝ) (hmax : ∀ j ∈ sp.S, sp.Ω j ≤ ωmax) :
    sp.A ≤ lamS sp * ωmax / 2 := by
  have h : ∀ j ∈ sp.S, sp.p j ≤ sp.p j / sp.Ω j * ωmax := by
    intro j hj
    have hΩ := sp.Ω_pos j hj
    have hq : 0 ≤ sp.p j / sp.Ω j := div_nonneg (sp.p_nonneg j hj) hΩ.le
    calc sp.p j = sp.p j / sp.Ω j * sp.Ω j := by field_simp
      _ ≤ sp.p j / sp.Ω j * ωmax := mul_le_mul_of_nonneg_left (hmax j hj) hq
  have h2 : sp.A ≤ (∑ j ∈ sp.S, sp.p j / sp.Ω j) * ωmax := by
    rw [sum_mul]; exact sum_le_sum h
  unfold lamS; linarith

/-- **L-MOM2 (squared form).** `A² ≤ λ m / 4` (Cauchy–Schwarz). -/
theorem lmom2_sq (sp : Spectrum ι) : sp.A ^ 2 ≤ lamS sp * momS sp / 4 := by
  have h := sum_sq_le_sum_mul_sum_of_sq_le_mul sp.S (r := sp.p)
    (f := fun j => sp.p j / sp.Ω j) (g := fun j => sp.p j * sp.Ω j)
    (fun j hj => div_nonneg (sp.p_nonneg j hj) (sp.Ω_pos j hj).le)
    (fun j hj => mul_nonneg (sp.p_nonneg j hj) (sp.Ω_pos j hj).le)
    (fun j hj => le_of_eq (by have := (sp.Ω_pos j hj).ne'; field_simp))
  unfold Spectrum.A lamS momS
  linarith

/-- **L-MOM2.** `A ≤ λ ω₂ / 2` for `λ > 0`. -/
theorem lmom2 (sp : Spectrum ι) (hl : 0 < lamS sp) : sp.A ≤ lamS sp * omega2 sp / 2 := by
  have hq : 0 ≤ momS sp / lamS sp := div_nonneg (momS_nonneg sp) hl.le
  have hs := Real.sq_sqrt hq
  have hX : (lamS sp * omega2 sp / 2) ^ 2 = lamS sp * momS sp / 4 := by
    unfold omega2
    have : (lamS sp * Real.sqrt (momS sp / lamS sp) / 2) ^ 2
        = lamS sp ^ 2 * (Real.sqrt (momS sp / lamS sp)) ^ 2 / 4 := by ring
    rw [this, hs]; field_simp
  have hX0 : 0 ≤ lamS sp * omega2 sp / 2 :=
    div_nonneg (mul_nonneg hl.le (Real.sqrt_nonneg _)) (by norm_num)
  have hA := lmom2_sq sp
  have hA0 := A_nonneg sp
  nlinarith [hX, hA, hA0, hX0]

/-- `λ > 0` forces `A > 0` (some atom has positive weight). -/
theorem A_pos_of_lamS_pos (sp : Spectrum ι) (hl : 0 < lamS sp) : 0 < sp.A := by
  by_contra h
  have hA0 : sp.A = 0 := le_antisymm (not_lt.mp h) (A_nonneg sp)
  have hz : ∀ j ∈ sp.S, sp.p j = 0 := (sum_eq_zero_iff_of_nonneg sp.p_nonneg).mp hA0
  have : lamS sp = 0 := by
    unfold lamS
    rw [sum_eq_zero (fun j hj => by rw [hz j hj, zero_div])]
    ring
  linarith

/-- **L-MOM3.** `A ≥ λ ω_log / 2` for `λ > 0` (weighted Jensen for `ln`, via `ln x ≤ x - 1`). -/
theorem lmom3 (sp : Spectrum ι) (hl : 0 < lamS sp) : lamS sp * omegaLog sp / 2 ≤ sp.A := by
  set s := ∑ j ∈ sp.S, sp.p j / sp.Ω j with hs_def
  have hls : lamS sp = 2 * s := rfl
  have hs : 0 < s := by linarith
  have hA := A_pos_of_lamS_pos sp hl
  set M := sp.A / s with hM_def
  have hM : 0 < M := div_pos hA hs
  -- termwise: (p/Ω) (ln Ω - ln M) ≤ (p/Ω)(Ω/M - 1) = p/M - p/Ω
  have hterm : ∀ j ∈ sp.S, sp.p j / sp.Ω j * (Real.log (sp.Ω j) - Real.log M)
      ≤ sp.p j / M - sp.p j / sp.Ω j := by
    intro j hj
    have hΩ := sp.Ω_pos j hj
    have hq : 0 ≤ sp.p j / sp.Ω j := div_nonneg (sp.p_nonneg j hj) hΩ.le
    have hlog : Real.log (sp.Ω j) - Real.log M ≤ sp.Ω j / M - 1 := by
      rw [← Real.log_div hΩ.ne' hM.ne']
      exact Real.log_le_sub_one_of_pos (div_pos hΩ hM)
    calc sp.p j / sp.Ω j * (Real.log (sp.Ω j) - Real.log M)
        ≤ sp.p j / sp.Ω j * (sp.Ω j / M - 1) := mul_le_mul_of_nonneg_left hlog hq
      _ = sp.p j / M - sp.p j / sp.Ω j := by field_simp
  have hsum := sum_le_sum hterm
  have hL : ∑ j ∈ sp.S, sp.p j / sp.Ω j * (Real.log (sp.Ω j) - Real.log M)
      = (∑ j ∈ sp.S, sp.p j / sp.Ω j * Real.log (sp.Ω j)) - s * Real.log M := by
    rw [hs_def, sum_mul, ← sum_sub_distrib]
    apply sum_congr rfl; intro j _; ring
  have hR : ∑ j ∈ sp.S, (sp.p j / M - sp.p j / sp.Ω j) = sp.A / M - s := by
    rw [sum_sub_distrib, ← sum_div]; rfl
  have hAM : sp.A / M = s := by rw [hM_def]; field_simp
  rw [hL, hR, hAM] at hsum
  -- exponent ≤ ln M
  have hexp : (2 / lamS sp) * ∑ j ∈ sp.S, sp.p j / sp.Ω j * Real.log (sp.Ω j) ≤ Real.log M := by
    rw [hls]
    have : (2 / (2 * s)) = 1 / s := by field_simp
    rw [this, one_div, inv_mul_le_iff₀ hs]
    linarith
  have hω : omegaLog sp ≤ M := by
    unfold omegaLog
    calc Real.exp _ ≤ Real.exp (Real.log M) := Real.exp_le_exp.mpr hexp
      _ = M := Real.exp_log hM
  calc lamS sp * omegaLog sp / 2 = s * omegaLog sp := by rw [hls]; ring
    _ ≤ s * M := mul_le_mul_of_nonneg_left hω hs.le
    _ = sp.A := by rw [hM_def]; field_simp

end

end TargetLogic
