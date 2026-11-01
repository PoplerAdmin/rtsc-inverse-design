/-
P3: monotonicity of `g_T(k, ω)` in `T`, the cell bounds on `λ(k; T)`, and the
relative-weight envelope.
-/
import PipelineLogic.Common

open Finset

namespace PipelineLogic

open TcR1

noncomputable section

/-- **P3.** `g_T(k, ω) = 2ω/(ω² + (κkT)²)` is antitone in `T ≥ 0` (for `ω > 0`; any `k`, any
`κ`): `0 ≤ T1 ≤ T2 ⇒ g_{T2}(k, ω) ≤ g_{T1}(k, ω)`. -/
theorem g_antitone (κ T1 T2 : ℝ) (k : ℕ) (w : ℝ) (hT1 : 0 ≤ T1) (h12 : T1 ≤ T2) (hw : 0 < w) :
    g κ T2 k w ≤ g κ T1 k w := by
  unfold g
  have hsq : T1 ^ 2 ≤ T2 ^ 2 := by nlinarith
  have hden : w ^ 2 + (κ * k * T1) ^ 2 ≤ w ^ 2 + (κ * k * T2) ^ 2 := by
    have e1 : (κ * k * T1) ^ 2 = (κ * k) ^ 2 * T1 ^ 2 := by ring
    have e2 : (κ * k * T2) ^ 2 = (κ * k) ^ 2 * T2 ^ 2 := by ring
    rw [e1, e2]
    nlinarith [sq_nonneg (κ * (k : ℝ))]
  apply div_le_div_of_nonneg_left (by positivity) (by positivity) hden

/-- `λ(k; T)` is antitone in `T ≥ 0`. -/
theorem lam_antitone {ι : Type*} (κ : ℝ) (sp : Spectrum ι) (T1 T2 : ℝ) (hT1 : 0 ≤ T1)
    (h12 : T1 ≤ T2) (k : ℕ) : lam κ sp T2 k ≤ lam κ sp T1 k := by
  unfold lam
  apply sum_le_sum; intro j hj
  exact mul_le_mul_of_nonneg_left (g_antitone κ T1 T2 k _ hT1 h12 (sp.Ω_pos j hj))
    (sp.p_nonneg j hj)

/-- **P3 (cell).** For `T ∈ [T1, T2]`, `T1 ≥ 0`: `λ(k; T2) ≤ λ(k; T) ≤ λ(k; T1)`. -/
theorem lam_mem_cell {ι : Type*} (κ : ℝ) (sp : Spectrum ι) (T1 T2 T : ℝ) (hT1 : 0 ≤ T1)
    (h1 : T1 ≤ T) (h2 : T ≤ T2) (k : ℕ) :
    lam κ sp T2 k ≤ lam κ sp T k ∧ lam κ sp T k ≤ lam κ sp T1 k :=
  ⟨lam_antitone κ sp T T2 (le_trans hT1 h1) h2 k, lam_antitone κ sp T1 T hT1 h1 k⟩

/-- `λ` of an envelope spectrum, written over the nominal atoms. -/
theorem lam_env_eq {ι : Type*} (κ : ℝ) (e : ℝ) (sp sp' : Spectrum ι) (h : InEnv e sp sp')
    (T : ℝ) (k : ℕ) : lam κ sp' T k = ∑ j ∈ sp.S, sp'.p j * g κ T k (sp.Ω j) := by
  obtain ⟨hS, hΩ, _⟩ := h
  unfold lam
  rw [hS]
  apply sum_congr rfl; intro j hj
  rw [hΩ j hj]

/-- **P3 (envelope scaling).** For `sp'` in the envelope `e` of `sp`:
`(1-e) λ_sp(k;T) ≤ λ_sp'(k;T) ≤ (1+e) λ_sp(k;T)`. -/
theorem lam_env {ι : Type*} (κ e : ℝ) (sp sp' : Spectrum ι) (h : InEnv e sp sp') (T : ℝ)
    (k : ℕ) :
    (1 - e) * lam κ sp T k ≤ lam κ sp' T k ∧ lam κ sp' T k ≤ (1 + e) * lam κ sp T k := by
  rw [lam_env_eq κ e sp sp' h]
  obtain ⟨_, _, hp⟩ := h
  unfold lam
  rw [mul_sum, mul_sum]
  constructor
  · apply sum_le_sum; intro j hj
    have hg := g_nonneg κ T k _ (sp.Ω_pos j hj)
    rw [← mul_assoc]
    exact mul_le_mul_of_nonneg_right (hp j hj).1 hg
  · apply sum_le_sum; intro j hj
    have hg := g_nonneg κ T k _ (sp.Ω_pos j hj)
    rw [← mul_assoc]
    exact mul_le_mul_of_nonneg_right (hp j hj).2 hg

/-- **P3 (cell + envelope), the R-CELLS box.** For `sp'` in the envelope of `sp`,
`T ∈ [T1, T2]`, `T1 ≥ 0`:
`(1-e) λ_sp(k; T2) ≤ λ_sp'(k; T) ≤ (1+e) λ_sp(k; T1)`. -/
theorem lam_box {ι : Type*} (κ e : ℝ) (sp sp' : Spectrum ι) (h : InEnv e sp sp')
    (T1 T2 T : ℝ) (hT1 : 0 ≤ T1) (h1 : T1 ≤ T) (h2 : T ≤ T2) (k : ℕ) :
    (1 - e) * lam κ sp T2 k ≤ lam κ sp' T k ∧ lam κ sp' T k ≤ (1 + e) * lam κ sp T1 k := by
  -- work atom by atom: p'_j g_T ∈ [(1-e) p_j g_{T2}, (1+e) p_j g_{T1}]
  rw [lam_env_eq κ e sp sp' h]
  obtain ⟨hS, _, hp⟩ := h
  unfold lam
  rw [mul_sum, mul_sum]
  have hT : 0 ≤ T := le_trans hT1 h1
  have hpj' : ∀ j ∈ sp.S, 0 ≤ sp'.p j := fun j hj => sp'.p_nonneg j (hS ▸ hj)
  constructor
  · apply sum_le_sum; intro j hj
    have hw := sp.Ω_pos j hj
    rw [← mul_assoc]
    calc (1 - e) * sp.p j * g κ T2 k (sp.Ω j) ≤ sp'.p j * g κ T2 k (sp.Ω j) :=
          mul_le_mul_of_nonneg_right (hp j hj).1 (g_nonneg κ T2 k _ hw)
      _ ≤ sp'.p j * g κ T k (sp.Ω j) :=
          mul_le_mul_of_nonneg_left (g_antitone κ T T2 k _ hT h2 hw) (hpj' j hj)
  · apply sum_le_sum; intro j hj
    have hw := sp.Ω_pos j hj
    rw [← mul_assoc]
    calc sp'.p j * g κ T k (sp.Ω j) ≤ sp'.p j * g κ T1 k (sp.Ω j) :=
          mul_le_mul_of_nonneg_left (g_antitone κ T1 T k _ hT1 h1 hw) (hpj' j hj)
      _ ≤ (1 + e) * sp.p j * g κ T1 k (sp.Ω j) :=
          mul_le_mul_of_nonneg_right (hp j hj).2 (g_nonneg κ T1 k _ hw)

/-- Total weight in the envelope: `A(sp') ≤ (1+e) A(sp)`. -/
theorem A_env {ι : Type*} (e : ℝ) (sp sp' : Spectrum ι) (h : InEnv e sp sp') :
    sp'.A ≤ (1 + e) * sp.A := by
  obtain ⟨hS, _, hp⟩ := h
  unfold Spectrum.A
  rw [hS, mul_sum]
  exact sum_le_sum (fun j hj => (hp j hj).2)

/-- The nominal spectrum is in its own envelope when `e ≥ 0`. -/
theorem InEnv_refl {ι : Type*} (e : ℝ) (he : 0 ≤ e) (sp : Spectrum ι) : InEnv e sp sp := by
  refine ⟨rfl, fun _ _ => rfl, fun j hj => ⟨?_, ?_⟩⟩
  · have := sp.p_nonneg j hj; nlinarith
  · have := sp.p_nonneg j hj; nlinarith

end

end PipelineLogic
