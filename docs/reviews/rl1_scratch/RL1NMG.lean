/-
RL1 scratch 2 (review of TargetLogic.NMG): `NMHyp` is not vacuous, and a weak NM rule is a
THEOREM given `H`.

* `stable_of_zero`  : a spectrum with all weights 0 is stable (needed because an antitone `g`
                      with `g > 0` forces `NMHyp` to say something about the zero spectrum).
* `nmhyp_of_H`      : from `H κ a0 η` (a0, η > 0) and `T* > 0`, `NMHyp κ gT T*` holds for
                      `gT(λ) = 4 (a0 T*)² / max(λ, ε)`, via TcR1.stable_of_ge + L-MOM2.
* `gT_antitone`     : `gT` is antitone on `[0, ∞)`.
* `nmg_rule_inst`   : all hypotheses of `nmg_rule_stable` instantiated (g := gT, up := id,
                      λ_cap = 100, g∞ = 0) with a concrete spectrum on which the rule fires.
-/
import TargetLogic

open TargetLogic TcR1 Finset

namespace RL1N

noncomputable section

theorem C_le_neg_nsq (N : ℕ) (x : ℕ → ℝ) : C N x ≤ -1 * nsq N x := by
  unfold C nsq
  have : ∑ n ∈ range N, x n ^ 2 ≤ ∑ n ∈ range N, (2 * (n : ℝ) + 1) * x n ^ 2 := by
    apply sum_le_sum
    intro n _
    have h1 : (1 : ℝ) ≤ 2 * (n : ℝ) + 1 := by have := Nat.cast_nonneg (α := ℝ) n; linarith
    nlinarith [sq_nonneg (x n)]
  linarith

theorem Kp_zero (n m : ℕ) : Kp (fun _ => (0 : ℝ)) n m = 0 := by
  unfold Kp Ssum; split_ifs <;> simp

theorem Km_zero (n m : ℕ) : Km (fun _ => (0 : ℝ)) n m = 0 := by
  unfold Km Ssum; split_ifs <;> simp

/-- The zero spectrum is stable at every `T`, every `μ* ≥ 0`, every cutoff. -/
theorem stable_of_zero {ι : Type*} (κ : ℝ) (sp : Spectrum ι) (h0 : ∀ j ∈ sp.S, sp.p j = 0)
    (μ : ℝ) (hμ : 0 ≤ μ) (Nc : WithTop ℕ) (T : ℝ) : PipelineLogic.Stable κ sp μ Nc T := by
  have hlam : lam κ sp T = fun _ => 0 := by
    funext k; unfold lam; exact sum_eq_zero (fun j hj => by simp [h0 j hj])
  constructor
  · apply not_unstable_of_bound _ 1 one_pos
    intro N x
    have hF : Fp (lam κ sp T) N x = 0 := by
      rw [hlam]; unfold Fp; simp [Kp_zero]
    have hM : 0 ≤ muTerm μ Nc N x := by unfold muTerm; positivity
    have := C_le_neg_nsq N x
    unfold Qp; rw [hF]; linarith
  · apply not_unstable_of_bound _ 1 one_pos
    intro N x
    have hF : Fm (lam κ sp T) N x = 0 := by
      rw [hlam]; unfold Fm; simp [Km_zero]
    have := C_le_neg_nsq N x
    unfold Qm; rw [hF]; linarith

/-- The provable moment floor `gT(λ) = 4 (a0 T*)² / max(λ, ε)`. -/
def gT (a0 Ts ε : ℝ) (l : ℝ) : ℝ := 4 * (a0 * Ts) ^ 2 / max l ε

theorem nmhyp_of_H (κ a0 η : ℝ) (hH : H κ a0 η) (ha0 : 0 < a0) (hη : 0 < η)
    (Ts : ℝ) (hTs : 0 < Ts) (ε : ℝ) (hε : 0 < ε) : NMHyp κ (gT a0 Ts ε) Ts := by
  intro ι sp μ Nc T hμ hT hns
  have hA : a0 * T < sp.A := by
    by_contra hle
    replace hle := not_lt.mp hle
    rcases (A_nonneg sp).lt_or_eq with hpos | hz
    · exact hns (stable_of_ge κ a0 η hH ha0 hη sp hpos μ hμ Nc T
        (by rw [div_le_iff₀ ha0]; linarith))
    · apply hns
      apply stable_of_zero κ sp _ μ hμ
      exact (sum_eq_zero_iff_of_nonneg sp.p_nonneg).mp hz.symm
  have hc : a0 * Ts ≤ a0 * T := mul_le_mul_of_nonneg_left hT ha0.le
  have hc0 : 0 < a0 * Ts := mul_pos ha0 hTs
  have hm := lmom2_sq sp
  have hl0 := lamS_nonneg sp
  have hm0 := momS_nonneg sp
  have hmax : 0 < max (lamS sp) ε := lt_max_of_lt_right hε
  unfold gT
  rw [div_le_iff₀ hmax]
  have h1 : 4 * (a0 * Ts) ^ 2 < 4 * sp.A ^ 2 := by nlinarith
  have h2 : lamS sp * momS sp ≤ max (lamS sp) ε * momS sp :=
    mul_le_mul_of_nonneg_right (le_max_left _ _) hm0
  nlinarith

theorem gT_antitone (a0 Ts ε : ℝ) (hε : 0 < ε) : AntitoneOn (gT a0 Ts ε) (Set.Ici 0) := by
  intro a _ b _ hab
  unfold gT
  apply div_le_div_of_nonneg_left (by positivity) (lt_max_of_lt_right hε)
  exact max_le_max hab le_rfl

/-- One atom, p = 1 at Ω = 1: λ = 2, m = 2. -/
def sp1 : Spectrum Unit where
  S := {()}
  p := fun _ => 1
  Ω := fun _ => 1
  p_nonneg := by intro j _; norm_num
  Ω_pos := by intro j _; norm_num

theorem lam_sp1 : lamS sp1 = 2 := by simp [lamS, sp1]
theorem mom_sp1 : momS sp1 = 2 := by simp [momS, sp1]

/-- Every hypothesis of `nmg_rule_stable` except `H` instantiated, and the rule fires
(m_hi = 2 < g_lookup(2) = 2 (a0 T*)²). -/
theorem nmg_rule_inst (κ η : ℝ) (hH : H κ 4.33148145 η) (hη : 0 < η) :
    ∀ μ : ℝ, 0 ≤ μ → ∀ Nc : WithTop ℕ, ∀ T : ℝ, 300 ≤ T →
      PipelineLogic.Stable κ sp1 μ Nc T :=
  nmg_rule_stable κ (gT 4.33148145 300 1) 300
    (nmhyp_of_H κ _ η hH (by norm_num) hη 300 (by norm_num) 1 one_pos)
    (gT_antitone _ _ 1 one_pos) id (fun _ => le_rfl) 100 0
    (fun l _ => by unfold gT; positivity) sp1 2 2 (by rw [mom_sp1]) (by rw [lam_sp1])
    (by unfold gLookup gT; norm_num)

end

end RL1N

#print axioms RL1N.stable_of_zero
#print axioms RL1N.nmhyp_of_H
#print axioms RL1N.gT_antitone
#print axioms RL1N.nmg_rule_inst
