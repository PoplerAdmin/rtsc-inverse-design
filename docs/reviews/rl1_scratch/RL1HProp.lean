/-
RL1 scratch 3: the two numerical hypotheses `TcR1.H tp a0 η` and `MultiBand.PropN tp a0` are
jointly satisfiable with `a0 > 0`, `η > 0` (so `decide_rejectedProven_sound`, `n1_dos_stable`,
`n1_rowmax_stable` are not vacuous; together with RL1Decide.ex_sound_inst every hypothesis of
`decide_rejectedProven_sound` is met by one concrete instance).

Witness: tp = 1, a0 = 1/64, η = 1/2. This says nothing about the certified value
a0 = 4.33148145 at tp = 2π (that is the interval-arithmetic certificate); it only rules out
that the hypothesis set is contradictory.

Route: |Fp(ge 1 v)(x)| ≤ 5 S |x|², with S = 2π a uniform (in v and N) bound on the prefix sums
∑_{1≤k<K} 2v/(v²+k²) (telescoping arctan, mean value theorem), and the row-sum lemmas of
MultiBand.Odd (`sum_shift`, `sum_adist_le`, `sum_adist_ge`).
-/
import TargetLogic
import Mathlib.Analysis.SpecialFunctions.Trigonometric.ArctanDeriv
import Mathlib.Analysis.Calculus.Deriv.MeanValue
import Mathlib.Analysis.Real.Pi.Bounds

open TcR1 MultiBand Finset

namespace RL1H

noncomputable section

/-- Per-term bound: `2v/(v²+k²) ≤ 2 arctan(k/v) − 2 arctan((k−1)/v)` for `k ≥ 1`. -/
theorem ge_le_arctan (v : ℝ) (hv : 0 < v) (k : ℕ) (hk : 1 ≤ k) :
    ge 1 v k ≤ 2 * Real.arctan ((k : ℝ) / v) - 2 * Real.arctan (((k : ℝ) - 1) / v) := by
  have hk' : (1 : ℝ) ≤ k := by exact_mod_cast hk
  obtain ⟨c, ⟨hc1, hc2⟩, hc⟩ := exists_hasDerivAt_eq_slope
    (fun t => Real.arctan (t / v)) (fun t => 1 / (1 + (t / v) ^ 2) * (1 / v))
    (show (k : ℝ) - 1 < k by linarith)
    ((Real.continuous_arctan.comp (continuous_id.div_const v)).continuousOn)
    (fun t _ => ((hasDerivAt_id t).div_const v).arctan)
  have hsub : (k : ℝ) - ((k : ℝ) - 1) = 1 := by ring
  rw [hsub, div_one] at hc
  have hc0 : 0 ≤ c := by linarith
  have hck : c ^ 2 ≤ (k : ℝ) ^ 2 := by nlinarith
  have key : v / (v ^ 2 + (k : ℝ) ^ 2) ≤ 1 / (1 + (c / v) ^ 2) * (1 / v) := by
    have e : 1 / (1 + (c / v) ^ 2) * (1 / v) = v / (v ^ 2 + c ^ 2) := by
      field_simp
    rw [e]
    apply div_le_div_of_nonneg_left hv.le (by positivity) (by linarith)
  have hge : ge 1 v k = 2 * (v / (v ^ 2 + (k : ℝ) ^ 2)) := by
    unfold ge; ring
  rw [hge]
  linarith

/-- `l' k = 0` at `k = 0`, `ge 1 v k` otherwise. -/
def lp (v : ℝ) (k : ℕ) : ℝ := if k = 0 then 0 else ge 1 v k

theorem lp_nonneg (v : ℝ) (hv : 0 < v) (k : ℕ) : 0 ≤ lp v k := by
  unfold lp; split_ifs
  · exact le_rfl
  · unfold ge; positivity

theorem pre_lp_le (v : ℝ) (hv : 0 < v) (K : ℕ) : pre (lp v) K ≤ 2 * Real.pi := by
  set F : ℕ → ℝ := fun j => 2 * Real.arctan (((j : ℝ) - 1) / v) with hF
  have hstep : ∀ k ∈ range K, lp v k ≤ F (k + 1) - F k := by
    intro k _
    simp only [hF]
    push_cast
    rcases Nat.eq_zero_or_pos k with h0 | hpos
    · subst h0
      simp only [lp, ite_true]
      have : Real.arctan ((0 - 1) / v) < 0 := by
        rw [Real.arctan_lt_zero]; apply div_neg_of_neg_of_pos (by norm_num) hv
      simp; linarith
    · have := ge_le_arctan v hv k hpos
      simp only [lp, Nat.pos_iff_ne_zero.mp hpos, ite_false]
      have e : (k : ℝ) + 1 - 1 = k := by ring
      rw [e]; linarith
  unfold pre
  calc ∑ k ∈ range K, lp v k ≤ ∑ k ∈ range K, (F (k + 1) - F k) := sum_le_sum hstep
    _ = F K - F 0 := sum_range_sub F K
    _ ≤ 2 * Real.pi := by
      simp only [hF]
      have h1 := Real.arctan_lt_pi_div_two (((K : ℝ) - 1) / v)
      have h2 := Real.neg_pi_div_two_lt_arctan ((((0 : ℕ) : ℝ) - 1) / v)
      linarith

/-- Row sums of `P_nm = l(|n−m|) + l(n+m+1)` for `l ≥ 0`, `l 0 = 0`, prefix sums `≤ S`. -/
theorem rowP_le (l : ℕ → ℝ) (hl : ∀ k, 0 ≤ l k) (hl0 : l 0 = 0) (S : ℝ)
    (hS : ∀ K, pre l K ≤ S) (n N : ℕ) :
    ∑ m ∈ range N, (l (adist n m) + l (n + m + 1)) ≤ 3 * S := by
  rw [sum_add_distrib, sum_shift]
  have hS0 : 0 ≤ S := (pre_nonneg l hl 0).trans (hS 0)
  have ha : ∑ m ∈ range N, l (adist n m) ≤ 2 * S := by
    rcases le_total N (n + 1) with h | h
    · rw [sum_adist_le l n N h]
      linarith [pre_nonneg l hl (n + 1 - N), hS (n + 1)]
    · rw [sum_adist_ge l n N h, hl0]
      linarith [hS (n + 1), hS (N - n)]
  linarith [pre_nonneg l hl (n + 1), hS (n + N + 1)]

/-- `Kp L = P − diag(2 Ssum L n)` with `P` built from `l' = L` off `0`, `l' 0 = 0`. -/
theorem Kp_split (L : ℕ → ℝ) (n m : ℕ) :
    Kp L n m = ((fun k => if k = 0 then 0 else L k) (adist n m) +
                (fun k => if k = 0 then 0 else L k) (n + m + 1)) -
               (if n = m then 2 * Ssum L n else 0) := by
  unfold Kp
  by_cases h : n = m
  · subst h
    have : adist n n = 0 := by unfold adist; omega
    simp [this]
    have : 2 * n + 1 = n + n + 1 := by ring
    rw [this]
  · have h1 : adist n m ≠ 0 := by unfold adist; omega
    simp [h, h1]

theorem Ssum_eq_lp (v : ℝ) (n : ℕ) : Ssum (ge 1 v) n = pre (lp v) (n + 1) := by
  rw [show pre (lp v) (n + 1) = Ssum (lp v) n + lp v 0 by rw [Ssum_eq_pre]; ring]
  simp [Ssum, lp]

/-- Symmetric quadratic-form bound: `∑∑ P_nm |x_n||x_m| ≤ ∑_n x_n² ∑_m P_nm`. -/
theorem quad_le (P : ℕ → ℕ → ℝ) (hP : ∀ n m, 0 ≤ P n m) (hs : ∀ n m, P n m = P m n)
    (N : ℕ) (x : ℕ → ℝ) :
    ∑ n ∈ range N, ∑ m ∈ range N, P n m * |x n| * |x m| ≤
      ∑ n ∈ range N, x n ^ 2 * ∑ m ∈ range N, P n m := by
  have h1 : ∀ n m, P n m * |x n| * |x m| ≤ P n m * x n ^ 2 / 2 + P n m * x m ^ 2 / 2 := by
    intro n m
    have : |x n| * |x m| ≤ x n ^ 2 / 2 + x m ^ 2 / 2 := by
      nlinarith [sq_nonneg (|x n| - |x m|), sq_abs (x n), sq_abs (x m)]
    nlinarith [mul_le_mul_of_nonneg_left this (hP n m)]
  calc ∑ n ∈ range N, ∑ m ∈ range N, P n m * |x n| * |x m|
      ≤ ∑ n ∈ range N, ∑ m ∈ range N, (P n m * x n ^ 2 / 2 + P n m * x m ^ 2 / 2) :=
        sum_le_sum (fun n _ => sum_le_sum (fun m _ => h1 n m))
    _ = ∑ n ∈ range N, ∑ m ∈ range N, P n m * x n ^ 2 / 2 +
        ∑ n ∈ range N, ∑ m ∈ range N, P n m * x m ^ 2 / 2 := by
        rw [← sum_add_distrib]; apply sum_congr rfl; intro n _; rw [sum_add_distrib]
    _ = ∑ n ∈ range N, ∑ m ∈ range N, P n m * x n ^ 2 / 2 +
        ∑ n ∈ range N, ∑ m ∈ range N, P n m * x n ^ 2 / 2 := by
        congr 1
        rw [sum_comm]
        apply sum_congr rfl; intro n _; apply sum_congr rfl; intro m _; rw [hs]
    _ = ∑ n ∈ range N, x n ^ 2 * ∑ m ∈ range N, P n m := by
        rw [← two_mul, mul_sum]
        apply sum_congr rfl; intro n _
        rw [mul_sum, mul_sum]
        apply sum_congr rfl; intro m _; ring

/-- `|Fp(ge 1 v)(x)| ≤ 5·(2π)·|x|²`. -/
theorem Fp_bound (v : ℝ) (hv : 0 < v) (N : ℕ) (x : ℕ → ℝ) :
    |Fp (ge 1 v) N x| ≤ 10 * Real.pi * nsq N x := by
  set S := 2 * Real.pi
  set P : ℕ → ℕ → ℝ := fun n m => lp v (adist n m) + lp v (n + m + 1) with hPdef
  have hP : ∀ n m, 0 ≤ P n m := fun n m => add_nonneg (lp_nonneg v hv _) (lp_nonneg v hv _)
  have hs : ∀ n m, P n m = P m n := by
    intro n m; simp only [hPdef]; rw [adist_comm' n m, Nat.add_comm n m]
  have hrow : ∀ n, ∑ m ∈ range N, P n m ≤ 3 * S := fun n =>
    rowP_le (lp v) (lp_nonneg v hv) (by simp [lp]) S (pre_lp_le v hv) n N
  have hdiag : ∀ n, 0 ≤ Ssum (ge 1 v) n ∧ Ssum (ge 1 v) n ≤ S := by
    intro n; rw [Ssum_eq_lp]
    exact ⟨pre_nonneg _ (lp_nonneg v hv) _, pre_lp_le v hv _⟩
  -- Fp = Σ P x x − Σ 2 Ssum x²
  have hsplit : Fp (ge 1 v) N x = ∑ n ∈ range N, ∑ m ∈ range N, P n m * x n * x m -
      ∑ n ∈ range N, 2 * Ssum (ge 1 v) n * x n ^ 2 := by
    unfold Fp
    rw [← sum_sub_distrib]
    apply sum_congr rfl; intro n hn
    rw [show 2 * Ssum (ge 1 v) n * x n ^ 2 =
        ∑ m ∈ range N, (if n = m then 2 * Ssum (ge 1 v) n else 0) * x n * x m by
          rw [sum_eq_single n]
          · simp; ring
          · intro m _ hmn; simp [Ne.symm hmn]
          · intro h; exact absurd hn h]
    rw [← sum_sub_distrib]
    apply sum_congr rfl; intro m _
    rw [Kp_split]; simp only [hPdef, lp]; ring
  have hq := quad_le P hP hs N x
  have hPx : |∑ n ∈ range N, ∑ m ∈ range N, P n m * x n * x m| ≤
      ∑ n ∈ range N, ∑ m ∈ range N, P n m * |x n| * |x m| := by
    refine (abs_sum_le_sum_abs _ _).trans (sum_le_sum (fun n _ => ?_))
    refine (abs_sum_le_sum_abs _ _).trans (le_of_eq (sum_congr rfl (fun m _ => ?_)))
    rw [abs_mul, abs_mul, abs_of_nonneg (hP n m)]
  have hR : ∑ n ∈ range N, x n ^ 2 * ∑ m ∈ range N, P n m ≤ 3 * S * nsq N x := by
    unfold nsq; rw [mul_sum]
    apply sum_le_sum; intro n _
    nlinarith [hrow n, sq_nonneg (x n)]
  have hD0 : 0 ≤ ∑ n ∈ range N, 2 * Ssum (ge 1 v) n * x n ^ 2 :=
    sum_nonneg (fun n _ => by nlinarith [(hdiag n).1, sq_nonneg (x n)])
  have hD1 : ∑ n ∈ range N, 2 * Ssum (ge 1 v) n * x n ^ 2 ≤ 2 * S * nsq N x := by
    unfold nsq; rw [mul_sum]
    apply sum_le_sum; intro n _
    nlinarith [(hdiag n).2, sq_nonneg (x n)]
  rw [hsplit, abs_le]
  have hab := abs_le.mp hPx
  constructor <;> nlinarith [hab.1, hab.2]

theorem H_witness : H 1 (1 / 64) (1 / 2) := by
  intro v hv N x
  rw [Qe_linear]
  have hb := abs_le.mp (Fp_bound v hv N x)
  have hC := C_le N x
  have hn := nsq_nonneg N x
  have hpi := Real.pi_lt_d2
  nlinarith [hb.1, hb.2]

theorem PropN_witness : PropN 1 (1 / 64) := by
  intro w hw N y
  rw [Qe_linear]
  have hb := abs_le.mp (Fp_bound w hw N y)
  have hC := C_le N y
  have hn := nsq_nonneg N y
  have hpi := Real.pi_lt_d2
  nlinarith [hb.1, hb.2]

theorem H_and_PropN_consistent :
    ∃ tp a0 η : ℝ, 0 < a0 ∧ 0 < η ∧ H tp a0 η ∧ PropN tp a0 :=
  ⟨1, 1 / 64, 1 / 2, by norm_num, by norm_num, H_witness, PropN_witness⟩

end

end RL1H

#print axioms RL1H.Fp_bound
#print axioms RL1H.H_and_PropN_consistent

/-! ### A fully hypothesis-free instance of `decide_rejectedProven_sound`
(witness constants tp = 1, a0 = 1/64; T* = 64000 so that floor = a0·T* = 1000; the input and
model are those of RL1Decide.x0 / sp1, duplicated here because scratch files cannot import
each other). -/

namespace RL1H
open TargetLogic

noncomputable section

def cfgW : Config := ⟨1 / 64, 64000, fun _ => 0⟩

def valW : Validity :=
  { EF := some 100000, omegaMax := 1000, maxImagModeK := 0, containsHLiBe := false,
    anharmonicCheck := false, anharmonicShift := 0, alpha2FRecomputedAnharmonic := false,
    muStarNegativeSuspected := false, admittedSmearings := 3, qMeshes := 2, nSheets := 1,
    bandResolvedAvailable := false, dosFractionsKnown := false, kResolvedCheckDone := false }

def xW : Inputs := ⟨.fullDFPT, valW, ⟨500, 1, 1e6, 0.05, 0.05, 0.05, 0.10, 0.02⟩,
  ⟨none, none, none⟩, []⟩

theorem decide_xW : decideT cfgW xW = .rejectedProven := by
  have hG : GValid xW := by
    refine ⟨⟨100000, rfl, by norm_num, ?_⟩, ?_⟩
    · simp only [xW, valW, Envelope.lamHi]; norm_num
    · simp only [xW, valW, Envelope.WF]; norm_num
  have hF : FiresCertified cfgW xW := by
    refine Or.inr (Or.inr (Or.inr ⟨rfl, ?_, Or.inl ?_⟩)) <;>
      simp only [xW, cfgW, Config.floor, Envelope.Ahi] <;> norm_num
  have hB : ¬ BlockedRej xW := by simp [BlockedRej, xW]
  have h1 : xW.tier.proofGrade = true := rfl
  unfold decideT applyPolicy decideCore
  simp [hG, hF, hB, h1]

def spW : MSpec 1 Unit where
  S := {()}
  Ω := fun _ => 1000
  W := fun _ _ _ => 500
  Ω_pos := by intro s _; norm_num
  W_nonneg := by intro s _ i j; norm_num
  W_symm := by intro s _ i j; rfl

theorem sound_xW : SoundMB xW spW (fun _ => 1) where
  kappa := by intro k hk; simp [xW] at hk
  rowmax := by intro r hr; simp [xW] at hr
  dosw := by intro w hw; simp [xW] at hw
  iso := by
    intro _ i
    simp only [rowA, Amat, spW, xW, Envelope.Ahi]
    simp
    norm_num

/-- No hypotheses left: every premise of `decide_rejectedProven_sound` is discharged. -/
theorem fully_instantiated :
    ∀ μ : Fin 1 → Fin 1 → ℝ, CoulPSD (fun _ => (1 : ℝ)) μ → ∀ Nc : WithTop ℕ, ∀ T : ℝ,
      cfgW.Tstar ≤ T → MStable 1 spW (fun _ => 1) μ Nc T := by
  have hH : H 1 (1 / 64) (1 / 2) ∧ PropN 1 (1 / 64) := ⟨H_witness, PropN_witness⟩
  exact decide_rejectedProven_sound cfgW xW 1 (1 / 2) hH.1 (by norm_num) (by simp [cfgW])
    hH.2 (by simp [cfgW]) spW (fun _ => 1) (fun _ => one_pos) sound_xW decide_xW

end

end RL1H

#print axioms RL1H.fully_instantiated

/-! ### `n1_dos_stable` with every hypothesis discharged (two bands, w_min = 1/2) -/
namespace RL1H
open TargetLogic
noncomputable section

def sp2 : MSpec 2 Unit where
  S := {()}
  Ω := fun _ => 1000
  W := fun _ _ _ => 100
  Ω_pos := by intro s _; norm_num
  W_nonneg := by intro s _ i j; norm_num
  W_symm := by intro s _ i j; rfl

theorem n1_dos_inst :
    ∀ μ : Fin 2 → Fin 2 → ℝ, CoulPSD (fun _ => (1 : ℝ)) μ → ∀ Nc : WithTop ℕ, ∀ T : ℝ,
      64000 ≤ T → MStable 1 sp2 (fun _ => 1) μ Nc T :=
  (n1_dos_stable 1 (1 / 64) (1 / 2) H_witness (by norm_num) (by norm_num) PropN_witness sp2
    (fun _ => 1) (fun _ => one_pos) 200 (1 / 2) 64000 (by norm_num)
    (fun i => by simp [wfrac])
    (by simp [Aiso, sp2]; norm_num) (by norm_num) (by norm_num)).2

end
end RL1H

#print axioms RL1H.n1_dos_inst
