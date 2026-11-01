/-
TargetLogic.Monotone: worse data never upgrades (RL1 findings 2, 6, 7c).

`Worse x x'`: `x'` has the same validity data, tags, nominal values and provenance as `x`, but
larger spreads δ_A, δ_λ, δ_m, δ_𝒜, δ_w and larger ε_sys, ε_ω (any q-spread), and possibly
drops the band-resolved matrix or the sheet-DOS data. Then, if `x` passes both gates, the rank
of `decideT` (REJECTED_PROVEN 2 > screenedOut 1 > UNDECIDED / PASS 0) does not increase.
-/
import TargetLogic.Decide

open Finset

set_option linter.unusedSectionVars false

namespace TargetLogic

section

variable {K : Type*} [Field K] [LinearOrder K] [IsStrictOrderedRing K]

def MatWorse : Option (MatData K) → Option (MatData K) → Prop
  | _, none => True
  | some md, some md' => ∃ d, md' = { md with dCalA := d } ∧ md.dCalA ≤ d
  | none, some _ => False

def DosWorse : Option (DosData K) → Option (DosData K) → Prop
  | _, none => True
  | some dd, some dd' => ∃ d, dd' = { dd with dW := d } ∧ dd.dW ≤ d
  | none, some _ => False

/-- `x'` is `x` with wider envelopes and possibly less multiband data. -/
structure Worse (x x' : Inputs K) : Prop where
  val : x'.val = x.val
  tags : x'.tags = x.tags
  iso : ∃ a b c, x'.iso = { x.iso with dA := a, dLam := b, dM := c } ∧
    x.iso.dA ≤ a ∧ x.iso.dLam ≤ b ∧ x.iso.dM ≤ c
  eps : ∃ s o q, x'.eps = { x.eps with epsSys := s, epsOmega := o, qSpread := q } ∧
    x.eps.epsSys ≤ s ∧ x.eps.epsOmega ≤ o
  mat : MatWorse x.mat x'.mat
  dos : DosWorse x.dos x'.dos

theorem frac_mono {n n' d d' : K} (hn : 0 ≤ n) (hnn : n ≤ n') (hd' : 0 < d') (hdd : d' ≤ d) :
    n / d ≤ n' / d' :=
  (div_le_div_of_nonneg_left hn hd' hdd).trans (div_le_div_of_nonneg_right hnn hd'.le)

theorem prod3_mono {a u u' s s' : K} (ha : 0 ≤ a) (hu : 0 ≤ u) (huu : u ≤ u') (hs : 0 ≤ s)
    (hss : s ≤ s') : a * (1 + u) * (1 + s) ≤ a * (1 + u') * (1 + s') := by
  have h1 : a * (1 + u) ≤ a * (1 + u') := mul_le_mul_of_nonneg_left (by linarith) ha
  exact mul_le_mul h1 (by linarith) (by linarith) (mul_nonneg ha (by linarith))

/-- Facts used repeatedly: from `wf` of both inputs. -/
theorem env_mono {x x' : Inputs K} (hW : Worse x x') (hx : wf x = true) (hx' : wf x' = true) :
    Ahi x ≤ Ahi x' ∧ lamHi x ≤ lamHi x' ∧ mHi x ≤ mHi x' := by
  obtain ⟨a, b, c, hi, ha, hb, hc⟩ := hW.iso
  obtain ⟨s, o, q, he, hs, ho⟩ := hW.eps
  obtain ⟨hA, hl, hm, hdA, hdl, hdm, hes, heo0, heo1⟩ := wf_spec hx
  obtain ⟨-, -, -, -, -, -, hes', heo0', heo1'⟩ := wf_spec hx'
  have ho1 : o < 1 := by have h := heo1'; rw [he] at h; exact h
  refine ⟨?_, ?_, ?_⟩
  · unfold Ahi; rw [hi, he]; dsimp only
    exact frac_mono (by positivity) (prod3_mono hA.le hdA ha (by linarith) hs)
      (by linarith) (by linarith)
  · unfold lamHi; rw [hi, he]; dsimp only
    have h0 : 0 < 1 - o := by linarith
    exact frac_mono (by positivity) (prod3_mono hl.le hdl hb (by linarith) hs)
      (pow_pos h0 2) (by nlinarith)
  · unfold mHi; rw [hi, he]; dsimp only
    exact prod3_mono hm.le hdm hc (by linarith) hs

theorem kResolved_mono {x x' : Inputs K} (hW : Worse x x') (h : kResolved x' = true) :
    kResolved x = true := by
  have hm := hW.mat
  unfold kResolved at *
  cases hx : x.mat <;> cases hx' : x'.mat <;> simp_all [MatWorse]

/-- Each floor test is monotone: if it holds on the worse data it holds on `x`. -/
theorem floor_mono (cfg : Config K)
    (hg : ∀ a b, 0 ≤ a → a ≤ b → cfg.gUp b ≤ cfg.gUp a)
    {x x' : Inputs K} (hW : Worse x x') (hx : wf x = true) (hx' : wf x' = true) (r : Rule)
    (h : floorR cfg x' r = true) : floorR cfg x r = true := by
  obtain ⟨hAm, hlm, hmm⟩ := env_mono hW hx hx'
  have hval := hW.val
  obtain ⟨a, b, c, hi, ha, hb, hc⟩ := hW.iso
  obtain ⟨s, o, q, he, hs, ho⟩ := hW.eps
  obtain ⟨hA, hl, hm, hdA, hdl, hdm, hes, heo0, heo1⟩ := wf_spec hx
  obtain ⟨-, -, -, -, -, -, hes', heo0', heo1'⟩ := wf_spec hx'
  have ho1 : o < 1 := by have h := heo1'; rw [he] at h; exact h
  have hkr := kResolved_mono hW
  cases r with
  | kappa =>
    have hmw := hW.mat
    simp only [floorR] at h ⊢
    cases hxm : x.mat with
    | none => rw [hxm] at hmw; cases hxm' : x'.mat <;> simp_all [MatWorse]
    | some md =>
      cases hxm' : x'.mat with
      | none => simp_all
      | some md' =>
        rw [hxm, hxm'] at hmw
        obtain ⟨d, rfl, hd⟩ := hmw
        rw [hxm'] at h
        simp only [decide_eq_true_eq] at h ⊢
        have hs0 := (matWF_spec' (hxm ▸ (wf_mat' hx)))
        refine lt_of_le_of_lt ?_ h
        apply finMax_mono
        intro i
        apply div_le_div_of_nonneg_right _ (hs0.2.2.2.1 i).le
        apply sum_le_sum
        intro j _
        apply mul_le_mul_of_nonneg_right _ (hs0.2.2.2.1 j).le
        unfold calAhi
        rw [he]
        dsimp only
        exact frac_mono (by have := hs0.2.1 i j; have := hs0.2.2.2.2.1; positivity)
          (prod3_mono (hs0.2.1 i j) hs0.2.2.2.2.1 hd (by linarith) hs)
          (by linarith) (by linarith)
  | rowmax =>
    have hmw := hW.mat
    simp only [floorR] at h ⊢
    cases hxm : x.mat with
    | none => rw [hxm] at hmw; cases hxm' : x'.mat <;> simp_all [MatWorse]
    | some md =>
      cases hxm' : x'.mat with
      | none => simp_all
      | some md' =>
        rw [hxm, hxm'] at hmw
        obtain ⟨d, rfl, hd⟩ := hmw
        rw [hxm'] at h
        simp only [decide_eq_true_eq] at h ⊢
        have hs0 := (matWF_spec' (hxm ▸ (wf_mat' hx)))
        refine lt_of_le_of_lt ?_ h
        apply finMax_mono
        intro i
        apply sum_le_sum
        intro j _
        unfold calAhi
        rw [he]
        dsimp only
        exact frac_mono (by have := hs0.2.1 i j; have := hs0.2.2.2.2.1; positivity)
          (prod3_mono (hs0.2.1 i j) hs0.2.2.2.2.1 hd (by linarith) hs)
          (by linarith) (by linarith)
  | dos =>
    have hdw := hW.dos
    simp only [floorR] at h ⊢
    cases hxd : x.dos with
    | none => rw [hxd] at hdw; cases hxd' : x'.dos <;> simp_all [DosWorse]
    | some dd =>
      cases hxd' : x'.dos with
      | none => simp_all
      | some dd' =>
        rw [hxd, hxd'] at hdw
        obtain ⟨d, rfl, hd⟩ := hdw
        rw [hxd'] at h
        simp only [Bool.and_eq_true, decide_eq_true_eq] at h ⊢
        obtain ⟨hw0', hlt⟩ := h
        have hs0 := dosWF_spec' (hxd ▸ wf_dos' hx)
        have hwle : wminLo { dd with dW := d } ≤ wminLo dd := by
          apply finMin1_mono
          intro i
          exact mul_le_mul_of_nonneg_left (by linarith) (hs0.2.1 i).1.le
        refine ⟨lt_of_lt_of_le hw0' hwle, lt_of_le_of_lt ?_ hlt⟩
        exact frac_mono ((wf_spec hx).1.le.trans (Ahi_ge x hx)) hAm hw0' hwle
  | iso =>
    simp only [floorR, Bool.and_eq_true, Bool.or_eq_true, decide_eq_true_eq] at h ⊢
    obtain ⟨⟨h1, h2⟩, h3⟩ := h
    refine ⟨⟨by rw [← hval]; exact h1, lt_of_le_of_lt hAm h2⟩, ?_⟩
    rcases h3 with h3 | h3
    · exact Or.inl (lt_of_le_of_lt hAm h3)
    · exact Or.inr (hkr h3)
  | nmCert =>
    simp only [floorR, Bool.and_eq_true, decide_eq_true_eq] at h ⊢
    obtain ⟨h1, h2⟩ := h
    refine ⟨by rw [← hval]; exact h1, lt_of_le_of_lt hmm (lt_of_lt_of_le h2 ?_)⟩
    exact div_le_div_of_nonneg_left (by positivity) (lamHi_pos x hx) hlm
  | nmG =>
    simp only [floorR, Bool.and_eq_true, Bool.or_eq_true, decide_eq_true_eq] at h ⊢
    obtain ⟨⟨h1, h2⟩, h3⟩ := h
    have hgg : cfg.gUp (lamHi x') ≤ cfg.gUp (lamHi x) := hg _ _ (lamHi_pos x hx).le hlm
    refine ⟨⟨by rw [← hval]; exact h1, by linarith⟩, ?_⟩
    rcases h3 with h3 | h3
    · exact Or.inl (by nlinarith)
    · exact Or.inr (hkr h3)

theorem provOf_worse {x x' : Inputs K} (hW : Worse x x') (n : InputName) :
    provOf x' n = none ∨ provOf x' n = provOf x n := by
  obtain ⟨a, b, c, hi, -⟩ := hW.iso
  obtain ⟨s, o, q, he, -⟩ := hW.eps
  have hmw := hW.mat
  have hdw := hW.dos
  cases n with
  | isoA => right; simp [provOf, hi]
  | isoLam => right; simp [provOf, hi]
  | isoM => right; simp [provOf, hi]
  | eps => right; simp [provOf, he]
  | mat =>
    simp only [provOf]
    cases hxm : x.mat <;> cases hxm' : x'.mat <;> rw [hxm, hxm'] at hmw <;> simp_all [MatWorse]
    obtain ⟨d, rfl, -⟩ := hmw; rfl
  | dos =>
    simp only [provOf]
    cases hxd : x.dos <;> cases hxd' : x'.dos <;> rw [hxd, hxd'] at hdw <;> simp_all [DosWorse]
    obtain ⟨d, rfl, -⟩ := hdw; rfl

theorem provOK_mono {x x' : Inputs K} (hW : Worse x x') (r : Rule) (h : provOK x' r = true) :
    provOK x r = true := by
  unfold provOK at *
  obtain ⟨p, -, hall⟩ := fullSameRun_spec _ h
  have hl : r.uses.map (provOf x') = r.uses.map (provOf x) := by
    apply List.map_congr_left
    intro n hn
    have hs := hall _ (List.mem_map_of_mem hn)
    rcases provOf_worse hW n with h0 | h0
    · rw [h0] at hs; cases hs
    · exact h0
  rw [← hl]; exact h

/-- **Worse data never upgrades.** If `x` passes WELL-FORMED and G-VALID, `cfg.gUp` is
antitone on `[0, ∞)` (true for a table with nonincreasing g values and `g(∞)` at most the last
one), and `x'` is `x` with wider envelopes / less multiband data, then
`rank (decideT x') ≤ rank (decideT x)`. -/
theorem decide_worse (cfg : Config K) (hg : ∀ a b, 0 ≤ a → a ≤ b → cfg.gUp b ≤ cfg.gUp a)
    {x x' : Inputs K} (hW : Worse x x') (hx : (wf x && gvalid x) = true) :
    (decideT cfg x').rank ≤ (decideT cfg x).rank := by
  have hbl : blocked x' = blocked x := by simp [blocked, hW.tags]
  have hxw : wf x = true := by simp only [Bool.and_eq_true] at hx; exact hx.1
  unfold decideT
  rw [hbl]
  apply policy_rank_mono
  have hcases : (decideCore cfg x').rank = 0 ∨ (decideCore cfg x').rank = 1 ∨
      (decideCore cfg x').rank = 2 := by
    rcases decideCore_cases cfg x' with ⟨r, h⟩ | h | h | h <;> rw [h] <;> simp [Outcome.rank]
  rcases hcases with h0 | h1 | h2
  · rw [h0]; exact Nat.zero_le _
  · obtain ⟨hg', r, hr⟩ := core_of_rank_one (cfg := cfg) (x := x') (by omega)
    have hx'w : wf x' = true := by simp only [Bool.and_eq_true] at hg'; exact hg'.1
    have := core_rank_one hx r (floor_mono cfg hg hW hxw hx'w r hr)
    omega
  · obtain ⟨hg', r, hr⟩ := core_of_rank_two (cfg := cfg) (x := x') h2
    have hx'w : wf x' = true := by simp only [Bool.and_eq_true] at hg'; exact hg'.1
    simp only [fires, Bool.and_eq_true] at hr
    have hfx : fires cfg x r = true := by
      simp only [fires, Bool.and_eq_true]
      exact ⟨floor_mono cfg hg hW hxw hx'w r hr.1, provOK_mono hW r hr.2⟩
    rw [h2, core_rank_two hx r hfx]

end

end TargetLogic
