/-
TargetLogic.NMG: rule NM-G of `01_eliashberg_targets.md` §0.4 with the numerical content
isolated as hypotheses (in the style of `TcR1.H`).

Model: single isotropic α²F (recorded assumption M-ISO), a finite atomic nonnegative spectrum
`TcR1.Spectrum`; stability is `PipelineLogic.Stable` (neither sector of `TcR1` unstable).

Hypotheses (numerical / derived, NOT proved here):
* `NMHyp κ g T*` — "every spectrum unstable at some `T ≥ T*` (μ* ≥ 0, any cutoff) satisfies
  `m ≥ g(λ)`". This bundles the extreme-point reduction, NUM-EINSTEIN-MIN and the μ* = 0 /
  T-supremum steps of §0.3.
* `AntitoneOn g (Ici 0)` — NUM-G-MONOTONE.
* `gInf ≤ g λ` for all `λ ≥ 0` — the `g(∞)` row (KAY III asymptotics + monotonicity).

Proved: the table lookup at `λ_up ≥ λ_hi` (or `g(∞)`) is conservative, and
`m_hi < g_lookup(λ_hi)`, `m ≤ m_hi`, `λ ≤ λ_hi` ⇒ stable at every `T ≥ T*`.
-/
import PipelineLogic.Common
import TargetLogic.Moments

open Finset

namespace TargetLogic

open TcR1

noncomputable section

/-- The numerical hypothesis of rule NM: every single-band spectrum that is unstable at some
`T ≥ T*` (for some `μ* ≥ 0` and cutoff) has `m ≥ g(λ)`. -/
def NMHyp (κ : ℝ) (g : ℝ → ℝ) (Tstar : ℝ) : Prop :=
  ∀ {ι : Type} (sp : Spectrum ι) (μ : ℝ) (Nc : WithTop ℕ) (T : ℝ),
    0 ≤ μ → Tstar ≤ T → ¬ PipelineLogic.Stable κ sp μ Nc T → g (lamS sp) ≤ momS sp

/-- The lookup actually used by the rule: `g(λ_up)` with `λ_up = up λ_hi ≥ λ_hi` (smallest
tabulated λ ≥ λ_hi) if `λ_hi ≤ λ_cap` (= 100), else `g(∞)`. -/
def gLookup (g : ℝ → ℝ) (up : ℝ → ℝ) (lamCap gInf : ℝ) (lamHi : ℝ) : ℝ :=
  if lamHi ≤ lamCap then g (up lamHi) else gInf

/-- The lookup is conservative: for every true `λ ∈ [0, λ_hi]`, `g_lookup(λ_hi) ≤ g(λ)`. -/
theorem gLookup_le (g : ℝ → ℝ) (hg : AntitoneOn g (Set.Ici 0)) (up : ℝ → ℝ)
    (hup : ∀ l, l ≤ up l) (lamCap gInf : ℝ) (hInf : ∀ l, 0 ≤ l → gInf ≤ g l)
    (lamHi lam : ℝ) (hl0 : 0 ≤ lam) (hl : lam ≤ lamHi) :
    gLookup g up lamCap gInf lamHi ≤ g lam := by
  unfold gLookup
  split_ifs
  · have h1 : lam ≤ up lamHi := hl.trans (hup lamHi)
    exact hg (Set.mem_Ici.mpr hl0) (Set.mem_Ici.mpr (hl0.trans h1)) h1
  · exact hInf lam hl0

/-- **NM-G (logic).** Under `NMHyp`: if `m ≤ m_hi`, and `m_hi < g_up` where `g_up ≤ g(λ)`
(e.g. by `gLookup_le`), the spectrum is stable at every `T ≥ T*`, every `μ* ≥ 0`, every
cutoff. -/
theorem nmg_stable (κ : ℝ) (g : ℝ → ℝ) (Tstar : ℝ) (hNM : NMHyp κ g Tstar)
    {ι : Type} (sp : Spectrum ι) (mHi gUp : ℝ) (hm : momS sp ≤ mHi)
    (hgUp : gUp ≤ g (lamS sp)) (hfire : mHi < gUp) :
    ∀ μ : ℝ, 0 ≤ μ → ∀ Nc : WithTop ℕ, ∀ T : ℝ, Tstar ≤ T → PipelineLogic.Stable κ sp μ Nc T := by
  intro μ hμ Nc T hT
  by_contra hns
  have := hNM sp μ Nc T hμ hT hns
  linarith

/-- **NM-G (rule as implemented).** `NMHyp`, g antitone on `[0,∞)`, `g(∞) ≤ g`, `λ ≤ λ_hi`,
`m ≤ m_hi`, `m_hi < g_lookup(λ_hi)` ⇒ stable at every `T ≥ T*`. -/
theorem nmg_rule_stable (κ : ℝ) (g : ℝ → ℝ) (Tstar : ℝ) (hNM : NMHyp κ g Tstar)
    (hg : AntitoneOn g (Set.Ici 0)) (up : ℝ → ℝ) (hup : ∀ l, l ≤ up l) (lamCap gInf : ℝ)
    (hInf : ∀ l, 0 ≤ l → gInf ≤ g l)
    {ι : Type} (sp : Spectrum ι) (mHi lamHi : ℝ) (hm : momS sp ≤ mHi) (hl : lamS sp ≤ lamHi)
    (hfire : mHi < gLookup g up lamCap gInf lamHi) :
    ∀ μ : ℝ, 0 ≤ μ → ∀ Nc : WithTop ℕ, ∀ T : ℝ, Tstar ≤ T → PipelineLogic.Stable κ sp μ Nc T :=
  nmg_stable κ g Tstar hNM sp mHi _ hm
    (gLookup_le g hg up hup lamCap gInf hInf lamHi (lamS sp) (lamS_nonneg sp) hl) hfire

/-- **N1-ISO (single band, via `TcR1.stable_of_ge`).** `H(a0, η)` (`η > 0`, `a0 > 0`),
`0 < A ≤ A_hi < a0 T*`: stable at every `T ≥ T*`, every `μ* ≥ 0`, every cutoff.
(`0 < A` is the hypothesis of `TcR1.stable_of_ge`.) -/
theorem n1_iso_stable (κ a0 η : ℝ) (hH : H κ a0 η) (ha0 : 0 < a0) (hη : 0 < η)
    {ι : Type*} (sp : Spectrum ι) (hpos : 0 < sp.A) (Ahi Tstar : ℝ) (hA : sp.A ≤ Ahi)
    (hU : Ahi < a0 * Tstar) :
    ∀ μ : ℝ, 0 ≤ μ → ∀ Nc : WithTop ℕ, ∀ T : ℝ, Tstar ≤ T → PipelineLogic.Stable κ sp μ Nc T := by
  intro μ hμ Nc T hT
  apply stable_of_ge κ a0 η hH ha0 hη sp hpos μ hμ Nc T
  rw [div_le_iff₀ ha0]
  have : a0 * Tstar ≤ a0 * T := mul_le_mul_of_nonneg_left hT ha0.le
  linarith
/-! ## Zero spectrum, A ≥ 0 versions, NM-CERT, table lookup -/

theorem C_le_neg_nsq (N : ℕ) (x : ℕ → ℝ) : C N x ≤ -1 * nsq N x := by
  have := C_le N x; linarith

theorem Kp_zero (n m : ℕ) : Kp (fun _ => (0 : ℝ)) n m = 0 := by
  unfold Kp Ssum; split_ifs <;> simp

theorem Km_zero (n m : ℕ) : Km (fun _ => (0 : ℝ)) n m = 0 := by
  unfold Km Ssum; split_ifs <;> simp

/-- The zero spectrum is stable at every `T`, every `μ* ≥ 0`, every cutoff.
(Proof taken from the RL1 review scratch file `docs/reviews/rl1_scratch/RL1NMG.lean`.) -/
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

/-- **N1-ISO, any `A ≥ 0`.** `H(a0, η)` (`η > 0`, `a0 > 0`), `A ≤ A_hi < a0 T*`, `T* > 0`:
stable at every `T ≥ T*`, every `μ* ≥ 0`, every cutoff. -/
theorem n1_iso_stable' (κ a0 η : ℝ) (hH : H κ a0 η) (ha0 : 0 < a0) (hη : 0 < η)
    {ι : Type*} (sp : Spectrum ι) (Ahi Tstar : ℝ) (hA : sp.A ≤ Ahi) (hU : Ahi < a0 * Tstar) :
    ∀ μ : ℝ, 0 ≤ μ → ∀ Nc : WithTop ℕ, ∀ T : ℝ, Tstar ≤ T → PipelineLogic.Stable κ sp μ Nc T := by
  intro μ hμ Nc T hT
  rcases lt_or_eq_of_le (A_nonneg sp) with hpos | hzero
  · exact n1_iso_stable κ a0 η hH ha0 hη sp hpos Ahi Tstar hA hU μ hμ Nc T hT
  · have hz : ∀ j ∈ sp.S, sp.p j = 0 :=
      (sum_eq_zero_iff_of_nonneg sp.p_nonneg).mp hzero.symm
    exact stable_of_zero κ sp hz μ hμ Nc T

/-- **NM-CERT (single band).** `H(a0, η)`, `m ≤ m_hi`, `0 < λ ≤ λ_hi` … more precisely
`λ ≤ λ_hi`, `0 < λ_hi`, and `m_hi < 4 (a0 T*)² / λ_hi` ⇒ stable at every `T ≥ T*`.
Proof: `A² ≤ λ m / 4 ≤ λ_hi m_hi / 4 < (a0 T*)²` (L-MOM2), then N1-ISO. -/
theorem nmcert_stable (κ a0 η : ℝ) (hH : H κ a0 η) (ha0 : 0 < a0) (hη : 0 < η)
    {ι : Type*} (sp : Spectrum ι) (mHi lamHi Tstar : ℝ) (hT0 : 0 < Tstar)
    (hm : momS sp ≤ mHi) (hl : lamS sp ≤ lamHi) (hl0 : 0 < lamHi)
    (hfire : mHi < 4 * (a0 * Tstar) ^ 2 / lamHi) :
    ∀ μ : ℝ, 0 ≤ μ → ∀ Nc : WithTop ℕ, ∀ T : ℝ, Tstar ≤ T → PipelineLogic.Stable κ sp μ Nc T := by
  have hF : 0 < a0 * Tstar := mul_pos ha0 hT0
  have h1 : mHi * lamHi < 4 * (a0 * Tstar) ^ 2 := by
    rw [lt_div_iff₀ hl0] at hfire; linarith
  have hA2 := lmom2_sq sp
  have hL := lamS_nonneg sp
  have hM := momS_nonneg sp
  have h2 : lamS sp * momS sp ≤ lamHi * mHi :=
    mul_le_mul hl hm hM (hL.trans hl)
  have hA0 := A_nonneg sp
  have hsq : sp.A ^ 2 < (a0 * Tstar) ^ 2 := by nlinarith
  have hlt : sp.A < a0 * Tstar := by nlinarith
  exact n1_iso_stable' κ a0 η hH ha0 hη sp sp.A Tstar le_rfl hlt

/-- Table lookup used by NM-G: the first entry `(λ_t, g_t)` with `λ_hi ≤ λ_t`, else `g(∞)`.
(The §0.3 table is sorted by λ and ends at λ = 100, so `λ_hi > 100` gives `g(∞)`.) -/
def tableLookup {K : Type*} [LinearOrder K] (tbl : List (K × K)) (gInf lh : K) : K :=
  match tbl.find? (fun p => decide (lh ≤ p.1)) with
  | some p => p.2
  | none => gInf

/-- The table lookup is conservative when every table value is `≤ g` at its λ, `g` is
antitone on `[0, ∞)` and `g(∞) ≤ g`. -/
theorem tableLookup_le (g : ℝ → ℝ) (hg : AntitoneOn g (Set.Ici 0)) (tbl : List (ℝ × ℝ))
    (htbl : ∀ p ∈ tbl, p.2 ≤ g p.1) (gInf : ℝ) (hInf : ∀ l, 0 ≤ l → gInf ≤ g l)
    (lh l : ℝ) (hl0 : 0 ≤ l) (hl : l ≤ lh) : tableLookup tbl gInf lh ≤ g l := by
  unfold tableLookup
  cases hf : tbl.find? (fun p => decide (lh ≤ p.1)) with
  | none => exact hInf l hl0
  | some p =>
    have hmem : p ∈ tbl := List.mem_of_find?_eq_some hf
    have hp : lh ≤ p.1 := by simpa using List.find?_some hf
    have h1 : l ≤ p.1 := hl.trans hp
    exact (htbl p hmem).trans (hg (Set.mem_Ici.mpr hl0) (Set.mem_Ici.mpr (hl0.trans h1)) h1)

end

end TargetLogic
