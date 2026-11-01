/-
RL1 scratch 1 (review of TargetLogic.Decide): concrete inputs.

* `ex_proven`        : a concrete, doc-conformant input reaches `rejectedProven` (N1-ISO).
* `ex_sound_inst`    : the same input together with a concrete one-band model satisfies every
                       hypothesis of `decide_rejectedProven_sound` except the two numerical
                       hypotheses `H` and `PropN` (kept abstract) -> non-vacuity modulo H, PropN.
* `ex_unstable`      : the conclusion `MStable` is not a tautology (an unstable one-band model at
                       T = 300 for every `tp`).
* counterexamples    : inputs that pass `GValid` and reach `rejectedProven` although the doc
                       says they must not (negative A_nom; eps_sys = eps_omega = 0; negative
                       imaginary-mode frequency; negative anharmonic shift).
* `ex_B0`            : the soundness hypothesis is satisfiable by the empty (B = 0) model for an
                       input whose `SoundMB` is unsatisfiable at B = 1.
-/
import TargetLogic

open TargetLogic MultiBand TcR1 Finset

namespace RL1

noncomputable section

/-- Doc constants: a0 = 4.33148145, T* = 300 K, g lookup irrelevant here (g(inf)). -/
def cfg : Config := ⟨4.33148145, 300, fun _ => 2.6955e6⟩

/-- Generic: proof grade, G-VALID, a certified floor and no blocking tag give rejectedProven. -/
theorem proven_of (c : Config) (x : Inputs) (h1 : x.tier.proofGrade = true) (hG : GValid x)
    (hF : FiresCertified c x) (hB : ¬ BlockedRej x) : decideT c x = .rejectedProven := by
  unfold decideT applyPolicy decideCore
  simp [hG, hF, hB, h1]

def val0 : Validity :=
  { EF := some 100000, omegaMax := 1000, maxImagModeK := 0, containsHLiBe := false,
    anharmonicCheck := false, anharmonicShift := 0, alpha2FRecomputedAnharmonic := false,
    muStarNegativeSuspected := false, admittedSmearings := 3, qMeshes := 2, nSheets := 1,
    bandResolvedAvailable := false, dosFractionsKnown := false, kResolvedCheckDone := false }

/-- Doc-conformant envelope: eps_sys = 0.10, eps_omega = 0.02. A_hi = 500*1.05*1.1/0.98 ≈ 589. -/
def env0 : Envelope := ⟨500, 1, 1e6, 0.05, 0.05, 0.05, 0.10, 0.02⟩

def x0 : Inputs := ⟨.fullDFPT, val0, env0, ⟨none, none, none⟩, []⟩

theorem gvalid_x0 : GValid x0 := by
  refine ⟨⟨100000, rfl, by norm_num, ?_⟩, ?_⟩
  · simp only [x0, env0, val0, Envelope.lamHi]; norm_num
  · simp only [x0, env0, val0, Envelope.WF]; norm_num

theorem ex_proven : decideT cfg x0 = .rejectedProven := by
  have hF : FiresCertified cfg x0 := by
    refine Or.inr (Or.inr (Or.inr ⟨rfl, ?_, Or.inl ?_⟩)) <;>
      simp only [x0, env0, cfg, Config.floor, Envelope.Ahi] <;> norm_num
  have hB : ¬ BlockedRej x0 := by simp [BlockedRej, x0]
  exact proven_of cfg x0 rfl gvalid_x0 hF hB

/-- One band, one atom: W = 500 at Ω = 1000, N_0 = 1, so 𝒜 = (500) and A_iso = 500. -/
def sp1 : MSpec 1 Unit where
  S := {()}
  Ω := fun _ => 1000
  W := fun _ _ _ => 500
  Ω_pos := by intro s _; norm_num
  W_nonneg := by intro s _ i j; norm_num
  W_symm := by intro s _ i j; rfl

def dos1 : Fin 1 → ℝ := fun _ => 1

theorem sound_x0 : SoundMB x0 sp1 dos1 where
  kappa := by intro k hk; simp [x0] at hk
  rowmax := by intro r hr; simp [x0] at hr
  dosw := by intro w hw; simp [x0] at hw
  iso := by
    intro _ i
    simp only [rowA, Amat, sp1, dos1, x0, env0, Envelope.Ahi]
    simp
    norm_num

/-- Every hypothesis of `decide_rejectedProven_sound` other than `H`, `PropN` holds for a
concrete input/model pair; the theorem therefore has content as soon as `H ∧ PropN` is
consistent (see `RL1HProp.lean`). -/
theorem ex_sound_inst (tp η : ℝ) (hH : H tp cfg.a0 η) (hη : 0 ≤ η) (hN : PropN tp cfg.a0) :
    ∀ μ : Fin 1 → Fin 1 → ℝ, CoulPSD dos1 μ → ∀ Nc : WithTop ℕ, ∀ T : ℝ, cfg.Tstar ≤ T →
      MStable tp sp1 dos1 μ Nc T :=
  decide_rejectedProven_sound cfg x0 tp η hH hη (by simp [cfg]; norm_num) hN
    (by simp [cfg]) sp1 dos1 (fun _ => one_pos) sound_x0 ex_proven

/-! ### The conclusion is not a tautology -/

/-- A strongly coupled one-band model (W = 300(1+tp²) at Ω = 300) is unstable at T = 300
(μ = 0), for every `tp`. So `MStable … T` is a real restriction. -/
def spU (tp : ℝ) : MSpec 1 Unit where
  S := {()}
  Ω := fun _ => 300
  W := fun _ _ _ => 300 * (1 + tp ^ 2)
  Ω_pos := by intro s _; norm_num
  W_nonneg := by intro s _ i j; positivity
  W_symm := by intro s _ i j; rfl

theorem ex_unstable (tp : ℝ) :
    ¬ MStable tp (spU tp) dos1 (fun _ _ => 0) ⊤ 300 := by
  intro ⟨hP, _⟩
  apply hP
  intro ε hε
  refine ⟨1, fun _ _ => 1, ⟨0, 0, by simp, one_ne_zero⟩, ?_⟩
  have hq : Qplus tp (spU tp) dos1 (fun _ _ => 0) ⊤ 300 1 (fun _ _ => 1) = 1 := by
    have h1 : (0:ℝ) < 1 + tp ^ 2 := by positivity
    simp [Qplus, spU, dos1, Cw, Bp, Dq, Dn, Ssum, coul, lw, ge, adist]
    field_simp
    ring
  have hn : nrm2 dos1 1 (fun _ _ => 1) = 1 := by simp [nrm2, dos1, nsq]
  rw [hq, hn]; linarith

/-! ### Counterexamples: rejectedProven on inputs the doc excludes -/

/-- (C1) Negative nominal A (e.g. a signed, non-admitted α²F): G-VALID has no sign check. -/
def xNeg : Inputs := ⟨.fullDFPT, val0, ⟨-1, 1, 1e6, 0.05, 0.05, 0.05, 0.10, 0.02⟩,
  ⟨none, none, none⟩, []⟩

theorem ce_negA : decideT cfg xNeg = .rejectedProven := by
  have hG : GValid xNeg := by
    refine ⟨⟨100000, rfl, by norm_num, ?_⟩, ?_⟩
    · simp only [xNeg, val0, Envelope.lamHi]; norm_num
    · simp only [xNeg, val0, Envelope.WF]; norm_num
  have hF : FiresCertified cfg xNeg := by
    refine Or.inr (Or.inr (Or.inr ⟨rfl, ?_, Or.inl ?_⟩)) <;>
      simp only [xNeg, cfg, Config.floor, Envelope.Ahi] <;> norm_num
  have hB : ¬ BlockedRej xNeg := by simp [BlockedRej, xNeg]
  exact proven_of cfg xNeg rfl hG hF hB

/-- …and no B = 1 model is sound for it, so `decide_rejectedProven_sound` says nothing. -/
theorem ce_negA_no_model (sp : MSpec 1 Unit) (dos : Fin 1 → ℝ) (hdos : ∀ i, 0 < dos i) :
    ¬ SoundMB xNeg sp dos := by
  intro hS
  have h1 := hS.iso rfl 0
  have h2 := rowA_nonneg sp dos hdos 0
  have h3 : xNeg.env.Ahi < 0 := by simp only [xNeg, Envelope.Ahi]; norm_num
  linarith

/-- …but the empty model (B = 0) is "sound" for it, and every B = 0 model is stable. -/
theorem ex_B0_sound (sp : MSpec 0 Unit) (dos : Fin 0 → ℝ) : SoundMB xNeg sp dos where
  kappa := by intro k hk; simp [xNeg] at hk
  rowmax := by intro r hr; simp [xNeg] at hr
  dosw := by intro w hw; simp [xNeg] at hw
  iso := by intro _ i; exact i.elim0

theorem ex_B0_stable (tp : ℝ) (sp : MSpec 0 Unit) (dos : Fin 0 → ℝ) μ Nc T :
    MStable tp sp dos μ Nc T := by
  constructor <;>
  · intro hU
    obtain ⟨N, φ, ⟨i, _⟩, _⟩ := hU 1 one_pos
    exact i.elim0

/-- (C2) The doc fixes eps_sys = 0.10 and eps_omega ≥ 0.02; `Envelope.WF` accepts 0 for both.
With A_nom = 1250, δ_A = 0.03: A_hi(Lean, eps = 0) = 1287.5 < 1299.44 → rejectedProven, while
the doc envelope gives A_hi = 1250·1.03·1.1/0.98 ≈ 1445 > 1299.44 → pass. -/
def valK : Validity := { val0 with kResolvedCheckDone := true }
def xEps0 : Inputs := ⟨.fullDFPT, valK, ⟨1250, 1, 1e6, 0.03, 0.03, 0.03, 0, 0⟩,
  ⟨none, none, none⟩, []⟩
def xEpsDoc : Inputs := ⟨.fullDFPT, valK, ⟨1250, 1, 1e6, 0.03, 0.03, 0.03, 0.10, 0.02⟩,
  ⟨none, none, none⟩, []⟩

theorem ce_eps0 : decideT cfg xEps0 = .rejectedProven := by
  have hG : GValid xEps0 := by
    refine ⟨⟨100000, rfl, by norm_num, ?_⟩, ?_⟩
    · simp only [xEps0, valK, val0, Envelope.lamHi]; norm_num
    · simp only [xEps0, valK, val0, Envelope.WF]; norm_num
  have hF : FiresCertified cfg xEps0 := by
    refine Or.inr (Or.inr (Or.inr ⟨rfl, ?_, Or.inr rfl⟩))
    simp only [xEps0, cfg, Config.floor, Envelope.Ahi]; norm_num
  have hB : ¬ BlockedRej xEps0 := by simp [BlockedRej, xEps0]
  exact proven_of cfg xEps0 rfl hG hF hB

theorem ce_epsDoc_not_fired : ¬ FiresCertified cfg xEpsDoc := by
  rintro (⟨k, hk, -⟩ | ⟨r, hr, -⟩ | ⟨hd, -⟩ | ⟨-, h, -⟩)
  · simp [xEpsDoc] at hk
  · simp [xEpsDoc] at hr
  · simp [xEpsDoc, valK, val0] at hd
  · simp only [xEpsDoc, cfg, Config.floor, Envelope.Ahi] at h; norm_num at h

/-- (C3) Sign conventions: QE prints imaginary modes as NEGATIVE frequencies; G-VALID tests
`maxImagModeK ≤ 30`, so −500 K passes. A negative (softening) anharmonic shift of −30 % with
no α²F recomputation also passes (the doc uses |shift|). -/
def valSign : Validity :=
  { val0 with maxImagModeK := -500, containsHLiBe := true, anharmonicCheck := true,
              anharmonicShift := -0.3, alpha2FRecomputedAnharmonic := false }
def xSign : Inputs := ⟨.fullDFPT, valSign, env0, ⟨none, none, none⟩, []⟩

theorem ce_sign : decideT cfg xSign = .rejectedProven := by
  have hG : GValid xSign := by
    refine ⟨⟨100000, rfl, by norm_num, ?_⟩, ?_⟩
    · simp only [xSign, env0, valSign, val0, Envelope.lamHi]; norm_num
    · simp only [xSign, env0, valSign, val0, Envelope.WF]; norm_num
  have hF : FiresCertified cfg xSign := by
    refine Or.inr (Or.inr (Or.inr ⟨rfl, ?_, Or.inl ?_⟩)) <;>
      simp only [xSign, env0, cfg, Config.floor, Envelope.Ahi] <;> norm_num
  have hB : ¬ BlockedRej xSign := by simp [BlockedRej, xSign]
  exact proven_of cfg xSign rfl hG hF hB

/-- (C4) Mixed provenance: one `tier` covers every number. A κ_hi of any origin (ML, coarse)
attached to a full-DFPT isotropic record yields rejectedProven via N1-KAPPA. -/
def xMix : Inputs := ⟨.fullDFPT, { val0 with nSheets := 2, bandResolvedAvailable := true },
  ⟨5000, 1, 1e6, 0.05, 0.05, 0.05, 0.10, 0.02⟩, ⟨some 1000, none, none⟩, []⟩

theorem ce_mix : decideT cfg xMix = .rejectedProven := by
  have hG : GValid xMix := by
    refine ⟨⟨100000, rfl, by norm_num, ?_⟩, ?_⟩
    · simp only [xMix, val0, Envelope.lamHi]; norm_num
    · simp only [xMix, val0, Envelope.WF]; norm_num
  have hF : FiresCertified cfg xMix := by
    refine Or.inl ⟨1000, rfl, ?_⟩
    simp only [cfg, Config.floor]; norm_num
  have hB : ¬ BlockedRej xMix := by simp [BlockedRej, xMix]
  exact proven_of cfg xMix rfl hG hF hB

end

end RL1

#print axioms RL1.ex_proven
#print axioms RL1.ex_sound_inst
#print axioms RL1.ex_unstable
#print axioms RL1.ce_negA
#print axioms RL1.ce_negA_no_model
#print axioms RL1.ex_B0_sound
#print axioms RL1.ex_B0_stable
#print axioms RL1.ce_eps0
#print axioms RL1.ce_epsDoc_not_fired
#print axioms RL1.ce_sign
#print axioms RL1.ce_mix
