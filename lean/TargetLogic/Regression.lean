/-
TargetLogic.Regression: the RL1 counterexamples (`docs/reviews/rl1_scratch/RL1Decide.lean`,
`ce_*`) re-encoded for the revised input format, each now proved NOT to give REJECTED_PROVEN,
plus a positive control (a doc-conformant input that is rejected by N1-ISO), so the gates are
not vacuous.

Config: `Config.ofCert a0Cert 300` (certified NM-G table).
-/
import TargetLogic.Sound

namespace TargetLogic.Regression

open TargetLogic

noncomputable section

def cfgDoc : Config ℝ := Config.ofCert a0Cert 300

def pf (r : ℕ) : Prov := ⟨.fullDFPT, r⟩

def val0 : Validity ℝ :=
  { EF := some 100000, omegaMax := 1000, imagAbsMaxK := 0, containsHLiBe := false,
    anharmonicCheck := false, anharmShiftAbs := 0, alpha2FRecomputedAnharmonic := false,
    muStarNegativeSuspected := false, admittedSmearings := 3, qMeshes := 2, nSheets := 1 }

def iso0 (A : ℝ) : IsoData ℝ := ⟨A, 1, 1000000, 1/20, 1/20, 1/20, pf 1, pf 1, pf 1⟩

/-- Doc envelope: ε_sys = 0.10, ε_ω = 0.02. -/
def eps0 : EpsData ℝ := ⟨1/10, 1/50, 0, pf 1⟩

/-- Positive control: A_nom = 500, A_hi = 500·1.05·1.1/0.98 ≈ 589 < 0.75·1299.44. -/
def x0 : Inputs ℝ := ⟨val0, iso0 500, eps0, none, none, []⟩

theorem wf_x0 : wf x0 = true := by
  norm_num [wf, x0, val0, iso0, eps0, in0half, matWF, dosWF]

theorem gvalid_x0 : gvalid x0 = true := by
  norm_num [gvalid, x0, val0, iso0, eps0, lamHi, mbPresent]

theorem ex_proven : decideT cfgDoc x0 = .rejectedProven .iso := by
  have hiso : fires cfgDoc x0 .iso = true := by
    norm_num [fires, floorR, provOK, fullSameRun, provOf, Rule.uses, x0, val0, iso0, eps0, Ahi,
      cfgDoc, Config.ofCert, Config.floor, a0Cert, kResolved, pf]
  have hk : fires cfgDoc x0 .kappa = false := by simp [fires, floorR, x0]
  have hr : fires cfgDoc x0 .rowmax = false := by simp [fires, floorR, x0]
  have hd : fires cfgDoc x0 .dos = false := by simp [fires, floorR, x0]
  have hc : (wf x0 && gvalid x0) = true := by rw [wf_x0, gvalid_x0]; rfl
  have hf : firstRule cfgDoc x0 = some .iso := by
    simp only [firstRule, allRules, List.find?_cons, hk, hr, hd, hiso]
  have hb : blocked x0 = false := rfl
  unfold decideT decideCore
  rw [if_pos hc, hf, hb]
  rfl

/-- (C1, RL1 `ce_negA`) Negative nominal A: WELL-FORMED fails ⇒ UNDECIDED. -/
def xNeg : Inputs ℝ := ⟨val0, iso0 (-1), eps0, none, none, []⟩

theorem ce_negA_fixed : decideT cfgDoc xNeg = .undecided :=
  decide_of_not_wf _ _ (by norm_num [wf, xNeg, iso0])

/-- (C3, RL1 `ce_sign`) Imaginary mode entered with its sign (−500 K): WELL-FORMED fails. -/
def xSignImag : Inputs ℝ := ⟨{ val0 with imagAbsMaxK := -500 }, iso0 500, eps0, none, none, []⟩

theorem ce_signImag_fixed : decideT cfgDoc xSignImag = .undecided :=
  decide_of_not_wf _ _ (by norm_num [wf, xSignImag, val0, iso0, eps0])

/-- …and entered correctly as a magnitude (500 K): G-VALID fails. -/
def xImag : Inputs ℝ := ⟨{ val0 with imagAbsMaxK := 500 }, iso0 500, eps0, none, none, []⟩

theorem ce_imag_fixed : decideT cfgDoc xImag = .undecided :=
  decide_of_not_gvalid _ _ (by norm_num [gvalid, xImag, val0, iso0, eps0, lamHi])

/-- (C3) Anharmonic shift −30 % with H present and α²F not recomputed: WELL-FORMED fails. -/
def valShift (s : ℝ) : Validity ℝ :=
  { val0 with containsHLiBe := true, anharmonicCheck := true, anharmShiftAbs := s,
              alpha2FRecomputedAnharmonic := false }
def xSignShift : Inputs ℝ := ⟨valShift (-3/10), iso0 500, eps0, none, none, []⟩

theorem ce_signShift_fixed : decideT cfgDoc xSignShift = .undecided :=
  decide_of_not_wf _ _ (by norm_num [wf, xSignShift, valShift, val0, iso0, eps0])

/-- …and as a magnitude (0.30): G-VALID fails (no recomputation; also ε_ω < |shift|). -/
def xShift : Inputs ℝ := ⟨valShift (3/10), iso0 500, eps0, none, none, []⟩

theorem ce_shift_fixed : decideT cfgDoc xShift = .undecided :=
  decide_of_not_gvalid _ _ (by norm_num [gvalid, xShift, valShift, val0, iso0, eps0, lamHi])

/-- (C2, RL1 `ce_eps0`) ε_sys = ε_ω = 0: WELL-FORMED fails (ε_sys < 0.10). -/
def xEps0 : Inputs ℝ := ⟨val0, iso0 1250, ⟨0, 0, 0, pf 1⟩, none, none, []⟩

theorem ce_eps0_fixed : decideT cfgDoc xEps0 = .undecided :=
  decide_of_not_wf _ _ (by norm_num [wf, xEps0, val0, iso0, in0half])

/-- ε_sys = 0.10 but ε_ω = 0 < 0.02: G-VALID fails. -/
def xEpsOm0 : Inputs ℝ := ⟨val0, iso0 1250, ⟨1/10, 0, 0, pf 1⟩, none, none, []⟩

theorem ce_epsOmega0_fixed : decideT cfgDoc xEpsOm0 = .undecided :=
  decide_of_not_gvalid _ _ (by norm_num [gvalid, xEpsOm0, val0, iso0, lamHi])

/-- A 1×1 coupling matrix (𝒜 = 500 K) with given provenance. -/
def mat1 (p : Prov) : MatData ℝ :=
  ⟨1, fun _ _ => 500, fun _ => 1, fun _ => 1, 1/20, p⟩

/-- (C4, RL1 `ce_mix`) Two sheets, full-DFPT isotropic data, κ input from an ML surrogate:
no rule may give REJECTED_PROVEN. -/
def xMix : Inputs ℝ :=
  ⟨{ val0 with nSheets := 2 }, iso0 5000, eps0, some (mat1 ⟨.ml, 1⟩), none, []⟩

/-- Same, but the matrix is full-DFPT from a DIFFERENT run (run 2 vs run 1). -/
def xMixRun : Inputs ℝ :=
  ⟨{ val0 with nSheets := 2 }, iso0 5000, eps0, some (mat1 (pf 2)), none, []⟩

theorem not_proven_of_floor_false {x : Inputs ℝ} (r : Rule) (h : floorR cfgDoc x r = false) :
    decideT cfgDoc x ≠ .rejectedProven r := by
  intro hd
  rw [(decide_rejectedProven hd).2.2.2.1] at h
  cases h

theorem xMix_mat_low : ∀ p, provOf xMix .mat = some p → p.tier ≠ .fullDFPT := by
  intro p hp
  simp only [provOf, xMix, mat1, Option.map_some, Option.some.injEq] at hp
  subst hp; simp

theorem ce_mix_fixed (r : Rule) : decideT cfgDoc xMix ≠ .rejectedProven r := by
  cases r
  · exact decide_not_proven_of_input (cfg := cfgDoc) (x := xMix) .kappa .mat
      (by simp [Rule.uses]) xMix_mat_low
  · exact decide_not_proven_of_input (cfg := cfgDoc) (x := xMix) .rowmax .mat
      (by simp [Rule.uses]) xMix_mat_low
  · exact not_proven_of_floor_false _ (by simp [floorR, xMix])
  · exact not_proven_of_floor_false _ (by simp [floorR, xMix, val0])
  · exact not_proven_of_floor_false _ (by simp [floorR, xMix, val0])
  · exact not_proven_of_floor_false _ (by simp [floorR, xMix, val0])

theorem xMixRun_mat : provOf xMixRun .mat = some (pf 2) := rfl
theorem xMixRun_eps : provOf xMixRun .eps = some (pf 1) := rfl

theorem ce_mixRun_fixed (r : Rule) : decideT cfgDoc xMixRun ≠ .rejectedProven r := by
  cases r
  · exact decide_not_proven_of_mixed_runs (cfg := cfgDoc) .kappa .mat .eps
      (by simp [Rule.uses]) (by simp [Rule.uses]) (pf 2) (pf 1) xMixRun_mat xMixRun_eps
      (by simp [pf])
  · exact decide_not_proven_of_mixed_runs (cfg := cfgDoc) .rowmax .mat .eps
      (by simp [Rule.uses]) (by simp [Rule.uses]) (pf 2) (pf 1) xMixRun_mat xMixRun_eps
      (by simp [pf])
  · exact not_proven_of_floor_false _ (by simp [floorR, xMixRun])
  · exact not_proven_of_floor_false _ (by simp [floorR, xMixRun, val0])
  · exact not_proven_of_floor_false _ (by simp [floorR, xMixRun, val0])
  · exact not_proven_of_floor_false _ (by simp [floorR, xMixRun, val0])

/-- Fail-closed tags: an unknown tag on the positive control downgrades it to UNDECIDED. -/
theorem unknown_tag_blocks :
    decideT cfgDoc { x0 with tags := [.unknown "A-FILE-INCOMPLET"] } = .undecided := by
  have h := decide_downgrade cfgDoc x0 [] [.unknown "A-FILE-INCOMPLET"] (by simp)
  have h0 : decideT cfgDoc { x0 with tags := [] } = .rejectedProven .iso := ex_proven
  rcases h with h | h
  · exfalso
    have hb : blocked { x0 with tags := [InputTag.unknown "A-FILE-INCOMPLET"] } = true := by
      simp [blocked, InputTag.blocksRejection]
    rcases decide_of_blocked cfgDoc _ hb with h1 | h1 <;> rw [h, h0] at h1 <;> cases h1
  · exact h

end

end TargetLogic.Regression
