/-
P7: the decision function `classify` (ARCH §4, combined across alternative spectra by S5)
and its soundness.

`classify` consumes only booleans. Each boolean is the result of a rigorous (Arb) check in
`rtt.core.certify`; its semantic meaning is the certificate predicate named in `EvidenceSound`.
-/
import PipelineLogic.Assembly
import PipelineLogic.Rayleigh

open Finset

namespace PipelineLogic

open TcR1

/-- Mirrors `rtt.core.status.Status`. Only the first three are produced by `classify`. -/
inductive Status
  | rejectedProven
  | certifiedCandidate
  | undecided
  | excludedAssumption
  | deferred
  | error
  deriving DecidableEq, Repr

/-- The three rigorous checks for one alternative spectrum (run parameters `T*`, `e`, `μ_hi`,
`a0`, cutoff rule fixed).
* `univFires`: R-UNIV, `A_hi < a0 · T*`.
* `cellsBelowTarget`: R-CELLS succeeded and its bound satisfies `T_up < T*`.
* `rayleighAtOrAboveTarget`: C-RAYLEIGH holds at some `T_low ≥ T*`. -/
structure AltEvidence where
  univFires : Bool
  cellsBelowTarget : Bool
  rayleighAtOrAboveTarget : Bool
  deriving DecidableEq, Repr

/-- ARCH §4 `classify` for one spectrum: R-UNIV, else R-CELLS, else C-RAYLEIGH, else
UNDECIDED. -/
def classifyAlt (ev : AltEvidence) : Status :=
  if ev.univFires then .rejectedProven
  else if ev.cellsBelowTarget then .rejectedProven
  else if ev.rayleighAtOrAboveTarget then .certifiedCandidate
  else .undecided

/-- ARCH §4 + S5: combine the alternatives by conjunction. An empty list of alternatives is
UNDECIDED (a vacuous conjunction must never reject or certify). -/
def classify (evs : List AltEvidence) : Status :=
  if evs.isEmpty then .undecided
  else if evs.all (fun ev => classifyAlt ev == .rejectedProven) then .rejectedProven
  else if evs.all (fun ev => classifyAlt ev == .certifiedCandidate) then .certifiedCandidate
  else .undecided

example : classify [] = .undecided := rfl
example : classify [⟨true, false, false⟩, ⟨false, true, true⟩] = .rejectedProven := rfl
example : classify [⟨false, false, true⟩, ⟨false, false, true⟩] = .certifiedCandidate := rfl
example : classify [⟨true, false, false⟩, ⟨false, false, true⟩] = .undecided := rfl
example : classify [⟨false, false, true⟩, ⟨false, false, false⟩] = .undecided := rfl

/-! ## Combinatorics of `classify` -/

theorem classifyAlt_cases (ev : AltEvidence) :
    classifyAlt ev = .rejectedProven ∨ classifyAlt ev = .certifiedCandidate ∨
      classifyAlt ev = .undecided := by
  unfold classifyAlt; split_ifs <;> simp

/-- **P7 (exhaustive).** `classify` returns REJECTED_PROVEN, CERTIFIED_CANDIDATE or
UNDECIDED (never EXCLUDED_ASSUMPTION, DEFERRED or ERROR). -/
theorem classify_cases (evs : List AltEvidence) :
    classify evs = .rejectedProven ∨ classify evs = .certifiedCandidate ∨
      classify evs = .undecided := by
  unfold classify; split_ifs <;> simp

/-- **P7 (exclusive).** The statuses are pairwise distinct, so `classify` (a function) yields
exactly one of the three. -/
theorem status_exclusive :
    Status.rejectedProven ≠ Status.certifiedCandidate ∧
    Status.rejectedProven ≠ Status.undecided ∧
    Status.certifiedCandidate ≠ Status.undecided := by
  decide

theorem classify_eq_rejected_iff (evs : List AltEvidence) :
    classify evs = .rejectedProven ↔
      evs ≠ [] ∧ ∀ ev ∈ evs, classifyAlt ev = .rejectedProven := by
  unfold classify
  split_ifs with h1 h2 h3
  · simp only [List.isEmpty_iff] at h1
    simp [h1]
  · simp only [List.isEmpty_iff] at h1
    simp only [List.all_eq_true, beq_iff_eq] at h2
    simp only [true_iff]; exact ⟨h1, h2⟩
  · simp only [List.all_eq_true, beq_iff_eq] at h2
    simp only [false_iff, not_and]
    intro _; exact h2
  · simp only [List.all_eq_true, beq_iff_eq] at h2
    simp only [false_iff, not_and]
    intro _; exact h2

theorem classify_eq_certified_iff (evs : List AltEvidence) :
    classify evs = .certifiedCandidate ↔
      evs ≠ [] ∧ ∀ ev ∈ evs, classifyAlt ev = .certifiedCandidate := by
  unfold classify
  split_ifs with h1 h2 h3
  · simp only [List.isEmpty_iff] at h1
    simp [h1]
  · simp only [List.isEmpty_iff] at h1
    simp only [List.all_eq_true, beq_iff_eq] at h2
    simp only [false_iff, not_and]
    intro _ hall
    obtain ⟨ev, hev⟩ := List.exists_mem_of_ne_nil evs h1
    have := (h2 ev hev).symm.trans (hall ev hev)
    exact absurd this (by decide)
  · simp only [List.isEmpty_iff] at h1
    simp only [List.all_eq_true, beq_iff_eq] at h3
    simp only [true_iff]; exact ⟨h1, h3⟩
  · simp only [List.all_eq_true, beq_iff_eq] at h3
    simp only [false_iff, not_and]
    intro _; exact h3

/-- With a single alternative, `classify` is `classifyAlt`. -/
theorem classify_singleton (ev : AltEvidence) : classify [ev] = classifyAlt ev := by
  obtain ⟨a, b, c⟩ := ev
  cases a <;> cases b <;> cases c <;> rfl

/-! ## Semantics of the booleans (certificate predicates) -/

section Semantics

variable {ι : Type*}

/-- R-UNIV certificate: an upper bound `A_hi ≥ (1+e) A` with `A_hi < a0 T*`. -/
def UnivCert (e a0 Tstar : ℝ) (sp : Spectrum ι) : Prop :=
  ∃ Ahi : ℝ, (1 + e) * sp.A ≤ Ahi ∧ Ahi < a0 * Tstar

/-- R-CELLS certificate with `T_up < T*`: cells tiling `[T_up, t_end]` exactly, each certified
(`CellOK`), and `(1+e) A ≤ A_hi ≤ a0 · t_end` (so R-UNIV covers `T ≥ t_end`). -/
def CellsCert (κ e a0 Tstar : ℝ) (sp : Spectrum ι) : Prop :=
  ∃ Tup tend Ahi : ℝ, ∃ cells : List (ℝ × ℝ), Tup < Tstar ∧ Tiles Tup tend cells ∧
    (∀ c ∈ cells, CellOK κ sp e c.1 c.2) ∧ (1 + e) * sp.A ≤ Ahi ∧ Ahi ≤ a0 * tend

/-- C-RAYLEIGH certificate at some `T_low ≥ T*`: a nonzero finitely supported `x` whose
envelope box minimum of `Q⁺_{T_low}(x)` at `μ* = μ_hi`, cutoff `Ncf T_low`, is `≥ 0`. -/
noncomputable def RayleighCert (κ e μhi Tstar : ℝ) (Ncf : ℝ → WithTop ℕ) (sp : Spectrum ι) :
    Prop :=
  ∃ T : ℝ, Tstar ≤ T ∧ ∃ (N : ℕ) (x : ℕ → ℝ), (∃ n ∈ range N, x n ≠ 0) ∧
    0 ≤ boxMinQ κ sp e μhi (Ncf T) T N x

/-- Each `true` boolean is backed by its certificate (this is what the Arb checks establish). -/
def EvidenceSound (κ a0 e μhi Tstar : ℝ) (Ncf : ℝ → WithTop ℕ) (sp : Spectrum ι)
    (ev : AltEvidence) : Prop :=
  (ev.univFires = true → UnivCert e a0 Tstar sp) ∧
  (ev.cellsBelowTarget = true → CellsCert κ e a0 Tstar sp) ∧
  (ev.rayleighAtOrAboveTarget = true → RayleighCert κ e μhi Tstar Ncf sp)

/-! ## Soundness for one alternative -/

/-- One alternative classified REJECTED_PROVEN: there is `T_b < T*` such that every spectrum
of the envelope is stable at every `T ≥ T_b`, for every `μ* ≥ 0` and every cutoff
(so `Tc_model ≤ T_b < T*`). -/
theorem classifyAlt_rejected_sound (κ a0 η : ℝ) (hH : H κ a0 η) (ha0 : 0 < a0) (hη : 0 < η)
    (e μhi Tstar : ℝ) (Ncf : ℝ → WithTop ℕ) (sp : Spectrum ι) (ev : AltEvidence)
    (hs : EvidenceSound κ a0 e μhi Tstar Ncf sp ev) (h : classifyAlt ev = .rejectedProven) :
    ∃ Tb : ℝ, Tb < Tstar ∧ ∀ sp' : Spectrum ι, InEnv e sp sp' → ∀ μ : ℝ, 0 ≤ μ →
      ∀ Nc : WithTop ℕ, ∀ T : ℝ, Tb ≤ T → Stable κ sp' μ Nc T := by
  obtain ⟨hu, hc, _⟩ := hs
  unfold classifyAlt at h
  split_ifs at h with h1 h2 h3
  · obtain ⟨Ahi, hAhi, hlt⟩ := hu h1
    exact ⟨Ahi / a0, R_UNIV κ a0 η hH ha0 hη e sp Ahi Tstar hAhi hlt⟩
  · obtain ⟨Tup, tend, Ahi, cells, hbelow, htile, hcells, hAhi, hend⟩ := hc h2
    exact ⟨Tup, hbelow, fun sp' hin μ hμ Nc T hT =>
      cells_stable κ a0 η hH ha0 hη e sp Tup tend Ahi cells htile hcells hAhi hend sp' hin μ hμ
        Nc T hT⟩

/-- One alternative classified CERTIFIED_CANDIDATE: there is `T_low ≥ T*` at which, for every
spectrum of the envelope and every `μ* ≤ μ_hi`, the even sector (cutoff `Ncf T_low`) is
unstable; hence `T_low` is not stable and `Tc_model ≥ T_low ≥ T*`. -/
theorem classifyAlt_certified_sound (κ a0 e μhi Tstar : ℝ) (Ncf : ℝ → WithTop ℕ)
    (sp : Spectrum ι) (ev : AltEvidence) (hs : EvidenceSound κ a0 e μhi Tstar Ncf sp ev)
    (h : classifyAlt ev = .certifiedCandidate) :
    ∃ T : ℝ, Tstar ≤ T ∧ ∀ sp' : Spectrum ι, InEnv e sp sp' → ∀ μ : ℝ, μ ≤ μhi →
      Unstable (Qp κ sp' μ (Ncf T) T) ∧ ¬ Stable κ sp' μ (Ncf T) T := by
  obtain ⟨_, _, hr⟩ := hs
  unfold classifyAlt at h
  split_ifs at h with h1 h2 h3
  obtain ⟨T, hT, N, x, hx, h0⟩ := hr h3
  refine ⟨T, hT, fun sp' hin μ hμ => ?_⟩
  have hU := rayleigh_cert κ e μhi sp (Ncf T) T N x hx h0 sp' hin μ hμ
  exact ⟨hU, not_stable_of_unstable κ sp' μ (Ncf T) T hU⟩

/-- **Consistency.** With `e ≥ 0` and `μ_hi ≥ 0`, sound evidence can never have both a
rejection rule and C-RAYLEIGH true. So the priority order inside `classifyAlt` is immaterial. -/
theorem evidence_consistent (κ a0 η : ℝ) (hH : H κ a0 η) (ha0 : 0 < a0) (hη : 0 < η)
    (e μhi Tstar : ℝ) (he : 0 ≤ e) (hμhi : 0 ≤ μhi) (Ncf : ℝ → WithTop ℕ) (sp : Spectrum ι)
    (ev : AltEvidence) (hs : EvidenceSound κ a0 e μhi Tstar Ncf sp ev) :
    ¬ ((ev.univFires = true ∨ ev.cellsBelowTarget = true) ∧
        ev.rayleighAtOrAboveTarget = true) := by
  rintro ⟨hrej, hray⟩
  have hcl : classifyAlt ev = .rejectedProven := by
    unfold classifyAlt
    rcases hrej with h | h <;> simp [h]
  obtain ⟨Tb, hTb, hst⟩ := classifyAlt_rejected_sound κ a0 η hH ha0 hη e μhi Tstar Ncf sp ev hs hcl
  obtain ⟨T, hT, N, x, hx, h0⟩ := hs.2.2 hray
  have hU := rayleigh_cert κ e μhi sp (Ncf T) T N x hx h0 sp (InEnv_refl e he sp) μhi le_rfl
  exact (hst sp (InEnv_refl e he sp) μhi hμhi (Ncf T) T (by linarith)).1 hU

/-! ## Soundness of `classify` (S5 conjunction over alternatives) -/

/-- **P7 (REJECTED_PROVEN soundness).** Under `H(a0, η)` (`a0, η > 0`), if every alternative's
booleans are backed by their certificates and `classify` returns REJECTED_PROVEN, then the
list of alternatives is nonempty and for every alternative spectrum there is `T_b < T*` such
that every spectrum in its envelope is stable at every `T ≥ T_b`, for every `μ* ≥ 0` and
every cutoff. -/
theorem classify_rejected_sound (κ a0 η : ℝ) (hH : H κ a0 η) (ha0 : 0 < a0) (hη : 0 < η)
    (e μhi Tstar : ℝ) (Ncf : ℝ → WithTop ℕ) (alts : List (Spectrum ι × AltEvidence))
    (hs : ∀ a ∈ alts, EvidenceSound κ a0 e μhi Tstar Ncf a.1 a.2)
    (h : classify (alts.map Prod.snd) = .rejectedProven) :
    alts ≠ [] ∧ ∀ a ∈ alts, ∃ Tb : ℝ, Tb < Tstar ∧ ∀ sp' : Spectrum ι, InEnv e a.1 sp' →
      ∀ μ : ℝ, 0 ≤ μ → ∀ Nc : WithTop ℕ, ∀ T : ℝ, Tb ≤ T → Stable κ sp' μ Nc T := by
  rw [classify_eq_rejected_iff] at h
  obtain ⟨hne, hall⟩ := h
  refine ⟨fun h0 => hne (by simp [h0]), fun a ha => ?_⟩
  exact classifyAlt_rejected_sound κ a0 η hH ha0 hη e μhi Tstar Ncf a.1 a.2 (hs a ha)
    (hall a.2 (List.mem_map_of_mem ha))

/-- **P7 (REJECTED_PROVEN, as stated in ARCH S5).** … in particular: for every alternative,
every spectrum in its envelope, every `μ* ≥ 0`, every cutoff and every `T ≥ T*`: Stable(T). -/
theorem classify_rejected_sound' (κ a0 η : ℝ) (hH : H κ a0 η) (ha0 : 0 < a0) (hη : 0 < η)
    (e μhi Tstar : ℝ) (Ncf : ℝ → WithTop ℕ) (alts : List (Spectrum ι × AltEvidence))
    (hs : ∀ a ∈ alts, EvidenceSound κ a0 e μhi Tstar Ncf a.1 a.2)
    (h : classify (alts.map Prod.snd) = .rejectedProven) :
    ∀ a ∈ alts, ∀ sp' : Spectrum ι, InEnv e a.1 sp' → ∀ μ : ℝ, 0 ≤ μ →
      ∀ Nc : WithTop ℕ, ∀ T : ℝ, Tstar ≤ T → Stable κ sp' μ Nc T := by
  intro a ha sp' hin μ hμ Nc T hT
  obtain ⟨Tb, hTb, hst⟩ :=
    (classify_rejected_sound κ a0 η hH ha0 hη e μhi Tstar Ncf alts hs h).2 a ha
  exact hst sp' hin μ hμ Nc T (by linarith)

/-- **P7 (CERTIFIED_CANDIDATE soundness).** If every alternative's booleans are backed by their
certificates and `classify` returns CERTIFIED_CANDIDATE, then the list is nonempty and for
every alternative there is `T_low ≥ T*` such that for every spectrum in its envelope and every
`μ* ≤ μ_hi` the even sector at `T_low` (cutoff `Ncf T_low`) is unstable, so `T_low` is not
stable. (No hypothesis `H` is needed.) -/
theorem classify_certified_sound (κ a0 e μhi Tstar : ℝ) (Ncf : ℝ → WithTop ℕ)
    (alts : List (Spectrum ι × AltEvidence))
    (hs : ∀ a ∈ alts, EvidenceSound κ a0 e μhi Tstar Ncf a.1 a.2)
    (h : classify (alts.map Prod.snd) = .certifiedCandidate) :
    alts ≠ [] ∧ ∀ a ∈ alts, ∃ T : ℝ, Tstar ≤ T ∧ ∀ sp' : Spectrum ι, InEnv e a.1 sp' →
      ∀ μ : ℝ, μ ≤ μhi →
        Unstable (Qp κ sp' μ (Ncf T) T) ∧ ¬ Stable κ sp' μ (Ncf T) T := by
  rw [classify_eq_certified_iff] at h
  obtain ⟨hne, hall⟩ := h
  refine ⟨fun h0 => hne (by simp [h0]), fun a ha => ?_⟩
  exact classifyAlt_certified_sound κ a0 e μhi Tstar Ncf a.1 a.2 (hs a ha)
    (hall a.2 (List.mem_map_of_mem ha))

end Semantics

end PipelineLogic
