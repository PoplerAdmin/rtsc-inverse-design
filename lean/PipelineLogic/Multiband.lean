/-
PipelineLogic.Multiband: Lean mirror of `rtt.core.certify.classify.combine_multiband`
(ARCH §4c, `rejection_model = "isotropic" | "multiband"`) and its soundness.

Python (verbatim control flow):

    def combine_multiband(status_sb, mb_fires, rejection_model, n_spectra=None):
        if rejection_model == "isotropic":
            return status_sb
        if n_spectra is not None and n_spectra == 0:
            return Status.UNDECIDED
        mb_rej = len(mb_fires) > 0 and all(mb_fires)
        if mb_rej and status_sb != Status.CERTIFIED_CANDIDATE:
            return Status.REJECTED_PROVEN
        if status_sb == Status.CERTIFIED_CANDIDATE and not mb_rej:
            return Status.CERTIFIED_CANDIDATE
        return Status.UNDECIDED

and `classify` returns
`combine_multiband(classify_evidence(evs), mb_fires, model, n_spectra=len(spectra))`, with one
`AltEvidence` per spectrum (so `n_spectra = len(evs)`). `n_spectra : int | None` is
`nSpectra : Option ℕ` (`None` ↦ `none`). The assumption `BAND-DATA-NOT-PER-ALTERNATIVE`
(added when the numbers of band couplings and spectra differ) is an annotation only; it does
not change the decision.
`rejection_model` is a string in Python: `"isotropic"` ↦ `.isotropic`, every other string
(only `"multiband"` passes `classify`'s validation) ↦ `.multiband`.
-/
import PipelineLogic.Classify
import MultiBand

open Finset

namespace PipelineLogic

open TcR1 MultiBand

/-- `CertifyParams.rejection_model`. -/
inductive RejectionModel
  | isotropic
  | multiband
  deriving DecidableEq, Repr

/-- Mirror of `combine_multiband`. `sb` is the single-band status (`classify evs`),
`mbFires` the R-MULTI booleans, one per band-coupling alternative, `nSpectra` the number of
single-band spectra (`n_spectra`, default `None`). -/
def combineMultiband (sb : Status) (mbFires : List Bool) (model : RejectionModel)
    (nSpectra : Option ℕ := none) : Status :=
  match model with
  | .isotropic => sb
  | .multiband =>
    if nSpectra == some 0 then .undecided
    else
    let mbRej : Bool := !mbFires.isEmpty && mbFires.all id
    if mbRej && sb != .certifiedCandidate then .rejectedProven
    else if sb == .certifiedCandidate && !mbRej then .certifiedCandidate
    else .undecided

/-- The whole decision of `rtt.core.certify.classify`: single-band `classify` (S5 over the
alternative spectra), then `combine_multiband`. -/
def classifyFull (evs : List AltEvidence) (mbFires : List Bool) (model : RejectionModel) :
    Status :=
  combineMultiband (classify evs) mbFires model (some evs.length)

example : combineMultiband .rejectedProven [] .multiband = .undecided := rfl
example : combineMultiband .undecided [true, true] .multiband = .rejectedProven := rfl
example : combineMultiband .certifiedCandidate [true] .multiband = .undecided := rfl
example : combineMultiband .certifiedCandidate [true, false] .multiband = .certifiedCandidate :=
  rfl
example : combineMultiband .rejectedProven [] .isotropic = .rejectedProven := rfl
example : combineMultiband .undecided [true] .multiband (some 0) = .undecided := rfl
example : combineMultiband .undecided [true] .multiband (some 2) = .rejectedProven := rfl
example : combineMultiband .rejectedProven [true] .isotropic (some 0) = .rejectedProven := rfl

/-! ## Combinatorics -/

/-- R-MULTI rejection flag: nonempty and every band-coupling alternative fires. -/
def MBRej (mbFires : List Bool) : Prop := mbFires ≠ [] ∧ ∀ f ∈ mbFires, f = true

theorem mbRej_bool (mbFires : List Bool) :
    (!mbFires.isEmpty && mbFires.all id) = true ↔ MBRej mbFires := by
  unfold MBRej
  simp [List.all_eq_true]

theorem combineMultiband_isotropic (sb : Status) (mbFires : List Bool) (n : Option ℕ) :
    combineMultiband sb mbFires .isotropic n = sb := rfl

/-- **Zero spectra (multiband mode).** `n_spectra = 0` gives UNDECIDED. -/
theorem combine_mb_zero_spectra (sb : Status) (mbFires : List Bool) :
    combineMultiband sb mbFires .multiband (some 0) = .undecided := rfl

/-- Multiband mode: REJECTED iff `n_spectra ≠ 0`, R-MULTI fires on a nonempty list of
band-coupling alternatives and the single-band status is not CERTIFIED. -/
theorem combine_mb_rejected_iff (sb : Status) (mbFires : List Bool) (n : Option ℕ) :
    combineMultiband sb mbFires .multiband n = .rejectedProven ↔
      n ≠ some 0 ∧ MBRej mbFires ∧ sb ≠ .certifiedCandidate := by
  rw [← mbRej_bool]
  unfold combineMultiband
  by_cases hn : n = some 0
  · subst hn; simp
  · have : (n == some 0) = false := by simpa using hn
    simp only [this, Bool.false_eq_true, ↓reduceIte, ne_eq, hn, not_false_eq_true, true_and]
    cases h : (!mbFires.isEmpty && mbFires.all id) <;> cases sb <;> simp

/-- Multiband mode: CERTIFIED iff `n_spectra ≠ 0`, the single-band status is CERTIFIED and
R-MULTI does not reject. -/
theorem combine_mb_certified_iff (sb : Status) (mbFires : List Bool) (n : Option ℕ) :
    combineMultiband sb mbFires .multiband n = .certifiedCandidate ↔
      n ≠ some 0 ∧ sb = .certifiedCandidate ∧ ¬ MBRej mbFires := by
  rw [← mbRej_bool]
  unfold combineMultiband
  by_cases hn : n = some 0
  · subst hn; simp
  · have : (n == some 0) = false := by simpa using hn
    simp only [this, Bool.false_eq_true, ↓reduceIte, ne_eq, hn, not_false_eq_true, true_and]
    cases h : (!mbFires.isEmpty && mbFires.all id) <;> cases sb <;> simp

/-- **Exhaustive (multiband mode).** For every input status (even one `classify` never
produces), the multiband result is REJECTED, CERTIFIED or UNDECIDED. -/
theorem combine_mb_cases (sb : Status) (mbFires : List Bool) (n : Option ℕ) :
    combineMultiband sb mbFires .multiband n = .rejectedProven ∨
    combineMultiband sb mbFires .multiband n = .certifiedCandidate ∨
    combineMultiband sb mbFires .multiband n = .undecided := by
  unfold combineMultiband
  dsimp only
  split_ifs <;> simp

/-- **Exhaustive (full decision, both modes).** -/
theorem classifyFull_cases (evs : List AltEvidence) (mbFires : List Bool)
    (model : RejectionModel) :
    classifyFull evs mbFires model = .rejectedProven ∨
    classifyFull evs mbFires model = .certifiedCandidate ∨
    classifyFull evs mbFires model = .undecided := by
  cases model
  · exact classify_cases evs
  · exact combine_mb_cases _ _ _

/-- **Exclusive.** The three statuses are pairwise distinct (so `classifyFull`, a function,
yields exactly one). -/
theorem classifyFull_exclusive (evs : List AltEvidence) (mbFires : List Bool)
    (model : RejectionModel) :
    ¬ (classifyFull evs mbFires model = .rejectedProven ∧
        classifyFull evs mbFires model = .certifiedCandidate) ∧
    ¬ (classifyFull evs mbFires model = .rejectedProven ∧
        classifyFull evs mbFires model = .undecided) ∧
    ¬ (classifyFull evs mbFires model = .certifiedCandidate ∧
        classifyFull evs mbFires model = .undecided) := by
  refine ⟨fun ⟨h1, h2⟩ => ?_, fun ⟨h1, h2⟩ => ?_, fun ⟨h1, h2⟩ => ?_⟩ <;>
    rw [h1] at h2 <;> exact absurd h2 (by decide)

/-- **Empty band-coupling list.** In multiband mode an empty list of band-coupling
alternatives never yields REJECTED (whatever the single-band status and `n_spectra`). -/
theorem combine_mb_empty_not_rejected (sb : Status) (n : Option ℕ) :
    combineMultiband sb [] .multiband n ≠ .rejectedProven := by
  rw [Ne, combine_mb_rejected_iff]
  rintro ⟨_, ⟨h, _⟩, _⟩
  exact h rfl

theorem classifyFull_empty_not_rejected (evs : List AltEvidence) :
    classifyFull evs [] .multiband ≠ .rejectedProven :=
  combine_mb_empty_not_rejected _ _

/-- **Zero spectra.** With no single-band spectrum the full decision is UNDECIDED in both
modes (multiband: the `n_spectra = 0` guard; isotropic: `classify [] = undecided`), whatever
the band-coupling evidence; in particular never REJECTED. -/
theorem classifyFull_zero_spectra (mbFires : List Bool) (model : RejectionModel) :
    classifyFull [] mbFires model = .undecided := by
  cases model <;> rfl

theorem classifyFull_zero_spectra_not_rejected (mbFires : List Bool) (model : RejectionModel) :
    classifyFull [] mbFires model ≠ .rejectedProven := by
  rw [classifyFull_zero_spectra]; decide

/-- In multiband mode a single-band rejection alone never yields REJECTED. -/
theorem combine_mb_rejected_needs_mb (sb : Status) (mbFires : List Bool) (n : Option ℕ)
    (h : combineMultiband sb mbFires .multiband n = .rejectedProven) : MBRej mbFires :=
  ((combine_mb_rejected_iff sb mbFires n).mp h).2.1

/-- CERTIFIED (either mode) only passes through a single-band CERTIFIED. -/
theorem combine_certified_imp (sb : Status) (mbFires : List Bool) (model : RejectionModel)
    (n : Option ℕ) (h : combineMultiband sb mbFires model n = .certifiedCandidate) :
    sb = .certifiedCandidate := by
  cases model
  · exact h
  · exact ((combine_mb_certified_iff sb mbFires n).mp h).2.1

/-! ## Semantics of an R-MULTI boolean -/

/-- One band-coupling alternative: the (nominal) multiband model it stands for
(`MSpec` with total weights `W^tot_ij = N_i 𝒜_ij`, DOS `N_i`), its relative envelope `e`, and
the R-MULTI boolean computed by `rmulti_evidence`. -/
structure MBAlt (α : Type*) where
  B : ℕ
  sp : MSpec B α
  dos : Fin B → ℝ
  e : ℝ
  fires : Bool

/-- `sp'` is in the multiband envelope of the alternative: `𝒜'_ij ≤ (1+e) 𝒜_ij` entrywise
(same DOS; any atoms, any frequencies). -/
def MBEnv {α β : Type*} (a : MBAlt α) (sp' : MSpec a.B β) : Prop :=
  ∀ i j, Amat sp' a.dos i j ≤ (1 + a.e) * Amat a.sp a.dos i j

/-- Collatz–Wielandt core of every route: `v > 0`, `𝒜 v ≤ κ v`, `(1+e) κ ≤ κ_hi`. -/
def CWCore {B : ℕ} {α : Type*} (sp : MSpec B α) (dos : Fin B → ℝ) (e κhi : ℝ) : Prop :=
  ∃ (v : Fin B → ℝ) (κ : ℝ), (∀ i, 0 < v i) ∧ (∀ i, ∑ j, Amat sp dos i j * v j ≤ κ * v i) ∧
    (1 + e) * κ ≤ κhi

/-- Semantic content of `fires = True` (`rmulti_evidence`): `T* > 0`, `0 ≤ e` (checked
`0 ≤ e < 1`), `N_i > 0` (model), and a route-certified `κ_hi` with `κ_hi < a0 T*`. -/
def RMultiCert {α : Type*} (a0 Tstar : ℝ) (a : MBAlt α) : Prop :=
  0 < Tstar ∧ 0 ≤ a.e ∧ (∀ i, 0 < a.dos i) ∧
    ∃ κhi : ℝ, κhi < a0 * Tstar ∧ CWCore a.sp a.dos a.e κhi

section Routes

variable {B : ℕ} {α : Type*}

/-- Route `"cw"`: `(1+e)(𝒜 v)_i ≤ κ_hi v_i` with `v > 0` (exact rationals). -/
theorem route_cw (sp : MSpec B α) (dos : Fin B → ℝ) (e κhi : ℝ) (he : 0 < 1 + e)
    (v : Fin B → ℝ) (hv : ∀ i, 0 < v i)
    (h : ∀ i, (1 + e) * ∑ j, Amat sp dos i j * v j ≤ κhi * v i) : CWCore sp dos e κhi := by
  refine ⟨v, κhi / (1 + e), hv, fun i => ?_, by rw [mul_div_cancel₀ _ he.ne']⟩
  rw [div_mul_eq_mul_div, le_div_iff₀ he, mul_comm]
  exact h i

/-- Route `"rowsum"` (source `A_K`): `(1+e) ∑_j 𝒜_ij ≤ κ_hi` for every `i`. -/
theorem route_rowsum (sp : MSpec B α) (dos : Fin B → ℝ) (e κhi : ℝ) (he : 0 < 1 + e)
    (h : ∀ i, (1 + e) * ∑ j, Amat sp dos i j ≤ κhi) : CWCore sp dos e κhi :=
  route_cw sp dos e κhi he (fun _ => 1) (fun _ => one_pos) (fun i => by simpa using h i)

/-- Route `"rowsum"` (source `row_sum_upper_K`): declared bounds `∑_j 𝒜_ij ≤ R_i` (a data
assumption, `DATA BOUND`) and `(1+e) R_i ≤ κ_hi`. -/
theorem route_rowsum_declared (sp : MSpec B α) (dos : Fin B → ℝ) (e κhi : ℝ)
    (he : 0 < 1 + e) (R : Fin B → ℝ) (hR : ∀ i, ∑ j, Amat sp dos i j ≤ R i)
    (h : ∀ i, (1 + e) * R i ≤ κhi) : CWCore sp dos e κhi :=
  route_rowsum sp dos e κhi he (fun i => le_trans (mul_le_mul_of_nonneg_left (hR i) he.le) (h i))

/-- Route `"eig"`: `r_i = √N_i` (`r_i > 0`, `r_i² = N_i`), `u > 0` and
`(1+e) ∑_j W^tot_ij/(r_i r_j) u_j ≤ κ_hi u_i` (Arb), where `W^tot` is the model's total
weight matrix (the symmetrized `W̄` when the data were not exactly symmetric). -/
theorem route_eig (sp : MSpec B α) (dos r u : Fin B → ℝ) (hr : ∀ i, 0 < r i)
    (hr2 : ∀ i, r i ^ 2 = dos i) (hu : ∀ i, 0 < u i) (e κhi : ℝ) (he : 0 < 1 + e)
    (h : ∀ i, (1 + e) * ∑ j, (∑ s ∈ sp.S, sp.W s i j) / (r i * r j) * u j ≤ κhi * u i) :
    CWCore sp dos e κhi := by
  refine ⟨fun i => u i / r i, κhi / (1 + e), fun i => div_pos (hu i) (hr i),
    sub_of_K sp dos r u hr hr2 (κhi / (1 + e)) (fun i => ?_), by rw [mul_div_cancel₀ _ he.ne']⟩
  rw [div_mul_eq_mul_div, le_div_iff₀ he, mul_comm]
  exact h i

end Routes

/-! ## Soundness -/

/-- **Multiband REJECTED soundness.** Hypotheses: `H(a0, η)` (`η ≥ 0`, `a0 > 0`), Prop. N, and
every `fires = True` backed by `RMultiCert`. If the full decision in multiband mode is
REJECTED_PROVEN, then the list of single-band spectra and the band-coupling list are nonempty
and for every band-coupling alternative there is `κ_hi` with `κ_hi / a0 < T*` such that every
envelope member `sp'` (`𝒜' ≤ (1+e)𝒜`, any atoms/frequencies), every PSD Coulomb matrix
`N_i μ*_ij`, every cutoff and every `T > 0` with `T > κ_hi / a0` is stable in the multiband
model; in particular every `T ≥ T*` (see `classifyFull_mb_rejected_sound'`). -/
theorem classifyFull_mb_rejected_sound {α β : Type*} (tp a0 η : ℝ) (hH : H tp a0 η)
    (hη : 0 ≤ η) (ha0 : 0 < a0) (hN : PropN tp a0) (Tstar : ℝ) (evs : List AltEvidence)
    (alts : List (MBAlt α)) (hs : ∀ a ∈ alts, a.fires = true → RMultiCert a0 Tstar a)
    (h : classifyFull evs (alts.map MBAlt.fires) .multiband = .rejectedProven) :
    evs ≠ [] ∧ alts ≠ [] ∧ ∀ a ∈ alts, ∃ κhi : ℝ, κhi / a0 < Tstar ∧
      ∀ sp' : MSpec a.B β, MBEnv a sp' → ∀ μ : Fin a.B → Fin a.B → ℝ, CoulPSD a.dos μ →
        ∀ Nc : WithTop ℕ, ∀ T : ℝ, 0 < T → κhi / a0 < T → MStable tp sp' a.dos μ Nc T := by
  obtain ⟨hn0, ⟨hne, hall⟩, _⟩ := (combine_mb_rejected_iff _ _ _).mp h
  refine ⟨fun h0 => hn0 (by simp [h0]), fun h0 => hne (by simp [h0]), fun a ha => ?_⟩
  obtain ⟨hT0, he, hdos, κhi, hlt, v, κ, hv, hsub, hκ⟩ :=
    hs a ha (hall a.fires (List.mem_map_of_mem ha))
  refine ⟨κhi, by rw [div_lt_iff₀ ha0]; linarith, fun sp' henv μ hμ Nc T hT hTb => ?_⟩
  apply mstable_of_gt tp a0 η hH hη ha0 hN sp' a.dos hdos v hv κhi
    (sub_env a.sp sp' a.dos v hv κ a.e κhi (by linarith) hsub henv hκ) μ hμ Nc T hT
  rw [div_lt_iff₀ ha0] at hTb; linarith

/-- **Multiband REJECTED soundness, ARCH wording.** … for every band-coupling alternative,
every envelope member, every PSD Coulomb matrix, every cutoff and every `T ≥ T*`: stable. -/
theorem classifyFull_mb_rejected_sound' {α β : Type*} (tp a0 η : ℝ) (hH : H tp a0 η)
    (hη : 0 ≤ η) (ha0 : 0 < a0) (hN : PropN tp a0) (Tstar : ℝ) (evs : List AltEvidence)
    (alts : List (MBAlt α)) (hs : ∀ a ∈ alts, a.fires = true → RMultiCert a0 Tstar a)
    (h : classifyFull evs (alts.map MBAlt.fires) .multiband = .rejectedProven) :
    ∀ a ∈ alts, ∀ sp' : MSpec a.B β, MBEnv a sp' → ∀ μ : Fin a.B → Fin a.B → ℝ,
      CoulPSD a.dos μ → ∀ Nc : WithTop ℕ, ∀ T : ℝ, Tstar ≤ T → MStable tp sp' a.dos μ Nc T := by
  intro a ha sp' henv μ hμ Nc T hT
  obtain ⟨_, _, hall⟩ := classifyFull_mb_rejected_sound (β := β) tp a0 η hH hη ha0 hN Tstar evs
    alts hs h
  obtain ⟨κhi, hκ, hst⟩ := hall a ha
  have hT0 : 0 < Tstar := (hs a ha (((combine_mb_rejected_iff _ _ _).mp h).2.1.2 a.fires
    (List.mem_map_of_mem ha))).1
  exact hst sp' henv μ hμ Nc T (by linarith) (by linarith)

/-- **Isotropic-mode REJECTED** is exactly the single-band REJECTED of `classify` (Lean
`classify_rejected_sound`): a statement about the single-band (DOS-averaged) model only
(Python: assumption `M-ISO-SINGLEBAND`, flag `MULTIBAND_RISK`). -/
theorem classifyFull_iso_rejected_sound {ι : Type*} (κ a0 η : ℝ) (hH : H κ a0 η)
    (ha0 : 0 < a0) (hη : 0 < η) (e μhi Tstar : ℝ) (Ncf : ℝ → WithTop ℕ)
    (alts : List (Spectrum ι × AltEvidence))
    (hs : ∀ a ∈ alts, EvidenceSound κ a0 e μhi Tstar Ncf a.1 a.2) (mbFires : List Bool)
    (h : classifyFull (alts.map Prod.snd) mbFires .isotropic = .rejectedProven) :
    alts ≠ [] ∧ ∀ a ∈ alts, ∃ Tb : ℝ, Tb < Tstar ∧ ∀ sp' : Spectrum ι, InEnv e a.1 sp' →
      ∀ μ : ℝ, 0 ≤ μ → ∀ Nc : WithTop ℕ, ∀ T : ℝ, Tb ≤ T → Stable κ sp' μ Nc T :=
  classify_rejected_sound κ a0 η hH ha0 hη e μhi Tstar Ncf alts hs h

/-- **CERTIFIED (either mode), exactly what is proved.** A CERTIFIED result of the full
decision implies the single-band CERTIFIED, hence (`classify_certified_sound`): the list of
single-band alternatives is nonempty and for every alternative spectrum there is `T_low ≥ T*`
at which every spectrum of its envelope, for every `μ* ≤ μ_hi` and the run's cutoff
`Ncf T_low`, is unstable in the **single-band (DOS-averaged, isotropic) model**. Nothing is
claimed about the multiband model (Python assumption `M-ISO-SINGLEBAND`); see
`certified_transfers_to_multiband` for the extra hypotheses under which it transfers. -/
theorem classifyFull_certified_sound {ι : Type*} (κ a0 e μhi Tstar : ℝ)
    (Ncf : ℝ → WithTop ℕ) (alts : List (Spectrum ι × AltEvidence))
    (hs : ∀ a ∈ alts, EvidenceSound κ a0 e μhi Tstar Ncf a.1 a.2) (mbFires : List Bool)
    (model : RejectionModel)
    (h : classifyFull (alts.map Prod.snd) mbFires model = .certifiedCandidate) :
    alts ≠ [] ∧ ∀ a ∈ alts, ∃ T : ℝ, Tstar ≤ T ∧ ∀ sp' : Spectrum ι, InEnv e a.1 sp' →
      ∀ μ : ℝ, μ ≤ μhi →
        Unstable (Qp κ sp' μ (Ncf T) T) ∧ ¬ Stable κ sp' μ (Ncf T) T :=
  classify_certified_sound κ a0 e μhi Tstar Ncf alts hs (combine_certified_imp _ _ _ _ h)

/-! ## Optional: when a single-band instability transfers to the multiband model

Constant-gap trial vector `φ_i = y` for all bands. This is NOT used by the pipeline (it does
not check the hypotheses); it states precisely what would lift `M-ISO-SINGLEBAND`. -/

section Transfer

variable {B : ℕ} {α : Type*}

/-- `N_tot = ∑_i N_i`. -/
def Ntot (dos : Fin B → ℝ) : ℝ := ∑ i, dos i

/-- The DOS-averaged isotropic spectrum of a multiband spectrum:
`α²F_iso = ∑_ij W_ij / N_tot` (same atoms and frequencies). -/
noncomputable def isoAvg (sp : MSpec B α) (dos : Fin B → ℝ) (hN : 0 < Ntot dos) : Spectrum α :=
  { S := sp.S
    p := fun s => (∑ i, ∑ j, sp.W s i j) / Ntot dos
    Ω := sp.Ω
    p_nonneg := fun s hs => div_nonneg
      (sum_nonneg (fun i _ => sum_nonneg (fun j _ => sp.W_nonneg s hs i j))) hN.le
    Ω_pos := sp.Ω_pos }

/-- DOS-averaged Coulomb pseudopotential `μ*_iso = ∑_ij N_i μ*_ij / N_tot`. -/
noncomputable def muAvg (dos : Fin B → ℝ) (μ : Fin B → Fin B → ℝ) : ℝ :=
  (∑ i, ∑ j, dos i * μ i j) / Ntot dos

/-- Constant-gap identity: `Q⁺_multi(φ_i ≡ y) = N_tot · Q⁺_iso(y)` with the averaged spectrum
and averaged `μ*` (`T > 0`). -/
theorem Qplus_const (tp : ℝ) (sp : MSpec B α) (dos : Fin B → ℝ) (hN : 0 < Ntot dos)
    (μ : Fin B → Fin B → ℝ) (Nc : WithTop ℕ) (T : ℝ) (hT : 0 < T) (N : ℕ) (y : ℕ → ℝ) :
    Qplus tp sp dos μ Nc T N (fun _ => y)
      = Ntot dos * Qp tp (isoAvg sp dos hN) (muAvg dos μ) Nc T N y := by
  have hNne : Ntot dos ≠ 0 := hN.ne'
  rw [L1]
  -- per-atom coefficient: Fp(g_T(·, Ω)) = (1/T) Fp(l_{Ω/T})
  have hatom : ∀ s ∈ sp.S, Fp (fun k => g tp T k (sp.Ω s)) N y
      = (1 / T) * Fp (lw tp (sp.Ω s / T)) N y := by
    intro s hs
    rw [← Fp_smul]
    congr 1; funext k
    have := L4_pointwise tp 1 T (sp.Ω s) k hT (sp.Ω_pos s hs)
    simp only [one_mul] at this
    rw [this]
  have hpair : ∀ s ∈ sp.S, ∑ i, ∑ j, sp.W s i j *
      (Bp tp (sp.Ω s / T) N y y - Dq tp (sp.Ω s / T) N y y)
      = (∑ i, ∑ j, sp.W s i j) * Fp (lw tp (sp.Ω s / T)) N y := by
    intro s _
    have hG : Bp tp (sp.Ω s / T) N y y - Dq tp (sp.Ω s / T) N y y
        = Fp (lw tp (sp.Ω s / T)) N y := Gp_self tp (sp.Ω s / T) N y
    rw [hG, sum_mul]
    apply sum_congr rfl; intro i _; rw [sum_mul]
  have hcoul : coul dos μ Nc N (fun _ => y)
      = Ntot dos * muTerm (muAvg dos μ) Nc N y := by
    unfold coul muTerm muAvg Sc
    have e : ∑ i, ∑ j, dos i * μ i j * (∑ n ∈ (range N).filter
        (fun n : ℕ => ((n : ℕ) : WithTop ℕ) < Nc), y n) * (∑ n ∈ (range N).filter
        (fun n : ℕ => ((n : ℕ) : WithTop ℕ) < Nc), y n)
        = (∑ i, ∑ j, dos i * μ i j) * (∑ n ∈ (range N).filter
        (fun n : ℕ => ((n : ℕ) : WithTop ℕ) < Nc), y n) ^ 2 := by
      rw [sum_mul]; apply sum_congr rfl; intro i _
      rw [sum_mul]; apply sum_congr rfl; intro j _; ring
    rw [e]
    field_simp
  have hC : ∑ i, dos i * Cw N y = Ntot dos * (-C N y) := by
    unfold Ntot; rw [← sum_mul, Cw_eq]
  have hsum : ∑ s ∈ (isoAvg sp dos hN).S, (isoAvg sp dos hN).p s
        * Fp (fun k => g tp T k ((isoAvg sp dos hN).Ω s)) N y
      = ∑ s ∈ sp.S, (∑ i, ∑ j, sp.W s i j) / Ntot dos
        * (1 / T * Fp (lw tp (sp.Ω s / T)) N y) :=
    sum_congr rfl (fun s hs => by
      show (∑ i, ∑ j, sp.W s i j) / Ntot dos * Fp (fun k => g tp T k (sp.Ω s)) N y = _
      rw [hatom s hs])
  rw [hsum]
  unfold Qplus
  rw [hcoul, hC, sum_congr rfl hpair]
  have : ∀ s ∈ sp.S, Ntot dos * ((∑ i, ∑ j, sp.W s i j) / Ntot dos *
      (1 / T * Fp (lw tp (sp.Ω s / T)) N y))
      = 1 / T * ((∑ i, ∑ j, sp.W s i j) * Fp (lw tp (sp.Ω s / T)) N y) := by
    intro s _; field_simp
  have hR : Ntot dos * ∑ s ∈ sp.S, (∑ i, ∑ j, sp.W s i j) / Ntot dos *
      (1 / T * Fp (lw tp (sp.Ω s / T)) N y)
      = 1 / T * ∑ s ∈ sp.S, (∑ i, ∑ j, sp.W s i j) * Fp (lw tp (sp.Ω s / T)) N y := by
    rw [mul_sum, mul_sum]; exact sum_congr rfl this
  rw [mul_add, mul_sub, hR]
  ring

theorem nrm2_const (dos : Fin B → ℝ) (N : ℕ) (y : ℕ → ℝ) :
    nrm2 dos N (fun _ => y) = Ntot dos * nsq N y := by
  unfold nrm2 Ntot; rw [sum_mul]

/-- Single-band instability of the averaged spectrum implies multiband instability. -/
theorem munstable_of_iso (tp : ℝ) (sp : MSpec B α) (dos : Fin B → ℝ) (hN : 0 < Ntot dos)
    (μ : Fin B → Fin B → ℝ) (Nc : WithTop ℕ) (T : ℝ) (hT : 0 < T)
    (hU : Unstable (Qp tp (isoAvg sp dos hN) (muAvg dos μ) Nc T)) :
    MUnstable dos (Qplus tp sp dos μ Nc T) := by
  obtain ⟨i0⟩ : Nonempty (Fin B) := by
    rcases Nat.eq_zero_or_pos B with hB | hB
    · subst hB; simp [Ntot] at hN
    · exact ⟨⟨0, hB⟩⟩
  intro ε hε
  obtain ⟨N, y, ⟨n, hn, hyn⟩, hlt⟩ := hU ε hε
  refine ⟨N, fun _ => y, ⟨i0, n, hn, hyn⟩, ?_⟩
  rw [Qplus_const tp sp dos hN μ Nc T hT, nrm2_const]
  nlinarith

/-- **Transfer of a CERTIFIED verdict to the multiband model (conditional).** Suppose the
single-band C-RAYLEIGH certificate holds for the nominal spectrum `sp0` with envelope `e`
at `T_low > 0` (`RayleighCert` data `N, x`, `boxMinQ(μ_hi) ≥ 0`). If the true multiband
spectrum `sp` has its DOS average `isoAvg sp dos` in the envelope of `sp0`, and the averaged
Coulomb parameter `∑_ij N_i μ*_ij / N_tot ≤ μ_hi`, then the multiband model is not stable at
`T_low` (same cutoff). The pipeline does NOT check these two hypotheses. -/
theorem certified_transfers_to_multiband (tp e μhi : ℝ) (sp0 : Spectrum α) (Nc : WithTop ℕ)
    (T : ℝ) (hT : 0 < T) (N : ℕ) (x : ℕ → ℝ) (hx : ∃ n ∈ range N, x n ≠ 0)
    (h0 : 0 ≤ boxMinQ tp sp0 e μhi Nc T N x)
    (sp : MSpec B α) (dos : Fin B → ℝ) (hN : 0 < Ntot dos) (henv : InEnv e sp0 (isoAvg sp dos hN))
    (μ : Fin B → Fin B → ℝ) (hμ : muAvg dos μ ≤ μhi) :
    ¬ MStable tp sp dos μ Nc T := fun hs =>
  hs.1 (munstable_of_iso tp sp dos hN μ Nc T hT
    (rayleigh_cert tp e μhi sp0 Nc T N x hx h0 _ henv _ hμ))

end Transfer

end PipelineLogic
