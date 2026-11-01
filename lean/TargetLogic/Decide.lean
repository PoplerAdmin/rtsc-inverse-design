/-
TargetLogic.Decide: the decision layer of `docs/theory/01_eliashberg_targets.md` §0.4 (YAML,
revised after review RL1) and §0.6 (input domains, multiband envelopes, per-input provenance),
followed by the downgrade-only input policy of the core pipeline (`apply_input_policy` /
`blocking_tags`, review R12 finding 1).

Everything here is generic over an ordered field `K` and uses only `Bool`-valued tests, so the
SAME definitions run over `ℚ` (exported truth table, `ExportCases.lean`) and are reasoned about
over `ℝ` (soundness, `TargetLogic.Sound`).

Order (YAML): WELL-FORMED, G-VALID, then the reject rules N1-KAPPA, N1-ROWMAX, N1-DOS, N1-ISO,
NM-CERT, NM-G. A rule returns REJECTED_PROVEN only if its floor test holds AND every input it
uses is `fullDFPT` with one common run id (§0.6c). A floor that holds on lower-tier inputs gives
`screenedOut` (doc status UNDECIDED, "would fire, ranking only"). RANK-AIM never decides and is
not modelled.
-/
import TargetLogic.NMG
import TargetLogic.NMTable
import Mathlib.Algebra.Order.BigOperators.Group.Finset
import Mathlib.Data.Finset.Fold
import Mathlib.Order.Lattice

open Finset

set_option linter.unusedSectionVars false

namespace TargetLogic

/-! ## Provenance, rules, assumptions, outcomes (field independent) -/

/-- Tier of one input (§0.6c). -/
inductive Tier
  | fullDFPT
  | coarseDFPT
  | ml
  | external
  deriving DecidableEq, Repr

/-- Provenance of one input: its tier and the run / structure id it came from. -/
structure Prov where
  tier : Tier
  run : ℕ
  deriving DecidableEq, Repr

/-- The input groups that carry their own provenance (envelope parameters count as inputs). -/
inductive InputName
  | isoA      -- A_nom with δ_A
  | isoLam    -- λ_nom with δ_λ
  | isoM      -- m_nom with δ_m
  | eps       -- ε_sys, ε_ω, q-spread of ω_log
  | mat       -- 𝒜_ij (band/k resolved) with δ_𝒜, DOS weights, test vector
  | dos       -- sheet DOS fractions w_i with δ_w
  deriving DecidableEq, Repr

inductive Rule
  | kappa
  | rowmax
  | dos
  | iso
  | nmCert
  | nmG
  deriving DecidableEq, Repr

/-- Recorded assumption tags (YAML `assumptions`, `common_assumptions`). -/
inductive Assumption
  | meModel
  | certH
  | envSound
  | mDiscrete
  | mSheetIso
  | sheetsCorrect
  | mIso
  | certNM
  deriving DecidableEq, Repr

def commonAssumptions : List Assumption := [.meModel, .certH, .envSound]

/-- Per-rule assumptions (YAML). -/
def Rule.assumptions : Rule → List Assumption
  | .kappa => [.mDiscrete]
  | .rowmax => [.mDiscrete]
  | .dos => [.mSheetIso, .sheetsCorrect]
  | .iso => [.mIso]
  | .nmCert => [.mIso]
  | .nmG => [.mIso, .certNM]

/-- Inputs used by each rule (§0.6c table). -/
def Rule.uses : Rule → List InputName
  | .kappa => [.mat, .eps]
  | .rowmax => [.mat, .eps]
  | .dos => [.isoA, .eps, .dos]
  | .iso => [.isoA, .eps]
  | .nmCert => [.isoM, .isoLam, .eps]
  | .nmG => [.isoM, .isoLam, .eps]

def allRules : List Rule := [.kappa, .rowmax, .dos, .iso, .nmCert, .nmG]

theorem mem_allRules (r : Rule) : r ∈ allRules := by cases r <;> simp [allRules]

/-- Outcome. `rejectedProven r` = doc REJECTED_PROVEN by rule `r`; `screenedOut` = doc
UNDECIDED with "a rule would fire on lower-tier inputs (ranking only)"; `undecided`, `pass`. -/
inductive Outcome
  | rejectedProven (r : Rule)
  | screenedOut
  | undecided
  | pass
  deriving DecidableEq, Repr

/-- Recorded assumptions of an outcome (§0.6d): common list plus the fired rule's list. -/
def Outcome.assumptions : Outcome → List Assumption
  | .rejectedProven r => commonAssumptions ++ r.assumptions
  | _ => []

/-- Inputs used by the outcome's fired rule. -/
def Outcome.inputsUsed : Outcome → List InputName
  | .rejectedProven r => r.uses
  | _ => []

def Outcome.rank : Outcome → ℕ
  | .rejectedProven _ => 2
  | .screenedOut => 1
  | .undecided => 0
  | .pass => 0

def Outcome.isRemoval (o : Outcome) : Bool := decide (1 ≤ o.rank)

/-! ## Inputs (generic field `K`) -/

/-- Isotropic α²F moments with their ladder spreads and provenance. -/
structure IsoData (K : Type*) where
  A : K
  lam : K
  m : K
  dA : K
  dLam : K
  dM : K
  provA : Prov
  provLam : Prov
  provM : Prov

/-- Envelope parameters. -/
structure EpsData (K : Type*) where
  epsSys : K
  epsOmega : K
  qSpread : K
  prov : Prov

/-- Band- or k-resolved coupling-weight matrix `𝒜_ij` (nominal), DOS weights `N_i` (for the
symmetry check), a Collatz–Wielandt test vector `v > 0`, and the spread `δ_𝒜`. -/
structure MatData (K : Type*) where
  B : ℕ
  calA : Fin B → Fin B → K
  n : Fin B → K
  v : Fin B → K
  dCalA : K
  prov : Prov

/-- Sheet DOS fractions `w_i = N_i/N` with their spread `δ_w`. -/
structure DosData (K : Type*) where
  n : ℕ
  w : Fin n → K
  dW : K
  prov : Prov

/-- G-VALID inputs. `imagAbsMaxK` and `anharmShiftAbs` are MAGNITUDES (§0.6a); the mapper takes
the absolute value, and WELL-FORMED rejects negative entries. -/
structure Validity (K : Type*) where
  EF : Option K
  omegaMax : K
  imagAbsMaxK : K
  containsHLiBe : Bool
  anharmonicCheck : Bool
  anharmShiftAbs : K
  alpha2FRecomputedAnharmonic : Bool
  muStarNegativeSuspected : Bool
  admittedSmearings : ℕ
  qMeshes : ℕ
  nSheets : ℕ

/-- Input-quality tags. `unknown` blocks (fail-closed, RL1 finding 7b); `informational` does
not. -/
inductive InputTag (K : Type*)
  | fileIncomplete
  | fileCompleteUnverified
  | dropNonpos (massK : K)
  | negClip (massK : K)
  | sigmaExcluded (blocksRejection blocksCertification : Bool)
  | informational (name : String)
  | unknown (name : String)

structure Inputs (K : Type*) where
  val : Validity K
  iso : IsoData K
  eps : EpsData K
  mat : Option (MatData K)
  dos : Option (DosData K)
  tags : List (InputTag K)

/-- Run configuration: `a0`, `T* = Tc_target`, the §0.3 table `(λ_t, g_t)` (300 K values,
sorted by λ, last λ = 100) and `g(∞)`. -/
structure Config (K : Type*) where
  a0 : K
  Tstar : K
  tbl : List (K × K)
  gInf : K

section Generic

variable {K : Type*} [Field K] [LinearOrder K] [IsStrictOrderedRing K]

/-! ## Derived quantities -/

/-- `max(0, max_i f i)` (computable). -/
def finMax {B : ℕ} (f : Fin B → K) : K := (univ : Finset (Fin B)).fold max 0 f

/-- `min(1, min_i f i)` (computable). -/
def finMin1 {B : ℕ} (f : Fin B → K) : K := (univ : Finset (Fin B)).fold min 1 f

theorem le_finMax {B : ℕ} (f : Fin B → K) (i : Fin B) : f i ≤ finMax f :=
  (Finset.le_fold_max _).mpr (Or.inr ⟨i, mem_univ _, le_rfl⟩)

theorem finMax_nonneg {B : ℕ} (f : Fin B → K) : 0 ≤ finMax f :=
  (Finset.le_fold_max _).mpr (Or.inl le_rfl)

theorem finMax_mono {B : ℕ} (f g : Fin B → K) (h : ∀ i, f i ≤ g i) : finMax f ≤ finMax g :=
  (Finset.fold_max_le _).mpr ⟨finMax_nonneg g, fun i _ => (h i).trans (le_finMax g i)⟩

theorem finMin1_le {B : ℕ} (f : Fin B → K) (i : Fin B) : finMin1 f ≤ f i :=
  (Finset.fold_min_le _).mpr (Or.inr ⟨i, mem_univ _, le_rfl⟩)

theorem finMin1_le_one {B : ℕ} (f : Fin B → K) : finMin1 f ≤ 1 :=
  (Finset.fold_min_le _).mpr (Or.inl le_rfl)

theorem finMin1_mono {B : ℕ} (f g : Fin B → K) (h : ∀ i, f i ≤ g i) : finMin1 f ≤ finMin1 g :=
  (Finset.le_fold_min _).mpr ⟨finMin1_le_one f, fun i _ => (finMin1_le f i).trans (h i)⟩

/-- `A_hi = A_nom (1+δ_A)(1+ε_sys)/(1−ε_ω)`. -/
def Ahi (x : Inputs K) : K :=
  x.iso.A * (1 + x.iso.dA) * (1 + x.eps.epsSys) / (1 - x.eps.epsOmega)
/-- `λ_hi = λ_nom (1+δ_λ)(1+ε_sys)/(1−ε_ω)²`. -/
def lamHi (x : Inputs K) : K :=
  x.iso.lam * (1 + x.iso.dLam) * (1 + x.eps.epsSys) / (1 - x.eps.epsOmega) ^ 2
/-- `m_hi = m_nom (1+δ_m)(1+ε_sys)`. -/
def mHi (x : Inputs K) : K := x.iso.m * (1 + x.iso.dM) * (1 + x.eps.epsSys)

/-- `𝒜_hi,ij = 𝒜_ij (1+δ_𝒜)(1+ε_sys)/(1−ε_ω)` (A-type scaling, §0.6b). -/
def calAhi (md : MatData K) (e : EpsData K) (i j : Fin md.B) : K :=
  md.calA i j * (1 + md.dCalA) * (1 + e.epsSys) / (1 - e.epsOmega)

/-- `κ_hi = max_i (𝒜_hi v)_i / v_i` (Collatz–Wielandt, §0.6b). -/
def kappaHi (md : MatData K) (e : EpsData K) : K :=
  finMax (fun i => (∑ j, calAhi md e i j * md.v j) / md.v i)

/-- `rowA_hi,i = Σ_j 𝒜_hi,ij`. -/
def rowAhi (md : MatData K) (e : EpsData K) (i : Fin md.B) : K := ∑ j, calAhi md e i j

def rowMaxHi (md : MatData K) (e : EpsData K) : K := finMax (rowAhi md e)

/-- `w_min,lo = min_i w_i (1 − δ_w)`. -/
def wminLo (dd : DosData K) : K := finMin1 (fun i => dd.w i * (1 - dd.dW))

def Config.floor (c : Config K) : K := c.a0 * c.Tstar

/-- `g(λ_up)·s2`, `s2 = (T*/300)²`, `λ_up` = first table λ ≥ λ_hi, else `g(∞)`. -/
def Config.gUp (c : Config K) (lh : K) : K := (c.Tstar / 300) ^ 2 * tableLookup c.tbl c.gInf lh

/-- The certified NM-G table `(λ_k, g_lo(λ_k) = 90000·m0_k)` (doc 01 §0.3, `NMTable`). -/
def certTable : List (K × K) := nmM0Table.map (fun p => (p.1, 90000 * p.2))

/-- The doc configuration: given `a0`, `T*`, the certified table and `g_lo(∞) = 90000·m0_∞`. -/
def Config.ofCert (a0 Tstar : K) : Config K := ⟨a0, Tstar, certTable, 90000 * nmM0Inf⟩

/-- Derived flags (§0.6a): never set independently. -/
def mbPresent (x : Inputs K) : Bool := x.mat.isSome || x.dos.isSome
def kResolved (x : Inputs K) : Bool := x.mat.isSome

/-! ## Gates -/

def in0half (d : K) : Bool := decide (0 ≤ d) && decide (d ≤ 1 / 2)

def matWF : Option (MatData K) → Bool
  | none => true
  | some md =>
    decide (1 ≤ md.B) && decide (∀ i j, 0 ≤ md.calA i j) && decide (∀ i, 0 < md.n i) &&
    decide (∀ i, 0 < md.v i) && in0half md.dCalA &&
    decide (∀ i j, |md.n i * md.calA i j - md.n j * md.calA j i| ≤
      1 / 100 * max (md.n i * md.calA i j) (md.n j * md.calA j i))

def dosWF (nS : ℕ) : Option (DosData K) → Bool
  | none => true
  | some dd =>
    decide (dd.n = nS) && decide (∀ i, 0 < dd.w i ∧ dd.w i ≤ 1) &&
    decide (99 / 100 ≤ ∑ i, dd.w i ∧ ∑ i, dd.w i ≤ 101 / 100) && in0half dd.dW &&
    decide (∑ i, dd.w i * (1 - dd.dW) ≤ 1)

/-- **WELL-FORMED** (§0.6a). -/
def wf (x : Inputs K) : Bool :=
  decide (0 < x.iso.A) && decide (0 < x.iso.lam) && decide (0 < x.iso.m) &&
  decide (0 < x.val.omegaMax) &&
  (match x.val.EF with | some e => decide (0 < e) | none => true) &&
  decide (0 ≤ x.val.imagAbsMaxK) && decide (0 ≤ x.val.anharmShiftAbs) &&
  decide (1 ≤ x.val.nSheets) &&
  in0half x.iso.dA && in0half x.iso.dLam && in0half x.iso.dM && in0half x.eps.qSpread &&
  decide (1 / 10 ≤ x.eps.epsSys) && decide (0 ≤ x.eps.epsOmega) && decide (x.eps.epsOmega < 1) &&
  decide (x.iso.A ^ 2 ≤ x.iso.lam * x.iso.m / 4 * (1 + 1 / 100)) &&
  matWF x.mat && dosWF x.val.nSheets x.dos

/-- **G-VALID** (YAML, revised). -/
def gvalid (x : Inputs K) : Bool :=
  (match x.val.EF with
    | some e => decide (0 < e) && decide (lamHi x * x.val.omegaMax ≤ e / 10)
    | none => false) &&
  decide (x.val.imagAbsMaxK ≤ 30) &&
  (!x.val.containsHLiBe || x.val.anharmonicCheck) &&
  (!decide (1 / 10 < x.val.anharmShiftAbs) || x.val.alpha2FRecomputedAnharmonic) &&
  !x.val.muStarNegativeSuspected &&
  decide (3 ≤ x.val.admittedSmearings) && decide (2 ≤ x.val.qMeshes) &&
  decide (1 / 10 ≤ x.eps.epsSys) &&
  decide (max (1 / 50) (max x.val.anharmShiftAbs x.eps.qSpread) ≤ x.eps.epsOmega) &&
  (!decide (1 < x.val.nSheets) || mbPresent x)

/-! ## Rules -/

/-- The floor test of each rule (including its own `pre` conditions). -/
def floorR (cfg : Config K) (x : Inputs K) : Rule → Bool
  | .kappa => match x.mat with
    | some md => decide (kappaHi md x.eps < cfg.floor)
    | none => false
  | .rowmax => match x.mat with
    | some md => decide (rowMaxHi md x.eps < cfg.floor)
    | none => false
  | .dos => match x.dos with
    | some dd => decide (0 < wminLo dd) && decide (Ahi x / wminLo dd < cfg.floor)
    | none => false
  | .iso => decide (x.val.nSheets = 1) && decide (Ahi x < cfg.floor) &&
      (decide (Ahi x < 3 / 4 * cfg.floor) || kResolved x)
  | .nmCert => decide (x.val.nSheets = 1) && decide (mHi x < 4 * cfg.floor ^ 2 / lamHi x)
  | .nmG => decide (x.val.nSheets = 1) && decide (mHi x < cfg.gUp (lamHi x)) &&
      (decide (mHi x < 3 / 4 * cfg.gUp (lamHi x)) || kResolved x)

def provOf (x : Inputs K) : InputName → Option Prov
  | .isoA => some x.iso.provA
  | .isoLam => some x.iso.provLam
  | .isoM => some x.iso.provM
  | .eps => some x.eps.prov
  | .mat => x.mat.map (·.prov)
  | .dos => x.dos.map (·.prov)

end Generic

/-- All entries present, the first is `fullDFPT`, and all equal it (same tier and run). -/
def fullSameRun : List (Option Prov) → Bool
  | [] => false
  | none :: _ => false
  | some p :: rest => (p.tier == .fullDFPT) && rest.all (fun q => q == some p)

theorem fullSameRun_spec (ps : List (Option Prov)) (h : fullSameRun ps = true) :
    ∃ p : Prov, p.tier = .fullDFPT ∧ ∀ q ∈ ps, q = some p := by
  match ps, h with
  | some p :: rest, h =>
    simp only [fullSameRun, Bool.and_eq_true, beq_iff_eq, List.all_eq_true] at h
    refine ⟨p, h.1, ?_⟩
    intro q hq
    rcases List.mem_cons.mp hq with rfl | hq
    · rfl
    · exact h.2 q hq

section Generic2

variable {K : Type*} [Field K] [LinearOrder K] [IsStrictOrderedRing K]

/-- Every input the rule uses is `fullDFPT` from one run. -/
def provOK (x : Inputs K) (r : Rule) : Bool := fullSameRun (r.uses.map (provOf x))

def fires (cfg : Config K) (x : Inputs K) (r : Rule) : Bool := floorR cfg x r && provOK x r

def firstRule (cfg : Config K) (x : Inputs K) : Option Rule := allRules.find? (fires cfg x)

def wouldFire (cfg : Config K) (x : Inputs K) : Bool := allRules.any (floorR cfg x)

def decideCore (cfg : Config K) (x : Inputs K) : Outcome :=
  if wf x && gvalid x then
    match firstRule cfg x with
    | some r => .rejectedProven r
    | none => if wouldFire cfg x then .screenedOut else .pass
  else .undecided

def InputTag.blocksRejection : InputTag K → Bool
  | .fileIncomplete => true
  | .fileCompleteUnverified => true
  | .dropNonpos m => decide (0 < m)
  | .negClip _ => false
  | .sigmaExcluded r _ => r
  | .informational _ => false
  | .unknown _ => true

def blocked (x : Inputs K) : Bool := x.tags.any InputTag.blocksRejection

/-- Downgrade-only policy: a removal with a blocking tag becomes UNDECIDED. -/
def applyPolicy (b : Bool) (o : Outcome) : Outcome :=
  if b && o.isRemoval then .undecided else o

/-- **The decision.** -/
def decideT (cfg : Config K) (x : Inputs K) : Outcome := applyPolicy (blocked x) (decideCore cfg x)

/-! ## Structural theorems -/

theorem firstRule_some {cfg : Config K} {x : Inputs K} {r : Rule}
    (h : firstRule cfg x = some r) : fires cfg x r = true := List.find?_some h

theorem firstRule_none {cfg : Config K} {x : Inputs K} (h : firstRule cfg x = none) (r : Rule) :
    fires cfg x r = false := by
  have := List.find?_eq_none.mp h r (mem_allRules r)
  simpa using this

theorem decideCore_cases (cfg : Config K) (x : Inputs K) :
    (∃ r, decideCore cfg x = .rejectedProven r) ∨ decideCore cfg x = .screenedOut ∨
    decideCore cfg x = .undecided ∨ decideCore cfg x = .pass := by
  unfold decideCore
  by_cases hg : (wf x && gvalid x) = true
  · rw [if_pos hg]
    cases firstRule cfg x with
    | some r => exact Or.inl ⟨r, rfl⟩
    | none =>
      by_cases hw : wouldFire cfg x = true
      · simp [hw]
      · simp [hw]
  · rw [if_neg hg]; simp

theorem applyPolicy_eq (b : Bool) (o : Outcome) :
    applyPolicy b o = o ∨ applyPolicy b o = .undecided := by
  unfold applyPolicy; split_ifs <;> simp

/-- Unpack `decideT = rejectedProven r`. -/
theorem decide_rejectedProven_core {cfg : Config K} {x : Inputs K} {r : Rule}
    (h : decideT cfg x = .rejectedProven r) :
    blocked x = false ∧ (wf x && gvalid x) = true ∧ firstRule cfg x = some r := by
  unfold decideT applyPolicy at h
  by_cases hb : (blocked x && (decideCore cfg x).isRemoval) = true
  · rw [if_pos hb] at h; cases h
  · rw [if_neg hb] at h
    have hrem : (decideCore cfg x).isRemoval = true := by rw [h]; rfl
    have hbl : blocked x = false := by
      cases hx : blocked x
      · rfl
      · rw [hx, hrem] at hb; exact absurd rfl hb
    refine ⟨hbl, ?_⟩
    unfold decideCore at h
    by_cases hg : (wf x && gvalid x) = true
    · rw [if_pos hg] at h
      refine ⟨hg, ?_⟩
      cases hf : firstRule cfg x with
      | some r' => rw [hf] at h; simp at h; rw [h]
      | none => rw [hf] at h; simp only at h; split_ifs at h
    · rw [if_neg hg] at h; cases h

/-- **(i) REJECTED_PROVEN, combinatorial.** The gates passed, no blocking tag, the rule's floor
test (with its `pre` clauses) held, and every input it uses is `fullDFPT` from one run. -/
theorem decide_rejectedProven {cfg : Config K} {x : Inputs K} {r : Rule}
    (h : decideT cfg x = .rejectedProven r) :
    wf x = true ∧ gvalid x = true ∧ blocked x = false ∧ floorR cfg x r = true ∧
      ∃ p : Prov, p.tier = .fullDFPT ∧ ∀ n ∈ r.uses, provOf x n = some p := by
  obtain ⟨hb, hg, hf⟩ := decide_rejectedProven_core h
  have hfi := firstRule_some hf
  simp only [Bool.and_eq_true] at hg
  simp only [fires, Bool.and_eq_true] at hfi
  obtain ⟨p, hp, hall⟩ := fullSameRun_spec _ hfi.2
  refine ⟨hg.1, hg.2, hb, hfi.1, p, hp, fun n hn => hall _ (List.mem_map_of_mem hn)⟩

/-- The outcome records the common assumptions and the fired rule's assumptions. -/
theorem rejectedProven_assumptions (r : Rule) :
    (Outcome.rejectedProven r).assumptions = commonAssumptions ++ r.assumptions := rfl

theorem iso_records_mIso : Assumption.mIso ∈ (Outcome.rejectedProven .iso).assumptions := by
  decide

theorem nmG_records_certNM :
    Assumption.certNM ∈ (Outcome.rejectedProven .nmG).assumptions ∧
    Assumption.mIso ∈ (Outcome.rejectedProven .nmG).assumptions := by decide

/-- **(ii) WELL-FORMED fails ⇒ UNDECIDED.** -/
theorem decide_of_not_wf (cfg : Config K) (x : Inputs K) (h : wf x = false) :
    decideT cfg x = .undecided := by
  simp [decideT, applyPolicy, decideCore, h, Outcome.isRemoval, Outcome.rank]

/-- **(ii) G-VALID fails ⇒ UNDECIDED.** -/
theorem decide_of_not_gvalid (cfg : Config K) (x : Inputs K) (h : gvalid x = false) :
    decideT cfg x = .undecided := by
  simp [decideT, applyPolicy, decideCore, h, Outcome.isRemoval, Outcome.rank]

/-- **(ii) Provenance.** If some input used by rule `r` is absent or not `fullDFPT`, `r` is not
the recorded rule of a REJECTED_PROVEN. -/
theorem decide_not_proven_of_input {cfg : Config K} {x : Inputs K} (r : Rule) (n : InputName)
    (hn : n ∈ r.uses) (hlow : ∀ p, provOf x n = some p → p.tier ≠ .fullDFPT) :
    decideT cfg x ≠ .rejectedProven r := by
  intro h
  obtain ⟨-, -, -, -, p, hp, hall⟩ := decide_rejectedProven h
  exact hlow p (hall n hn) hp

/-- **(ii) No `fullDFPT` input at all (ML / coarse / external data) ⇒ never REJECTED_PROVEN.** -/
theorem decide_not_proven_of_lowTier (cfg : Config K) (x : Inputs K)
    (hlow : ∀ n p, provOf x n = some p → p.tier ≠ .fullDFPT) (r : Rule) :
    decideT cfg x ≠ .rejectedProven r := by
  exact decide_not_proven_of_input r .eps (by cases r <;> simp [Rule.uses]) (hlow _)

/-- Mixed runs: two used inputs from different runs ⇒ not REJECTED_PROVEN by that rule. -/
theorem decide_not_proven_of_mixed_runs {cfg : Config K} {x : Inputs K} (r : Rule)
    (n₁ n₂ : InputName) (h₁ : n₁ ∈ r.uses) (h₂ : n₂ ∈ r.uses) (p₁ p₂ : Prov)
    (hp₁ : provOf x n₁ = some p₁) (hp₂ : provOf x n₂ = some p₂) (hne : p₁.run ≠ p₂.run) :
    decideT cfg x ≠ .rejectedProven r := by
  intro h
  obtain ⟨-, -, -, -, p, -, hall⟩ := decide_rejectedProven h
  have e₁ := (hall n₁ h₁).symm.trans hp₁
  have e₂ := (hall n₂ h₂).symm.trans hp₂
  simp only [Option.some.injEq] at e₁ e₂
  exact hne (by rw [← e₁, ← e₂])

/-- **(ii) Blocking tag ⇒ UNDECIDED or PASS.** -/
theorem decide_of_blocked (cfg : Config K) (x : Inputs K) (h : blocked x = true) :
    decideT cfg x = .undecided ∨ decideT cfg x = .pass := by
  unfold decideT applyPolicy
  rcases decideCore_cases cfg x with ⟨r, h1⟩ | h1 | h1 | h1 <;>
    simp [h1, h, Outcome.isRemoval, Outcome.rank]

/-- `decideCore` ignores the tags. -/
theorem decideCore_tags (cfg : Config K) (x : Inputs K) (ts : List (InputTag K)) :
    decideCore cfg { x with tags := ts } = decideCore cfg x := rfl

theorem blocked_mono (x : Inputs K) (ts ts' : List (InputTag K)) (hsub : ∀ t ∈ ts, t ∈ ts')
    (h : blocked { x with tags := ts } = true) : blocked { x with tags := ts' } = true := by
  simp only [blocked, List.any_eq_true] at *
  obtain ⟨t, ht, hb⟩ := h
  exact ⟨t, hsub t ht, hb⟩

/-- **(iii) Downgrade only (tags).** -/
theorem decide_downgrade (cfg : Config K) (x : Inputs K) (ts ts' : List (InputTag K))
    (hsub : ∀ t ∈ ts, t ∈ ts') :
    decideT cfg { x with tags := ts' } = decideT cfg { x with tags := ts } ∨
    decideT cfg { x with tags := ts' } = .undecided := by
  have e1 : decideCore cfg { x with tags := ts' } = decideCore cfg x := rfl
  have e2 : decideCore cfg { x with tags := ts } = decideCore cfg x := rfl
  unfold decideT
  rw [e1, e2]
  unfold applyPolicy
  by_cases h' : (blocked { x with tags := ts' } && (decideCore cfg x).isRemoval) = true
  · right; rw [if_pos h']
  · left
    have h : ¬ (blocked { x with tags := ts } && (decideCore cfg x).isRemoval) = true := by
      intro hb
      simp only [Bool.and_eq_true] at hb h'
      exact h' ⟨blocked_mono x ts ts' hsub hb.1, hb.2⟩
    rw [if_neg h', if_neg h]

theorem decide_downgrade_rank (cfg : Config K) (x : Inputs K) (ts ts' : List (InputTag K))
    (hsub : ∀ t ∈ ts, t ∈ ts') :
    (decideT cfg { x with tags := ts' }).rank ≤ (decideT cfg { x with tags := ts }).rank := by
  rcases decide_downgrade cfg x ts ts' hsub with h | h
  · rw [h]
  · rw [h]; simp [Outcome.rank]

/-! ## Envelope direction -/

theorem wf_spec {x : Inputs K} (h : wf x = true) :
    0 < x.iso.A ∧ 0 < x.iso.lam ∧ 0 < x.iso.m ∧ 0 ≤ x.iso.dA ∧ 0 ≤ x.iso.dLam ∧ 0 ≤ x.iso.dM ∧
    1 / 10 ≤ x.eps.epsSys ∧ 0 ≤ x.eps.epsOmega ∧ x.eps.epsOmega < 1 := by
  simp only [wf, in0half, Bool.and_eq_true, decide_eq_true_eq] at h
  tauto

theorem Ahi_ge (x : Inputs K) (h : wf x = true) : x.iso.A ≤ Ahi x := by
  obtain ⟨hA, -, -, hdA, -, -, hs, hw0, hw1⟩ := wf_spec h
  unfold Ahi
  rw [le_div_iff₀ (by linarith)]
  have h1 : 1 ≤ (1 + x.iso.dA) * (1 + x.eps.epsSys) := by nlinarith
  nlinarith [mul_nonneg hA.le hw0, mul_nonneg hA.le (sub_nonneg.mpr h1)]

theorem mHi_ge (x : Inputs K) (h : wf x = true) : x.iso.m ≤ mHi x := by
  obtain ⟨-, -, hm, -, -, hdm, hs, -, -⟩ := wf_spec h
  unfold mHi
  have h1 : 1 ≤ (1 + x.iso.dM) * (1 + x.eps.epsSys) := by nlinarith
  nlinarith [mul_nonneg hm.le (sub_nonneg.mpr h1)]

theorem lamHi_ge (x : Inputs K) (h : wf x = true) : x.iso.lam ≤ lamHi x := by
  obtain ⟨-, hl, -, -, hdl, -, hs, hw0, hw1⟩ := wf_spec h
  unfold lamHi
  have h0 : 0 < 1 - x.eps.epsOmega := by linarith
  have hp : 0 < (1 - x.eps.epsOmega) ^ 2 := pow_pos h0 2
  rw [le_div_iff₀ hp]
  have hq : (1 - x.eps.epsOmega) ^ 2 ≤ 1 := by nlinarith
  have h1 : 1 ≤ (1 + x.iso.dLam) * (1 + x.eps.epsSys) := by nlinarith
  nlinarith [mul_nonneg hl.le (sub_nonneg.mpr h1), mul_nonneg hl.le (sub_nonneg.mpr hq)]

theorem lamHi_pos (x : Inputs K) (h : wf x = true) : 0 < lamHi x :=
  lt_of_lt_of_le (wf_spec h).2.1 (lamHi_ge x h)

theorem wf_mat' {x : Inputs K} (h : wf x = true) : matWF x.mat = true := by
  simp only [wf, Bool.and_eq_true] at h; tauto

theorem wf_dos' {x : Inputs K} (h : wf x = true) : dosWF x.val.nSheets x.dos = true := by
  simp only [wf, Bool.and_eq_true] at h; tauto

theorem matWF_spec' {md : MatData K} (h : matWF (some md) = true) :
    1 ≤ md.B ∧ (∀ i j, 0 ≤ md.calA i j) ∧ (∀ i, 0 < md.n i) ∧ (∀ i, 0 < md.v i) ∧
    0 ≤ md.dCalA ∧ md.dCalA ≤ 1 / 2 := by
  simp only [matWF, in0half, Bool.and_eq_true, decide_eq_true_eq] at h; tauto

theorem dosWF_spec' {nS : ℕ} {dd : DosData K} (h : dosWF nS (some dd) = true) :
    dd.n = nS ∧ (∀ i, 0 < dd.w i ∧ dd.w i ≤ 1) ∧ 0 ≤ dd.dW ∧ dd.dW ≤ 1 / 2 ∧
    99 / 100 ≤ ∑ i, dd.w i ∧ ∑ i, dd.w i * (1 - dd.dW) ≤ 1 := by
  simp only [dosWF, in0half, Bool.and_eq_true, decide_eq_true_eq] at h; tauto

/-! ## Rank lemmas for the core -/

theorem core_rank_two {cfg : Config K} {x : Inputs K} (hg : (wf x && gvalid x) = true) (r : Rule)
    (hr : fires cfg x r = true) : (decideCore cfg x).rank = 2 := by
  unfold decideCore
  rw [if_pos hg]
  cases hf : firstRule cfg x with
  | some r' => rfl
  | none => rw [firstRule_none hf r] at hr; cases hr

theorem core_rank_one {cfg : Config K} {x : Inputs K} (hg : (wf x && gvalid x) = true) (r : Rule)
    (hr : floorR cfg x r = true) : 1 ≤ (decideCore cfg x).rank := by
  have hw : wouldFire cfg x = true := List.any_eq_true.mpr ⟨r, mem_allRules r, hr⟩
  unfold decideCore
  rw [if_pos hg]
  cases firstRule cfg x with
  | some r' => simp [Outcome.rank]
  | none => simp [hw, Outcome.rank]

theorem core_of_rank_two {cfg : Config K} {x : Inputs K} (h : (decideCore cfg x).rank = 2) :
    (wf x && gvalid x) = true ∧ ∃ r, fires cfg x r = true := by
  unfold decideCore at h
  by_cases hg : (wf x && gvalid x) = true
  · rw [if_pos hg] at h
    refine ⟨hg, ?_⟩
    cases hf : firstRule cfg x with
    | some r => exact ⟨r, firstRule_some hf⟩
    | none =>
      rw [hf] at h
      by_cases hw : wouldFire cfg x = true <;> simp [hw, Outcome.rank] at h
  · rw [if_neg hg] at h; simp [Outcome.rank] at h

theorem core_of_rank_one {cfg : Config K} {x : Inputs K} (h : 1 ≤ (decideCore cfg x).rank) :
    (wf x && gvalid x) = true ∧ ∃ r, floorR cfg x r = true := by
  unfold decideCore at h
  by_cases hg : (wf x && gvalid x) = true
  · rw [if_pos hg] at h
    refine ⟨hg, ?_⟩
    cases hf : firstRule cfg x with
    | some r =>
      have := firstRule_some hf
      simp only [fires, Bool.and_eq_true] at this
      exact ⟨r, this.1⟩
    | none =>
      rw [hf] at h
      by_cases hw : wouldFire cfg x = true
      · obtain ⟨r, -, hr⟩ := List.any_eq_true.mp hw; exact ⟨r, hr⟩
      · simp [hw, Outcome.rank] at h
  · rw [if_neg hg] at h; simp [Outcome.rank] at h

theorem policy_rank_mono (b : Bool) (o o' : Outcome) (h : o'.rank ≤ o.rank) :
    (applyPolicy b o').rank ≤ (applyPolicy b o).rank := by
  unfold applyPolicy Outcome.isRemoval
  cases b <;> simp only [Bool.false_and, Bool.true_and, decide_eq_true_eq] <;>
    split_ifs <;> simp_all [Outcome.rank]

end Generic2

end TargetLogic
