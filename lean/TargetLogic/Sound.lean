/-
TargetLogic.Sound: semantic soundness of every REJECTED_PROVEN rule of `decideT` (over ℝ),
with the certified constants pinned (RL1 finding 9) and one theorem per rule that uses only the
hypotheses that rule needs (RL1 finding 3).

Constants: `tpPhys = 2π` (kernel `2πkT`), `a0Cert = 4.33148145`, `etaCert = 1e-9`.
`CertH` := `H tpPhys a0Cert etaCert ∧ PropN tpPhys a0Cert` — the core interval-arithmetic
certificates (Prop. 6, wMB Prop. N). They are hypotheses here, never proved in Lean.

Each `*_sound` theorem: `CertH`, `cfg.a0 = a0Cert`, `T* > 0`, the rule's ENV-SOUND / model
hypotheses, and `decideT cfg x = rejectedProven r` ⇒ stability at every `T ≥ T*`
(Tc ≤ T*). The multiband rules also give the strict form `∃ T' < T*` (stable on `[T', ∞)`).
-/
import TargetLogic.Decide
import TargetLogic.Matrix
import Mathlib.Analysis.SpecialFunctions.Trigonometric.Basic

open Finset

namespace TargetLogic

open MultiBand TcR1

noncomputable section

/-- Kernel constant `2π` of `(2πkT)²`. -/
def tpPhys : ℝ := 2 * Real.pi
/-- Certified `a0`. -/
def a0Cert : ℝ := 433148145 / 100000000
/-- Certified `η`. -/
def etaCert : ℝ := 1 / 1000000000

/-- CERT-H: the two certified numerical hypotheses at the pinned constants. -/
def CertH : Prop := H tpPhys a0Cert etaCert ∧ PropN tpPhys a0Cert

theorem a0Cert_pos : 0 < a0Cert := by norm_num [a0Cert]
theorem etaCert_pos : 0 < etaCert := by norm_num [etaCert]

/-- Between `b/a0` and `T*`: a positive `T' < T*` with `b < a0 T'`. -/
theorem exists_T' (a0 Ts b : ℝ) (ha0 : 0 < a0) (hT : 0 < Ts) (hb : b < a0 * Ts) :
    ∃ T', T' < Ts ∧ 0 < T' ∧ b < a0 * T' := by
  have hb' : b / a0 < Ts := by rw [div_lt_iff₀ ha0]; linarith
  refine ⟨max (Ts / 2) ((b / a0 + Ts) / 2), max_lt (by linarith) (by linarith),
    lt_max_of_lt_left (by linarith), ?_⟩
  have h1 : b / a0 < max (Ts / 2) ((b / a0 + Ts) / 2) :=
    lt_max_of_lt_right (by linarith)
  rw [div_lt_iff₀ ha0] at h1
  linarith

/-- Collatz–Wielandt route with the pinned constants: `𝒜 v ≤ b v` (`v > 0`), `b < a0 T*`
⇒ stable on `[T*, ∞)` and on `[T', ∞)` for some `0 < T' < T*`. -/
theorem stable_of_cw (hC : CertH) {B : ℕ} {α : Type*} (sp : MSpec B α) (dos : Fin B → ℝ)
    (hdos : ∀ i, 0 < dos i) (v : Fin B → ℝ) (hv : ∀ i, 0 < v i) (b : ℝ)
    (hsub : ∀ i, ∑ j, Amat sp dos i j * v j ≤ b * v i) (Ts : ℝ) (hT : 0 < Ts)
    (hb : b < a0Cert * Ts) :
    (∀ μ : Fin B → Fin B → ℝ, CoulPSD dos μ → ∀ Nc : WithTop ℕ, ∀ T : ℝ, Ts ≤ T →
      MStable tpPhys sp dos μ Nc T) ∧
    ∃ T', T' < Ts ∧ 0 < T' ∧ ∀ μ : Fin B → Fin B → ℝ, CoulPSD dos μ → ∀ Nc : WithTop ℕ,
      ∀ T : ℝ, T' ≤ T → MStable tpPhys sp dos μ Nc T := by
  have key : ∀ T0, 0 < T0 → b < a0Cert * T0 → ∀ μ : Fin B → Fin B → ℝ, CoulPSD dos μ →
      ∀ Nc : WithTop ℕ, ∀ T : ℝ, T0 ≤ T → MStable tpPhys sp dos μ Nc T := by
    intro T0 hT0 hb0 μ hμ Nc T hT
    apply mstable_of_gt tpPhys a0Cert etaCert hC.1 etaCert_pos.le a0Cert_pos hC.2 sp dos hdos v hv
      b hsub μ hμ Nc T (by linarith)
    have : a0Cert * T0 ≤ a0Cert * T := mul_le_mul_of_nonneg_left hT a0Cert_pos.le
    linarith
  obtain ⟨T', h1, h2, h3⟩ := exists_T' a0Cert Ts b a0Cert_pos hT hb
  exact ⟨key Ts hT hb, T', h1, h2, key T' h2 h3⟩

/-! ## Extracting facts from the gates -/

theorem wf_mat {x : Inputs ℝ} (h : wf x = true) : matWF x.mat = true := by
  simp only [wf, Bool.and_eq_true] at h; tauto

theorem wf_dos {x : Inputs ℝ} (h : wf x = true) : dosWF x.val.nSheets x.dos = true := by
  simp only [wf, Bool.and_eq_true] at h; tauto

theorem matWF_spec {md : MatData ℝ} (h : matWF (some md) = true) :
    1 ≤ md.B ∧ (∀ i j, 0 ≤ md.calA i j) ∧ (∀ i, 0 < md.n i) ∧ (∀ i, 0 < md.v i) ∧
    0 ≤ md.dCalA ∧ md.dCalA ≤ 1 / 2 := by
  simp only [matWF, in0half, Bool.and_eq_true, decide_eq_true_eq] at h; tauto

theorem dosWF_spec {nS : ℕ} {dd : DosData ℝ} (h : dosWF nS (some dd) = true) :
    dd.n = nS ∧ (∀ i, 0 < dd.w i ∧ dd.w i ≤ 1) ∧ 0 ≤ dd.dW ∧ dd.dW ≤ 1 / 2 := by
  simp only [dosWF, in0half, Bool.and_eq_true, decide_eq_true_eq] at h; tauto

theorem calAhi_nonneg {x : Inputs ℝ} (hw : wf x = true) {md : MatData ℝ} (hm : x.mat = some md)
    (i j : Fin md.B) : 0 ≤ calAhi md x.eps i j := by
  have hs := matWF_spec (hm ▸ wf_mat hw)
  obtain ⟨-, -, -, -, -, -, hes, hw0, hw1⟩ := wf_spec hw
  unfold calAhi
  apply div_nonneg _ (by linarith)
  have := hs.2.1 i j
  have : 0 ≤ 1 + md.dCalA := by linarith [hs.2.2.2.2.1]
  have : 0 ≤ 1 + x.eps.epsSys := by linarith
  positivity

/-! ## Per-rule soundness -/

variable (cfg : Config ℝ) (x : Inputs ℝ)

/-- **N1-KAPPA.** Model: band/k-discretised `MSpec md.B α` (M-DISCRETE); ENV-SOUND:
`𝒜_true ≤ 𝒜_hi` entrywise. Then REJECTED_PROVEN(kappa) ⇒ stable for every `T ≥ T*`
(and on `[T', ∞)` for some `T' < T*`), every PSD Coulomb matrix and cutoff. -/
theorem kappa_sound (hC : CertH) (ha0 : cfg.a0 = a0Cert) (hT : 0 < cfg.Tstar)
    (h : decideT cfg x = .rejectedProven .kappa) (md : MatData ℝ) (hm : x.mat = some md)
    {α : Type*} (sp : MSpec md.B α) (dos : Fin md.B → ℝ) (hdos : ∀ i, 0 < dos i)
    (henv : ∀ i j, Amat sp dos i j ≤ calAhi md x.eps i j) :
    (∀ μ, CoulPSD dos μ → ∀ Nc : WithTop ℕ, ∀ T : ℝ, cfg.Tstar ≤ T →
      MStable tpPhys sp dos μ Nc T) ∧
    ∃ T', T' < cfg.Tstar ∧ 0 < T' ∧ ∀ μ, CoulPSD dos μ → ∀ Nc : WithTop ℕ, ∀ T : ℝ, T' ≤ T →
      MStable tpPhys sp dos μ Nc T := by
  obtain ⟨hw, -, -, hf, -⟩ := decide_rejectedProven h
  have hv := (matWF_spec (hm ▸ wf_mat hw)).2.2.2.1
  simp only [floorR, hm, decide_eq_true_eq, Config.floor, ha0] at hf
  apply stable_of_cw hC sp dos hdos md.v hv (kappaHi md x.eps) _ cfg.Tstar hT hf
  intro i
  have h1 : ∑ j, Amat sp dos i j * md.v j ≤ ∑ j, calAhi md x.eps i j * md.v j :=
    sum_le_sum (fun j _ => mul_le_mul_of_nonneg_right (henv i j) (hv j).le)
  have h2 : (∑ j, calAhi md x.eps i j * md.v j) / md.v i ≤ kappaHi md x.eps :=
    le_finMax (fun i => (∑ j, calAhi md x.eps i j * md.v j) / md.v i) i
  rw [div_le_iff₀ (hv i)] at h2
  linarith

/-- **N1-ROWMAX.** Same model and ENV-SOUND as N1-KAPPA. -/
theorem rowmax_sound (hC : CertH) (ha0 : cfg.a0 = a0Cert) (hT : 0 < cfg.Tstar)
    (h : decideT cfg x = .rejectedProven .rowmax) (md : MatData ℝ) (hm : x.mat = some md)
    {α : Type*} (sp : MSpec md.B α) (dos : Fin md.B → ℝ) (hdos : ∀ i, 0 < dos i)
    (henv : ∀ i j, Amat sp dos i j ≤ calAhi md x.eps i j) :
    (∀ μ, CoulPSD dos μ → ∀ Nc : WithTop ℕ, ∀ T : ℝ, cfg.Tstar ≤ T →
      MStable tpPhys sp dos μ Nc T) ∧
    ∃ T', T' < cfg.Tstar ∧ 0 < T' ∧ ∀ μ, CoulPSD dos μ → ∀ Nc : WithTop ℕ, ∀ T : ℝ, T' ≤ T →
      MStable tpPhys sp dos μ Nc T := by
  obtain ⟨hw, -, -, hf, -⟩ := decide_rejectedProven h
  simp only [floorR, hm, decide_eq_true_eq, Config.floor, ha0] at hf
  apply stable_of_cw hC sp dos hdos (fun _ => 1) (fun _ => one_pos) (rowMaxHi md x.eps) _
    cfg.Tstar hT hf
  intro i
  have h1 : ∑ j, Amat sp dos i j * 1 ≤ rowAhi md x.eps i := by
    simp only [mul_one, rowAhi]; exact sum_le_sum (fun j _ => henv i j)
  have h2 := le_finMax (rowAhi md x.eps) i
  simp only [rowMaxHi] at *
  linarith

/-- **N1-DOS.** Model: sheet model `MSpec dd.n α` (M-SHEET-ISO, SHEETS-CORRECT); ENV-SOUND:
`w_i (1−δ_w) ≤` true DOS fraction for every sheet, and `A_iso,true ≤ A_hi`. -/
theorem dos_sound (hC : CertH) (ha0 : cfg.a0 = a0Cert) (hT : 0 < cfg.Tstar)
    (h : decideT cfg x = .rejectedProven .dos) (dd : DosData ℝ) (hd : x.dos = some dd)
    {α : Type*} (sp : MSpec dd.n α) (dos : Fin dd.n → ℝ) (hdos : ∀ i, 0 < dos i)
    (hw : ∀ i, dd.w i * (1 - dd.dW) ≤ wfrac dos i) (hA : Aiso sp dos ≤ Ahi x) :
    (∀ μ, CoulPSD dos μ → ∀ Nc : WithTop ℕ, ∀ T : ℝ, cfg.Tstar ≤ T →
      MStable tpPhys sp dos μ Nc T) ∧
    ∃ T', T' < cfg.Tstar ∧ 0 < T' ∧ ∀ μ, CoulPSD dos μ → ∀ Nc : WithTop ℕ, ∀ T : ℝ, T' ≤ T →
      MStable tpPhys sp dos μ Nc T := by
  obtain ⟨-, -, -, hf, -⟩ := decide_rejectedProven h
  simp only [floorR, hd, Bool.and_eq_true, decide_eq_true_eq, Config.floor, ha0] at hf
  obtain ⟨hw0, hlt⟩ := hf
  have hwi : ∀ i, wminLo dd ≤ wfrac dos i := fun i =>
    (finMin1_le (fun i => dd.w i * (1 - dd.dW)) i).trans (hw i)
  apply stable_of_cw hC sp dos hdos (fun _ => 1) (fun _ => one_pos) (Ahi x / wminLo dd) _
    cfg.Tstar hT hlt
  intro i
  have h1 := rowA_le_dos sp dos hdos (wminLo dd) hw0 hwi i
  have h2 : Aiso sp dos / wminLo dd ≤ Ahi x / wminLo dd := div_le_div_of_nonneg_right hA hw0.le
  simp only [mul_one]
  exact h1.trans h2

/-- **N1-ISO.** Any model `MSpec B α`; M-ISO + ENV-SOUND: every true row sum `≤ A_hi`
(for a one-band model this is `A_iso ≤ A_hi`, `rowA_single`). -/
theorem iso_sound (hC : CertH) (ha0 : cfg.a0 = a0Cert) (hT : 0 < cfg.Tstar)
    (h : decideT cfg x = .rejectedProven .iso)
    {B : ℕ} {α : Type*} (sp : MSpec B α) (dos : Fin B → ℝ) (hdos : ∀ i, 0 < dos i)
    (hiso : ∀ i, rowA sp dos i ≤ Ahi x) :
    (∀ μ, CoulPSD dos μ → ∀ Nc : WithTop ℕ, ∀ T : ℝ, cfg.Tstar ≤ T →
      MStable tpPhys sp dos μ Nc T) ∧
    ∃ T', T' < cfg.Tstar ∧ 0 < T' ∧ ∀ μ, CoulPSD dos μ → ∀ Nc : WithTop ℕ, ∀ T : ℝ, T' ≤ T →
      MStable tpPhys sp dos μ Nc T := by
  obtain ⟨-, -, -, hf, -⟩ := decide_rejectedProven h
  simp only [floorR, Bool.and_eq_true, decide_eq_true_eq, Config.floor, ha0] at hf
  apply stable_of_cw hC sp dos hdos (fun _ => 1) (fun _ => one_pos) (Ahi x) _ cfg.Tstar hT
    hf.1.2
  intro i
  simpa [rowA] using hiso i

/-- **NM-CERT** (single isotropic band, M-ISO). ENV-SOUND: `m ≤ m_hi`, `λ ≤ λ_hi`. -/
theorem nmCert_sound (hC : CertH) (ha0 : cfg.a0 = a0Cert) (hT : 0 < cfg.Tstar)
    (h : decideT cfg x = .rejectedProven .nmCert) {ι : Type*} (sp : Spectrum ι)
    (hm : momS sp ≤ mHi x) (hl : lamS sp ≤ lamHi x) :
    ∀ μ : ℝ, 0 ≤ μ → ∀ Nc : WithTop ℕ, ∀ T : ℝ, cfg.Tstar ≤ T →
      PipelineLogic.Stable tpPhys sp μ Nc T := by
  obtain ⟨hw, -, -, hf, -⟩ := decide_rejectedProven h
  simp only [floorR, Bool.and_eq_true, decide_eq_true_eq, Config.floor, ha0] at hf
  exact nmcert_stable tpPhys a0Cert etaCert hC.1 a0Cert_pos etaCert_pos sp (mHi x) (lamHi x)
    cfg.Tstar hT hm hl (lamHi_pos x hw) hf.2

/-- CERT-NM for one certificate: Theorem NM(λ_k, m0) of doc 03 §2.2, "every α²F with
`λ ≤ λ_k` and `m ≤ m0·T²` is stable at `T`" (every `μ* ≥ 0`, every cutoff). A hypothesis
(Arb certificate + paper reduction), like `CertH`. -/
def CertNM (lk m0 : ℝ) : Prop :=
  ∀ {ι : Type} (sp : Spectrum ι), lamS sp ≤ lk → ∀ T : ℝ, 0 < T → momS sp ≤ m0 * T ^ 2 →
    ∀ μ : ℝ, 0 ≤ μ → ∀ Nc : WithTop ℕ, PipelineLogic.Stable tpPhys sp μ Nc T

/-- CERT-NM for the whole certificate set `docs/theory/certificates/nm_*.json`
(`NMTable`): each tabulated `(λ_k, m0_k)`, and the λ = ∞ certificate (any λ). -/
def CertNMSet : Prop :=
  (∀ p ∈ (nmM0Table : List (ℝ × ℝ)), CertNM p.1 p.2) ∧ ∀ lk, CertNM lk (nmM0Inf : ℝ)

/-- The lookup used by the rule lands on a certificate: some `(λ_k, m0)` with `λ_hi ≤ λ_k`,
CERT-NM for it, and `gUp(λ_hi) = m0·T*²`. -/
theorem gUp_cert (hNM : CertNMSet) (hcfg : cfg.tbl = certTable ∧ cfg.gInf = 90000 * nmM0Inf)
    (lh : ℝ) : ∃ lk m0, CertNM lk m0 ∧ lh ≤ lk ∧ cfg.gUp lh = m0 * cfg.Tstar ^ 2 := by
  unfold Config.gUp tableLookup
  rw [hcfg.1, hcfg.2]
  cases hf : (certTable : List (ℝ × ℝ)).find? (fun p => decide (lh ≤ p.1)) with
  | none => exact ⟨lh, nmM0Inf, hNM.2 lh, le_rfl, by ring⟩
  | some p =>
    have hmem : p ∈ (certTable : List (ℝ × ℝ)) := List.mem_of_find?_eq_some hf
    have hp : lh ≤ p.1 := by simpa using List.find?_some hf
    obtain ⟨q, hq, rfl⟩ := List.mem_map.mp hmem
    exact ⟨q.1, q.2, hNM.1 q hq, hp, by simp only; ring⟩

/-- **NM-G** ([THM-C]; single isotropic band, M-ISO). Hypotheses: CERT-NM for the certificate
set, the doc configuration (`cfg.tbl`, `cfg.gInf` = the certified table), `T* > 0`, and
ENV-SOUND `m ≤ m_hi`, `λ ≤ λ_hi`. No monotonicity of g and no Einstein-minimiser numerics are
needed (each certificate covers every `λ ≤ λ_k`). -/
theorem nmG_sound (hT : 0 < cfg.Tstar) (hNM : CertNMSet)
    (hcfg : cfg.tbl = certTable ∧ cfg.gInf = 90000 * nmM0Inf)
    (h : decideT cfg x = .rejectedProven .nmG)
    {ι : Type} (sp : Spectrum ι) (hm : momS sp ≤ mHi x) (hl : lamS sp ≤ lamHi x) :
    ∀ μ : ℝ, 0 ≤ μ → ∀ Nc : WithTop ℕ, ∀ T : ℝ, cfg.Tstar ≤ T →
      PipelineLogic.Stable tpPhys sp μ Nc T := by
  obtain ⟨-, -, -, hf, -⟩ := decide_rejectedProven h
  simp only [floorR, Bool.and_eq_true, decide_eq_true_eq] at hf
  obtain ⟨lk, m0, hc, hlk, hg⟩ := gUp_cert cfg hNM hcfg (lamHi x)
  have hlt : momS sp < m0 * cfg.Tstar ^ 2 := by
    have := hf.1.2; rw [hg] at this; linarith
  have hm0 : 0 ≤ m0 := by
    by_contra hneg
    have : m0 * cfg.Tstar ^ 2 < 0 := mul_neg_of_neg_of_pos (not_le.mp hneg) (by positivity)
    linarith [momS_nonneg sp]
  intro μ hμ Nc T hTT
  apply hc sp (hl.trans hlk) T (by linarith) _ μ hμ Nc
  have hTsq : cfg.Tstar ^ 2 ≤ T ^ 2 := by nlinarith
  nlinarith [mul_le_mul_of_nonneg_left hTsq hm0]

/-! ## The ENV-SOUND hypotheses are satisfiable on every well-formed input (non-vacuity,
RL1 findings 1 and 8): the zero spectrum on `B ≥ 1` bands meets them. -/

/-- The zero `MSpec` on `B` bands. -/
def zeroMSpec (B : ℕ) : MSpec B Unit where
  S := ∅
  Ω := fun _ => 1
  W := fun _ _ _ => 0
  Ω_pos := by intro s hs; simp at hs
  W_nonneg := by intro s hs; simp at hs
  W_symm := by intro s hs; simp at hs

theorem Amat_zero (B : ℕ) (dos : Fin B → ℝ) (i j : Fin B) : Amat (zeroMSpec B) dos i j = 0 := by
  simp [Amat, zeroMSpec]

theorem kappa_env_satisfiable (hw : wf x = true) (md : MatData ℝ) (hm : x.mat = some md) :
    1 ≤ md.B ∧ ∀ i j, Amat (zeroMSpec md.B) (fun _ => 1) i j ≤ calAhi md x.eps i j :=
  ⟨(matWF_spec (hm ▸ wf_mat hw)).1, fun i j => by
    rw [Amat_zero]; exact calAhi_nonneg hw hm i j⟩

theorem iso_env_satisfiable (hw : wf x = true) :
    ∀ i, rowA (zeroMSpec 1) (fun _ => 1) i ≤ Ahi x := by
  intro i
  have h1 : rowA (zeroMSpec 1) (fun _ => 1) i = 0 := by simp [rowA, Amat_zero]
  rw [h1]
  exact (wf_spec hw).1.le.trans (Ahi_ge x hw)

/-- **N1-DOS non-vacuity** (uses the §0.6a gate `Σ_i w_i(1−δ_w) ≤ 1`): on a well-formed input
with DOS data, the sheet model with `N_i = w_i` and the zero spectrum satisfies the N1-DOS
ENV-SOUND hypotheses of `dos_sound`. -/
theorem dos_env_satisfiable (hw : wf x = true) (dd : DosData ℝ) (hd : x.dos = some dd) :
    (∀ i, 0 < dd.w i) ∧ (∀ i, dd.w i * (1 - dd.dW) ≤ wfrac dd.w i) ∧
    Aiso (zeroMSpec dd.n) dd.w ≤ Ahi x := by
  obtain ⟨-, hwi, -, -, hS, hsum⟩ := dosWF_spec' (hd ▸ wf_dos' hw)
  have hS0 : 0 < ∑ i, dd.w i := by linarith
  refine ⟨fun i => (hwi i).1, fun i => ?_, ?_⟩
  · unfold wfrac
    rw [le_div_iff₀ hS0]
    have e : dd.w i * (1 - dd.dW) * ∑ k, dd.w k = dd.w i * ∑ k, dd.w k * (1 - dd.dW) := by
      rw [← sum_mul]; ring
    rw [e]
    have := mul_le_mul_of_nonneg_left hsum (hwi i).1.le
    linarith
  · have h0 : Aiso (zeroMSpec dd.n) dd.w = 0 := by simp [Aiso, zeroMSpec]
    rw [h0]
    exact (wf_spec hw).1.le.trans (Ahi_ge x hw)

/-- The zero single-band spectrum meets the NM-CERT / NM-G ENV-SOUND hypotheses. -/
def zeroSpectrum : Spectrum Unit where
  S := ∅
  p := fun _ => 0
  Ω := fun _ => 1
  p_nonneg := by intro j hj; simp at hj
  Ω_pos := by intro j hj; simp at hj

theorem nm_env_satisfiable (hw : wf x = true) :
    momS zeroSpectrum ≤ mHi x ∧ lamS zeroSpectrum ≤ lamHi x := by
  have h1 : momS zeroSpectrum = 0 := by simp [momS, zeroSpectrum]
  have h2 : lamS zeroSpectrum = 0 := by simp [lamS, zeroSpectrum]
  rw [h1, h2]
  exact ⟨(wf_spec hw).2.2.1.le.trans (mHi_ge x hw), (lamHi_pos x hw).le⟩

end

end TargetLogic
