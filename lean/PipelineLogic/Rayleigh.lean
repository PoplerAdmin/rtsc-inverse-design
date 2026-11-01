/-
P5: Rayleigh-quotient instability certificate, box minimum over the weight envelope,
monotonicity in μ*.
-/
import PipelineLogic.Monotone

open Finset

namespace PipelineLogic

open TcR1

noncomputable section

/-- **P5 (Rayleigh).** If some finitely supported `x ≠ 0` has `Q(x) ≥ 0`, then `Q` is
`TcR1.Unstable` (i.e. `ρ ≥ 0`). -/
theorem rayleigh_unstable (Q : ℕ → (ℕ → ℝ) → ℝ) (N : ℕ) (x : ℕ → ℝ)
    (hx : ∃ n ∈ range N, x n ≠ 0) (hQ : 0 ≤ Q N x) : Unstable Q := by
  intro ε hε
  refine ⟨N, x, hx, ?_⟩
  have := nsq_pos_of_ne hx
  nlinarith

/-- An unstable even sector means `T` is not stable. -/
theorem not_stable_of_unstable {ι : Type*} (κ : ℝ) (sp : Spectrum ι) (μ : ℝ)
    (Nc : WithTop ℕ) (T : ℝ) (h : Unstable (Qp κ sp μ Nc T)) : ¬ Stable κ sp μ Nc T :=
  fun hs => hs.1 h

/-- **P5 (monotone in μ*).** `μ ≤ μ_hi ⇒ Q⁺_{μ_hi}(x) ≤ Q⁺_{μ}(x)` (same cutoff). -/
theorem Qp_anti_mu {ι : Type*} (κ : ℝ) (sp : Spectrum ι) (μ μhi : ℝ) (h : μ ≤ μhi)
    (Nc : WithTop ℕ) (T : ℝ) (N : ℕ) (x : ℕ → ℝ) :
    Qp κ sp μhi Nc T N x ≤ Qp κ sp μ Nc T N x := by
  unfold Qp muTerm
  have : 0 ≤ (∑ n ∈ (range N).filter (fun n : ℕ => ((n : ℕ) : WithTop ℕ) < Nc), x n) ^ 2 :=
    sq_nonneg _
  nlinarith

/-- Sign-chosen endpoint of the box `[lo, hi]` for the linear functional with coefficient `F`:
the lower end where `F ≥ 0`, the upper end where `F < 0`. -/
def sel (F lo hi : ℝ) : ℝ := if 0 ≤ F then lo else hi

theorem sel_mul_le (F lo hi c : ℝ) (hlo : lo ≤ c) (hhi : c ≤ hi) : sel F lo hi * F ≤ c * F := by
  unfold sel
  split_ifs with h
  · exact mul_le_mul_of_nonneg_right hlo h
  · exact mul_le_mul_of_nonpos_right hhi (le_of_lt (not_le.mp h))

theorem sel_mem (F lo hi : ℝ) (h : lo ≤ hi) : lo ≤ sel F lo hi ∧ sel F lo hi ≤ hi := by
  unfold sel; split_ifs <;> constructor <;> linarith

/-- The sign-chosen term is the smaller endpoint term: `sel F lo hi · F = min(lo·F, hi·F)` for
`lo ≤ hi`. (So an Arb lower bound of `∑_j min((1-e)p_j F_j, (1+e)p_j F_j)` is a lower bound of
the box minimum even when the sign of a ball `F_j` is undetermined.) -/
theorem sel_mul_eq_min (F lo hi : ℝ) (h : lo ≤ hi) : sel F lo hi * F = min (lo * F) (hi * F) := by
  unfold sel
  split_ifs with hF
  · exact (min_eq_left (mul_le_mul_of_nonneg_right h hF)).symm
  · exact (min_eq_right (mul_le_mul_of_nonpos_right h (le_of_lt (not_le.mp hF)))).symm

/-- **P5 (box minimum, abstract).** For a linear functional `c ↦ ∑ c_j F_j` on the box
`∏ [lo_j, hi_j]`, the sign-chosen endpoint gives a lower bound … -/
theorem affine_box_min {ι : Type*} (S : Finset ι) (F lo hi c : ι → ℝ)
    (hc : ∀ j ∈ S, lo j ≤ c j ∧ c j ≤ hi j) :
    ∑ j ∈ S, sel (F j) (lo j) (hi j) * F j ≤ ∑ j ∈ S, c j * F j :=
  sum_le_sum (fun j hj => sel_mul_le _ _ _ _ (hc j hj).1 (hc j hj).2)

/-- … which is attained, since the sign-chosen endpoint lies in the box (`lo ≤ hi`). -/
theorem affine_box_min_mem {ι : Type*} (S : Finset ι) (F lo hi : ι → ℝ)
    (h : ∀ j ∈ S, lo j ≤ hi j) :
    ∀ j ∈ S, lo j ≤ sel (F j) (lo j) (hi j) ∧ sel (F j) (lo j) (hi j) ≤ hi j :=
  fun j hj => sel_mem _ _ _ (h j hj)

/-- Per-atom coefficient of the even form: `F_T(ω; x) = F⁺(g_T(·, ω); x)`. -/
def atomF (κ T w : ℝ) (N : ℕ) (x : ℕ → ℝ) : ℝ := Fp (fun k => g κ T k w) N x

/-- Sign-chosen envelope weight of atom `j`. -/
def boxWeight {ι : Type*} (κ : ℝ) (sp : Spectrum ι) (e T : ℝ) (N : ℕ) (x : ℕ → ℝ)
    (j : ι) : ℝ :=
  sel (atomF κ T (sp.Ω j) N x) ((1 - e) * sp.p j) ((1 + e) * sp.p j)

/-- The box minimum of `Q⁺_T(x)` over the weight envelope `e` of `sp` (ARCH C-RAYLEIGH):
`C(x) - 2μ*(∑_{n<N_c} x_n)² + ∑_j w_j F_T(Ω_j; x)` with `w_j = (1-e)p_j` if `F_T(Ω_j;x) ≥ 0`
and `w_j = (1+e)p_j` otherwise. -/
def boxMinQ {ι : Type*} (κ : ℝ) (sp : Spectrum ι) (e μ : ℝ) (Nc : WithTop ℕ) (T : ℝ)
    (N : ℕ) (x : ℕ → ℝ) : ℝ :=
  C N x - muTerm μ Nc N x
    + ∑ j ∈ sp.S, boxWeight κ sp e T N x j * atomF κ T (sp.Ω j) N x

/-- `Q⁺` of an envelope spectrum written over the nominal atoms (TcR1.L1). -/
theorem Qp_env_eq {ι : Type*} (κ e : ℝ) (sp sp' : Spectrum ι) (h : InEnv e sp sp') (μ : ℝ)
    (Nc : WithTop ℕ) (T : ℝ) (N : ℕ) (x : ℕ → ℝ) :
    Qp κ sp' μ Nc T N x = C N x - muTerm μ Nc N x
      + ∑ j ∈ sp.S, sp'.p j * atomF κ T (sp.Ω j) N x := by
  rw [L1]
  obtain ⟨hS, hΩ, _⟩ := h
  rw [hS]
  congr 1
  apply sum_congr rfl; intro j hj
  rw [hΩ j hj]; rfl

/-- **P5 (box-min formula, lower bound).** `boxMinQ ≤ Q⁺_T(x)` for every spectrum in the
envelope. -/
theorem boxMin_le {ι : Type*} (κ e : ℝ) (sp sp' : Spectrum ι) (h : InEnv e sp sp') (μ : ℝ)
    (Nc : WithTop ℕ) (T : ℝ) (N : ℕ) (x : ℕ → ℝ) :
    boxMinQ κ sp e μ Nc T N x ≤ Qp κ sp' μ Nc T N x := by
  rw [Qp_env_eq κ e sp sp' h]
  unfold boxMinQ boxWeight
  have := affine_box_min sp.S (fun j => atomF κ T (sp.Ω j) N x) (fun j => (1 - e) * sp.p j)
    (fun j => (1 + e) * sp.p j) sp'.p h.2.2
  linarith

/-- **P5 (box-min formula, attained).** For `0 ≤ e ≤ 1` the sign-chosen endpoint is a
spectrum of the envelope, and its `Q⁺_T(x)` equals `boxMinQ`. So `boxMinQ` is exactly the
minimum of `Q⁺_T(x)` over the envelope. -/
theorem boxMin_attained {ι : Type*} (κ e : ℝ) (he0 : 0 ≤ e) (he1 : e ≤ 1) (sp : Spectrum ι)
    (μ : ℝ) (Nc : WithTop ℕ) (T : ℝ) (N : ℕ) (x : ℕ → ℝ) :
    ∃ sp' : Spectrum ι, InEnv e sp sp' ∧ Qp κ sp' μ Nc T N x = boxMinQ κ sp e μ Nc T N x := by
  have hnn : ∀ j ∈ sp.S, 0 ≤ boxWeight κ sp e T N x j := by
    intro j hj
    have := sp.p_nonneg j hj
    unfold boxWeight sel
    split_ifs <;> nlinarith
  let sp' : Spectrum ι :=
    { S := sp.S, p := boxWeight κ sp e T N x, Ω := sp.Ω, p_nonneg := hnn,
      Ω_pos := sp.Ω_pos }
  have hin : InEnv e sp sp' := by
    refine ⟨rfl, fun _ _ => rfl, fun j hj => ?_⟩
    have := sp.p_nonneg j hj
    exact sel_mem _ _ _ (by nlinarith)
  refine ⟨sp', hin, ?_⟩
  rw [Qp_env_eq κ e sp sp' hin]
  rfl

/-- **P5 (C-RAYLEIGH soundness).** If `x ≠ 0` and the envelope box minimum of `Q⁺_T(x)` at
`μ* = μ_hi` (cutoff `N_c`) is `≥ 0`, then for every spectrum in the envelope and every
`μ* ≤ μ_hi`, the even sector at `T` (cutoff `N_c`) is unstable. -/
theorem rayleigh_cert {ι : Type*} (κ e μhi : ℝ) (sp : Spectrum ι) (Nc : WithTop ℕ) (T : ℝ)
    (N : ℕ) (x : ℕ → ℝ) (hx : ∃ n ∈ range N, x n ≠ 0)
    (h0 : 0 ≤ boxMinQ κ sp e μhi Nc T N x) :
    ∀ sp' : Spectrum ι, InEnv e sp sp' → ∀ μ : ℝ, μ ≤ μhi → Unstable (Qp κ sp' μ Nc T) := by
  intro sp' hin μ hμ
  apply rayleigh_unstable _ N x hx
  calc (0 : ℝ) ≤ boxMinQ κ sp e μhi Nc T N x := h0
    _ ≤ Qp κ sp' μhi Nc T N x := boxMin_le κ e sp sp' hin μhi Nc T N x
    _ ≤ Qp κ sp' μ Nc T N x := Qp_anti_mu κ sp' μ μhi hμ Nc T N x

end

end PipelineLogic
