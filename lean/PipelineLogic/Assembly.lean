/-
P6: assembly of the stability rules R-UNIV and R-CELLS.
-/
import PipelineLogic.Schur
import PipelineLogic.Monotone
import PipelineLogic.Cover

open Finset

namespace PipelineLogic

open TcR1

noncomputable section

/-- A spectrum of total weight `0` is stable at every `T`, `μ* ≥ 0`, cutoff. -/
theorem stable_of_A_zero {ι : Type*} (κ : ℝ) (sp : Spectrum ι) (hA : sp.A = 0) (μ : ℝ)
    (hμ : 0 ≤ μ) (Nc : WithTop ℕ) (T : ℝ) : Stable κ sp μ Nc T := by
  have hp : ∀ j ∈ sp.S, sp.p j = 0 := (sum_eq_zero_iff_of_nonneg sp.p_nonneg).mp hA
  have hlam : lam κ sp T = fun k => (0 : ℝ) * lam κ sp T k := by
    funext k
    unfold lam
    rw [zero_mul]
    exact sum_eq_zero (fun j hj => by rw [hp j hj, zero_mul])
  apply stable_of_even_bound κ sp T 1 one_pos _ μ hμ Nc
  intro N x
  rw [hlam, Fp_smul, zero_mul, add_zero]
  have := C_le N x
  linarith

/-- MAIN in the form used by the pipeline: under `H(a0, η)` (`a0, η > 0`), any spectrum with
`A ≤ a0 T` is stable at `T` for every `μ* ≥ 0` and every cutoff (TcR1.stable_of_ge, plus the
trivial case `A = 0`). -/
theorem stable_of_A_le {ι : Type*} (κ a0 η : ℝ) (hH : H κ a0 η) (ha0 : 0 < a0) (hη : 0 < η)
    (sp : Spectrum ι) (μ : ℝ) (hμ : 0 ≤ μ) (Nc : WithTop ℕ) (T : ℝ) (hAT : sp.A ≤ a0 * T) :
    Stable κ sp μ Nc T := by
  have hA0 : 0 ≤ sp.A := sum_nonneg sp.p_nonneg
  rcases hA0.lt_or_eq with hA | hA
  · apply stable_of_ge κ a0 η hH ha0 hη sp hA μ hμ Nc T
    rw [div_le_iff₀ ha0]; linarith
  · exact stable_of_A_zero κ sp hA.symm μ hμ Nc T

/-- **P6 (universal, envelope form).** Under `H(a0, η)`: if `(1+e)·A(sp) ≤ A_hi`, every
spectrum in the envelope `e` of `sp` is stable at every `T` with `A_hi ≤ a0 T`, for every
`μ* ≥ 0` and every cutoff. -/
theorem univ_stable {ι : Type*} (κ a0 η : ℝ) (hH : H κ a0 η) (ha0 : 0 < a0) (hη : 0 < η)
    (e : ℝ) (sp : Spectrum ι) (Ahi : ℝ) (hAhi : (1 + e) * sp.A ≤ Ahi)
    (sp' : Spectrum ι) (hin : InEnv e sp sp') (μ : ℝ) (hμ : 0 ≤ μ) (Nc : WithTop ℕ) (T : ℝ)
    (hT : Ahi ≤ a0 * T) : Stable κ sp' μ Nc T :=
  stable_of_A_le κ a0 η hH ha0 hη sp' μ hμ Nc T
    (le_trans (A_env e sp sp' hin) (le_trans hAhi hT))

/-- **P6 / R-UNIV.** Under `H(a0, η)`: `A_hi ≥ (1+e) A` and `A_hi < a0 T*` imply, for every
spectrum in the envelope, every `μ* ≥ 0`, every cutoff and every `T ≥ T*`: Stable(T).
Moreover the threshold `T_b = A_hi / a0` is `< T*` (so `Tc_model ≤ A_hi/a0 < T*`). -/
theorem R_UNIV {ι : Type*} (κ a0 η : ℝ) (hH : H κ a0 η) (ha0 : 0 < a0) (hη : 0 < η)
    (e : ℝ) (sp : Spectrum ι) (Ahi Tstar : ℝ) (hAhi : (1 + e) * sp.A ≤ Ahi)
    (hU : Ahi < a0 * Tstar) :
    Ahi / a0 < Tstar ∧
    ∀ sp' : Spectrum ι, InEnv e sp sp' → ∀ μ : ℝ, 0 ≤ μ → ∀ Nc : WithTop ℕ, ∀ T : ℝ,
      Ahi / a0 ≤ T → Stable κ sp' μ Nc T := by
  refine ⟨by rw [div_lt_iff₀ ha0]; linarith, ?_⟩
  intro sp' hin μ hμ Nc T hT
  apply univ_stable κ a0 η hH ha0 hη e sp Ahi hAhi sp' hin μ hμ Nc T
  rw [div_le_iff₀ ha0] at hT; linarith

/-- Semantic content of one certified R-CELLS cell `[T1, T2]` for the nominal spectrum `sp`
with envelope `e`: `T1 ≥ 0`, and a margin `η > 0` and a positive test vector `v` (head + tail)
such that every truncated row of the interval envelope
`Mhat(hi, lo)`, `hi(k) = (1+e) λ_sp(k; T1)`, `lo(k) = (1-e) λ_sp(k; T2)`,
satisfies `∑_{m<N} Mhat_nm v_m ≤ -η v_n` (`n < N`). By `rows_of_eventually` it suffices to
certify the full infinite row sums. -/
def CellOK {ι : Type*} (κ : ℝ) (sp : Spectrum ι) (e T1 T2 : ℝ) : Prop :=
  0 ≤ T1 ∧ ∃ η : ℝ, 0 < η ∧ ∃ v : ℕ → ℝ, (∀ n, 0 < v n) ∧
    ∀ N, ∀ n < N,
      rowSum (Mhat (fun k => (1 + e) * lam κ sp T1 k) (fun k => (1 - e) * lam κ sp T2 k))
        v N n ≤ -η * v n

/-- **R-CELLS, one cell.** A certified cell `[T1, T2]` makes every spectrum of the envelope
stable at every `T ∈ [T1, T2]`, for every `μ* ≥ 0` and every cutoff. -/
theorem cell_stable {ι : Type*} (κ e T1 T2 : ℝ) (sp : Spectrum ι) (hc : CellOK κ sp e T1 T2)
    (sp' : Spectrum ι) (hin : InEnv e sp sp') (T : ℝ) (h1 : T1 ≤ T) (h2 : T ≤ T2)
    (μ : ℝ) (hμ : 0 ≤ μ) (Nc : WithTop ℕ) : Stable κ sp' μ Nc T := by
  obtain ⟨hT1, η, hη, v, hv, hrow⟩ := hc
  apply stable_of_even_bound κ sp' T η hη _ μ hμ Nc
  intro N x
  exact envelope_form_bound (lam κ sp' T) _ _ (fun k _ => lam_nonneg κ sp' T k)
    (fun k => (lam_box κ e sp sp' hin T1 T2 T hT1 h1 h2 k).1)
    (fun k => (lam_box κ e sp sp' hin T1 T2 T hT1 h1 h2 k).2) v hv η hrow N x

/-- **P6 (cells).** Under `H(a0, η)`: if the certified cells tile `[T_up, t_end]` contiguously
and `(1+e) A ≤ A_hi ≤ a0 · t_end`, then every spectrum of the envelope is stable at every
`T ≥ T_up`, for every `μ* ≥ 0` and every cutoff. -/
theorem cells_stable {ι : Type*} (κ a0 η : ℝ) (hH : H κ a0 η) (ha0 : 0 < a0) (hη : 0 < η)
    (e : ℝ) (sp : Spectrum ι) (Tup tend Ahi : ℝ) (cells : List (ℝ × ℝ))
    (htile : Tiles Tup tend cells) (hcells : ∀ c ∈ cells, CellOK κ sp e c.1 c.2)
    (hAhi : (1 + e) * sp.A ≤ Ahi) (hend : Ahi ≤ a0 * tend)
    (sp' : Spectrum ι) (hin : InEnv e sp sp') (μ : ℝ) (hμ : 0 ≤ μ) (Nc : WithTop ℕ)
    (T : ℝ) (hT : Tup ≤ T) : Stable κ sp' μ Nc T := by
  by_cases h : T ≤ tend
  · obtain ⟨c, hc, h1, h2⟩ := cover cells Tup tend T htile hT h
    exact cell_stable κ e c.1 c.2 sp (hcells c hc) sp' hin T h1 h2 μ hμ Nc
  · apply univ_stable κ a0 η hH ha0 hη e sp Ahi hAhi sp' hin μ hμ Nc T
    have : a0 * tend ≤ a0 * T := mul_le_mul_of_nonneg_left (le_of_lt (not_le.mp h)) ha0.le
    linarith

/-- **R-CELLS.** With `T_up < T*`: `T_b = T_up < T*` and every `T ≥ T_up` is stable
(so `Tc_model ≤ T_up < T*`); in particular every `T ≥ T*`. -/
theorem R_CELLS {ι : Type*} (κ a0 η : ℝ) (hH : H κ a0 η) (ha0 : 0 < a0) (hη : 0 < η)
    (e : ℝ) (sp : Spectrum ι) (Tup tend Ahi Tstar : ℝ) (cells : List (ℝ × ℝ))
    (htile : Tiles Tup tend cells) (hcells : ∀ c ∈ cells, CellOK κ sp e c.1 c.2)
    (hAhi : (1 + e) * sp.A ≤ Ahi) (hend : Ahi ≤ a0 * tend) (hbelow : Tup < Tstar) :
    ∀ sp' : Spectrum ι, InEnv e sp sp' → ∀ μ : ℝ, 0 ≤ μ → ∀ Nc : WithTop ℕ, ∀ T : ℝ,
      Tstar ≤ T → Stable κ sp' μ Nc T :=
  fun sp' hin μ hμ Nc T hT =>
    cells_stable κ a0 η hH ha0 hη e sp Tup tend Ahi cells htile hcells hAhi hend sp' hin μ hμ
      Nc T (le_trans hbelow.le hT)

end

end PipelineLogic
