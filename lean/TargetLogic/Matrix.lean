/-
TargetLogic.Matrix: rules N1-ROWMAX and N1-DOS of `docs/theory/01_eliashberg_targets.md` §0.4.

Part 1 (pure linear algebra, any finite index `Fin B`):
* `le_rowSum_of_super`: a nonnegative matrix `M`, a nonnegative nonzero `v` with
  `κ v ≤ M v` (in particular an eigenvector, `M v = κ v`) ⇒ `κ ≤ rowSum M i` for some `i`.
* `eig_le_of_rowSum_le`: hence `κ ≤ r` whenever every row sum is `≤ r` (κ ≤ max row sum).
* `rowSum_le_weighted_div`: nonnegative row sums `R_i`, weights `w_i ≥ w_min > 0`
  ⇒ `R_i ≤ (∑_k w_k R_k) / w_min`.
* `eig_le_iso_div_wmin`: hence `κ ≤ A_iso / w_min` with `A_iso = ∑_i w_i R_i`.

Part 2 (connection to `MultiBand`): with `𝒜 = MultiBand.Amat sp dos` (the matrix whose
Collatz–Wielandt / Perron root is the `κ` of `MultiBand.mstable_of_gt`), row sums
`rowA i = ∑_j 𝒜_ij`, DOS fractions `wfrac i = N_i / ∑_k N_k` and
`Aiso = (∑_ij W^tot_ij) / ∑_k N_k` (the `A` of the DOS-averaged α²F_iso = ∑_ij (N_i/N) α²F_ij):
* `Aiso_eq`: `Aiso = ∑_i wfrac_i · rowA_i`;
* `n1_rowmax_stable`, `n1_dos_stable`: `bound < a0·T*` ⇒ stable at every `T ≥ T*`;
* `kappa_le_rowmax`, `kappa_le_dos`: the Perron root is below both bounds;
* `tc_le_rowmax`, `tc_le_dos`: an unstable `T > 0` satisfies `T ≤ bound / a0`.
-/
import MultiBand
import Mathlib.Algebra.BigOperators.Field
import Mathlib.Algebra.BigOperators.Fin

open Finset

namespace TargetLogic

noncomputable section

/-! ## Part 1: linear algebra -/

section LinAlg

variable {B : ℕ}

/-- Row sum `∑_j M_ij`. -/
def rowSum (M : Fin B → Fin B → ℝ) (i : Fin B) : ℝ := ∑ j, M i j

/-- **κ ≤ some row sum.** `M ≥ 0` entrywise, `v ≥ 0`, `v ≠ 0`, `κ v_i ≤ (M v)_i` for every
`i` (e.g. `M v = κ v`). Then `κ ≤ ∑_j M_ij` at an index `i` where `v` is maximal. -/
theorem le_rowSum_of_super (M : Fin B → Fin B → ℝ) (hM : ∀ i j, 0 ≤ M i j)
    (v : Fin B → ℝ) (hv : ∀ i, 0 ≤ v i) (hv0 : ∃ i, v i ≠ 0) (κ : ℝ)
    (hκ : ∀ i, κ * v i ≤ ∑ j, M i j * v j) :
    ∃ i, κ ≤ rowSum M i := by
  obtain ⟨i0, hi0⟩ := hv0
  have hne : (univ : Finset (Fin B)).Nonempty := ⟨i0, mem_univ _⟩
  obtain ⟨i, -, hmax⟩ := exists_max_image univ v hne
  have hvi0 : 0 < v i0 := lt_of_le_of_ne (hv i0) (Ne.symm hi0)
  have hvi : 0 < v i := lt_of_lt_of_le hvi0 (hmax i0 (mem_univ _))
  refine ⟨i, ?_⟩
  have h1 : ∑ j, M i j * v j ≤ ∑ j, M i j * v i :=
    sum_le_sum (fun j _ => mul_le_mul_of_nonneg_left (hmax j (mem_univ _)) (hM i j))
  have h2 : ∑ j, M i j * v i = rowSum M i * v i := by rw [rowSum, sum_mul]
  have h3 : κ * v i ≤ rowSum M i * v i := by linarith [hκ i]
  by_contra h
  have := mul_lt_mul_of_pos_right (not_le.mp h) hvi
  linarith

/-- **N1-ROWMAX (algebra).** An eigenvalue `κ` of a nonnegative matrix with a nonnegative
nonzero eigenvector is at most any common upper bound `r` of the row sums. -/
theorem eig_le_of_rowSum_le (M : Fin B → Fin B → ℝ) (hM : ∀ i j, 0 ≤ M i j)
    (v : Fin B → ℝ) (hv : ∀ i, 0 ≤ v i) (hv0 : ∃ i, v i ≠ 0) (κ : ℝ)
    (heig : ∀ i, ∑ j, M i j * v j = κ * v i) (r : ℝ) (hrow : ∀ i, rowSum M i ≤ r) :
    κ ≤ r := by
  obtain ⟨i, hi⟩ := le_rowSum_of_super M hM v hv hv0 κ (fun i => (heig i).ge)
  exact hi.trans (hrow i)

/-- **N1-DOS (algebra).** Nonnegative `R`, weights `w_i ≥ w_min > 0`:
`R_i ≤ (∑_k w_k R_k) / w_min` for every `i`. (No normalisation `∑ w = 1` is needed.) -/
theorem rowSum_le_weighted_div (R w : Fin B → ℝ) (hR : ∀ i, 0 ≤ R i) (wmin : ℝ)
    (hw0 : 0 < wmin) (hw : ∀ i, wmin ≤ w i) (i : Fin B) :
    R i ≤ (∑ k, w k * R k) / wmin := by
  rw [le_div_iff₀ hw0]
  have h1 : R i * wmin ≤ w i * R i := by nlinarith [hR i, hw i]
  have h2 : w i * R i ≤ ∑ k, w k * R k :=
    single_le_sum (f := fun k => w k * R k)
      (fun k _ => mul_nonneg (hw0.le.trans (hw k)) (hR k)) (mem_univ i)
  linarith

/-- **N1-DOS (κ form).** `M ≥ 0`, nonnegative nonzero eigenvector for `κ`, weights
`w_i ≥ w_min > 0`: `κ ≤ (∑_i w_i rowSum_i) / w_min`. -/
theorem eig_le_iso_div_wmin (M : Fin B → Fin B → ℝ) (hM : ∀ i j, 0 ≤ M i j)
    (v : Fin B → ℝ) (hv : ∀ i, 0 ≤ v i) (hv0 : ∃ i, v i ≠ 0) (κ : ℝ)
    (heig : ∀ i, ∑ j, M i j * v j = κ * v i) (w : Fin B → ℝ) (wmin : ℝ) (hw0 : 0 < wmin)
    (hw : ∀ i, wmin ≤ w i) :
    κ ≤ (∑ k, w k * rowSum M k) / wmin :=
  eig_le_of_rowSum_le M hM v hv hv0 κ heig _
    (rowSum_le_weighted_div (rowSum M) w (fun i => sum_nonneg (fun j _ => hM i j)) wmin hw0 hw)

end LinAlg

/-! ## Part 2: connection to `MultiBand` -/

section MB

open MultiBand TcR1

variable {B : ℕ} {α : Type*}

/-- Row sum of the coupling-weight matrix: `A_i = ∑_j 𝒜_ij` (the YAML input of N1-ROWMAX). -/
def rowA (sp : MSpec B α) (dos : Fin B → ℝ) (i : Fin B) : ℝ := ∑ j, Amat sp dos i j

/-- DOS fraction `w_i = N_i / N`, `N = ∑_k N_k`. -/
def wfrac (dos : Fin B → ℝ) (i : Fin B) : ℝ := dos i / ∑ k, dos k

/-- `A_iso = ∫ α²F_iso`, `α²F_iso = ∑_ij (N_i/N) α²F_ij = (1/N) ∑_ij W_ij`. -/
def Aiso (sp : MSpec B α) (dos : Fin B → ℝ) : ℝ :=
  (∑ i, ∑ j, ∑ s ∈ sp.S, sp.W s i j) / ∑ k, dos k

theorem Amat_nonneg (sp : MSpec B α) (dos : Fin B → ℝ) (hdos : ∀ i, 0 < dos i) (i j : Fin B) :
    0 ≤ Amat sp dos i j :=
  div_nonneg (sum_nonneg (fun s hs => sp.W_nonneg s hs i j)) (hdos i).le

theorem rowA_nonneg (sp : MSpec B α) (dos : Fin B → ℝ) (hdos : ∀ i, 0 < dos i) (i : Fin B) :
    0 ≤ rowA sp dos i :=
  sum_nonneg (fun j _ => Amat_nonneg sp dos hdos i j)

/-- `A_iso = ∑_i w_i A_i`: the DOS-weighted mean of the row sums of `𝒜`. -/
theorem Aiso_eq (sp : MSpec B α) (dos : Fin B → ℝ) (hdos : ∀ i, 0 < dos i) :
    Aiso sp dos = ∑ i, wfrac dos i * rowA sp dos i := by
  unfold Aiso wfrac rowA
  rw [sum_div]
  apply sum_congr rfl
  intro i _
  rw [sum_div, mul_sum]
  apply sum_congr rfl
  intro j _
  unfold Amat
  rw [div_mul_div_comm, mul_comm (dos i), mul_div_mul_right _ _ (ne_of_gt (hdos i))]

/-- With a single band (`B = 1`) the row sum is `A_iso` (this is the content of the recorded
assumption M-ISO for a one-sheet material). -/
theorem rowA_single (sp : MSpec 1 α) (dos : Fin 1 → ℝ) (i : Fin 1) :
    rowA sp dos i = Aiso sp dos := by
  have hi : i = 0 := Subsingleton.elim _ _
  subst hi
  simp [rowA, Aiso, Amat]

/-- Every row sum is `≤ A_iso / w_min` when every DOS fraction is `≥ w_min > 0`. -/
theorem rowA_le_dos (sp : MSpec B α) (dos : Fin B → ℝ) (hdos : ∀ i, 0 < dos i)
    (wmin : ℝ) (hw0 : 0 < wmin) (hw : ∀ i, wmin ≤ wfrac dos i) (i : Fin B) :
    rowA sp dos i ≤ Aiso sp dos / wmin := by
  rw [Aiso_eq sp dos hdos]
  exact rowSum_le_weighted_div (rowA sp dos) (wfrac dos) (rowA_nonneg sp dos hdos) wmin hw0 hw i

/-- **κ ≤ max row sum for `𝒜`.** A nonnegative nonzero eigenvector `𝒜 v = κ v` gives
`κ ≤ r` for every common bound `r` on the row sums `A_i`. -/
theorem kappa_le_rowmax (sp : MSpec B α) (dos : Fin B → ℝ) (hdos : ∀ i, 0 < dos i)
    (v : Fin B → ℝ) (hv : ∀ i, 0 ≤ v i) (hv0 : ∃ i, v i ≠ 0) (κ : ℝ)
    (heig : ∀ i, ∑ j, Amat sp dos i j * v j = κ * v i) (r : ℝ) (hrow : ∀ i, rowA sp dos i ≤ r) :
    κ ≤ r :=
  eig_le_of_rowSum_le (Amat sp dos) (Amat_nonneg sp dos hdos) v hv hv0 κ heig r hrow

/-- **κ ≤ A_iso / w_min for `𝒜`** (proof of §5.1 (3), by the row-sum route). -/
theorem kappa_le_dos (sp : MSpec B α) (dos : Fin B → ℝ) (hdos : ∀ i, 0 < dos i)
    (v : Fin B → ℝ) (hv : ∀ i, 0 ≤ v i) (hv0 : ∃ i, v i ≠ 0) (κ : ℝ)
    (heig : ∀ i, ∑ j, Amat sp dos i j * v j = κ * v i)
    (wmin : ℝ) (hw0 : 0 < wmin) (hw : ∀ i, wmin ≤ wfrac dos i) :
    κ ≤ Aiso sp dos / wmin :=
  kappa_le_rowmax sp dos hdos v hv hv0 κ heig _ (rowA_le_dos sp dos hdos wmin hw0 hw)

/-- **N1-ROWMAX (composite).** Under `H(a0, η)`, Prop. N: if every row sum `A_i ≤ A_rowmax,hi`
and `A_rowmax,hi < a0 T*` (`T* > 0`), then `A_rowmax,hi / a0 < T*` and the material is stable
at every `T ≥ T*`, for every PSD Coulomb matrix and every cutoff. -/
theorem n1_rowmax_stable (tp a0 η : ℝ) (hH : H tp a0 η) (hη : 0 ≤ η) (ha0 : 0 < a0)
    (hN : PropN tp a0) (sp : MSpec B α) (dos : Fin B → ℝ) (hdos : ∀ i, 0 < dos i)
    (rhi Tstar : ℝ) (hrow : ∀ i, rowA sp dos i ≤ rhi) (hU : rhi < a0 * Tstar) (hT0 : 0 < Tstar) :
    rhi / a0 < Tstar ∧
    ∀ μ : Fin B → Fin B → ℝ, CoulPSD dos μ → ∀ Nc : WithTop ℕ, ∀ T : ℝ, Tstar ≤ T →
      MStable tp sp dos μ Nc T :=
  R_UNIV_MB_rowsum tp a0 η hH hη ha0 hN sp dos hdos rhi Tstar hrow hU hT0

/-- **N1-DOS (composite).** Under `H(a0, η)`, Prop. N: if `A_iso ≤ A_iso,hi`, every DOS fraction
`w_i ≥ w_min > 0` and `A_iso,hi / w_min < a0 T*` (`T* > 0`), the material is stable at every
`T ≥ T*`, for every PSD Coulomb matrix and every cutoff. -/
theorem n1_dos_stable (tp a0 η : ℝ) (hH : H tp a0 η) (hη : 0 ≤ η) (ha0 : 0 < a0)
    (hN : PropN tp a0) (sp : MSpec B α) (dos : Fin B → ℝ) (hdos : ∀ i, 0 < dos i)
    (AisoHi wmin Tstar : ℝ) (hw0 : 0 < wmin) (hw : ∀ i, wmin ≤ wfrac dos i)
    (hA : Aiso sp dos ≤ AisoHi) (hU : AisoHi / wmin < a0 * Tstar) (hT0 : 0 < Tstar) :
    AisoHi / wmin / a0 < Tstar ∧
    ∀ μ : Fin B → Fin B → ℝ, CoulPSD dos μ → ∀ Nc : WithTop ℕ, ∀ T : ℝ, Tstar ≤ T →
      MStable tp sp dos μ Nc T := by
  have hdiv : Aiso sp dos / wmin ≤ AisoHi / wmin := div_le_div_of_nonneg_right hA hw0.le
  exact n1_rowmax_stable tp a0 η hH hη ha0 hN sp dos hdos (AisoHi / wmin) Tstar
    (fun i => (rowA_le_dos sp dos hdos wmin hw0 hw i).trans hdiv) hU hT0

/-- **Composite κ statement.** Positive Perron vector `𝒜 v = κ v`, any bound `κ ≤ κ_bound`
(e.g. from `kappa_le_rowmax` / `kappa_le_dos`) with `κ_bound < a0 T*` ⇒ stable at every
`T ≥ T*`. -/
theorem stable_of_kappa_bound (tp a0 η : ℝ) (hH : H tp a0 η) (hη : 0 ≤ η) (ha0 : 0 < a0)
    (hN : PropN tp a0) (sp : MSpec B α) (dos : Fin B → ℝ) (hdos : ∀ i, 0 < dos i)
    (v : Fin B → ℝ) (hv : ∀ i, 0 < v i) (κ : ℝ)
    (heig : ∀ i, ∑ j, Amat sp dos i j * v j = κ * v i)
    (κb Tstar : ℝ) (hκb : κ ≤ κb) (hU : κb < a0 * Tstar) (hT0 : 0 < Tstar) :
    ∀ μ : Fin B → Fin B → ℝ, CoulPSD dos μ → ∀ Nc : WithTop ℕ, ∀ T : ℝ, Tstar ≤ T →
      MStable tp sp dos μ Nc T :=
  (R_UNIV_MB_eig tp a0 η hH hη ha0 hN sp dos hdos v hv κ heig κb Tstar hκb hU hT0).2

/-- **Tc ≤ A_rowmax / a0.** Any `T > 0` that is unstable satisfies `T ≤ r / a0` when every row
sum is `≤ r`. -/
theorem tc_le_rowmax (tp a0 η : ℝ) (hH : H tp a0 η) (hη : 0 ≤ η) (ha0 : 0 < a0)
    (hN : PropN tp a0) (sp : MSpec B α) (dos : Fin B → ℝ) (hdos : ∀ i, 0 < dos i)
    (r : ℝ) (hrow : ∀ i, rowA sp dos i ≤ r)
    (μ : Fin B → Fin B → ℝ) (hU : CoulPSD dos μ) (Nc : WithTop ℕ) (T : ℝ) (hT : 0 < T)
    (hUn : ¬ MStable tp sp dos μ Nc T) : T ≤ r / a0 :=
  unstable_imp_le tp a0 η hH hη ha0 hN sp dos hdos (fun _ => 1) (fun _ => one_pos) r
    (fun i => by simpa [rowA] using hrow i) μ hU Nc T hT hUn

/-- **Tc ≤ A_iso / (w_min a0).** -/
theorem tc_le_dos (tp a0 η : ℝ) (hH : H tp a0 η) (hη : 0 ≤ η) (ha0 : 0 < a0)
    (hN : PropN tp a0) (sp : MSpec B α) (dos : Fin B → ℝ) (hdos : ∀ i, 0 < dos i)
    (wmin : ℝ) (hw0 : 0 < wmin) (hw : ∀ i, wmin ≤ wfrac dos i)
    (μ : Fin B → Fin B → ℝ) (hU : CoulPSD dos μ) (Nc : WithTop ℕ) (T : ℝ) (hT : 0 < T)
    (hUn : ¬ MStable tp sp dos μ Nc T) : T ≤ Aiso sp dos / wmin / a0 :=
  tc_le_rowmax tp a0 η hH hη ha0 hN sp dos hdos _ (rowA_le_dos sp dos hdos wmin hw0 hw)
    μ hU Nc T hT hUn

end MB

end

end TargetLogic
