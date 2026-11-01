/-
MultiBand.Pipeline: the multiband screening rule of wMB §7 in the style of
`PipelineLogic.R_UNIV`: `κ_hi < a0 · T*` ⇒ stable at every `T ≥ T*`.
-/
import MultiBand.Main

open Finset

namespace MultiBand

open TcR1

noncomputable section

variable {B : ℕ} {α β : Type*}

/-- Envelope transfer of a Collatz–Wielandt certificate: if `𝒜 v ≤ κ v` (`v > 0`) for the
nominal spectrum and the true spectrum satisfies `𝒜'_ij ≤ (1+e) 𝒜_ij` entrywise (any atoms,
any frequencies), then `𝒜' v ≤ κ_hi v` whenever `(1+e) κ ≤ κ_hi`, `1 + e ≥ 0`. -/
theorem sub_env (sp : MSpec B α) (sp' : MSpec B β) (dos v : Fin B → ℝ) (hv : ∀ i, 0 < v i)
    (κ e κhi : ℝ) (he : 0 ≤ 1 + e) (hsub : ∀ i, ∑ j, Amat sp dos i j * v j ≤ κ * v i)
    (henv : ∀ i j, Amat sp' dos i j ≤ (1 + e) * Amat sp dos i j) (hκ : (1 + e) * κ ≤ κhi) :
    ∀ i, ∑ j, Amat sp' dos i j * v j ≤ κhi * v i := by
  intro i
  calc ∑ j, Amat sp' dos i j * v j ≤ ∑ j, (1 + e) * Amat sp dos i j * v j :=
        sum_le_sum (fun j _ => mul_le_mul_of_nonneg_right (henv i j) (hv j).le)
    _ = (1 + e) * ∑ j, Amat sp dos i j * v j := by
        rw [mul_sum]; apply sum_congr rfl; intro j _; ring
    _ ≤ (1 + e) * (κ * v i) := mul_le_mul_of_nonneg_left (hsub i) he
    _ = ((1 + e) * κ) * v i := by ring
    _ ≤ κhi * v i := mul_le_mul_of_nonneg_right hκ (hv i).le

/-- Symmetric-gauge certificate (what the pipeline computes): with `r_i = √N_i`
(`r_i > 0`, `r_i² = N_i`), `K_ij = W^tot_ij / (r_i r_j)` (`= √(N_i/N_j) 𝒜_ij`), a positive `u`
with `K u ≤ κ u` gives the positive vector `v_i = u_i / r_i` with `𝒜 v ≤ κ v`. -/
theorem sub_of_K (sp : MSpec B α) (dos r u : Fin B → ℝ) (hr : ∀ i, 0 < r i)
    (hr2 : ∀ i, r i ^ 2 = dos i) (κ : ℝ)
    (hK : ∀ i, ∑ j, (∑ s ∈ sp.S, sp.W s i j) / (r i * r j) * u j ≤ κ * u i) :
    ∀ i, ∑ j, Amat sp dos i j * (u j / r j) ≤ κ * (u i / r i) := by
  intro i
  have hri := hr i
  have e : ∑ j, Amat sp dos i j * (u j / r j)
      = (∑ j, (∑ s ∈ sp.S, sp.W s i j) / (r i * r j) * u j) / r i := by
    rw [div_eq_mul_inv, sum_mul]; apply sum_congr rfl; intro j _
    have hrj := hr j
    unfold Amat; rw [← hr2 i]; field_simp
  rw [e, div_le_iff₀ hri]
  have : κ * (u i / r i) * r i = κ * u i := by field_simp
  linarith [hK i]

/-- **Multiband R-UNIV.** Under `H(a0, η)` (`η ≥ 0`, `a0 > 0`) and Prop. N: let `v > 0` with
`𝒜 v ≤ κ v` for the nominal spectrum, `(1+e) κ ≤ κ_hi` (`1+e ≥ 0`), `κ_hi < a0 T*`, `T* > 0`.
Then `κ_hi / a0 < T*`, and every spectrum `sp'` with `𝒜' ≤ (1+e) 𝒜` entrywise is stable at
every `T ≥ T*`, for every Coulomb matrix with `N_i μ*_ij` PSD and every cutoff. -/
theorem R_UNIV_MB (tp a0 η : ℝ) (hH : H tp a0 η) (hη : 0 ≤ η) (ha0 : 0 < a0)
    (hN : PropN tp a0) (sp : MSpec B α) (dos : Fin B → ℝ) (hdos : ∀ i, 0 < dos i)
    (v : Fin B → ℝ) (hv : ∀ i, 0 < v i) (κ : ℝ)
    (hsub : ∀ i, ∑ j, Amat sp dos i j * v j ≤ κ * v i)
    (e κhi Tstar : ℝ) (he : 0 ≤ 1 + e) (hκhi : (1 + e) * κ ≤ κhi)
    (hU : κhi < a0 * Tstar) (hT0 : 0 < Tstar) :
    κhi / a0 < Tstar ∧
    ∀ sp' : MSpec B β, (∀ i j, Amat sp' dos i j ≤ (1 + e) * Amat sp dos i j) →
      ∀ μ : Fin B → Fin B → ℝ, CoulPSD dos μ → ∀ Nc : WithTop ℕ, ∀ T : ℝ, Tstar ≤ T →
        MStable tp sp' dos μ Nc T := by
  refine ⟨by rw [div_lt_iff₀ ha0]; linarith, ?_⟩
  intro sp' henv μ hμ Nc T hT
  apply mstable_of_gt tp a0 η hH hη ha0 hN sp' dos hdos v hv κhi
    (sub_env sp sp' dos v hv κ e κhi he hsub henv hκhi) μ hμ Nc T (by linarith)
  have : a0 * Tstar ≤ a0 * T := mul_le_mul_of_nonneg_left hT ha0.le
  linarith

/-- **Multiband R-UNIV, Perron-vector form** (no envelope): `𝒜 v = κ v`, `v > 0`,
`κ ≤ κ_hi < a0 T*`, `T* > 0` ⇒ stable at every `T ≥ T*`. -/
theorem R_UNIV_MB_eig (tp a0 η : ℝ) (hH : H tp a0 η) (hη : 0 ≤ η) (ha0 : 0 < a0)
    (hN : PropN tp a0) (sp : MSpec B α) (dos : Fin B → ℝ) (hdos : ∀ i, 0 < dos i)
    (v : Fin B → ℝ) (hv : ∀ i, 0 < v i) (κ : ℝ)
    (heig : ∀ i, ∑ j, Amat sp dos i j * v j = κ * v i)
    (κhi Tstar : ℝ) (hκhi : κ ≤ κhi) (hU : κhi < a0 * Tstar) (hT0 : 0 < Tstar) :
    κhi / a0 < Tstar ∧
    ∀ μ : Fin B → Fin B → ℝ, CoulPSD dos μ → ∀ Nc : WithTop ℕ, ∀ T : ℝ, Tstar ≤ T →
      MStable tp sp dos μ Nc T := by
  obtain ⟨h1, h2⟩ := R_UNIV_MB (β := α) tp a0 η hH hη ha0 hN sp dos hdos v hv κ
    (sub_of_eig sp dos v κ heig) 0 κhi Tstar (by norm_num) (by linarith) hU hT0
  exact ⟨h1, fun μ hμ Nc T hT => h2 sp (fun i j => by linarith) μ hμ Nc T hT⟩

/-- **Fallback rule (wMB §7.3): max row sum.** If every row sum `∑_j 𝒜_ij ≤ κ_hi` and
`κ_hi < a0 T*` (`T* > 0`), the spectrum is stable at every `T ≥ T*` (Collatz–Wielandt with
`v = 1`; no Perron vector needed). -/
theorem R_UNIV_MB_rowsum (tp a0 η : ℝ) (hH : H tp a0 η) (hη : 0 ≤ η) (ha0 : 0 < a0)
    (hN : PropN tp a0) (sp : MSpec B α) (dos : Fin B → ℝ) (hdos : ∀ i, 0 < dos i)
    (κhi Tstar : ℝ) (hrow : ∀ i, ∑ j, Amat sp dos i j ≤ κhi)
    (hU : κhi < a0 * Tstar) (hT0 : 0 < Tstar) :
    κhi / a0 < Tstar ∧
    ∀ μ : Fin B → Fin B → ℝ, CoulPSD dos μ → ∀ Nc : WithTop ℕ, ∀ T : ℝ, Tstar ≤ T →
      MStable tp sp dos μ Nc T := by
  obtain ⟨h1, h2⟩ := R_UNIV_MB (β := α) tp a0 η hH hη ha0 hN sp dos hdos (fun _ => 1)
    (fun _ => one_pos) κhi (fun i => by simpa using hrow i) 0 κhi Tstar (by norm_num)
    (by linarith) hU hT0
  exact ⟨h1, fun μ hμ Nc T hT => h2 sp (fun i j => by linarith) μ hμ Nc T hT⟩

end

end MultiBand
