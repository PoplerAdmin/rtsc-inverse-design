/-
MultiBand.Main: Theorem 1 of `wMB/PROOF.md` (even sector, `k_B T_c ≤ κ / a0`) and the
two-sector stability statement.

Hypotheses (as in wMB §4.3):
* `H tp a0 η` (TcR1; Prop. 6 of final/PROOF.md), `η ≥ 0`, `a0 > 0`;
* `PropN tp a0` (Prop. N of wMB §4.2): `Q^E_{-a0,w}(y) ≤ 0` for all `w > 0`, `y`;
* `CoulPSD dos μ`: `U_ij = N_i μ*_ij` PSD;
* a positive vector `v` with `𝒜 v ≤ κ v` (Collatz–Wielandt form; the Perron vector with
  `𝒜 v = κ v` is the special case `sub_of_eig`).
-/
import MultiBand.Odd

open Finset

namespace MultiBand

open TcR1

noncomputable section

variable {B : ℕ} {α : Type*}

/-- The Perron eigen-equation implies the Collatz–Wielandt inequality. -/
theorem sub_of_eig (sp : MSpec B α) (dos v : Fin B → ℝ) (κ : ℝ)
    (heig : ∀ i, ∑ j, Amat sp dos i j * v j = κ * v i) :
    ∀ i, ∑ j, Amat sp dos i j * v j ≤ κ * v i := fun i => (heig i).le

/-- wMB §4.3, first bullet: if `Q^E_{±a0,w}(y) ≤ 0` and `0 ≤ a ≤ a0`, then
`Q^E_{±a,w}(y) ≤ -(1 - a/a0) C(y)`. -/
theorem QE_pm (tp a0 a w : ℝ) (ha0 : 0 < a0) (ha : 0 ≤ a) (N : ℕ) (y : ℕ → ℝ)
    (hp : Qe tp a0 w N y ≤ 0) (hm : Qe tp (-a0) w N y ≤ 0) :
    Qe tp a w N y ≤ -(1 - a / a0) * Cw N y ∧ Qe tp (-a) w N y ≤ -(1 - a / a0) * Cw N y := by
  rw [Qe_eq] at hp hm ⊢
  rw [Qe_eq]
  set G := Gp tp w N y y
  set c := Cw N y
  have hq : 0 ≤ a / a0 := div_nonneg ha ha0.le
  have e : a / a0 * a0 = a := by field_simp
  have h1 := mul_le_mul_of_nonneg_left (show a0 * G ≤ c by linarith) hq
  have h2 := mul_le_mul_of_nonneg_left (show -a0 * G ≤ c by linarith) hq
  constructor
  · have : a * G ≤ a / a0 * c := by
      calc a * G = a / a0 * (a0 * G) := by rw [← mul_assoc, e]
        _ ≤ a / a0 * c := h1
    linarith
  · have : -a * G ≤ a / a0 * c := by
      calc -a * G = -(a / a0 * a0) * G := by rw [e]
        _ = a / a0 * (-a0 * G) := by ring
        _ ≤ a / a0 * c := h2
    linarith

/-- `∑_i N_i C(φ_i) ≥ |x|²`. -/
theorem nrm2_le_X (dos : Fin B → ℝ) (hdos : ∀ i, 0 ≤ dos i) (N : ℕ) (φ : Fin B → ℕ → ℝ) :
    nrm2 dos N φ ≤ ∑ i, dos i * Cw N (φ i) :=
  sum_le_sum (fun i _ => mul_le_mul_of_nonneg_left (Cw_ge_nsq N _) (hdos i))

/-- **Theorem 1 (even sector), Collatz–Wielandt form.**
Under `H(a0, η)` (`η ≥ 0`, `a0 > 0`), Prop. N, `N_i > 0`, `v > 0` with `𝒜 v ≤ κ v`,
`N_i μ*_ij` PSD, any cutoff `N_c`, `T > 0` and `κ ≤ a0 T`:
`Q⁺(φ) ≤ -(1 - κ/(a0 T)) |x|²` for every finitely supported `φ`. -/
theorem thm1_even (tp a0 η : ℝ) (hH : H tp a0 η) (hη : 0 ≤ η) (ha0 : 0 < a0)
    (hN : PropN tp a0) (sp : MSpec B α) (dos : Fin B → ℝ) (hdos : ∀ i, 0 < dos i)
    (v : Fin B → ℝ) (hv : ∀ i, 0 < v i) (κ : ℝ)
    (hsub : ∀ i, ∑ j, Amat sp dos i j * v j ≤ κ * v i)
    (μ : Fin B → Fin B → ℝ) (hU : CoulPSD dos μ) (Nc : WithTop ℕ) (T : ℝ) (hT : 0 < T)
    (hκT : κ ≤ a0 * T) (N : ℕ) (φ : Fin B → ℕ → ℝ) :
    Qplus tp sp dos μ Nc T N φ ≤ -(1 - κ / (a0 * T)) * nrm2 dos N φ := by
  rw [lemma11_even]
  have hv' : ∀ i, v i ≠ 0 := fun i => (hv i).ne'
  set X := ∑ i, dos i * Cw N (φ i) with hX
  set Y := ∑ s ∈ sp.S, ∑ i, ∑ j, sp.W s i j * Gp tp (sp.Ω s / T) N (φ i) (φ j) with hYdef
  have hG : Y = ∑ s ∈ sp.S, ∑ i, ∑ j, mw sp v s i j *
            (Gp tp (sp.Ω s / T) N (psi v φ i j) (psi v φ i j)
              - Gp tp (sp.Ω s / T) N (chi v φ i j) (chi v φ i j)) := by
    apply sum_congr rfl; intro s _; apply sum_congr rfl; intro i _
    apply sum_congr rfl; intro j _
    exact gauge_G tp _ sp v hv' N φ s i j
  have hb : ∀ s ∈ sp.S, ∀ i j, a0 * (mw sp v s i j *
            (Gp tp (sp.Ω s / T) N (psi v φ i j) (psi v φ i j)
              - Gp tp (sp.Ω s / T) N (chi v φ i j) (chi v φ i j)))
        ≤ mw sp v s i j * (Cw N (psi v φ i j) + Cw N (chi v φ i j)) := by
    intro s hs i j
    have hw : 0 < sp.Ω s / T := div_pos (sp.Ω_pos s hs) hT
    have hψ := hH _ hw N (psi v φ i j)
    have hχ := hN _ hw N (chi v φ i j)
    rw [Qe_eq] at hψ hχ
    have hη' : 0 ≤ η * nsq N (psi v φ i j) := mul_nonneg hη (nsq_nonneg N _)
    have key : a0 * (Gp tp (sp.Ω s / T) N (psi v φ i j) (psi v φ i j)
              - Gp tp (sp.Ω s / T) N (chi v φ i j) (chi v φ i j))
        ≤ Cw N (psi v φ i j) + Cw N (chi v φ i j) := by linarith
    have := mul_le_mul_of_nonneg_left key (mw_nonneg sp v hv s hs i j)
    linarith
  have hY : a0 * Y ≤ κ * X := by
    rw [hG, mul_sum]
    calc ∑ s ∈ sp.S, a0 * ∑ i, ∑ j, mw sp v s i j *
            (Gp tp (sp.Ω s / T) N (psi v φ i j) (psi v φ i j)
              - Gp tp (sp.Ω s / T) N (chi v φ i j) (chi v φ i j))
        ≤ ∑ s ∈ sp.S, ∑ i, ∑ j, mw sp v s i j * (Cw N (psi v φ i j) + Cw N (chi v φ i j)) := by
          apply sum_le_sum; intro s hs
          rw [mul_sum]; apply sum_le_sum; intro i _
          rw [mul_sum]; apply sum_le_sum; intro j _
          exact hb s hs i j
      _ ≤ κ * X := gauge_C_sub sp dos hdos v hv κ hsub N φ
  have hYT : (1 / T) * Y ≤ κ / (a0 * T) * X := by
    have h1 : Y ≤ κ * X / a0 := by rw [le_div_iff₀ ha0]; linarith
    calc (1 / T) * Y ≤ (1 / T) * (κ * X / a0) :=
          mul_le_mul_of_nonneg_left h1 (by positivity)
      _ = κ / (a0 * T) * X := by field_simp
  have hP := Pterm_nonneg tp sp T hT N φ
  have hC := coul_nonneg dos μ hU Nc N φ
  have hc : 0 ≤ 1 - κ / (a0 * T) := by
    have : κ / (a0 * T) ≤ 1 := by rw [div_le_one (by positivity)]; exact hκT
    linarith
  have hXn := mul_le_mul_of_nonneg_left (nrm2_le_X dos (fun i => (hdos i).le) N φ) hc
  nlinarith

/-- **Theorem 1 (even sector), literally along wMB §4.3**: through the exact identity of
Lemma 4.1 (Perron vector `𝒜 v = κ v`, `κ > 0`) and `Q^E_{±a,w} ≤ -(1-a/a0) C`. -/
theorem thm1_even_via_identity (tp a0 η : ℝ) (hH : H tp a0 η) (hη : 0 ≤ η) (ha0 : 0 < a0)
    (hN : PropN tp a0) (sp : MSpec B α) (dos : Fin B → ℝ) (hdos : ∀ i, 0 < dos i)
    (v : Fin B → ℝ) (hv : ∀ i, 0 < v i) (κ : ℝ) (hκ : 0 < κ)
    (heig : ∀ i, ∑ j, Amat sp dos i j * v j = κ * v i)
    (μ : Fin B → Fin B → ℝ) (hU : CoulPSD dos μ) (Nc : WithTop ℕ) (T : ℝ) (hT : 0 < T)
    (hκT : κ ≤ a0 * T) (N : ℕ) (φ : Fin B → ℕ → ℝ) :
    Qplus tp sp dos μ Nc T N φ ≤ -(1 - κ / (a0 * T)) * nrm2 dos N φ := by
  have hv' : ∀ i, v i ≠ 0 := fun i => (hv i).ne'
  have hd' : ∀ i, dos i ≠ 0 := fun i => (hdos i).ne'
  rw [lemma41 tp sp dos μ Nc T N φ v κ hκ.ne' hT.ne' hv' hd' heig]
  set a := κ / T with ha_def
  have ha : 0 ≤ a := div_nonneg hκ.le hT.le
  have haa : a / a0 = κ / (a0 * T) := by rw [ha_def]; field_simp
  set c := 1 - κ / (a0 * T) with hc_def
  have hterm : ∀ s ∈ sp.S, ∀ i j, mw sp v s i j *
        (Qe tp a (sp.Ω s / T) N (psi v φ i j) + Qe tp (-a) (sp.Ω s / T) N (chi v φ i j))
      ≤ -c * (mw sp v s i j * (Cw N (psi v φ i j) + Cw N (chi v φ i j))) := by
    intro s hs i j
    have hw : 0 < sp.Ω s / T := div_pos (sp.Ω_pos s hs) hT
    have hp : ∀ y, Qe tp a0 (sp.Ω s / T) N y ≤ 0 := fun y =>
      le_trans (hH _ hw N y) (by nlinarith [nsq_nonneg N y])
    have h1 := (QE_pm tp a0 a _ ha0 ha N (psi v φ i j) (hp _) (hN _ hw N _)).1
    have h2 := (QE_pm tp a0 a _ ha0 ha N (chi v φ i j) (hp _) (hN _ hw N _)).2
    rw [haa] at h1 h2
    have := mul_le_mul_of_nonneg_left (add_le_add h1 h2) (mw_nonneg sp v hv s hs i j)
    linarith
  have hsum : ∑ s ∈ sp.S, ∑ i, ∑ j, mw sp v s i j *
        (Qe tp a (sp.Ω s / T) N (psi v φ i j) + Qe tp (-a) (sp.Ω s / T) N (chi v φ i j))
      ≤ -c * (κ * ∑ i, dos i * Cw N (φ i)) := by
    rw [← gauge_C_eig sp dos hd' v hv' κ heig N φ, mul_sum]
    apply sum_le_sum; intro s hs
    rw [mul_sum]; apply sum_le_sum; intro i _
    rw [mul_sum]; apply sum_le_sum; intro j _
    exact hterm s hs i j
  have hk : (1 / κ) * (-c * (κ * ∑ i, dos i * Cw N (φ i))) = -c * ∑ i, dos i * Cw N (φ i) := by
    field_simp
  have h3 := mul_le_mul_of_nonneg_left hsum (show 0 ≤ 1 / κ by positivity)
  rw [hk] at h3
  have hP := Pterm_nonneg tp sp T hT N φ
  have hC := coul_nonneg dos μ hU Nc N φ
  have hc : 0 ≤ c := by
    have : κ / (a0 * T) ≤ 1 := by rw [div_le_one (by positivity)]; exact hκT
    linarith
  have hXn := mul_le_mul_of_nonneg_left (nrm2_le_X dos (fun i => (hdos i).le) N φ) hc
  nlinarith

/-- A multiband form bounded by `-c|x|²` with `c > 0` is not unstable. -/
theorem not_munstable_of_bound (dos : Fin B → ℝ) (Q : ℕ → (Fin B → ℕ → ℝ) → ℝ) (c : ℝ)
    (hc : 0 < c) (hQ : ∀ N φ, Q N φ ≤ -c * nrm2 dos N φ) : ¬ MUnstable dos Q := by
  intro hU
  obtain ⟨N, φ, _, hφ⟩ := hU c hc
  linarith [hQ N φ]

/-- **Theorem 1 + Theorem 3.** Under the hypotheses of `thm1_even`, every `T > 0` with
`κ < a0 T` (i.e. `T > κ/a0`) is stable in both sectors: `Q⁺ ≤ -(1 - κ/(a0 T))|x|²` with
`1 - κ/(a0 T) > 0`, and `Q⁻ ≤ -|x|²`. Hence `T_c ≤ κ / a0`. -/
theorem mstable_of_gt (tp a0 η : ℝ) (hH : H tp a0 η) (hη : 0 ≤ η) (ha0 : 0 < a0)
    (hN : PropN tp a0) (sp : MSpec B α) (dos : Fin B → ℝ) (hdos : ∀ i, 0 < dos i)
    (v : Fin B → ℝ) (hv : ∀ i, 0 < v i) (κ : ℝ)
    (hsub : ∀ i, ∑ j, Amat sp dos i j * v j ≤ κ * v i)
    (μ : Fin B → Fin B → ℝ) (hU : CoulPSD dos μ) (Nc : WithTop ℕ) (T : ℝ) (hT : 0 < T)
    (hκT : κ < a0 * T) : MStable tp sp dos μ Nc T := by
  constructor
  · apply not_munstable_of_bound dos _ (1 - κ / (a0 * T))
    · have : κ / (a0 * T) < 1 := by rw [div_lt_one (by positivity)]; exact hκT
      linarith
    · exact thm1_even tp a0 η hH hη ha0 hN sp dos hdos v hv κ hsub μ hU Nc T hT hκT.le
  · apply not_munstable_of_bound dos _ 1 one_pos
    intro N φ
    have := thm3_odd tp sp dos (fun i => (hdos i).le) T hT N φ
    linarith

/-- **Corollary (`T_c ≤ κ/a0`).** If `T > 0` is unstable in either sector then `T ≤ κ/a0`. -/
theorem unstable_imp_le (tp a0 η : ℝ) (hH : H tp a0 η) (hη : 0 ≤ η) (ha0 : 0 < a0)
    (hN : PropN tp a0) (sp : MSpec B α) (dos : Fin B → ℝ) (hdos : ∀ i, 0 < dos i)
    (v : Fin B → ℝ) (hv : ∀ i, 0 < v i) (κ : ℝ)
    (hsub : ∀ i, ∑ j, Amat sp dos i j * v j ≤ κ * v i)
    (μ : Fin B → Fin B → ℝ) (hU : CoulPSD dos μ) (Nc : WithTop ℕ) (T : ℝ) (hT : 0 < T)
    (hUn : ¬ MStable tp sp dos μ Nc T) : T ≤ κ / a0 := by
  by_contra h
  apply hUn
  apply mstable_of_gt tp a0 η hH hη ha0 hN sp dos hdos v hv κ hsub μ hU Nc T hT
  have := not_le.mp h
  rw [div_lt_iff₀ ha0] at this; linarith

end

end MultiBand
