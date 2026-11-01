/-
Route R1: structural lemmas L1–L6 and MAIN (finite-atomic spectra).

Conventions
* A finitely supported real sequence is represented by `N : ℕ` and `x : ℕ → ℝ`;
  only `x 0, …, x (N-1)` enter (all sums run over `Finset.range N`).
* `κ : ℝ` is the constant `2π` appearing in `(2π k T)²`.  All statements below hold for
  every real `κ`, in particular for `κ = 2 * Real.pi`; keeping it abstract avoids importing
  the analytic definition of `π` (it enters the hypothesis `H` and the conclusions through
  the same symbol, so nothing is lost).
* The Matsubara cutoff `Nc : WithTop ℕ` (`⊤` = no cutoff).
-/
import Mathlib.Basic.Real.Basic
import Mathlib.Algebra.Order.BigOperators.Group.Finset
import Mathlib.Algebra.BigOperators.Ring.Finset
import Mathlib.Tactic.Linarith
import Mathlib.Tactic.Positivity
import Mathlib.Tactic.FieldSimp
import Mathlib.Tactic.Ring

open Finset

namespace TcR1

noncomputable section

/-! ## Definitions -/

/-- `|n - m|` on `ℕ`. -/
def adist (n m : ℕ) : ℕ := (n - m) + (m - n)

/-- `|x|² = ∑_{n<N} x_n²`. -/
def nsq (N : ℕ) (x : ℕ → ℝ) : ℝ := ∑ n ∈ range N, x n ^ 2

/-- `C(x) = -∑_{n<N} (2n+1) x_n²`. -/
def C (N : ℕ) (x : ℕ → ℝ) : ℝ := -∑ n ∈ range N, (2 * (n : ℝ) + 1) * x n ^ 2

/-- `∑_{k=1}^{n} L(k)`. -/
def Ssum (L : ℕ → ℝ) (n : ℕ) : ℝ := ∑ k ∈ range n, L (k + 1)

/-- Even-sector kernel, `λ(0)` already cancelled:
diagonal `L(2n+1) - 2∑_{k=1}^n L(k)`, off-diagonal `L(|n-m|) + L(n+m+1)`. -/
def Kp (L : ℕ → ℝ) (n m : ℕ) : ℝ :=
  if n = m then L (2 * n + 1) - 2 * Ssum L n else L (adist n m) + L (n + m + 1)

/-- Odd-sector kernel: diagonal `-L(2n+1) - 2∑_{k=1}^n L(k)`, off-diagonal
`L(|n-m|) - L(n+m+1)`. -/
def Km (L : ℕ → ℝ) (n m : ℕ) : ℝ :=
  if n = m then -L (2 * n + 1) - 2 * Ssum L n else L (adist n m) - L (n + m + 1)

/-- The `λ`-part of the even form for a coupling sequence `L`. -/
def Fp (L : ℕ → ℝ) (N : ℕ) (x : ℕ → ℝ) : ℝ :=
  ∑ n ∈ range N, ∑ m ∈ range N, Kp L n m * x n * x m

/-- The `λ`-part of the odd form for a coupling sequence `L`. -/
def Fm (L : ℕ → ℝ) (N : ℕ) (x : ℕ → ℝ) : ℝ :=
  ∑ n ∈ range N, ∑ m ∈ range N, Km L n m * x n * x m

/-- A finite atomic spectral measure `∑_j p_j δ_{Ω_j}` with `p_j ≥ 0`, `Ω_j > 0`. -/
structure Spectrum (ι : Type*) where
  S : Finset ι
  p : ι → ℝ
  Ω : ι → ℝ
  p_nonneg : ∀ j ∈ S, 0 ≤ p j
  Ω_pos : ∀ j ∈ S, 0 < Ω j

/-- Total weight `A = ∑_j p_j`. -/
def Spectrum.A {ι : Type*} (sp : Spectrum ι) : ℝ := ∑ j ∈ sp.S, sp.p j

/-- `g_T(k, w) = 2w / (w² + (κ k T)²)`. -/
def g (κ T : ℝ) (k : ℕ) (w : ℝ) : ℝ := 2 * w / (w ^ 2 + (κ * k * T) ^ 2)

/-- `λ_T(k) = ∫ g_T(k,w) dα²F(w) = ∑_j p_j g_T(k, Ω_j)`. -/
def lam {ι : Type*} (κ : ℝ) (sp : Spectrum ι) (T : ℝ) (k : ℕ) : ℝ :=
  ∑ j ∈ sp.S, sp.p j * g κ T k (sp.Ω j)

/-- `2 μ* (∑_{n<N, n<N_c} x_n)²`. -/
def muTerm (μ : ℝ) (Nc : WithTop ℕ) (N : ℕ) (x : ℕ → ℝ) : ℝ :=
  2 * μ * (∑ n ∈ (range N).filter (fun n : ℕ => ((n : ℕ) : WithTop ℕ) < Nc), x n) ^ 2

/-- Even form `Q⁺_T(x)`. -/
def Qp {ι : Type*} (κ : ℝ) (sp : Spectrum ι) (μ : ℝ) (Nc : WithTop ℕ) (T : ℝ)
    (N : ℕ) (x : ℕ → ℝ) : ℝ :=
  C N x + Fp (lam κ sp T) N x - muTerm μ Nc N x

/-- Odd form `Q⁻_T(x)`. -/
def Qm {ι : Type*} (κ : ℝ) (sp : Spectrum ι) (T : ℝ) (N : ℕ) (x : ℕ → ℝ) : ℝ :=
  C N x + Fm (lam κ sp T) N x

/-- Dimensionless Einstein coupling `2 v / (v² + (κ k)²)` (multiply by `a`). -/
def ge (κ v : ℝ) (k : ℕ) : ℝ := 2 * v / (v ^ 2 + (κ * k) ^ 2)

/-- Dimensionless Einstein form `Qe(a, v; x)`: the even form at `T = 1`, `μ* = 0`,
for the single atom of weight `a` at frequency `v`. -/
def Qe (κ a v : ℝ) (N : ℕ) (x : ℕ → ℝ) : ℝ :=
  C N x + Fp (fun k => a * ge κ v k) N x

/-- The numerical hypothesis `H(a0, η)` (certified separately). -/
def H (κ a0 η : ℝ) : Prop :=
  ∀ v : ℝ, 0 < v → ∀ (N : ℕ) (x : ℕ → ℝ), Qe κ a0 v N x ≤ -η * nsq N x

/-- `ρ(T) ≥ 0` for a quadratic form `Q` on finitely supported sequences:
for every `ε > 0` there is a nonzero `x` with `Q(x) > -ε |x|²`. -/
def Unstable (Q : ℕ → (ℕ → ℝ) → ℝ) : Prop :=
  ∀ ε : ℝ, 0 < ε → ∃ (N : ℕ) (x : ℕ → ℝ), (∃ n ∈ range N, x n ≠ 0) ∧ -ε * nsq N x < Q N x

/-! ## Faithfulness: the forms agree with the SPEC forms whenever `λ(0)` is finite -/

/-- The ROUTE_R1 display of `Fp`: off-diagonal sum plus diagonal sum. -/
theorem Fp_route_form (L : ℕ → ℝ) (N : ℕ) (x : ℕ → ℝ) :
    Fp L N x = ∑ n ∈ range N, ∑ m ∈ (range N).erase n,
        (L (adist n m) + L (n + m + 1)) * x n * x m
      + ∑ n ∈ range N, (L (2 * n + 1) - 2 * Ssum L n) * x n ^ 2 := by
  unfold Fp
  rw [← sum_add_distrib]
  apply sum_congr rfl
  intro n hn
  rw [← add_sum_erase _ _ hn]
  have : ∑ m ∈ (range N).erase n, Kp L n m * x n * x m
      = ∑ m ∈ (range N).erase n, (L (adist n m) + L (n + m + 1)) * x n * x m := by
    apply sum_congr rfl
    intro m hm
    have : n ≠ m := fun h => (ne_of_mem_erase hm) h.symm
    simp only [Kp, this, ite_false]
  rw [this]
  simp only [Kp, ite_true]
  ring

/-- The ROUTE_R1 display of `Fm`. -/
theorem Fm_route_form (L : ℕ → ℝ) (N : ℕ) (x : ℕ → ℝ) :
    Fm L N x = ∑ n ∈ range N, ∑ m ∈ (range N).erase n,
        (L (adist n m) - L (n + m + 1)) * x n * x m
      + ∑ n ∈ range N, (-L (2 * n + 1) - 2 * Ssum L n) * x n ^ 2 := by
  unfold Fm
  rw [← sum_add_distrib]
  apply sum_congr rfl
  intro n hn
  rw [← add_sum_erase _ _ hn]
  have : ∑ m ∈ (range N).erase n, Km L n m * x n * x m
      = ∑ m ∈ (range N).erase n, (L (adist n m) - L (n + m + 1)) * x n * x m := by
    apply sum_congr rfl
    intro m hm
    have : n ≠ m := fun h => (ne_of_mem_erase hm) h.symm
    simp only [Km, this, ite_false]
  rw [this]
  simp only [Km, ite_true]
  ring

/-- For any value `L0` given to `λ(0)`, the cancelled form `Fp` equals the SPEC expression
`∑_{n,m} [L(|n-m|) + L(n+m+1)] x_n x_m - ∑_n [L0 + 2∑_{k=1}^n L(k)] x_n²`. -/
theorem Fp_spec_form (L : ℕ → ℝ) (L0 : ℝ) (N : ℕ) (x : ℕ → ℝ) :
    Fp L N x = ∑ n ∈ range N, ∑ m ∈ range N,
        ((if adist n m = 0 then L0 else L (adist n m)) + L (n + m + 1)) * x n * x m
      - ∑ n ∈ range N, (L0 + 2 * Ssum L n) * x n ^ 2 := by
  rw [Fp_route_form, ← sum_add_distrib, ← sum_sub_distrib]
  apply sum_congr rfl
  intro n hn
  rw [← add_sum_erase _ _ hn]
  have h0 : adist n n = 0 := by unfold adist; omega
  have h1 : n + n + 1 = 2 * n + 1 := by omega
  have hs : ∑ m ∈ (range N).erase n,
      ((if adist n m = 0 then L0 else L (adist n m)) + L (n + m + 1)) * x n * x m
      = ∑ m ∈ (range N).erase n, (L (adist n m) + L (n + m + 1)) * x n * x m := by
    apply sum_congr rfl; intro m hm
    have hnm : n ≠ m := fun h => (ne_of_mem_erase hm) h.symm
    have h2 : adist n m ≠ 0 := by unfold adist; omega
    rw [ite_eq_right_iff.mpr (fun h => absurd h h2)]
  rw [hs, h0, h1]
  simp only [↓reduceIte]
  ring

/-- Same for the odd kernel:
`∑_{n,m} [L(|n-m|) - L(n+m+1)] x_n x_m - ∑_n [L0 + 2∑_{k=1}^n L(k)] x_n²`. -/
theorem Fm_spec_form (L : ℕ → ℝ) (L0 : ℝ) (N : ℕ) (x : ℕ → ℝ) :
    Fm L N x = ∑ n ∈ range N, ∑ m ∈ range N,
        ((if adist n m = 0 then L0 else L (adist n m)) - L (n + m + 1)) * x n * x m
      - ∑ n ∈ range N, (L0 + 2 * Ssum L n) * x n ^ 2 := by
  rw [Fm_route_form, ← sum_add_distrib, ← sum_sub_distrib]
  apply sum_congr rfl
  intro n hn
  rw [← add_sum_erase _ _ hn]
  have h0 : adist n n = 0 := by unfold adist; omega
  have h1 : n + n + 1 = 2 * n + 1 := by omega
  have hs : ∑ m ∈ (range N).erase n,
      ((if adist n m = 0 then L0 else L (adist n m)) - L (n + m + 1)) * x n * x m
      = ∑ m ∈ (range N).erase n, (L (adist n m) - L (n + m + 1)) * x n * x m := by
    apply sum_congr rfl; intro m hm
    have hnm : n ≠ m := fun h => (ne_of_mem_erase hm) h.symm
    have h2 : adist n m ≠ 0 := by unfold adist; omega
    rw [ite_eq_right_iff.mpr (fun h => absurd h h2)]
  rw [hs, h0, h1]
  simp only [↓reduceIte]
  ring

/-- `Q⁺_T` equals the SPEC even form for any finite value `L0` assigned to `λ(0)`. -/
theorem Qp_spec_form {ι : Type*} (κ : ℝ) (sp : Spectrum ι) (μ : ℝ) (Nc : WithTop ℕ)
    (T : ℝ) (L0 : ℝ) (N : ℕ) (x : ℕ → ℝ) :
    Qp κ sp μ Nc T N x =
      ∑ n ∈ range N, ∑ m ∈ range N,
        ((if adist n m = 0 then L0 else lam κ sp T (adist n m)) + lam κ sp T (n + m + 1))
          * x n * x m
      - ∑ n ∈ range N, ((2 * (n : ℝ) + 1) + L0 + 2 * Ssum (lam κ sp T) n) * x n ^ 2
      - 2 * μ * (∑ n ∈ (range N).filter (fun n : ℕ => ((n : ℕ) : WithTop ℕ) < Nc), x n) ^ 2 := by
  unfold Qp muTerm C
  rw [Fp_spec_form _ L0]
  have : ∑ n ∈ range N, ((2 * (n : ℝ) + 1) + L0 + 2 * Ssum (lam κ sp T) n) * x n ^ 2
      = ∑ n ∈ range N, (2 * (n : ℝ) + 1) * x n ^ 2
        + ∑ n ∈ range N, (L0 + 2 * Ssum (lam κ sp T) n) * x n ^ 2 := by
    rw [← sum_add_distrib]; apply sum_congr rfl; intro n _; ring
  rw [this]; ring

/-- `Q⁻_T` equals the SPEC odd form for any finite value `L0` assigned to `λ(0)`. -/
theorem Qm_spec_form {ι : Type*} (κ : ℝ) (sp : Spectrum ι) (T : ℝ) (L0 : ℝ) (N : ℕ)
    (x : ℕ → ℝ) :
    Qm κ sp T N x =
      ∑ n ∈ range N, ∑ m ∈ range N,
        ((if adist n m = 0 then L0 else lam κ sp T (adist n m)) - lam κ sp T (n + m + 1))
          * x n * x m
      - ∑ n ∈ range N, ((2 * (n : ℝ) + 1) + L0 + 2 * Ssum (lam κ sp T) n) * x n ^ 2 := by
  unfold Qm C
  rw [Fm_spec_form _ L0]
  have : ∑ n ∈ range N, ((2 * (n : ℝ) + 1) + L0 + 2 * Ssum (lam κ sp T) n) * x n ^ 2
      = ∑ n ∈ range N, (2 * (n : ℝ) + 1) * x n ^ 2
        + ∑ n ∈ range N, (L0 + 2 * Ssum (lam κ sp T) n) * x n ^ 2 := by
    rw [← sum_add_distrib]; apply sum_congr rfl; intro n _; ring
  rw [this]; ring

/-! ## Linearity of the kernels in the coupling sequence -/

theorem Kp_sum {ι : Type*} (S : Finset ι) (f : ι → ℕ → ℝ) (n m : ℕ) :
    Kp (fun k => ∑ j ∈ S, f j k) n m = ∑ j ∈ S, Kp (f j) n m := by
  unfold Kp Ssum
  split_ifs
  · rw [sum_comm, mul_sum, ← sum_sub_distrib]
  · rw [← sum_add_distrib]

theorem Kp_smul (c : ℝ) (L : ℕ → ℝ) (n m : ℕ) :
    Kp (fun k => c * L k) n m = c * Kp L n m := by
  unfold Kp Ssum
  split_ifs
  · rw [← mul_sum]; ring
  · ring

theorem Fp_sum {ι : Type*} (S : Finset ι) (f : ι → ℕ → ℝ) (N : ℕ) (x : ℕ → ℝ) :
    Fp (fun k => ∑ j ∈ S, f j k) N x = ∑ j ∈ S, Fp (f j) N x := by
  unfold Fp
  calc ∑ n ∈ range N, ∑ m ∈ range N, Kp (fun k => ∑ j ∈ S, f j k) n m * x n * x m
      = ∑ n ∈ range N, ∑ m ∈ range N, ∑ j ∈ S, Kp (f j) n m * x n * x m := by
        simp only [Kp_sum, sum_mul]
    _ = ∑ n ∈ range N, ∑ j ∈ S, ∑ m ∈ range N, Kp (f j) n m * x n * x m := by
        apply sum_congr rfl; intro n _; exact sum_comm
    _ = ∑ j ∈ S, ∑ n ∈ range N, ∑ m ∈ range N, Kp (f j) n m * x n * x m := sum_comm

theorem Fp_smul (c : ℝ) (L : ℕ → ℝ) (N : ℕ) (x : ℕ → ℝ) :
    Fp (fun k => c * L k) N x = c * Fp L N x := by
  unfold Fp
  rw [mul_sum]
  apply sum_congr rfl; intro n _
  rw [mul_sum]
  apply sum_congr rfl; intro m _
  rw [Kp_smul]; ring

/-! ## L1: affine decomposition in the spectrum -/

/-- **L1.** `Q⁺_T(x) = C(x) - 2μ*(…)² + ∑_j p_j F_T(Ω_j, x)` (the finite-atomic form of
`C(x) - 2μ*(…)² + ∫ F_T(w,x) dα²F(w)`). -/
theorem L1 {ι : Type*} (κ : ℝ) (sp : Spectrum ι) (μ : ℝ) (Nc : WithTop ℕ) (T : ℝ)
    (N : ℕ) (x : ℕ → ℝ) :
    Qp κ sp μ Nc T N x
      = C N x - muTerm μ Nc N x
        + ∑ j ∈ sp.S, sp.p j * Fp (fun k => g κ T k (sp.Ω j)) N x := by
  unfold Qp lam
  rw [Fp_sum sp.S (fun j k => sp.p j * g κ T k (sp.Ω j))]
  have : ∑ j ∈ sp.S, Fp (fun k => sp.p j * g κ T k (sp.Ω j)) N x
      = ∑ j ∈ sp.S, sp.p j * Fp (fun k => g κ T k (sp.Ω j)) N x := by
    apply sum_congr rfl; intro j _
    exact Fp_smul (sp.p j) (fun k => g κ T k (sp.Ω j)) N x
  rw [this]; ring

/-! ## L2: extreme-point bound (finite atomic) -/

/-- **L2**, abstract finite-atomic form: if `c + A F_j ≤ -η s` for every atom, `p_j ≥ 0`,
`∑ p_j = A > 0`, then `c + ∑ p_j F_j ≤ -η s`. -/
theorem L2_abstract {ι : Type*} (S : Finset ι) (p F : ι → ℝ) (c A η s : ℝ)
    (hp : ∀ j ∈ S, 0 ≤ p j) (hA : ∑ j ∈ S, p j = A) (hApos : 0 < A)
    (h : ∀ j ∈ S, c + A * F j ≤ -η * s) :
    c + ∑ j ∈ S, p j * F j ≤ -η * s := by
  have h1 : ∀ j ∈ S, p j * F j ≤ p j * ((-η * s - c) / A) := by
    intro j hj
    apply mul_le_mul_of_nonneg_left _ (hp j hj)
    rw [le_div_iff₀ hApos]
    linarith [h j hj, mul_comm (F j) A]
  have h2 := sum_le_sum h1
  rw [← sum_mul, hA] at h2
  have h3 : A * ((-η * s - c) / A) = -η * s - c := by field_simp
  linarith

/-- **L2** for the even form: if `C(x) + A F_T(w,x) ≤ -η|x|²` for every `w > 0`, then
`C(x) + ∑_j p_j F_T(Ω_j, x) ≤ -η|x|²` for every finite atomic spectrum of total weight
`A > 0`. -/
theorem L2 {ι : Type*} (κ : ℝ) (sp : Spectrum ι) (T η : ℝ) (N : ℕ) (x : ℕ → ℝ)
    (hA : 0 < sp.A)
    (h : ∀ w : ℝ, 0 < w → C N x + sp.A * Fp (fun k => g κ T k w) N x ≤ -η * nsq N x) :
    C N x + ∑ j ∈ sp.S, sp.p j * Fp (fun k => g κ T k (sp.Ω j)) N x ≤ -η * nsq N x :=
  L2_abstract sp.S sp.p (fun j => Fp (fun k => g κ T k (sp.Ω j)) N x) (C N x) sp.A η
    (nsq N x) sp.p_nonneg rfl hA (fun j hj => h (sp.Ω j) (sp.Ω_pos j hj))

/-! ## L3: the μ* term only lowers the even form -/

/-- **L3.** `-2μ*(∑_{n<N_c} x_n)² ≤ 0` for `μ* ≥ 0`. -/
theorem L3 (μ : ℝ) (hμ : 0 ≤ μ) (Nc : WithTop ℕ) (N : ℕ) (x : ℕ → ℝ) :
    -muTerm μ Nc N x ≤ 0 := by
  unfold muTerm
  have : 0 ≤ 2 * μ * (∑ n ∈ (range N).filter (fun n : ℕ => ((n : ℕ) : WithTop ℕ) < Nc), x n) ^ 2 := by
    positivity
  linarith

theorem muTerm_zero (Nc : WithTop ℕ) (N : ℕ) (x : ℕ → ℝ) : muTerm 0 Nc N x = 0 := by
  simp [muTerm]

/-! ## L4: scaling -/

/-- **L4**, pointwise: `A g_T(k,w) = 2 a v / (v² + (κ k)²)` with `a = A/T`, `v = w/T`. -/
theorem L4_pointwise (κ A T w : ℝ) (k : ℕ) (hT : 0 < T) (hw : 0 < w) :
    A * g κ T k w = (A / T) * ge κ (w / T) k := by
  unfold g ge
  have h1 : w ^ 2 + (κ * k * T) ^ 2 ≠ 0 := by positivity
  have h2 : (w / T) ^ 2 + (κ * k) ^ 2 ≠ 0 := by positivity
  have h3 : T ≠ 0 := hT.ne'
  have h4 : (w / T) ^ 2 + (κ * k) ^ 2 = (w ^ 2 + (κ * k * T) ^ 2) / T ^ 2 := by
    field_simp
  rw [h4]
  field_simp

/-- **L4.** `C(x) + A F_T(w,x) = Qe(A/T, w/T; x)` for `T > 0`, `w > 0`. -/
theorem L4 (κ A T w : ℝ) (hT : 0 < T) (hw : 0 < w) (N : ℕ) (x : ℕ → ℝ) :
    C N x + A * Fp (fun k => g κ T k w) N x = Qe κ (A / T) (w / T) N x := by
  unfold Qe
  rw [← Fp_smul]
  congr 2
  funext k
  exact L4_pointwise κ A T w k hT hw

/-! ## L5: convexity (linearity) in `a` -/

theorem C_le (N : ℕ) (x : ℕ → ℝ) : C N x ≤ -nsq N x := by
  unfold C nsq
  rw [neg_le_neg_iff]
  apply sum_le_sum
  intro n _
  have h0 : (0 : ℝ) ≤ n := Nat.cast_nonneg n
  nlinarith [mul_nonneg h0 (sq_nonneg (x n))]

theorem nsq_nonneg (N : ℕ) (x : ℕ → ℝ) : 0 ≤ nsq N x := by
  unfold nsq; positivity

/-- `Qe(a, v; x) = C(x) + a G(v; x)`. -/
theorem Qe_linear (κ a v : ℝ) (N : ℕ) (x : ℕ → ℝ) :
    Qe κ a v N x = C N x + a * Fp (ge κ v) N x := by
  unfold Qe
  rw [Fp_smul]

/-- **L5.** If `Qe(a0, v; x) ≤ -η|x|²` and `0 ≤ a ≤ a0` (`a0 > 0`) then
`Qe(a, v; x) ≤ -min(η,1)|x|²`. -/
theorem L5 (κ a0 a v η : ℝ) (N : ℕ) (x : ℕ → ℝ) (ha0 : 0 < a0) (ha : 0 ≤ a)
    (haa0 : a ≤ a0) (h : Qe κ a0 v N x ≤ -η * nsq N x) :
    Qe κ a v N x ≤ -min η 1 * nsq N x := by
  rw [Qe_linear] at h ⊢
  set G := Fp (ge κ v) N x
  set c := C N x
  set s := nsq N x
  set m := min η 1
  have hc : c ≤ -s := C_le N x
  have hs : 0 ≤ s := nsq_nonneg N x
  have hm1 : m ≤ η := min_le_left η 1
  have hm2 : m ≤ 1 := min_le_right η 1
  have h1 : a * (c + a0 * G) ≤ a * (-η * s) := mul_le_mul_of_nonneg_left h ha
  have h2 : (a0 - a) * c ≤ (a0 - a) * (-s) := mul_le_mul_of_nonneg_left hc (by linarith)
  have h3 : a * (-η * s) ≤ a * (-m * s) :=
    mul_le_mul_of_nonneg_left (by nlinarith [mul_le_mul_of_nonneg_right hm1 hs]) ha
  have h4 : (a0 - a) * (-s) ≤ (a0 - a) * (-m * s) :=
    mul_le_mul_of_nonneg_left (by nlinarith [mul_le_mul_of_nonneg_right hm2 hs])
      (by linarith)
  have key : a0 * (c + a * G) ≤ a0 * (-m * s) := by nlinarith
  exact le_of_mul_le_mul_left key ha0

/-! ## L6: odd sector is dominated by the even sector at `|x|` -/

theorem K_abs (L : ℕ → ℝ) (hL : ∀ k, 1 ≤ k → 0 ≤ L k) (n m : ℕ) (x : ℕ → ℝ) :
    Km L n m * x n * x m ≤ Kp L n m * |x n| * |x m| := by
  unfold Km Kp
  split_ifs with h
  · subst h
    have ha : 0 ≤ L (2 * n + 1) := hL _ (by omega)
    have hxx : |x n| * |x n| = x n * x n := abs_mul_abs_self (x n)
    have ht : 0 ≤ x n * x n := mul_self_nonneg (x n)
    have e : 0 ≤ L (2 * n + 1) * (x n * x n) := mul_nonneg ha ht
    calc (-L (2 * n + 1) - 2 * Ssum L n) * x n * x n
        = (L (2 * n + 1) - 2 * Ssum L n) * (x n * x n) - 2 * (L (2 * n + 1) * (x n * x n)) := by
          ring
      _ ≤ (L (2 * n + 1) - 2 * Ssum L n) * (x n * x n) := by linarith
      _ = (L (2 * n + 1) - 2 * Ssum L n) * |x n| * |x n| := by rw [mul_assoc, hxx]
  · have hd : 1 ≤ adist n m := by unfold adist; omega
    have ha : 0 ≤ L (adist n m) := hL _ hd
    have hb : 0 ≤ L (n + m + 1) := hL _ (by omega)
    have hyz : x n * x m ≤ |x n| * |x m| := by rw [← abs_mul]; exact le_abs_self _
    have hyz' : -(x n * x m) ≤ |x n| * |x m| := by rw [← abs_mul]; exact neg_le_abs _
    have e1 := mul_le_mul_of_nonneg_left hyz ha
    have e2 := mul_le_mul_of_nonneg_left hyz' hb
    calc (L (adist n m) - L (n + m + 1)) * x n * x m
        = L (adist n m) * (x n * x m) + L (n + m + 1) * (-(x n * x m)) := by ring
      _ ≤ L (adist n m) * (|x n| * |x m|) + L (n + m + 1) * (|x n| * |x m|) := by linarith
      _ = (L (adist n m) + L (n + m + 1)) * |x n| * |x m| := by ring

/-- `F⁻(L; x) ≤ F⁺(L; |x|)` when `L(k) ≥ 0` for `k ≥ 1`. -/
theorem Fm_le_Fp_abs (L : ℕ → ℝ) (hL : ∀ k, 1 ≤ k → 0 ≤ L k) (N : ℕ) (x : ℕ → ℝ) :
    Fm L N x ≤ Fp L N (fun n => |x n|) := by
  unfold Fm Fp
  apply sum_le_sum; intro n _
  apply sum_le_sum; intro m _
  exact K_abs L hL n m x

theorem C_abs (N : ℕ) (x : ℕ → ℝ) : C N (fun n => |x n|) = C N x := by
  simp [C, sq_abs]

theorem nsq_abs (N : ℕ) (x : ℕ → ℝ) : nsq N (fun n => |x n|) = nsq N x := by
  simp [nsq, sq_abs]

theorem g_nonneg (κ T : ℝ) (k : ℕ) (w : ℝ) (hw : 0 < w) : 0 ≤ g κ T k w := by
  unfold g; positivity

theorem lam_nonneg {ι : Type*} (κ : ℝ) (sp : Spectrum ι) (T : ℝ) (k : ℕ) :
    0 ≤ lam κ sp T k := by
  unfold lam
  apply sum_nonneg; intro j hj
  exact mul_nonneg (sp.p_nonneg j hj) (g_nonneg κ T k _ (sp.Ω_pos j hj))

/-- **L6.** `Q⁻_T(x) ≤ Q⁺_T(|x|)` at `μ* = 0` (any cutoff), and `| |x| | = |x|`. -/
theorem L6 {ι : Type*} (κ : ℝ) (sp : Spectrum ι) (Nc : WithTop ℕ) (T : ℝ) (N : ℕ)
    (x : ℕ → ℝ) :
    Qm κ sp T N x ≤ Qp κ sp 0 Nc T N (fun n => |x n|) := by
  unfold Qm Qp
  rw [muTerm_zero, C_abs]
  have := Fm_le_Fp_abs (lam κ sp T) (fun k _ => lam_nonneg κ sp T k) N x
  linarith

/-! ## MAIN -/

/-- **MAIN (even sector).** Under `H(a0, η)`, for every finite atomic spectrum of total
weight `A > 0`, every `μ* ≥ 0`, every cutoff `N_c`, and every `T > 0` with `A/T ≤ a0`:
`Q⁺_T(x) ≤ -min(η,1) |x|²` for all `N`, `x`. -/
theorem main_even {ι : Type*} (κ a0 η : ℝ) (hH : H κ a0 η) (sp : Spectrum ι)
    (hA : 0 < sp.A) (μ : ℝ) (hμ : 0 ≤ μ) (Nc : WithTop ℕ) (T : ℝ) (hT : 0 < T)
    (haT : sp.A / T ≤ a0) (N : ℕ) (x : ℕ → ℝ) :
    Qp κ sp μ Nc T N x ≤ -min η 1 * nsq N x := by
  have hAT : 0 < sp.A / T := div_pos hA hT
  have ha0 : 0 < a0 := lt_of_lt_of_le hAT haT
  rw [L1]
  have hpt : ∀ w : ℝ, 0 < w →
      C N x + sp.A * Fp (fun k => g κ T k w) N x ≤ -min η 1 * nsq N x := by
    intro w hw
    rw [L4 κ sp.A T w hT hw]
    exact L5 κ a0 (sp.A / T) (w / T) η N x ha0 hAT.le haT (hH (w / T) (div_pos hw hT) N x)
  have h2 := L2 κ sp T (min η 1) N x hA hpt
  have h3 := L3 μ hμ Nc N x
  linarith

/-- **MAIN (odd sector).** Same hypotheses (no `μ*` in the odd sector):
`Q⁻_T(x) ≤ -min(η,1) |x|²`. -/
theorem main_odd {ι : Type*} (κ a0 η : ℝ) (hH : H κ a0 η) (sp : Spectrum ι)
    (hA : 0 < sp.A) (T : ℝ) (hT : 0 < T) (haT : sp.A / T ≤ a0) (N : ℕ) (x : ℕ → ℝ) :
    Qm κ sp T N x ≤ -min η 1 * nsq N x := by
  calc Qm κ sp T N x ≤ Qp κ sp 0 ⊤ T N (fun n => |x n|) := L6 κ sp ⊤ T N x
    _ ≤ -min η 1 * nsq N (fun n => |x n|) :=
        main_even κ a0 η hH sp hA 0 le_rfl ⊤ T hT haT N _
    _ = -min η 1 * nsq N x := by rw [nsq_abs]

/-- **MAIN.** Both sectors. -/
theorem main {ι : Type*} (κ a0 η : ℝ) (hH : H κ a0 η) (sp : Spectrum ι)
    (hA : 0 < sp.A) (μ : ℝ) (hμ : 0 ≤ μ) (Nc : WithTop ℕ) (T : ℝ) (hT : 0 < T)
    (haT : sp.A / T ≤ a0) :
    (∀ (N : ℕ) (x : ℕ → ℝ), Qp κ sp μ Nc T N x ≤ -min η 1 * nsq N x) ∧
    (∀ (N : ℕ) (x : ℕ → ℝ), Qm κ sp T N x ≤ -min η 1 * nsq N x) :=
  ⟨main_even κ a0 η hH sp hA μ hμ Nc T hT haT, main_odd κ a0 η hH sp hA T hT haT⟩

/-- A form bounded by `-c|x|²` with `c > 0` is not unstable. -/
theorem not_unstable_of_bound (Q : ℕ → (ℕ → ℝ) → ℝ) (c : ℝ) (hc : 0 < c)
    (hQ : ∀ (N : ℕ) (x : ℕ → ℝ), Q N x ≤ -c * nsq N x) : ¬ Unstable Q := by
  intro hU
  obtain ⟨N, x, _, hx⟩ := hU c hc
  linarith [hQ N x]

/-- **Corollary.** Under `H(a0, η)` with `a0 > 0`, `η > 0`: if `T > 0` is unstable (even
sector with any `μ* ≥ 0`, `N_c`, or odd sector) then `A/T > a0`, i.e. `T < A/a0`.
Hence `Tc ≤ A / a0`. -/
theorem unstable_imp_lt {ι : Type*} (κ a0 η : ℝ) (hH : H κ a0 η) (ha0 : 0 < a0)
    (hη : 0 < η) (sp : Spectrum ι) (hA : 0 < sp.A) (μ : ℝ) (hμ : 0 ≤ μ)
    (Nc : WithTop ℕ) (T : ℝ) (hT : 0 < T)
    (hU : Unstable (Qp κ sp μ Nc T) ∨ Unstable (Qm κ sp T)) :
    T < sp.A / a0 := by
  by_contra hlt0
  have hlt := not_lt.mp hlt0
  have haT : sp.A / T ≤ a0 := by
    rw [div_le_iff₀ hT]
    rw [div_le_iff₀ ha0] at hlt
    linarith
  have hm : 0 < min η 1 := lt_min hη one_pos
  obtain ⟨h1, h2⟩ := main κ a0 η hH sp hA μ hμ Nc T hT haT
  rcases hU with hU | hU
  · exact not_unstable_of_bound _ _ hm h1 hU
  · exact not_unstable_of_bound _ _ hm h2 hU

/-- Equivalent form of the corollary: `T ≥ A / a0` is stable in both sectors. -/
theorem stable_of_ge {ι : Type*} (κ a0 η : ℝ) (hH : H κ a0 η) (ha0 : 0 < a0)
    (hη : 0 < η) (sp : Spectrum ι) (hA : 0 < sp.A) (μ : ℝ) (hμ : 0 ≤ μ)
    (Nc : WithTop ℕ) (T : ℝ) (hT : sp.A / a0 ≤ T) :
    ¬ Unstable (Qp κ sp μ Nc T) ∧ ¬ Unstable (Qm κ sp T) := by
  have hT0 : 0 < T := lt_of_lt_of_le (div_pos hA ha0) hT
  constructor
  · intro hU
    exact absurd hT (not_le.mpr
      (unstable_imp_lt κ a0 η hH ha0 hη sp hA μ hμ Nc T hT0 (Or.inl hU)))
  · intro hU
    exact absurd hT (not_le.mpr
      (unstable_imp_lt κ a0 η hH ha0 hη sp hA μ hμ Nc T hT0 (Or.inr hU)))

end

end TcR1
