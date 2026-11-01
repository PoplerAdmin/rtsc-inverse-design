/-
MultiBand.Defs: the multiband even/odd forms of `wMB/PROOF.md` §1 (finite atomic spectra).

Conventions
* Bands `i : Fin B`. A finitely supported band-resolved gap is `N : ℕ` and `φ : Fin B → ℕ → ℝ`
  (`φ i n = φ_{n,i}`, only `n < N` enters).
* `tp : ℝ` stands for `2π` (as `κ` in `TcR1`); every theorem holds for every real `tp`.
* `κ : ℝ` is the Perron root of the coupling-weight matrix `𝒜`.
* DOS weights `dos i = N_i`.  The symmetric variables are `x_{n,i} = √N_i φ_{n,i}`, so
  `|x|² = ∑_i N_i ∑_n φ_{n,i}² = nrm2 dos N φ` (see `nrm2_eq_x`).
* Finite atomic matrix-valued spectrum `W_ij = N_i α²F_ij = ∑_{s∈S} W s i j δ_{Ω_s}` with
  `W s i j ≥ 0`, `W s i j = W s j i`, `Ω_s > 0`.
-/
import TcR1.Basic

open Finset

namespace MultiBand

open TcR1

noncomputable section

/-- Finite atomic symmetric nonnegative matrix-valued spectrum
`W_ij(ω) = N_i α²F_ij(ω) = ∑_{s ∈ S} W s i j · δ(ω - Ω s)`. -/
structure MSpec (B : ℕ) (α : Type*) where
  S : Finset α
  Ω : α → ℝ
  W : α → Fin B → Fin B → ℝ
  Ω_pos : ∀ s ∈ S, 0 < Ω s
  W_nonneg : ∀ s ∈ S, ∀ i j, 0 ≤ W s i j
  W_symm : ∀ s ∈ S, ∀ i j, W s i j = W s j i

variable {B : ℕ} {α : Type*}

/-- Coupling-weight matrix `𝒜_ij = ∫ α²F_ij(ω) dω = (∑_s W s i j) / N_i`. -/
def Amat (sp : MSpec B α) (dos : Fin B → ℝ) (i j : Fin B) : ℝ :=
  (∑ s ∈ sp.S, sp.W s i j) / dos i

/-- `l_w(k) = 2w / (w² + (2πk)²)` (TcR1's `ge`). -/
abbrev lw (tp w : ℝ) (k : ℕ) : ℝ := ge tp w k

/-- `C(y) = ∑_n (2n+1) y_n²` (positive sign, as in wMB; `TcR1.C = -Cw`). -/
def Cw (N : ℕ) (y : ℕ → ℝ) : ℝ := ∑ n ∈ range N, (2 * (n : ℝ) + 1) * y n ^ 2

/-- `yᵀ B⁺_w z`, `(B⁺_w)_{nm} = l_w(|n-m|) + l_w(n+m+1)`. -/
def Bp (tp w : ℝ) (N : ℕ) (y z : ℕ → ℝ) : ℝ :=
  ∑ n ∈ range N, ∑ m ∈ range N, (lw tp w (adist n m) + lw tp w (n + m + 1)) * y n * z m

/-- `yᵀ B⁻_w z`, `(B⁻_w)_{nm} = l_w(|n-m|) - l_w(n+m+1)`. -/
def Bm (tp w : ℝ) (N : ℕ) (y z : ℕ → ℝ) : ℝ :=
  ∑ n ∈ range N, ∑ m ∈ range N, (lw tp w (adist n m) - lw tp w (n + m + 1)) * y n * z m

/-- `D_{w,n} = l_w(0) + 2 ∑_{k=1}^n l_w(k)`. -/
def Dn (tp w : ℝ) (n : ℕ) : ℝ := lw tp w 0 + 2 * Ssum (lw tp w) n

/-- `yᵀ D_w z` (diagonal). -/
def Dq (tp w : ℝ) (N : ℕ) (y z : ℕ → ℝ) : ℝ := ∑ n ∈ range N, Dn tp w n * y n * z n

/-- `G⁺_w(y,z) = yᵀ(B⁺_w - D_w)z`. -/
def Gp (tp w : ℝ) (N : ℕ) (y z : ℕ → ℝ) : ℝ := Bp tp w N y z - Dq tp w N y z

/-- `G⁻_w(y,z) = yᵀ(B⁻_w - D_w)z`. -/
def Gm (tp w : ℝ) (N : ℕ) (y z : ℕ → ℝ) : ℝ := Bm tp w N y z - Dq tp w N y z

/-- `|x|² = ∑_i N_i ∑_{n<N} φ_{n,i}²`. -/
def nrm2 (dos : Fin B → ℝ) (N : ℕ) (φ : Fin B → ℕ → ℝ) : ℝ := ∑ i, dos i * nsq N (φ i)

/-- `S_i = ∑_{n<N, n<N_c} φ_{n,i}`. -/
def Sc (Nc : WithTop ℕ) (N : ℕ) (y : ℕ → ℝ) : ℝ :=
  ∑ n ∈ (range N).filter (fun n : ℕ => ((n : ℕ) : WithTop ℕ) < Nc), y n

/-- Coulomb term `2 ∑_ij U_ij S_i S_j`, `U_ij = N_i μ*_ij`. -/
def coul (dos : Fin B → ℝ) (μ : Fin B → Fin B → ℝ) (Nc : WithTop ℕ) (N : ℕ)
    (φ : Fin B → ℕ → ℝ) : ℝ :=
  2 * ∑ i, ∑ j, (dos i * μ i j) * Sc Nc N (φ i) * Sc Nc N (φ j)

/-- **Even multiband form** (wMB §1):
`Q⁺(φ) = -∑_i N_i C(φ_i) + (1/T) ∑_ij ∫ dW_ij(ω) [φ_iᵀ B⁺_w φ_j - φ_iᵀ D_w φ_i]
          - 2 ∑_ij U_ij S_i S_j`,  `w = ω/T`.
The `D` (Z-renormalization) term of band `i` is weighted by the row sum `∑_j W_ij`,
the pairing term by the full matrix `W_ij`. -/
def Qplus (tp : ℝ) (sp : MSpec B α) (dos : Fin B → ℝ) (μ : Fin B → Fin B → ℝ)
    (Nc : WithTop ℕ) (T : ℝ) (N : ℕ) (φ : Fin B → ℕ → ℝ) : ℝ :=
  -∑ i, dos i * Cw N (φ i)
  + (1 / T) * ∑ s ∈ sp.S, ∑ i, ∑ j, sp.W s i j *
      (Bp tp (sp.Ω s / T) N (φ i) (φ j) - Dq tp (sp.Ω s / T) N (φ i) (φ i))
  - coul dos μ Nc N φ

/-- **Odd multiband form**: as `Qplus` with `B⁻` and no Coulomb term. -/
def Qminus (tp : ℝ) (sp : MSpec B α) (dos : Fin B → ℝ) (T : ℝ) (N : ℕ)
    (φ : Fin B → ℕ → ℝ) : ℝ :=
  -∑ i, dos i * Cw N (φ i)
  + (1 / T) * ∑ s ∈ sp.S, ∑ i, ∑ j, sp.W s i j *
      (Bm tp (sp.Ω s / T) N (φ i) (φ j) - Dq tp (sp.Ω s / T) N (φ i) (φ i))

/-- Pair-breaking term
`P(φ) = (1/2T) ∑_ij ∫ dW_ij ∑_n D_{w,n} (φ_{n,i} - φ_{n,j})²`. -/
def Pterm (tp : ℝ) (sp : MSpec B α) (T : ℝ) (N : ℕ) (φ : Fin B → ℕ → ℝ) : ℝ :=
  (1 / (2 * T)) * ∑ s ∈ sp.S, ∑ i, ∑ j, sp.W s i j *
      ∑ n ∈ range N, Dn tp (sp.Ω s / T) n * (φ i n - φ j n) ^ 2

/-- Perron gauge: `θ_i = φ_i / v_i`. -/
def theta (v : Fin B → ℝ) (φ : Fin B → ℕ → ℝ) (i : Fin B) (n : ℕ) : ℝ := φ i n / v i

/-- `ψ_ij = (θ_i + θ_j)/2`. -/
def psi (v : Fin B → ℝ) (φ : Fin B → ℕ → ℝ) (i j : Fin B) (n : ℕ) : ℝ :=
  (theta v φ i n + theta v φ j n) / 2

/-- `χ_ij = (θ_i - θ_j)/2`. -/
def chi (v : Fin B → ℝ) (φ : Fin B → ℕ → ℝ) (i j : Fin B) (n : ℕ) : ℝ :=
  (theta v φ i n - theta v φ j n) / 2

/-- Gauge weights `m_ij = v_i v_j W_ij` (per atom). -/
def mw (sp : MSpec B α) (v : Fin B → ℝ) (s : α) (i j : Fin B) : ℝ := v i * v j * sp.W s i j

/-- Proposition N of wMB §4.2 (a hypothesis here): for every `w > 0` and finitely supported `y`,
`Q^E_{-a0,w}(y) = -C(y) - a0 yᵀG⁺_w y ≤ 0`. -/
def PropN (tp a0 : ℝ) : Prop :=
  ∀ w : ℝ, 0 < w → ∀ (N : ℕ) (y : ℕ → ℝ), Qe tp (-a0) w N y ≤ 0

/-- `U_ij = N_i μ*_ij` is positive semidefinite (as a quadratic form). -/
def CoulPSD (dos : Fin B → ℝ) (μ : Fin B → Fin B → ℝ) : Prop :=
  ∀ S : Fin B → ℝ, 0 ≤ ∑ i, ∑ j, (dos i * μ i j) * S i * S j

/-- `ρ ≥ 0` for a multiband form: for every `ε > 0` there is a nonzero finitely supported `φ`
with `Q(φ) > -ε |x|²`. -/
def MUnstable (dos : Fin B → ℝ) (Q : ℕ → (Fin B → ℕ → ℝ) → ℝ) : Prop :=
  ∀ ε : ℝ, 0 < ε → ∃ (N : ℕ) (φ : Fin B → ℕ → ℝ),
    (∃ i, ∃ n ∈ range N, φ i n ≠ 0) ∧ -ε * nrm2 dos N φ < Q N φ

/-- Stable at `T`: neither sector is unstable. -/
def MStable (tp : ℝ) (sp : MSpec B α) (dos : Fin B → ℝ) (μ : Fin B → Fin B → ℝ)
    (Nc : WithTop ℕ) (T : ℝ) : Prop :=
  ¬ MUnstable dos (Qplus tp sp dos μ Nc T) ∧ ¬ MUnstable dos (Qminus tp sp dos T)

end

end MultiBand
