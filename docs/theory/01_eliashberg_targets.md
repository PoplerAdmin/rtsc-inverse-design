# 01 — Spectral (α²F-level) requirements for Tc ≥ 300 K in isotropic Migdal–Eliashberg theory

Status: theory deliverable T1. Every statement is tagged **[THM]** (theorem within the stated model,
with proof or certified proof referenced), **[DER]** (derived here; the proof is given, but a
numerical step was done in floating point and not certified), **[NUM]** (computed; script given),
**[CIT]** (cited and checked during this work), **[UNVERIFIED]** (recalled from memory, not
checked), or **[SOFT]** (empirical / shape-dependent).

Scripts: `src/theory/eliashberg.py` (solver), `validate_solver.py`, `scan_targets.py`,
`bounds_extra.py`, `make_tables.py`, `g_table.py`. Data: `docs/theory/data/*.csv`, `tables.md`, `validation.txt`.
Reproduce (WSL, ~1 min total on 2 threads):

```
export OMP_NUM_THREADS=2 OPENBLAS_NUM_THREADS=2 MKL_NUM_THREADS=2
cd src/theory
nice -n 19 python validate_solver.py > ../../docs/theory/data/validation.txt
nice -n 19 python scan_targets.py   ../../docs/theory/data
nice -n 19 python bounds_extra.py   ../../docs/theory/data
nice -n 19 python make_tables.py    ../../docs/theory/data
nice -n 19 python g_table.py        ../../docs/theory/data
```

---

## 0. Result: what the materials layer must meet

(Revised after review R01; the changes are listed in §7.)

Units: k_B = ħ = 1, so energies are in kelvin (1 meV = 11.6045 K). Definitions are in §1.

- A := ∫α²F dω and λ := 2∫α²F/ω dω.
- m := λ⟨ω²⟩ = 2∫ω α²F dω.
- ω_log := exp[(2/λ)∫ln ω α²F/ω dω].
- ⟨ω⟩ := 2A/λ and ω₂ := sqrt⟨ω²⟩.

All floors scale with the target: A-type floors ∝ Tc_target, and m-type floors ∝ Tc_target².
They are quoted for Tc_target = 300 K.

### 0.1 Necessary conditions (single isotropic α²F, inside the model's validity domain)

| # | condition | value for Tc = 300 K | status |
|---|---|---|---|
| N1 | A ≥ Tc · 4.3314814 | **A ≥ 1299.44 K (111.98 meV)** for any shape, any μ* ≥ 0, any cutoff | [THM] certified (`rtt.core`, §2.1) |
| N1-MB | κ ≥ Tc · 4.33148145 (κ = Perron root of the band/k coupling-weight matrix) | κ ≥ 1299.44 K. Weaker forms: max_i A_i ≥ 1299.44 K; A_iso/w_min ≥ 1299.44 K | [THM] wMB Thm 1 (`rtt.core`) and §5.1 here |
| N1μ | A ≥ Tc · a_E(μ*, ω_c) | Certified floors (03 §4.3): μ* = 0.10: A > 1671.31 / 1705.79 / 1703.17 / 1677.31 K at ω_c = 5000 / 10000 / 20000 / 50000 K; μ* = 0.13: A > 1781.43 / 1823.44 / 1811.15 / 1768.88 K. The μ* = 0.20 range 1954–2091 K is float only. | [THM-C] at the listed (μ*, ω_c), symmetric-block convention, including the T-supremum (03 §4); float [DER] elsewhere. Depends on a *lower* bound on μ*, so it is used for ranking only, not for rejection. |
| NM | m ≥ g(λ), g(λ) := λ·Ω_E(λ)² with Ω_E the μ* = 0 Einstein frequency giving Tc | Certified lower values (03 §0): g_lo(3) = 3.468968×10⁶ K², g_lo(2) = 4.014738×10⁶, g_lo(∞) = **2.695481×10⁶ K²** (table in §0.3) | [THM-C] computer-assisted at every tabulated λ (03 §2), symmetric-block convention for μ* > 0. This one rule replaces the old N2 and N3. |

The μ* = 0 versions are used because μ* ≥ 0 only lowers the even form (Lemma 1 / Lean `L3`).
They therefore hold for every μ* ≥ 0 and need no μ* or cutoff input. This is a quadratic-form
statement, so it holds in the symmetric Coulomb-block convention of §1. In the asymmetric
textbook convention it is not proved (03 §1 F1, R03 finding 3).

**Not a necessary condition: ω_log has no shape-independent bound, even at fixed λ.**
At λ = 2, two atoms with frequency ratio 128 give Tc/ω_log = 1.30 [NUM, Table 6; reproduced by
R01]. ω_log thresholds hold only for a given shape class [SOFT] (§3.4).

**Realistic requirements for hydride-like or covalent shapes** [NUM]/[SOFT]. These assume
μ* = 0.10–0.13 and ME with ω_c = 10 ω_max (symmetric block, §1), from Tables 1 and 5:

| λ | ω_log needed (K) | ω₂ needed (K) | A needed (K) | m needed (K²) |
|---|---|---|---|---|
| 2.0 | 1270–1830 | 1740–2440 | 1740–2020 | 6.1–11.9 ×10⁶ |
| 2.5 | 1010–1530 | 1460–1930 | 1820–2010 | 5.3–9.3 ×10⁶ |
| 3.0 | 850–1340 | 1280–1620 | 1890–2040 | 4.9–7.9 ×10⁶ |

- The lower end of each range is the LaH10-like broad shape; the upper end is Einstein- or
  MgB2-like.
- Allen–Dynes asks for about 5–10 % more ω_log (§4.3). For example, Einstein at λ = 3,
  μ* = 0.13 needs 1458 K instead of 1336 K.
- **Where the weight should sit** [NUM] (Bergmann–Rainer, §3.5): δTc/δα²F(ω) peaks at
  ω ≈ 7.2–7.6 Tc, which is 2150–2270 K (185–196 meV) for Tc = 300 K.

### 0.2 How a computed α²F is turned into a decision

**(a) Validity gate G-VALID.** A candidate gets **UNDECIDED (out-of-model)**, never REJECT, if
any of the following holds:

1. **Non-adiabatic:** λ_hi · ω_max / E_F > 0.1. E_F is measured from the bottom of the
   shallowest band that crosses the Fermi level (or half the bandwidth if smaller). If E_F is not
   reported, the gate fails.
2. **Unstable or anharmonic:**
   - the harmonic phonons have imaginary modes (|ω| > 30 K, i.e. outside numerical noise); or
   - the structure contains H, or Li/Be, and the α²F is harmonic only, with no anharmonic check
     (SSCHA, TDEP or MLFF-MD) showing frequency shifts ≤ 10 %; or
   - an anharmonic check exists and shows shifts > 10 %, but α²F was not recomputed on the
     renormalized phonons.
3. **Coulomb sign:** an attractive electronic contribution (μ* < 0, e.g. plasmonic or excitonic)
   is claimed or suspected.
4. **Missing data:** the convergence data needed for (b) are missing, or any δ_X there is > 0.5
   (unconverged).
5. **Multiband without data:** the material is multiband (n_sheets > 1, see (c)) and none of the
   N1-MB inputs are available. In that case **only** the N1-MB rules may reject.

**(b) Error envelopes**, per candidate and from its own convergence data (no global "1.2×").
The required ladder is:

- K ≥ 3 *admitted* smearings σ₁ < … < σ_K on the finest k-mesh. "Admitted" means mode-positive,
  the `rtt.core` criterion for QE α²F.
- The two finest q-meshes.

For each X ∈ {A, λ, m}:

- δ_X^σ = max_i |X(σ_i) − X_nom|/X_nom over the admitted smearings.
- δ_X^q = |X(q_fine) − X(q_coarse)|/X_nom.
- δ_X = δ_X^σ + δ_X^q.

Two further allowances:

- **Amplitude systematic** (functional and pseudopotential, not visible in the ladders):
  ε_sys = 0.10 [SOFT convention].
- **Frequency envelope:** ε_ω = max(q-ladder spread of ω_log, |anharmonic shift| if a check
  exists, 0.02).

Propagation uses the per-mode identity λ_ν ω_ν² = N_F⟨|dV|²⟩_ν / M_ν (Hopfield). At fixed
electronic matrix elements a frequency error rescales ω_ν, so λ_ν ∝ ω_ν⁻², A ∝ ω⁻¹ and m is
invariant [DER, §5.6]. The upper envelopes are:

    A_hi = A_nom (1+δ_A)(1+ε_sys) / (1−ε_ω)
    λ_hi = λ_nom (1+δ_λ)(1+ε_sys) / (1−ε_ω)²
    m_hi = m_nom (1+δ_m)(1+ε_sys)

Every floor uses only upper envelopes, so an error in either direction can only make a rejection
*less* likely.

**(c) Single-band justification.** n_sheets is the number of distinct Fermi-surface sheets (bands
crossing E_F, spin-degenerate pairs counted once), from the band structure.

- **"Single-band justified"** means n_sheets = 1. Even then the rejection carries the recorded
  assumption **M-ISO** (intra-sheet isotropy). Anisotropy within one sheet can in principle raise
  Tc above the isotropic value (§5.1).
- Isotropic-α²F rejections whose margin is < 25 % (X_hi > 0.75·floor) also need k-resolved data.
  Without it, the candidate becomes UNDECIDED.

### 0.3 Certified table g_lo(λ) for the moment rule NM-G (μ* = 0, Tc = 300 K)

Sources: `docs/theory/certificates/nm_<λ>.json` and `SUMMARY.json`, documented in
`03_certificates.md`.

- The certified value is exactly g_lo = 90000·m₀, where the decimal m₀ is stored in each JSON as
  `m0_over_T2`. Implementations should read m₀, not the float.
- The values below are floor-rounded to 7 significant digits.
- Ω_E is bracketed as [Ω_lo, Ω_hi], with the upper end from a Rayleigh instability certificate.

| λ_k | g_lo (10⁶ K²) | Ω_E ∈ [Ω_lo, Ω_hi] (K) |
|---|---|---|
| 0.8 | 9.452681 | [3437.41, 3437.79] |
| 0.9 | 7.871489 | [2957.38, 2957.70] |
| 1.0 | 6.850483 | [2617.34, 2617.62] |
| 1.1 | 6.147092 | [2363.95, 2364.20] |
| 1.2 | 5.637762 | [2167.51, 2167.75] |
| 1.3 | 5.254201 | [2010.39, 2010.61] |
| 1.4 | 4.956119 | [1881.51, 1881.71] |
| 1.5 | 4.718445 | [1773.59, 1773.78] |
| 1.6 | 4.524867 | [1681.67, 1681.86] |
| 1.7 | 4.364371 | [1602.27, 1602.45] |
| 1.8 | 4.229276 | [1532.83, 1533.01] |
| 1.9 | 4.114078 | [1471.49, 1471.66] |
| 2.0 | 4.014738 | [1416.81, 1416.97] |
| 2.25 | 3.817620 | [1302.58, 1302.72] |
| 2.5 | 3.671310 | [1211.82, 1211.96] |
| 2.75 | 3.558519 | [1137.54, 1137.67] |
| 3.0 | 3.468968 | [1075.32, 1075.44] |
| 3.25 | 3.396173 | [1022.24, 1022.35] |
| 3.5 | 3.335852 | [976.26, 976.38] |
| 3.75 | 3.285060 | [935.95, 936.06] |
| 4.0 | 3.241711 | [900.23, 900.34] |
| 4.5 | 3.171646 | [839.52, 839.62] |
| 5.0 | 3.117483 | [789.61, 789.71] |
| 6.0 | 3.039235 | [711.71, 711.80] |
| 7.0 | 2.985444 | [653.06, 653.14] |
| 8.0 | 2.946200 | [606.85, 606.93] |
| 10.0 | 2.892781 | [537.84, 537.91] |
| 15.0 | 2.824200 | [433.91, 433.96] |
| 20.0 | 2.790995 | [373.56, 373.61] |
| 50.0 | 2.732967 | [233.79, 233.82] |
| 100.0 | 2.714107 | [164.74, 164.77] |
| ∞ | **2.695481** | – |

**Use.** Look up g_lo at the smallest tabulated λ_k ≥ λ_hi. If λ_hi > 100 or unknown, use
g_lo(∞) = 2.695481×10⁶ K². g scales as (Tc_target/300 K)². Any sub-grid of this table may also be
used.

**Status [THM-C]** (03 §2). For each λ_k, an Arb certificate proves Theorem NM(λ_k, m₀):

> Every α²F with λ ≤ λ_k and m ≤ m₀T² is stable at T.

The proof works as follows:

- A Lagrange-multiplier reduction turns this into a one-parameter Einstein family.
- That family is certified on dyadic cells, with analytic ends.
- Two independent checkers re-check it (`src/certify/verify.py`, `src/certify/verify_indep.py`).

What the certificate makes unnecessary:

- The earlier numerical steps "Einstein minimiser" and "g decreasing", and the two-atom reduction.
  The lookup is sound without any monotonicity, because each certificate covers every λ ≤ λ_k
  (03 §2.5).
- Any separate T-supremum argument. That step is exact for μ* = 0, since m/T² decreases in T.

The float values in `data/g_table.csv` (and in the previous version of this table) lie 1×10⁻⁵
relative above g_lo. They are kept for reference only.

### 0.4 Rules as they should be implemented (machine-readable)

Revised for review RL1. Input domains, provenance and the multiband envelopes are defined in §0.6.

```yaml
# Units: kelvin. Tc_target default 300. s1 = Tc_target/300, s2 = s1**2.
# Constants: a0 = 4.33148145 (certified, eta = 1e-9, kernel 2*pi*k*T); F_A = a0*Tc_target = 1299.44*s1.
# Outcome = (status, fired_rule, assumptions[], inputs_used[{name, value, tier}]).
#   status: REJECTED_PROVEN | UNDECIDED | PASS.
#   A REJECTED_PROVEN outcome proves "stable at every T >= Tc_target" (i.e. Tc <= Tc_target),
#   relative to the listed assumptions.
# Order: WELL-FORMED, then G-VALID, then the reject rules. The first gate that fails gives
# UNDECIDED and no reject rule may fire.
# A rule fires only if every input it uses has tier == fullDFPT (Sec. 0.6c). Otherwise the
# rule is skipped (it may be logged as "would fire, ranking only").
common_assumptions: [ME-MODEL, CERT-H, ENV-SOUND]
  # ME-MODEL: linearized ME (isotropic or band/k-discretized), adiabatic, mu* >= 0,
  #           second-order transition.
  # CERT-H:   rtt.core Prop. 6 with a0 = 4.33148145.
  # ENV-SOUND: the *_hi / *_lo values bound the true values in the stated direction.

rules:
  - name: WELL-FORMED
    test: UNDECIDED unless every field present for the candidate is in its domain (Sec. 0.6a)
    status: policy

  - name: G-VALID
    inputs: [lambda_hi, omega_max, E_F, imag_abs_max_K, contains_H_Li_Be, anharmonic_check,
             anharm_shift_abs, alpha2F_recomputed_anharmonic, mu_star_sign, ladder, n_sheets,
             mb_inputs_present]
    test: >
      UNDECIDED if lambda_hi*omega_max/E_F > 0.1
      or imag_abs_max_K > 30
      or (contains_H_Li_Be and not anharmonic_check)
      or (anharm_shift_abs > 0.10 and not alpha2F_recomputed_anharmonic)
      or mu_star_sign == "negative_suspected"
      or ladder incomplete (K < 3 admitted smearings or < 2 q-meshes)
      or any delta > 0.5
      or eps_sys < 0.10
      or eps_omega < max(0.02, anharm_shift_abs, qspread_omega_log)
      or (n_sheets > 1 and not mb_inputs_present)
    status: policy (validity domain of the model)

  - name: ENVELOPE
    output: >
      A_hi   = A_nom*(1+d_A)*(1+eps_sys)/(1-eps_omega)
      lam_hi = lam_nom*(1+d_lam)*(1+eps_sys)/(1-eps_omega)**2
      m_hi   = m_nom*(1+d_m)*(1+eps_sys)
      calA_hi[i][j] = calA_nom[i][j]*(1+d_calA)*(1+eps_sys)/(1-eps_omega)
      kappa_hi = max_i (calA_hi v)_i / v_i   for any v > 0 (Collatz-Wielandt; v = nominal Perron vector)
      rowA_hi[i] = sum_j calA_hi[i][j]
      wmin_lo = min_i w_nom[i]*(1-d_w)
    status: "[DER] (Sec. 0.6b, 5.6); eps_sys convention [SOFT]"

  - name: N1-KAPPA
    uses: [kappa_hi]
    test: REJECT if kappa_hi < 1299.44*s1
    assumptions: [M-DISCRETE]
    status: "[THM] wMB Thm 1 (rtt.core)"

  - name: N1-ROWMAX
    uses: [rowA_hi]
    test: REJECT if max_i rowA_hi[i] < 1299.44*s1
    assumptions: [M-DISCRETE]
    status: "[THM] kappa <= max row sum"

  - name: N1-DOS
    uses: [A_hi, wmin_lo]
    test: REJECT if A_hi / wmin_lo < 1299.44*s1
    assumptions: [M-SHEET-ISO, SHEETS-CORRECT]
    status: "[THM] kappa <= A_iso / w_min (Sec. 5.1)"

  - name: N1-ISO
    uses: [A_hi]
    pre: n_sheets == 1 and (A_hi < 0.75*1299.44*s1 or kresolved_check_done)   # margin clause = policy
    test: REJECT if A_hi < 1299.44*s1
    assumptions: [M-ISO]
    status: "[THM] certified (rtt.core wB/wE1, Lean TcR1)"

  - name: NM-CERT
    uses: [m_hi, lam_hi]
    pre: n_sheets == 1
    test: REJECT if m_hi < 4*(1299.44*s1)**2 / lam_hi
    assumptions: [M-ISO]
    status: "[THM] N1 combined with Cauchy-Schwarz A^2 <= lam*m/4 (Sec. 0.6d)"

  - name: NM-G
    uses: [m_hi, lam_hi]
    pre: n_sheets == 1 and (m_hi < 0.75*g_lo(lam_up)*s2 or kresolved_check_done)
    test: >
      REJECT if m_hi < g_lo(lam_up)*s2, with lam_up = smallest tabulated lambda_k >= lam_hi,
      g_lo(lambda_k) = 90000*m0 (exact decimal m0 from certificates/nm_<lambda_k>.json)
      (use g_lo(inf) = 90000*29.94979021 = 2.695481e6 if lam_hi > 100 or unknown)
    assumptions: [M-ISO, CERT-NM]
    status: "[THM-C] computer-assisted (03_certificates.md Sec. 2); symmetric-block convention
             for the extension to mu* > 0; no monotonicity of g needed"

  - name: RANK-AIM
    use: ranking and generation conditioning only; never rejects
    target: "lambda in [2,3]; A >= 2060; omega_log >= 1210 (lam=3) to 2020 (lam=2);
             omega2 >= 1460 (lam=3) to 2010 (lam=2); m >= 6e6"
    status: "[SOFT] (Sec. 0.5)"
```

Tag meanings:

- **M-ISO**: one isotropic band. The true row sum (and κ) is ≤ A_iso.
- **M-SHEET-ISO**: coupling constant within each sheet.
- **SHEETS-CORRECT**: n_sheets and w_i come from the same structure and functional as α²F.
- **M-DISCRETE**: the band or k-discretized anisotropic ME problem represents the material.
- **CERT-NM**: Theorem NM(λ_k, m₀) of 03 §2.2 holds for the certificate files in
  `docs/theory/certificates/`.
  - It is computer-assisted (Arb ball arithmetic).
  - Two independent checkers re-check it.
  - The Lagrangian reduction and the end lemmas are proved on paper only (no Lean).
  - It replaces the former NUM-EINSTEIN-MIN and NUM-G-MONOTONE, which are dropped.

The 25 % margin clause is policy and is not part of any proof.

### 0.5 Aim box (ranking and generation conditioning, never rejection)

The box assumes μ* ≈ 0.10–0.13 and adds a +10 % Tc margin (aiming at 330 K). Frequencies and A
scale by 1.1, and m by 1.21.

| quantity | aim | reason |
|---|---|---|
| λ | 2.0 – 3.0 | The A-optimal λ for realistic shapes at μ* ≈ 0.1–0.13 is 1.8–3.2 (Table 1). Above ≈3–3.5, lattice-instability and ME-breakdown arguments apply (§5.3). |
| A = ∫α²F dω | ≥ 2060 K (≈ 178 meV) | 1.1 × 1875 K (H3S-like, μ* = 0.13, Table 1). The LaH10-like row gives 1.1 × 1969 = 2166 K. |
| ⟨ω⟩ = 2A/λ | ≥ 1370 K at λ = 3; ≥ 2060 K at λ = 2 | follows from A and λ |
| ω_log | ≥ 1210 K (λ = 3) to ≥ 2020 K (λ = 2) | 1.1 × Table 5b, from the H3S-like row (λ = 3, 1096 K) to the Einstein row (λ = 2, 1833 K) |
| ω₂ = sqrt⟨ω²⟩ | ≥ 1460 K (λ = 3) to ≥ 2010 K (λ = 2) | 1.1 × Ω_E(λ, μ* = 0.13) from Table 4 (1326 K and 1829 K) |
| m = λ⟨ω²⟩ | ≥ 6×10⁶ K² | Realistic shapes need 4.9–12 ×10⁶ K² at 300 K; × 1.21 gives about 5.9–14 ×10⁶ K² |
| where the weight sits | most of A between ≈ 1000 and 3500 K, centred near 2200 K | Bergmann–Rainer peak at 7.2–7.6 Tc (§3.5) |

The aim box is a statement about α²F. Whether any ambient-pressure, dynamically and
thermodynamically stable compound can have such an α²F is the open question for the materials
layer (§5.5).

### 0.6 Input domains, multiband envelopes, provenance (review RL1)

**(a) Domains and sign conventions.** These are checked by WELL-FORMED. Any violation gives
UNDECIDED.

| field | domain / convention |
|---|---|
| A_nom, λ_nom, m_nom | > 0. They come from the admitted (non-negative, mode-positive) α²F. A signed α²F is not admitted. |
| ω_max, E_F | > 0 (K). E_F is measured from the bottom of the shallowest crossing band. |
| imag_abs_max_K | ≥ 0: the **largest abs(ω) among negative (imaginary) eigenvalues**, 0 if there are none. QE and phonopy print these as negative numbers; the mapper must take the absolute value. |
| anharm_shift_abs | ≥ 0: **abs(ω_anh − ω_harm)/ω_harm**, maximized over the coupled modes. It is a magnitude, because softening raises λ and A. |
| n_sheets | integer ≥ 1 |
| δ_A, δ_λ, δ_m, δ_𝒜, δ_w, q-spread of ω_log | ∈ [0, 0.5] |
| ε_sys | ≥ 0.10 |
| ε_ω | ∈ [max(0.02, anharm_shift_abs, q-spread), 1) |
| w_i | ∈ (0, 1], with Σ_i w_i ∈ [0.99, 1.01] |
| w_i with δ_w (N1-DOS consistency, Lean review) | **Σ_i w_i,nom·(1 − δ_w) ≤ 1**. The true sheet fractions sum to exactly 1, so lower envelopes that sum to more than 1 cannot all be valid, and the N1-DOS envelope would be vacuous. Alternatively, normalise w_i to Σ w_i = 1 *before* applying δ_w. If neither holds: UNDECIDED. |
| moment consistency (data quality, Lean review) | **A_nom² ≤ λ_nom·m_nom/4 · (1 + 10⁻²)**. Cauchy–Schwarz (§0.6d) gives A² ≤ λm/4 for every non-negative α²F, so a violation beyond the relative tolerance 10⁻² (allowed for numerical integration) means the nominal moments are inconsistent: UNDECIDED. |
| 𝒜_ij | ≥ 0, with N_i 𝒜_ij = N_j 𝒜_ji within 1 % |
| κ_hi, rowA_hi | ≥ 0, and derived only via ENVELOPE (never entered directly) |

Flags such as mb_inputs_present and kresolved_check_done are *derived* from which fields are
present. They are never set independently.

**(b) Multiband envelopes from DFPT outputs.**

- 𝒜_ij = ∫α²F_ij dω = ½ Σ_{qν} w_q ω_{qν} λ_ij(qν). Here λ_ij(qν) is the electron–phonon
  coupling restricted to initial states on sheet i and final states on sheet j. It is obtained
  by partitioning the Fermi-surface double-delta sum by band index (EPW, or QE `elph` with
  band-resolved weights). A k-resolved version (EPW anisotropic) is the M-DISCRETE refinement.
- δ_𝒜 is the largest entrywise relative deviation of 𝒜 over the same smearing and q ladder as
  §0.2b. 𝒜_hi follows the A formula, because 𝒜 is A-type (∝ ω⁻¹ at fixed matrix elements).
- **κ_hi**: the Perron root is monotone in the entries of a non-negative matrix, so
  κ(𝒜_true) ≤ κ(𝒜_hi). For any v > 0, the Collatz–Wielandt bound κ(𝒜_hi) ≤ max_i (𝒜_hi v)_i/v_i
  makes κ_hi checkable. (The diagonal similarity to sqrt(N_i/N_j)𝒜_ij leaves the spectrum
  unchanged.)
- **rowA_hi[i]** = Σ_j 𝒜_hi,ij.
- **w_min,lo** = min_i w_i,nom (1 − δ_w), where w_i = N_i(E_F)/N(E_F) is the sheet-projected DOS
  at E_F. δ_w is its largest relative deviation over the smearing ladder and the two finest
  k-meshes.
- All three come from the same converged run (same structure hash) as the isotropic α²F.

**(c) Provenance per number.** Each input field carries its own tier:

- **fullDFPT**: from a DFPT (or EPW) run on this structure, with the full convergence ladder of
  §0.2b complete and δ ≤ 0.5;
- **coarseDFPT**: anything less than that;
- **ML**: a surrogate or MLFF;
- **external**: literature or hand entry.

A rule may return REJECTED_PROVEN only if **every input it uses** is fullDFPT and they share one
run/structure id. Envelope parameters (δ, ε) count as inputs.

| rule | inputs used |
|---|---|
| N1-ISO | A_nom, δ_A, ε |
| N1-DOS | the N1-ISO inputs plus w_i, δ_w |
| N1-KAPPA, N1-ROWMAX | 𝒜_ij, δ_𝒜, ε |
| NM-G, NM-CERT | m, λ and their δ, ε |

If any of them has a lower tier, the rule does not fire. The candidate stays UNDECIDED, or is
ranked only.

**(d) Recorded assumptions and NM-CERT.** Every outcome records its fired rule, the rule's
assumption list (YAML) and the common list [ME-MODEL, CERT-H, ENV-SOUND].

NM-CERT (suggested by RL1, finding 12) is a fully certified moment floor. By Cauchy–Schwarz,

  A² = (∫ sqrt(α²F/ω)·sqrt(ω α²F) dω)² ≤ (λ/2)(m/2),

so N1 gives λm ≥ 4(a0Tc)² = 6.754×10⁶ K² at 300 K. It is weaker than NM-G for λ ≳ 1, e.g.
2.25×10⁶ vs 3.47×10⁶ K² at λ = 3. It adds information only when the A envelope is looser than the
λ and m envelopes.

A REJECTED_PROVEN outcome means stability at every T ≥ Tc_target, i.e. **Tc ≤ Tc_target**.
"Tc < Tc_target" is not claimed.

---

## 1. Model and definitions (identical to `rtt.core`)

These are the model definitions of the certified core (`src/rtt/core/spectrum.py`), re-implemented
independently in `src/theory/eliashberg.py`.

- α²F is a finite positive Borel measure on (0, ∞). The definitions of A, λ, λ⟨ω²⟩ and ω_log are
  as in §0. Note that **A = λ⟨ω⟩/2**, where ⟨ω⟩ = (2/λ)∫α²F dω.
- Coupling function: λ_T(k) = ∫ 2ω α²F(ω) dω / (ω² + (2πkT)²) for k ≥ 1.
- Linearized isotropic imaginary-axis ME equations. With φ_n = Δ̃_n/|ω̃_n| and the exact
  cancellation of the λ(0) terms, the even-frequency sector is the quadratic form

  Q⁺_T(x) = Σ_{n≠m}[λ(|n−m|) + λ(n+m+1)] x_n x_m + Σ_n [λ(2n+1) − (2n+1) − 2Σ_{k=1}^{n} λ(k)] x_n² − 2μ* (Σ_{n<N_c} x_n)²,

  with the block cutoff N_c(T) = #{n ≥ 0 : (2n+1)πT ≤ ω_c}.
- **The μ* cutoff convention is a choice.** The Coulomb term is applied *symmetrically* to the
  N_c × N_c block, θ(ω_c − |ω_n|)θ(ω_c − |ω_m|). The textbook (Allen–Mitrović) form puts θ only on
  the summed index m. R01's independent solver (`docs/reviews/r01_scripts/spot_checks.py`, item 6)
  finds that the symmetric form gives the *higher* Tc:
  - by 0.13–0.25 % at ω_c = 10 Ω;
  - by 2.6 % at ω_c = 5000 K (N_c = 3).

  It is therefore the permissive choice, which is the right direction for necessary conditions.
  N1, N1-MB and NM do not use μ* at all, so they are unaffected.
- The odd-frequency sector is never unstable. This is the core's multiband proof (wMB)
  Theorem 3 with B = 1 [THM].
- ρ(T) = sup_x Q⁺_T(x)/|x|². T is *unstable* if ρ(T) ≥ 0, and **Tc := sup{T : T unstable}**. This
  Tc is ≥ the physical (highest-solution) Tc, so every upper bound on it is also an upper bound on
  the physical Tc.

**Conventions for μ\*.** μ* is always referred to a cutoff ω_c. Tables 1 and 5 use
ω_c = c·ω_max with c = 10, and the c = 5 sensitivity is reported. The extremal problems (§3) use
an absolute ω_c (5000–50000 K, reference 20000 K ≈ 1.7 eV ≈ 9 Ω_opt). The Tolmachev /
Morel–Anderson relation converts between cutoffs: 1/μ*(ω₂) = 1/μ*(ω₁) + ln(ω₁/ω₂)
(`eliashberg.mu_rescale`). A given number μ* is a *weaker* repulsion at a larger cutoff.

---

## 2. What the certified core proves (exact statements)

Sources:
- the route-B proof wB (main theorem);
- the wE1 certificate (Prop. 6);
- `lean/TcR1/Basic.lean` (Lean `main`, `unstable_imp_lt`, `stable_of_ge`);
- `src/rtt/core/certify/PROOFS.md` (rules R-UNIV, R-CELLS, C-RAYLEIGH, R-MULTI);
- the multiband proof wMB.

### 2.1 Universal bound (R-UNIV) [THM, computer-assisted, certified]

For every α²F (finite positive measure, mass A), every μ* ≥ 0 and every cutoff N_c ∈ {1, …, ∞}:

  **Tc ≤ A / 4.3314814 = 0.2308679 A.**

The constant is sharp to about 2×10⁻⁸. The best constant c* satisfies
1/4.3314815 ≤ c* ≤ 1/4.3314814. It is attained by an Einstein spectrum A δ(ω − Ω) with μ* = 0,
Ω/Tc = w* ≈ 7.57487 and λ* = 2a*/w* ≈ 1.14365, i.e. Ω = 1.74879 A. The float values are guidance.

How the proof goes:
1. **Sector and Coulomb reduction.** ρ⁻ ≤ ρ⁺_{μ*=0} and ρ⁺_{μ*} ≤ ρ⁺_{μ*=0}.
2. **Affinity.** Q_T[α²F](x) is affine in the measure, so ρ_T is convex in α²F. With
   ν = α²F/A this gives the extreme-point bound ρ_T(α²F) ≤ sup_ω ρ(A/T, ω/T) (Einstein family).
3. **Convexity in a = A/T** with ρ(0, w) = −1.
4. **Interval-arithmetic (Arb) certificates** for the Einstein family on w ∈ [1.3, 40], with
   analytic end-lemmas for w ≤ 1.3 and w ≥ 40.

The Lean theorem `TcR1.unstable_imp_lt` turns the hypothesis H(a0, η) (Prop. 6, a0 = 4.33148145,
η = 10⁻⁹) into "T unstable ⇒ T < A/a0" in both sectors.

**Corollary (core screening rule):** Tc ≥ 300 K ⇒ A ≥ 1299.44 K. The code uses
a0 = 4.33148145 exactly (`PROOFS.md` §6), i.e. A ≥ 1299.44 K.

Literature: this is the setting of Leavens' least-upper-bound problem. C. R. Leavens, "A least
upper bound on the superconducting transition temperature", Solid State Commun. **17**, 1499
(1975) [CIT: title, journal, volume, page and DOI 10.1016/0038-1098(75)90982-5 checked; the
paper's constant could not be retrieved]. From memory, Leavens' result is k_B Tc ≤ 0.2309 A for
μ* = 0 [UNVERIFIED]. If correct, the core certificate is a rigorous, sharp-to-10⁻⁸ version
of Leavens' bound (1/4.33148 = 0.230868), extended to μ* ≥ 0 and arbitrary cutoff.

### 2.2 Other certified rules
- **R-CELLS** [THM]: Arb certificates of stability on temperature cells [T₁, T₂] for a *given*
  atomic spectrum and envelope (1 ± e) on the weights. This gives a certified upper bound on
  Tc_model for that spectrum, for every μ* ≥ 0.
- **C-RAYLEIGH** [THM]: an Arb-checked test vector x with boxMinQ ≥ 0 at T. This certifies
  "T unstable" for every spectrum in the envelope, every μ* ∈ [0, μ_hi] and cutoff N_c(T), so
  Tc_model ≥ T. This lower bound depends on μ* and ω_c; the upper bounds do not.
- **R-MU / R-ODD** [THM]: μ* only lowers the even form; the odd sector is bounded by the even
  μ* = 0 form.
- **Multiband, R-MULTI** [THM, computer-assisted through Prop. 6]: Tc ≤ κ/4.33148145, where κ is
  the Perron root of K_ij = sqrt(N_i/N_j) 𝒜_ij and 𝒜_ij = ∫α²F_ij. The constant is sharp.
  The isotropic or band-averaged A is **not** a valid bound: Tc·a0/A_avg is unbounded, and an
  MgB2-like model beats it by 9 %. See §5.1.

---

## 3. New results: extremal problems for Tc = 300 K

### 3.0 Solver validation [NUM] (`data/validation.txt`)

| check | result |
|---|---|
| Einstein μ* = 0 optimum: min_w A/Tc | a* = 4.3314821 at w* = 7.57486, λ* = 1.14365. The certified enclosure is [4.3314814, 4.3314815]; the 1.4×10⁻⁷ relative excess is from float tolerance and the truncation sign (truncation only lowers ρ). |
| Large-λ limit Tc/(Ω√λ), Einstein, μ* = 0 | 0.1146 (λ = 1), 0.1764 (10), 0.1821 (100), 0.18266 (1000), increasing. The soft-kernel limit λ(k) = B/k² gives **c₂ = 0.18272625** (`bounds_extra.py`). This equals the rigorous asymptotic constant C∞ = 0.1827262477… of Kiessling–Altshuler–Yuzbashyan [CIT, KAY III]. |
| μ* = 0 vs Allen–Dynes, Pb-like (λ = 1.55) and Al-like (λ = 0.43) models | ratio ME/AD = 0.990 and 0.999 |
| μ* > 0 vs Allen–Dynes | ME (ω_c = 5–10 ω_max) gives 4–13 % *higher* Tc than AD at λ ≥ 1.5. The gap grows at weak coupling (×1.3–2.0 at λ ≈ 0.4–0.5). This is the known μ*-convention gap (§4.3). |
| Truncation (trunc_mult 30 → 60) | relative change ≤ 1.4×10⁻⁷ |
| Random μ* = 0 spectra (300 spectra, 1–5 atoms, 50–3000 K) | max Tc/A = 0.23073 (< 0.230868); max Tc/sqrt(λ⟨ω²⟩) = 0.18188 (< 0.18273) |

### 3.1 μ*-dependent version of the A bound (N1μ) [DER]

> **Update (03_certificates.md §4):** the floors at μ* ∈ {0.10, 0.13} and
> ω_c ∈ {5000, 10000, 20000, 50000} K, including the T-supremum step, are now [THM-C] (Arb
> certificates, symmetric-block convention). The float discussion below is kept for the other values.

**Claim.** Fix T > 0, μ* ≥ 0 and an *absolute* cutoff ω_c, so that N_c = N_c(T). Let
a_E(T; μ*, ω_c) := inf{A/T : ∃Ω > 0 with A δ_Ω unstable at T}. If α²F is unstable at T (in
particular, if it superconducts at T with a second-order transition), then A ≥ T·a_E(T; μ*, ω_c).

**Proof.** With ν = α²F/A,

Q⁺_T[α²F](x) = −Σ(2n+1)x_n² − 2μ*(Σ_{n<N_c} x_n)² + ∫ q_x dα²F = ∫ Q⁺_T[A δ_ω](x) dν(ω),

because the bare and μ* terms do not depend on the spectrum and ν is a probability measure.
Hence ρ_T(α²F) ≤ sup_ω ρ_T(A δ_ω).

Suppose A < T·a_E. Then every A δ_ω is stable at T. The sup is attained or approached only in
the interior, because ρ_T(A δ_ω) → −1 − O(μ*) both as ω → 0 (λ(k) → 0) and as ω → ∞ (λ(k) → 0).
By continuity the sup is < 0, so α²F is stable. The odd sector is never unstable (§1). ∎

**From "unstable at T" to "Tc ≥ 300 K".** Tc is defined as sup{T unstable}. So Tc ≥ 300 K only
gives instability at *some* T ≥ 300 K, or at points arbitrarily close to 300 K from above. The
requirement is therefore A ≥ inf_{T ≥ 300 K} T·a_E(T).

- For μ* = 0 the problem is scale-free, so T·a_E is linear in T and minimal at 300 K. [THM]
- With an absolute ω_c and μ* > 0, N_c(T) drops by one at T_j = ω_c/((2j+1)π), and a_E(T) jumps
  down there. R01's `reentrance_check.py` scanned T ∈ [300, 420] K on both sides of every jump
  (μ* = 0.10 and 0.13; ω_c = 5000, 10000 and 20000 K). Our `g_table.py` did the same scan for g
  and for the soft-mode requirement (`data/tsup_check.txt`). In every case the minimum is at
  T = 300 K. [NUM]

The proof is exact. The numbers a_E in Table 2 are float computations (`scan_targets.py`: ρ = 0
root in λ at T = 300 K, minimized over Ω/T). By the proof, the Einstein spectrum is the least
demanding shape at every μ* and absolute ω_c. R01's differential evolution over 2–3 atoms
confirmed this (ratio 1.00000). Table 1 uses a different cutoff convention (ω_c = 10 ω_max of each
shape), so a comparison across tables is only indicative. Within Table 1, Einstein needs the least
A in every column.

*Table 2 (excerpt, `data/einstein_abs.csv`). Minimal A for Tc = 300 K; minimizer Ω_opt ≈ 7.0–7.6 Tc:*

| μ* | ω_c = 5000 K | 10000 K | 20000 K | 50000 K | λ_opt (20000 K) |
|---|---|---|---|---|---|
| 0 | 1299.4 | 1299.4 | 1299.4 | 1299.4 | 1.144 |
| 0.08 | 1598 | 1627 | 1628 | 1612 | 1.49 |
| 0.10 | 1672 | 1706 | 1703 | 1678 | 1.58 |
| 0.13 | 1782 | 1824 | 1811 | 1769 | 1.70 |
| 0.16 | 1891 | 1940 | 1914 | 1853 | 1.82 |
| 0.20 | 2036 | 2091 | 2044 | 1954 | 1.97 |

The values are not monotone in ω_c at fixed μ*. At ω_c = 5000 K only N_c = 3 Matsubara
frequencies feel μ*, and for ω_c → ∞ the same number μ* means a weaker bare repulsion. The
dependence on ω_c is ≤ 7 %.

### 3.2 A bound in λ⟨ω²⟩: the λ → ∞ end of rule NM [DER]

> **Update (03 §2):** the λ → ∞ bound is now [THM-C]. It covers every α²F, not only Einstein
> spectra: Tc ≤ 0.18272717·sqrt(λ⟨ω²⟩) = C∞(1 + 5.0×10⁻⁶)·sqrt(λ⟨ω²⟩). This proves KAY III
> Conjecture 1 up to that factor, in non-strict form. The open step described below is closed.

This is g(∞) of §0.3. The fixed-λ version is in §3.3. Together they form the single rule NM.

**Derivation.** Let c := ∫ω dα²F, so λ⟨ω²⟩ = 2c. Set ν := ω α²F / c. The same affinity argument
gives Q⁺_T[α²F](x) = ∫ Q⁺_T[(c/ω) δ_ω](x) dν(ω). So it suffices to bound the Einstein family
with fixed λΩ² = 2c, i.e. weight c/Ω at Ω and λ = 2c/Ω².

As Ω → 0, λ_T(k) → B/k² with B = λ⟨ω²⟩/(4π²T²). This is the "soft-mode" or λ → ∞ limit. Let
B_c(μ*, ω_c, T) be the critical B of this kernel (`bounds_extra.py`, converged to 10⁻⁹ in N).
Then *if the Einstein family at fixed λΩ² is least stable in the limit Ω → 0*:

  **Tc ≤ c₂ · sqrt(λ⟨ω²⟩), with c₂ = 1/(2π sqrt(B_c)).**

The open step is this monotonicity. For μ* = 0 it is exactly Conjecture 1 of KAY III: Einstein
Tc ≤ C∞ Ω√λ for all λ, with C∞ = 0.1827262477…. That paper proves only the asymptotics and a
rigorous upper bound about 2.03 × larger [CIT]. Our evidence:
- Einstein Tc/(Ω√λ) increases monotonically in λ: 400-point dense grid over λ ∈ [0.8, 300],
  `g_table.py`;
- the 300 random spectra stay below 0.18273;
- for μ* > 0 the required λΩ²/Tc² decreases monotonically along Ω/Tc = 8, 4, 2, 1, 0.5, 0.25
  toward the soft-mode value (`einstein_abs.csv`).

*Table 3 (`data/soft_mode_constant.csv`, ω_c = 20000 K, T = 300 K):*

| μ* | c₂ | λ⟨ω²⟩ required | sqrt(λ⟨ω²⟩) required | ∫ω α²F dω required |
|---|---|---|---|---|
| 0 | 0.18273 | 2.696×10⁶ K² | 1642 K | 1.348×10⁶ K² |
| 0.08 | 0.16503 | 3.305×10⁶ | 1818 | 1.652×10⁶ |
| 0.10 | 0.16173 | 3.441×10⁶ | 1855 | 1.720×10⁶ |
| 0.13 | 0.15737 | 3.634×10⁶ | 1906 | 1.817×10⁶ |
| 0.16 | 0.15357 | 3.816×10⁶ | 1954 | 1.908×10⁶ |
| 0.20 | 0.14920 | 4.043×10⁶ | 2011 | 2.021×10⁶ |

Physical reading: A and λ⟨ω²⟩ = 2∫ωα²F are the two moments that cannot be "padded" by soft modes.
A mode with ω ≪ 2πTc adds to λ, and lowers ω_log, but changes λ_T(k) only by λ_mode ω²/(2πkT)²,
which is negligible.

### 3.3 Fixed λ: rule NM, m ≥ g(λ) [DER + NUM]

> **Update (03 §2):** superseded by Theorem NM(λ_k, m₀) [THM-C]. A Lagrange-multiplier reduction
> to a one-parameter Einstein family, certified in Arb at every tabulated λ_k, proves that
> every α²F with λ ≤ λ_k and m ≤ g_lo(λ_k)(T/300)² is stable at T. The two-atom reduction, the
> Einstein-minimiser numerics and the monotonicity of g below are no longer needed for the rule.

(This replaces the earlier N3, "ω₂ ≥ Ω_E(λ_max) if λ ≤ λ_max". That gate needed the *true* λ to be
≤ λ_max, which is the unsafe direction (R01 F3).)

Fixing λ = ∫(2/ω)dα²F adds a second linear constraint. The extreme points of
{α²F ≥ 0 : λ fixed, λ⟨ω²⟩ fixed} are measures with ≤ 2 atoms, so ρ_T is maximized over
two-atom spectra. `bounds_extra.py` searches two-atom spectra over a grid (ratio r ≤ 20, λ share
f = 0.1…0.9) with the scale solved for Tc = 300 K, ω_c = 20000 K. **In every case the minimizer
of λ⟨ω²⟩ was the Einstein spectrum (r = 1).** So at fixed λ, ω₂ ≥ Ω_E(λ, μ*), the Einstein
frequency that gives 300 K. Equivalently, m ≥ g(λ) := λ Ω_E(λ)².

Ω_E(λ) decreasing is the statement that Tc/Ω is increasing in λ. KAY III proves it for
λ ≳ 0.77: the critical surface is λ = 1/h(ϖ) with ϖ = Ω/2πTc, and h is invertible for ϖ ≤ √2
[CIT, abstract; direction checked numerically].

g decreasing is the stronger statement that Tc/(Ω√λ) is increasing. It is not proven, but it
holds on the dense grid. Given g decreasing, an *upper* envelope λ_hi ≥ λ_true gives
g(λ_true) ≥ g(λ_hi). Hence the direction-safe rule **reject if m_hi < g(λ_hi)**, with the table
in §0.3. With λ_hi = ∞ it reduces to the λ → ∞ bound of §3.2.

The rule is robust to a spurious soft mode that inflates λ: such a mode adds λ but almost no m,
and it lowers g, which loosens the test. The earlier N3 form could instead have been triggered
by the lowered ω₂.

*Table 4 (excerpt, `data/two_atom_fixed_lambda.csv`): minimal ω₂ = Ω_E (K) for Tc = 300 K*

| μ* \ λ | 1.5 | 2 | 2.5 | 3 | 4 |
|---|---|---|---|---|---|
| 0 | 1774 | 1417 | 1212 | 1075 | 900 |
| 0.10 | 2273 | 1737 | 1453 | 1273 | 1051 |
| 0.13 | 2431 | 1829 | 1520 | 1326 | 1090 |
| 0.16 | 2596 | 1921 | 1585 | 1377 | 1127 |

**A at fixed λ is not constrained by λ.** For λ above λ_opt, the A-minimizing two-atom spectrum
puts the excess λ into a very soft atom (r → ∞; the grid edge r = 20 was selected) and
approaches the unconstrained Einstein minimum of §3.1.

### 3.4 Shape-resolved feasibility [NUM/SOFT] (Table 1, Table 5; `data/t_vs_lambda.csv`, `min_A.csv`)

The shapes are scale-free. For each shape and μ*, the frequency axis is stretched so that
Tc = 300 K. The cutoff is ω_c = 10 ω_max (Table 1 also gives the 5 ω_max change, which is
+0.2 % to +3.5 %). The hydride-like shapes are **schematic three-Gaussian models, not fits to
any calculation**:
- H3S-like: λ shares 0.12 / 0.55 / 0.33 at 20 / 105 / 175 meV;
- LaH10-like: 0.15 / 0.45 / 0.40 at 12 / 60 / 140 meV;
- MgB2-like: 0.35 / 0.65 at 35 / 70 meV.

*Table 1: minimal A (K) for Tc = 300 K (λ at the optimum in parentheses)*

| shape | μ* = 0 | 0.08 | 0.10 | 0.13 | 0.16 | 0.20 |
|---|---|---|---|---|---|---|
| Einstein | 1299 (1.14) | 1627 (1.49) | 1701 (1.57) | 1811 (1.70) | 1914 (1.82) | 2036 (1.88) |
| 2δ, r = 3, f = 0.5 | 1410 (1.48) | 1768 (1.94) | 1849 (2.09) | 1961 (2.15) | 2072 (2.41) | 2205 (2.26) |
| 2δ, r = 8, f = 0.5 | 1414 (2.12) | 1772 (2.78) | 1853 (2.93) | 1969 (3.16) | 2079 (3.38) | 2217 (3.64) |
| 2δ, r = 8, f = 0.8 | 1329 (1.40) | 1664 (1.83) | 1740 (1.93) | 1849 (2.07) | 1953 (2.18) | 2083 (2.33) |
| box [0.1, 1] ω_max | 1458 (1.61) | 1822 (2.11) | 1900 (2.20) | 2011 (2.43) | 2111 (2.45) | 2237 (2.71) |
| Debye ω² | 1349 (1.29) | 1686 (1.70) | 1761 (1.79) | 1867 (1.91) | 1966 (1.89) | 2090 (2.14) |
| H3S-like | 1365 (1.38) | 1701 (1.80) | 1774 (1.90) | 1875 (2.02) | 1969 (2.13) | 2084 (2.25) |
| LaH10-like | 1437 (1.74) | 1790 (2.29) | 1865 (2.32) | 1969 (2.54) | 2068 (2.76) | 2182 (2.71) |
| MgB2-like | 1347 (1.28) | 1684 (1.68) | 1758 (1.77) | 1864 (1.90) | 1962 (1.89) | 2086 (2.14) |

Within ±1 %, the minimal A is nearly independent of the cutoff convention and only weakly
dependent on shape: +2.6 % (2δ r = 8 f = 0.8) to +12 % (box) over Einstein. **A is therefore the
most robust single spectral target.** The location of the optimum in λ depends more on the
shape: broad and bimodal spectra want larger λ. The λ grid has steps from the discrete N_c, so
λ_opt values are good to about ±0.1.

*Table 5b: required ω_log (K) at fixed λ, μ* = 0.13, ω_c = 10 ω_max (Table 5a, μ* = 0.10, is in `tables.md`)*

| shape | λ = 1.5 | 2 | 2.5 | 3 | 4 |
|---|---|---|---|---|---|
| Einstein | 2421 | 1833 | 1528 | 1336 | 1099 |
| H3S-like | 2141 | 1563 | 1273 | 1096 | 885 |
| LaH10-like | 1934 | 1349 | 1064 | 896 | 702 |
| MgB2-like | 2350 | 1740 | 1430 | 1236 | 1005 |
| box | 2288 | 1632 | 1301 | 1099 | 865 |
| 2δ, r = 8, f = 0.5 | 2087 | 1338 | 1006 | 826 | 633 |

The spread in ω_log at fixed λ is a factor of about 1.6 between shapes. Table 6 shows it is
unbounded in general.

*Table 6 (`data/unbounded_wlog_demo.csv`): two atoms, λ = 2 split equally, frequency ratio r, μ* = 0*

| r | 1 | 4 | 16 | 64 | 128 |
|---|---|---|---|---|---|
| Tc/ω_log | 0.212 | 0.260 | 0.463 | 0.918 | 1.297 |
| Tc/A | 0.212 | 0.208 | 0.218 | 0.226 | 0.228 |
| Tc/sqrt(λ⟨ω²⟩) | 0.150 | 0.126 | 0.116 | 0.115 | 0.115 |

### 3.5 Bergmann–Rainer functional derivative [NUM + CIT]

δTc/δα²F(ω) is computed by Hellmann–Feynman on the top eigenvector:
δρ = xᵀ(∂M/∂α²F(ω))x and δTc = −δρ/(∂ρ/∂T) (`eliashberg.functional_derivative`).

*Table 7: at Tc = 300 K (ω_c = 10 ω_max)*

| spectrum | λ | μ* | peak ω/Tc | peak ω | peak value | value at ω = 2Tc | value at 20Tc |
|---|---|---|---|---|---|---|---|
| Einstein optimum | 1.144 | 0 | 7.56 | 2269 K | 0.2309 | 0.119 | 0.170 |
| H3S-like | 2.0 | 0.13 | 7.16 | 2149 K | 0.171 | 0.090 | 0.115 |
| LaH10-like | 2.5 | 0.13 | 7.16 | 2149 K | 0.173 | 0.091 | 0.118 |
| box | 2.0 | 0.13 | 7.26 | 2179 K | 0.188 | 0.099 | 0.130 |

At the Einstein optimum the peak value equals Tc/A = 1/a* = 0.2309 and sits at Ω* = 7.57 Tc.
This is the KKT condition of the extremal problem: wB Lemma 11 [THM]. The result
agrees with Bergmann & Rainer (1973): δTc/δα²F > 0 at all ω, largest slightly above 2πTc [CIT].

---

## 4. Established analytic results (literature)

### 4.1 McMillan / Allen–Dynes
- W. L. McMillan, Phys. Rev. **167**, 331 (1968) [UNVERIFIED citation details]. Tc formula for
  λ ≲ 1.5.
- P. B. Allen & R. C. Dynes, "Transition temperature of strong-coupled superconductors
  reanalyzed", Phys. Rev. B **12**, 905 (1975) [CIT].
  - Tc = f₁f₂ (ω_log/1.2) exp[−1.04(1+λ)/(λ − μ*(1+0.62λ))].
  - f₁ = [1 + (λ/Λ₁)^{3/2}]^{1/3}, Λ₁ = 2.46(1+3.8μ*).
  - f₂ = 1 + (ω₂/ω_log − 1)λ²/(λ² + Λ₂²), Λ₂ = 1.82(1+6.3μ*)(ω₂/ω_log).

  The form of f₁ and f₂ was checked against secondary sources [CIT, secondary], implemented as
  `eliashberg.allen_dynes`. The abstract states the large-λ scaling, Tc ≈ 0.15 sqrt(λ⟨ω²⟩) at
  μ* = 0.1, and that McMillan's "λ = 2 limit" is spurious [CIT, abstract; R01].

  AD's 0.15 is 7 % below our c₂(μ* = 0.10) = 0.1617 at ω_c = 20000 K. This is the same μ*
  convention gap discussed in §4.3.
- Asymptotic constant: Tc → C∞ sqrt(λ⟨ω²⟩) with **C∞ = 0.1827262477…**. This was proved for the
  Einstein case in M. K.-H. Kiessling, B. L. Altshuler, E. A. Yuzbashyan, "Bounds on Tc in the
  Eliashberg theory of superconductivity. III: Einstein phonons", J. Stat. Phys. (2025),
  art. 93, DOI 10.1007/s10955-025-03469-y, arXiv:2409.02121 [CIT].
  - arXiv lists the volume as "92". R01 suggests it is probably 192, so the volume is
    [UNVERIFIED].
  - KAY III notes that AD's value was 0.182, obtained non-rigorously.
  - KAY III also shows that the map ϖ = Ω/2πT ↦ top eigenvalue is invertible for ϖ ≤ √2. For
    λ ≳ 0.77 this makes Tc = Ω f(λ) unique with f strictly monotone; this is used in §3.3.
  - Part II (dispersive phonons), J. Stat. Phys. (2025), art. 94, DOI 10.1007/s10955-025-03468-z,
    arXiv:2409.00532 [CIT], gives convergent lower bounds and an upper bound ≈ 2.03 × the sharp
    constant.

  Back to our own result: Our c₂(μ* = 0) = 0.18272625 reproduces C∞ to 9 digits. The
  lower bound Tc ≥ (Ω/2π) sqrt(λ − 1) for λ > 1 (KAY III, also in AD) can be used as a cheap
  certificate of instability.

### 4.2 Bergmann–Rainer, Carbotte, Leavens, Kresin
- G. Bergmann & D. Rainer, "The sensitivity of the transition temperature to changes in
  α²F(ω)", Z. Phys. **263**, 59 (1973) [CIT]. δTc/δα²F > 0 everywhere, and most effective
  slightly above 2πTc. Our 7.2–7.6 Tc ≈ 1.15–1.2 × 2πTc agrees.
- J. P. Carbotte, "Properties of boson-exchange superconductors", Rev. Mod. Phys. **62**, 1027
  (1990) [CIT: abstract]. Organizes strong-coupling corrections by Tc/ω_ln. The asymptotic limit
  Tc/ω_ln → ∞ is "of theoretical interest only" [CIT]. A recalled value: the largest Tc/ω_ln
  among known conventional superconductors is ≈ 0.2–0.25 [UNVERIFIED]. Tc = 300 K at λ = 3
  needs Tc/ω_log = 0.22–0.36, and at λ = 2 it needs 0.16–0.22 (Table 5b). That is at or above
  the strongest-coupled known materials.
- C. R. Leavens, Solid State Commun. **17**, 1499 (1975) [CIT bibliographic; constant
  UNVERIFIED]. See §2.1.
- V. Z. Kresin, "On the critical temperature for any strength of the electron-phonon coupling",
  Phys. Lett. A **122**, 434 (1987) [CIT bibliographic].
  - Formula as quoted in Kresin & Wolf, arXiv:0904.2038 (Rev. Mod. Phys. **81**, 481 (2009),
    journal details [UNVERIFIED]): Tc = 0.25 Ω̃ / (e^{2/λ_eff} − 1)^{1/2}, with μ* neglected
    [CIT, secondary].
  - Large λ: Tc → 0.25 Ω̃ sqrt(λ/2) = 0.177 Ω̃ √λ [DER]. With Ω̃ = sqrt⟨ω²⟩ [UNVERIFIED], this
    is 3 % below C∞ = 0.1827.
  - It is an interpolation formula, not a bound.

### 4.3 μ* and its renormalization
- P. Morel & P. W. Anderson, Phys. Rev. **125**, 1263 (1962) [CIT]:
  μ* = μ/(1 + μ ln(E_el/ω_ph)). This is the Bogoliubov–Tolmachev logarithm (Bogoliubov,
  Tolmachev & Shirkov 1958 [UNVERIFIED]).
- Consequence used here: μ* has meaning only together with ω_c, through
  1/μ*(ω₂) − 1/μ*(ω₁) = ln(ω₁/ω₂).
- The lowest-order Morel–Anderson result is uncontrolled for large μ. Bauer, Han & Gunnarsson
  (2012–13, arXiv:1202.5051) compute corrections [CIT: existence and scope from the abstract].
- **Convention gap** [NUM]: at the same number μ*, the ME solution with ω_c = 5–10 ω_max gives a
  4–13 % higher Tc than the AD formula for λ = 1.5–5 (Einstein scan in `validation.txt`). The gap
  grows to a factor 1.3–1.9 at λ ≈ 0.4–0.5.
  - Most of it is the Morel–Anderson renormalisation between our cutoff (ω_c = 10 ω_max) and the
    phonon scale implicit in AD. With ln 10 ≈ 2.3, μ* = 0.13 at 10 ω_max acts like about 0.10 at
    ω_max (R01 F9). Use `mu_rescale` to map a literature μ* onto our cutoff.
  - The μ*-dependent results (N1μ, Tables 1, 3–5) are stated in the permissive ME convention.
    N1, N1-MB and NM use μ* = 0 and are convention-free. Feeding AD-convention μ* values into
    the μ*-dependent tables is *not* conservative by up to ≈10 %.
  - For AD itself, the required ω_log at
  μ* = 0.13 is 1997 / 1663 / 1458 K at λ = 2 / 2.5 / 3 (Einstein), and the AD minimal A is
  1982 K (μ* = 0.13) [NUM, scratch computation with `allen_dynes`].
- Typical hydride practice: μ* = 0.10–0.16 [UNVERIFIED as a survey statement].

---

## 5. Hard vs soft, and what lies outside the model

### 5.1 Hard rules, and why anisotropic or multiband materials need different inputs

**Within a single isotropic α²F:**
- **N1** is certified for any μ*, cutoff and shape.
- **N1μ**: exact reduction to the Einstein family; float numerics for a_E.
- **NM** has three ingredients:
  - an exact reduction to ≤ 2 atoms;
  - the Einstein minimizer: grid search plus R01's differential evolution;
  - g decreasing: dense-grid numerics, stronger than KAY III's proven monotonicity of Tc/Ω.

**The isotropic α²F from DFPT is the band average α²F_iso = Σ_ij (N_i/N) α²F_ij, and none of N1, NM
is a necessary condition for it in a multiband or anisotropic material.**

*Decoupled two-band counterexample* (R01 F1, reproduced by its independent solver):

- Two sheets have equal DOS (w₁ = w₂ = ½).
- Band 1 has an Einstein α²F₁₁ with λ = 3 and Ω = 1075.3 K. Alone, it gives Tc = 300.0 K at μ* = 0
  (Table 4).
- Band 2 is uncoupled.

The average then has:

- λ = 1.5;
- ω₂ = 1075 K;
- m = 1.73×10⁶ K²;
- A = 806 K.

It fails N1 (806 < 1299.4), NM (1.73×10⁶ < g(1.5) = 4.72×10⁶) and the old N3. With a DOS
fraction w for the coupled band, A_iso and m_iso fall by the factor w, without limit. The
interband coupling is continuous, so weakly coupled versions also fail. Only κ = 𝒜₁₁ = 1613 K
gives the right answer.

**Valid N1 inputs in a multiband model** (constant 4.33148145 unchanged, wMB Thm 1 [THM]):

1. **κ.** The Perron root of K_ij = sqrt(N_i/N_j)𝒜_ij. This is the best input. It needs the band-
   or k-resolved α²F_ij. A k-point discretization is covered by wMB §6.
2. **max_i A_i**, with row sums A_i = Σ_j 𝒜_ij. Valid because κ ≤ Λ_A ≤ max_i A_i (wMB).
3. **A_iso / w_min**, where w_i = N_i/N are the sheet DOS fractions. [THM, derived here]

**Proof of (3).** Let S_ij := N_i 𝒜_ij. S is symmetric (by N_i α²F_ij = N_j α²F_ji) and has
non-negative entries, with Σ_ij S_ij = N·A_iso. Write K = D^{-1/2} S D^{-1/2} with D = diag(N_i).
For any x put y = D^{-1/2}x. Then

  xᵀKx = yᵀSy ≤ λ_max(S)|y|² ≤ λ_max(S)|x|²/N_min,

so κ ≤ λ_max(S)/N_min. For a non-negative symmetric matrix, λ_max(S) ≤ max row sum ≤ Σ_ij S_ij.
Hence

  **κ ≤ N·A_iso/N_min = A_iso/w_min.** ∎

Rule (3) needs only the isotropic α²F and the sheet-resolved DOS at E_F, which a band-structure
run gives cheaply. It assumes the coupling is constant within each sheet (recorded assumption
M-SHEET-ISO). It is weak but direction-safe: in the counterexample, 806/0.5 = 1612 ≥ 1299.4, so
there is no rejection.

**Intra-sheet anisotropy** (one sheet, k-dependent coupling) can still make κ > A_iso. There is no
bound without k-resolved data, which is why N1-ISO and NM carry the assumption M-ISO and use the
25 % margin clause (§0.2c).

**Multiband versions of NM do not exist yet** [open]. Multiband materials are therefore decided
only by N1-MB, or they are UNDECIDED.

### 5.2 Soft (empirical or shape-dependent)
- ω_log thresholds (no bound exists, Table 6).
- The "realistic" columns of Tables 1 and 5 (schematic shapes).
- Allen–Dynes accuracy: about 1 % at μ* = 0 in our tests, and a 5–13 % convention gap at μ* > 0.
- The λ window 2–3.
- ε_sys = 0.10 in the envelope.
- The thresholds in G-VALID (0.1, 30 K, 10 %).

### 5.3 Outside the model: lattice stability and the "λ ceiling"

**Lattice-stability arguments**

- **Cohen & Anderson (1972)** argued that lattice stability (positive static dielectric response)
  limits λ and caps phonon-mediated Tc at roughly 10–30 K. Reference: M. L. Cohen &
  P. W. Anderson, "Comments on the maximum superconducting transition temperature", in
  *Superconductivity in d- and f-Band Metals*, AIP Conf. Proc. **4**, 17 (1972) [UNVERIFIED
  page; content via secondary sources CIT: Sadovskii, Phys.-Usp. 65, 724 (2022),
  arXiv:2106.09948, and Weger & Birman, arXiv:cond-mat/0205279].
- Rebuttal:
  - Dolgov, Kirzhnits & Maksimov, Rev. Mod. Phys. **53**, 81 (1981) [CIT]: a negative static
    dielectric function at finite q is compatible with stability.
  - Empirically, Pb–Bi alloys and the hydrides have λ > 2.
  - H3S (Drozdov et al., Nature **525**, 73 (2015) [UNVERIFIED]) and LaH10 at about 250 K
    (Drozdov et al., Nature **569**, 528 (2019); Somayazulu et al., PRL **122**, 027001 (2019)
    [UNVERIFIED]) far exceed the 1972 estimate.
- **Moussa & Cohen**, Phys. Rev. B **74**, 094520 (2006) [CIT, abstract]: two bounds, (i) from
  the maximal doping-tuned coupling and (ii) from phonon softening to instability at strong
  coupling. MgB2 reaches the first.

**ME breakdown at large λ**

- Esterlis, Nosarzewski, Huang, Moritz, Devereaux, Scalapino & Kivelson, Phys. Rev. B **97**,
  140501(R) (2018) [CIT]. In the 2D Holstein model, ME agrees with DQMC only up to a *bare*
  dimensionless coupling λ₀ ≈ 0.4. This is not the renormalized λ used elsewhere in this
  document. CDW order and bipolarons take over beyond it.
- Esterlis, Kivelson & Scalapino, Phys. Rev. B **99**, 174516 (2019) [CIT, abstract per R01]:
  a metal → pseudogap crossover, with SC at small λ₀ and CDW at large λ₀. The abstract states no
  numerical threshold.
- The statement "Tc peaks near λ ∼ 1" is probably from Esterlis, Kivelson & Scalapino, npj Quantum
  Mater. **3**, 59 (2018) [UNVERIFIED].
- Yuzbashyan & Altshuler, PRB **106**, 054518 (2022) [CIT]: ME theory itself breaks down at
  λ_c ≈ 3.0–3.7, where the normal-state specific heat turns negative near 3.7.
- Reading: the Holstein numbers concern models with Ω comparable to the bandwidth and strong
  local coupling. Realistic broad-band hydrides with Ω/E_F ≈ 10⁻² operate at λ ≈ 2–3.5 with ME
  agreeing reasonably with experiment. We therefore adopt **λ ≤ 3 as the design ceiling** and
  treat 3–3.5 as a risk zone [SOFT]. This ceiling is used for *ranking* only. The reject rules
  do not assume it.

**Vertex corrections and non-adiabaticity**

- Pietronero, Strässler & Grimaldi, Phys. Rev. B **52**, 10516 (1995) [CIT].
- Corrections are of order λΩ/E_F. In that literature they are argued to be able to *raise* Tc.
  So a reject rule outside the adiabatic regime is not a theorem, which is the reason for G-VALID
  item 1.
- For hydrides (Ω ≈ 0.2 eV, E_F ≈ 10 eV, λ ≈ 2–3) the ratio is ≈ 0.05 [DER, order of magnitude].
  For flat-band or low-E_F candidates it is **not** small.

### 5.4 Anharmonicity and quantum nuclear effects
- In hydrides, harmonic DFPT phonons are unreliable. Errea et al., Nature **578**, 66 (2020)
  [CIT]: quantum fluctuations stabilize Fm-3m LaH10, with λ = 3.5 in the quantum structure.
- For H3S, anharmonicity changes λ by tens of percent (Errea et al., PRL **114**, 157004 (2015)
  [UNVERIFIED]).
- Consequence (G-VALID item 2): a harmonic-only α²F of an H/Li/Be compound, or of any structure
  with imaginary modes, cannot be hard-rejected. Its envelope is not controlled. If an anharmonic
  check exists, its frequency shift enters ε_ω (§0.2b).

### 5.5 The ambient-pressure context [CIT]

Gao, Cerqueira, Sanna, Fang, Dangić, Errea, Wang, Botti & Marques, "The maximum Tc of conventional
superconductors at ambient pressure", Nat. Commun. **16**, 8253 (2025),
DOI 10.1038/s41467-025-63702-w, arXiv:2502.18281. Their data cover >20,000 metals:
- ω_log rarely exceeds ≈ 1800 K.
- There is a λ–ω_log trade-off.
- The paper names ClB₂C₈ and Al₂OsH₇ (λ ≈ 3.3) as its largest-λ compounds. Yet it also gives
  λ ≈ 4 for Li₂AgH₆. The source is internally inconsistent on this point.
- Li₂AgH₆: Tc = 108.8 K (Eliashberg), 0.319 eV/atom above hull. Li₂AuH₆: 116.1 K, 0.172 eV/atom.
- The authors suggest hypothetical compounds reaching about 100–120 K can be designed. They say
  higher Tc is not forbidden by physical law, but room temperature at ambient pressure is
  "extremely unlikely". This is a judgement, not a bound.

Compared with our box: Tc = 300 K at λ = 3 needs ω_log ≈ 900–1340 K (shape-dependent), which on
its own lies inside their ω_log range. The hard part is the *joint* requirement A ≳ 1900–2050 K,
i.e. λ ≈ 3 *and* ⟨ω⟩ ≳ 1300 K at once. Their observed trade-off says this combination has not been
found in 20,000 metals. **A (or κ), not λ or ω_log separately, is the quantity to condition and
screen on.**

### 5.6 Error propagation under frequency errors [DER]

For a phonon mode ν, λ_ν = 2N_F⟨|g_ν|²⟩/ω_ν, and |g_ν|² = ⟨|∂V|²⟩/(2M_ν ω_ν). So

  λ_ν = N_F⟨|∂V|²⟩_ν / (M_ν ω_ν²)  and  λ_ν ω_ν² = N_F⟨|∂V|²⟩_ν/M_ν

(the Hopfield combination).

α²F = Σ_ν (λ_ν ω_ν/2) δ(ω − ω_ν), so A = Σ λ_ν ω_ν/2 and m = Σ λ_ν ω_ν². Suppose the frequencies
are wrong by a common factor s at fixed electronic matrix elements (anharmonic shift, mass or
force-constant error). Then:

- λ → λ/s²;
- A → A/s;
- m → m (unchanged).

An amplitude error (1+ε) multiplies λ, A and m alike. This gives the propagation formulas of
§0.2b. A uniform 20 % "1.2 × nominal" envelope has no meaning for ω₂ = sqrt(m/λ): ω₂ is invariant
under amplitude rescaling and changes as s under the frequency error. ω₂ is therefore no longer
used in any reject rule.

---

## 6. Open problems / what was not verified
1. ~~NM depends on two numerical facts~~ **Closed** (03 §2): NM-G is [THM-C] at every tabulated
   λ_k and at λ = ∞, re-checked by two independent checkers. What remains is that the new lemmas
   are proved on paper only (no Lean), and that the extension to μ* > 0 holds only in the
   symmetric-block convention.
2. N1μ is [THM-C] at μ* ∈ {0.10, 0.13} × ω_c ∈ {5000, 10000, 20000, 50000} K (03 §4). The other
   μ* > 0 tables are float computations, not interval-certified.
3. A multiband version of NM is not proved. Multiband materials can only be rejected via N1-MB.
4. A bound for intra-sheet anisotropy without k-resolved data does not exist (§5.1). The
   25 %-margin clause is a policy choice, not a theorem.
5. The hydride-like shapes are schematic. Replace them with computed α²F once available (e.g. the
   running H3S DFPT benchmark).
6. UNVERIFIED items:
   - the Leavens constant 0.2309;
   - McMillan, Drozdov (2015, 2019), Somayazulu and Errea 2015 citation details;
   - the Cohen–Anderson page range;
   - the Kresin–Wolf RMP volume/page and Ω̃ = sqrt⟨ω²⟩;
   - the Carbotte "Tc/ω_ln ≲ 0.25" statement;
   - the hydride μ* survey statement;
   - the KAY J. Stat. Phys. volume;
   - the Esterlis–Kivelson–Scalapino npj 2018 claim.
7. Superconductivity at 300 K is identified with "the normal state is unstable at some
   T ≥ 300 K". For a first-order transition this would not hold, and such a transition is
   outside the model.

---

## 7. Response to review R01

| finding | change made (or disagreement) |
|---|---|
| **F1** critical: floors are invalid for band-averaged α²F | Accepted, with the decoupled two-band example added (§5.1). Isotropic rules (N1-ISO, NM) now require n_sheets = 1, carry the assumption M-ISO, and need k-resolved data or a > 25 % margin. New N1-MB rules: κ; max row sum; and **κ ≤ A_iso/w_min**, a new short theorem (§5.1) that needs only the sheet DOS fractions. Multiband without any of these inputs gives UNDECIDED. I did not adopt a fixed δ_an: the counterexample shows the factor is unbounded (1/w), so no universal δ_an is justified, whereas A_iso/w_min is. |
| **F2** major: "1.2 × nominal" is ill-defined | Replaced by a per-candidate ladder (≥ 3 admitted smearings and 2 q-meshes) giving δ_A, δ_λ, δ_m, plus ε_sys = 0.10 and a frequency envelope ε_ω. Propagation (§0.2b) is derived in §5.6: A ∝ 1/s, λ ∝ 1/s², m invariant. A missing or unconverged ladder (δ > 0.5) gives UNDECIDED. ω₂ is no longer used in any reject rule. |
| **F3** major: N3's λ gate is in the unsafe direction | N2 and N3 are replaced by the single rule **NM: reject if m_hi < g(λ_hi)**, with g tabulated on a 31-point grid plus ∞ (§0.3, `g_table.py`, `data/g_table.csv`) and a lookup at the next larger tabulated λ. The proof status is split by step. KAY III is cited for the monotonicity of Tc/Ω for λ ≳ 0.77 (proven). Monotonicity of g is a separate numerical assumption: dense grid, all 399 differences > 0, max Tc/(Ω√λ) = 0.18252 < C∞. |
| **F4** major: no gate on model validity | Added **G-VALID** (§0.2a, YAML): non-adiabatic λΩ_max/E_F > 0.1, imaginary modes, H/Li/Be with harmonic-only phonons, unrecomputed anharmonic shifts > 10 %, μ* < 0 suspected, a missing ladder, or multiband without data → UNDECIDED. |
| **F5** Gao et al. | The article number is now 8253 and the DOI is added. The λ ≈ 3.3 vs 4 inconsistency is noted. "Ceiling" is softened to the source's wording. |
| **F6** Esterlis | PRB 97 is now quoted with the *bare* λ₀ ≈ 0.4. The PRB 99 statement is replaced by its abstract content. "Tc peaks near λ ∼ 1" is attributed to npj QM 3, 59 (2018) and marked UNVERIFIED. |
| **F7** KAY journal references | J. Stat. Phys. art. 93 (III) and art. 94 (II) with DOIs are added; the volume is UNVERIFIED (arXiv lists "92"). The proven monotonicity of Tc/Ω is cited in §3.3 and §4.1. |
| **F8** typo | The N1μ range at μ* = 0.10 is now 1672–1706 K. |
| **F9** μ* cutoff convention | §1 now states the symmetric-block choice and that it is permissive by 0.13–2.6 %. §4.3 explains the AD gap by Morel–Anderson renormalisation and points to `mu_rescale`. |
| **F10, F11, F13** | No change needed. R01's independent reproductions are cited where relevant. |
| **F12** T-supremum step | Added to §3.1. The μ* = 0 case is exact (scaling). For μ* > 0, R01's `reentrance_check.py` and our `g_table.py` → `data/tsup_check.txt` (g and the soft-mode requirement, both sides of every N_c jump, T ∈ [300, 420] K) put the minimum at T = 300 K in all cases. |
| **F14.1** rounding | All floors now use 1299.44 K, scaled by Tc_target/300. |
| **F14.2** aim-box arithmetic | A aim is now ≥ 2060 K (1.1 × 1875). The ω_log, ω₂ and m aims are recomputed with stated sources (§0.5). |
| **F14.3** mixed cutoff conventions | §3.1 now says that the cross-table comparison is indicative only. The claim rests on the proof and on R01's differential evolution. |
| **F14.4** AD 0.15 at μ* = 0.1 | Added to §4.1 together with the convention-gap remark. |

New or changed files:
- `src/theory/g_table.py` (new);
- `docs/theory/data/g_table.csv`, `g_table_mu.csv`, `tsup_check.txt` (new).

The solver code is unchanged: R01 found no errors in it.

### Response to RL1 (Lean review of §0.4)

| RL1 finding | change |
|---|---|
| 1 (domains / signs) | Added §0.6a and a WELL-FORMED gate that runs first. Imaginary modes are reported as the maximum abs(ω) over negative eigenvalues and the anharmonic shift as abs(shift). Positivity and ranges are stated for every field. Flags are derived from the fields present. |
| 2 (envelope convention not enforced) | G-VALID now requires ε_sys ≥ 0.10 and ε_ω ≥ max(0.02, abs(anharm shift), q-spread of ω_log). |
| 3 (assumptions not recorded) | Each outcome is now (status, fired_rule, assumptions, inputs_used with tiers). Per-rule assumption lists and common_assumptions are in the YAML. The margin clause is stated to be policy, not proof. |
| 4 (single tier) | §0.6c: one tier per input and the same run id. REJECTED_PROVEN requires every input used by the fired rule to be fullDFPT. |
| 5 (multiband envelopes) | §0.6b: 𝒜_hi (A-type scaling), κ_hi via monotonicity plus a Collatz–Wielandt certificate, rowA_hi, and w_min,lo = min w_i(1−δ_w). How each is obtained from DFPT/EPW output is stated. |
| 6 (Tc ≤ T*) | Outcomes now claim Tc ≤ Tc_target only. |
| 9 (constants) | a0 = 4.33148145, η = 10⁻⁹ and the kernel 2πkT are stated at the top of the YAML. |
| 12 (certified moment floor) | Added NM-CERT: m_hi < 4(a0Tc)²/λ_hi, [THM] under M-ISO, via A² ≤ λm/4. |
| 7, 8, 10, 11 | Lean/Python implementation points: fail-closed tags, B ≥ 1, wiring gUp to the table, derived flags. The doc now states derived flags and the g lookup rule. The rest is left to the Lean and Python authors. |

### Response to R03 (review of the certificates, 03_certificates.md)

R03 found no soundness defect. The changes made in this document (the full response is in 03 §7)
are:

| item | change in 01 |
|---|---|
| certified g | The §0.3 table is now the certified g_lo(λ_k) = 90000·m₀ on the full 31-point grid plus ∞, at δ = 10⁻⁵ below the float values, with Ω_E brackets. Implementations read the exact m₀ from the certificate JSON. |
| NM-G status | NM-G is [THM-C] with assumptions [M-ISO, CERT-NM]. NUM-EINSTEIN-MIN and NUM-G-MONOTONE are dropped. The g(∞) constant is 2.695481×10⁶ K² (certified), replacing 2.6955×10⁶ (float, 9×10⁻⁶ higher). |
| N1μ | Certified floors at the listed (μ*, ω_c) are in §0.1. It is still a ranking rule. |
| finding 3 (Coulomb convention) | §0.1 now states that "μ* ≥ 0 only lowers the even form" is proved in the symmetric-block convention only. |
| §3.1–3.3, §6 | Update notes point to 03. Open problems 1–2 are updated. |
| Lean review additions | §0.6a has two new WELL-FORMED checks: N1-DOS consistency Σ w_i(1−δ_w) ≤ 1 (or normalise before δ_w), and the moment-consistency check A_nom² ≤ λ_nom·m_nom/4·(1+10⁻²). |
