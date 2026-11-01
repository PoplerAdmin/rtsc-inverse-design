# lean/ — Lean 4 proofs for the targeted pipeline

Lean `v4.35.0`, Mathlib `v4.35.0` (commit `56a0058…`, pinned in `lake-manifest.json`).
Revision 3 (after review `docs/reviews/RL1_lean_review.md`, the post-RL1 update of
`docs/theory/01_eliashberg_targets.md` §0.4/§0.6, and the certified NM-G table [THM-C]).

| library | origin | content |
|---|---|---|
| `TcR1` | frozen core (sha256-pinned) | single-band bound `k_B Tc ≤ A/a0` (`TcR1.main`, `stable_of_ge`); certificate `H(a0, η)` is a hypothesis |
| `PipelineLogic` | frozen core | core `classify` logic, R-UNIV/R-CELLS/C-RAYLEIGH, `combine_multiband` mirror |
| `MultiBand` | frozen core | multiband bound `Tc ≤ κ/a0` (wMB Thm 1), `R_UNIV_MB*` |
| `TargetLogic` | **new** | the §0.4 rules: row-sum / DOS algebra, L-MOM1..3, NM-CERT, NM-G, the decision function `decideT` with per-input provenance and the downgrade-only tag policy, per-rule soundness, monotonicity, RL1 regressions |

The ledger state machine is not formalised in Lean; its invariants are SQL triggers.

## Status

* `lake build` (all four libraries): **success**, 2218 jobs (`logs/lake_build.log`).
* No `sorry`, `admit`, `axiom`, `native_decide`, `implemented_by`, `extern` in any `.lean` file
  (`scripts/check.sh` source scan).
* Axiom audit (`#print axioms`, `logs/*_axioms.txt`, `logs/summary.txt`):

  | check file | theorems | axioms used |
  |---|---|---|
  | `Check.lean` (TcR1) | 34 | `[propext, Classical.choice, Quot.sound]` |
  | `CheckPipeline.lean` | 82 | 75 use a subset of the three; 7 use none |
  | `CheckMultiBand.lean` | 67 | subset of the three |
  | `CheckTarget.lean` (TargetLogic) | 122 (every `theorem` in `TargetLogic/*.lean`) | subset of the three |

* Truth table for the Python mirror: `export/decide_cases.json` (563 cases, see below).

## Build / check / export

Build in an ext4 directory whose `.lake/packages` holds a prebuilt Mathlib `v4.35.0` (the
binary cache download is blocked here; a from-source Mathlib build takes hours):

```sh
# WSL. Pins lake to 2 CPUs (lake 5.0 has no -j flag), nice 19, LEAN_NUM_THREADS=2.
bash lean/scripts/build.sh "$PWD/lean" ~/lean_build            # rsync sources, lake build
bash lean/scripts/check.sh ~/lean_build ~/lean_build/axioms    # #print axioms + source scan
cd ~/lean_build && lake env lean ExportCases.lean              # writes export/decide_cases.json
```

---

## Hypotheses vs proved

**Never proved in Lean (hypotheses):**

| tag (YAML) | Lean | content |
|---|---|---|
| CERT-H | `TargetLogic.CertH` = `H tpPhys a0Cert etaCert ∧ PropN tpPhys a0Cert` | interval-arithmetic certificates of the core (`rtt.core.certify`) (Prop. 6, wMB Prop. N) at the **pinned** constants `tpPhys = 2π`, `a0Cert = 4.33148145`, `etaCert = 1e-9` |
| ENV-SOUND | explicit per-rule hypotheses (below) | the `*_hi` / `*_lo` numbers bound the true values in the stated direction |
| M-ISO, M-SHEET-ISO, SHEETS-CORRECT, M-DISCRETE | the choice of model in each soundness theorem (below) | which mathematical model represents the material |
| CERT-NM | `TargetLogic.CertNMSet` = `CertNM λ_k m0_k` for every entry of `NMTable` and `∀ λ, CertNM λ m0_∞` | Theorem NM(λ_k, m0) of doc 03 §2.2 (Arb certificates `docs/theory/certificates/nm_*.json`, paper Lagrangian reduction): every α²F with `λ ≤ λ_k` and `m ≤ m0·T²` is stable at `T` (μ* ≥ 0, any cutoff). Replaces NUM-EINSTEIN-MIN / NUM-G-MONOTONE / NMHyp. |

Everything else below is proved.

## A. `TargetLogic.Matrix` — linear algebra for N1-ROWMAX / N1-DOS

* `le_rowSum_of_super`, `eig_le_of_rowSum_le`: `M ≥ 0`, `v ≥ 0`, `v ≠ 0`, `M v = κ v` ⇒
  `κ ≤ max row sum`. `rowSum_le_weighted_div`, `eig_le_iso_div_wmin`: weights `w_i ≥ w_min > 0`
  ⇒ `κ ≤ (Σ_i w_i rowSum_i)/w_min`.
* Connection to `MultiBand` (`𝒜 = MultiBand.Amat`): `Aiso_eq` (`A_iso = Σ_i w_i A_i`),
  `rowA_single` (B = 1: row sum = A_iso), `rowA_le_dos`, `kappa_le_rowmax`, `kappa_le_dos`,
  composite `n1_rowmax_stable`, `n1_dos_stable`, `stable_of_kappa_bound`, `tc_le_rowmax`,
  `tc_le_dos` (for general `tp`, `a0`, with `H`, `PropN` as hypotheses).

## B. `TargetLogic.Moments` — L-MOM1..3 (documentation lemmas)

`lmom1`: `(∀ j, Ω_j ≤ ω_max) ⇒ A ≤ λ ω_max/2`; `lmom2_sq`: `A² ≤ λm/4`; `lmom2`: `λ > 0 ⇒
A ≤ λ ω₂/2`; `lmom3`: `λ > 0 ⇒ λ ω_log/2 ≤ A`. Nonnegative finite atomic spectra only.

## C. `TargetLogic.NMG` — single-band rules

* `stable_of_zero` (proof from the RL1 scratch file): the zero spectrum is stable.
* `n1_iso_stable'`: `H`, `A ≤ A_hi < a0 T*` ⇒ stable for `T ≥ T*` (any `A ≥ 0`).
* **`nmcert_stable` (NM-CERT):** `H`, `m ≤ m_hi`, `λ ≤ λ_hi`, `λ_hi > 0`,
  `m_hi < 4(a0T*)²/λ_hi` ⇒ stable for `T ≥ T*`. Proof: `A² ≤ λm/4 < (a0T*)²` (L-MOM2), then N1.
  Fully certified (CERT-H only).
* `tableLookup tbl g∞ λ_hi` = `g_t` of the first `(λ_t, g_t)` with `λ_hi ≤ λ_t`, else `g∞`
  (the §0.3 rule; the table ends at λ = 100). `tableLookup_le`: conservative (`≤ g(λ)` for every
  `0 ≤ λ ≤ λ_hi`) under the table hypotheses. `nmg_stable`, `nmg_rule_stable`, `gLookup_le`:
  the NM-G logic.

## D. `TargetLogic.Decide` — the decision function (generic ordered field `K`)

Definitions follow the post-RL1 YAML (§0.4) and §0.6 literally; all tests are `Bool`, so the
**same** `decideT` is evaluated over ℚ for the export and reasoned about over ℝ.

* **Inputs.** `Validity` (G-VALID fields; `imagAbsMaxK`, `anharmShiftAbs` are magnitudes),
  `IsoData` (A, λ, m nominal, δ_A, δ_λ, δ_m, one `Prov` each for A, λ, m), `EpsData`
  (ε_sys, ε_ω, q-spread of ω_log, `Prov`), optional `MatData` (B, 𝒜_ij, DOS weights N_i,
  Collatz–Wielandt test vector v, δ_𝒜, `Prov`), optional `DosData` (sheet fractions w_i, δ_w,
  `Prov`), tags. `Prov = (tier ∈ {fullDFPT, coarseDFPT, ML, external}, run id)`.
* **Derived (never independent).** `mbPresent = mat ∨ dos present`; `kResolved = mat present`;
  envelopes `Ahi`, `lamHi`, `mHi` (§0.2b), `calAhi = 𝒜(1+δ_𝒜)(1+ε_sys)/(1−ε_ω)`,
  `kappaHi = max(0, max_i (𝒜_hi v)_i/v_i)`, `rowAhi`, `rowMaxHi`, `wminLo = min(1, min_i w_i(1−δ_w))`,
  `Config.gUp λ = (T*/300)² · tableLookup` (RL1 finding 10: gUp IS the table lookup).
  `Config.ofCert a0 T*` uses the certified table `certTable = [(λ_k, 90000·m0_k)]` and
  `g_lo(∞) = 90000·m0_∞` from `TargetLogic/NMTable.lean`, GENERATED from the certificate JSONs
  (exact decimal `m0_over_T2`, sha256 of each file in its header; `scripts/gen_nm_table.py`).
* **WELL-FORMED `wf`** (§0.6a): A, λ, m, ω_max > 0; E_F > 0 if present; imag, shift ≥ 0;
  n_sheets ≥ 1; δ's and q-spread ∈ [0, 0.5]; ε_sys ≥ 0.10; 0 ≤ ε_ω < 1; matrix: B ≥ 1,
  𝒜 ≥ 0, N_i > 0, v_i > 0, δ_𝒜 ∈ [0, 0.5], `|N_i𝒜_ij − N_j𝒜_ji| ≤ 1 % · max`; DOS: n = n_sheets,
  w_i ∈ (0,1], Σw ∈ [0.99, 1.01], δ_w ∈ [0, 0.5], **Σ_i w_i(1−δ_w) ≤ 1**; and moment consistency
  **A_nom² ≤ λ_nom·m_nom/4·(1 + 10⁻²)** (§0.6a, post-certification revision).
* **G-VALID `gvalid`**: YAML test literally (adiabatic with E_F present, imag ≤ 30 K, H/Li/Be ⇒
  check, shift > 0.10 ⇒ recomputed, μ* not negative, ≥ 3 smearings, ≥ 2 q-meshes,
  ε_sys ≥ 0.10, ε_ω ≥ max(0.02, shift, q-spread), n_sheets > 1 ⇒ mbPresent).
* **Rules `floorR`** in order N1-KAPPA (`kappaHi < a0T*`), N1-ROWMAX (`rowMaxHi < a0T*`), N1-DOS
  (`0 < wminLo`, `A_hi/wminLo < a0T*`), N1-ISO (`n_sheets = 1`, `A_hi < a0T*`, margin or
  kResolved), NM-CERT (`n_sheets = 1`, `m_hi < 4(a0T*)²/λ_hi`), NM-G (`n_sheets = 1`,
  `m_hi < gUp(λ_hi)`, margin or kResolved). `Rule.uses` = the §0.6c table;
  `provOK` = every used input present, `fullDFPT`, one common run id.
* **Outcome** `rejectedProven r` | `screenedOut` | `undecided` | `pass`, with
  `Outcome.assumptions = [ME-MODEL, CERT-H, ENV-SOUND] ++ r.assumptions` and
  `Outcome.inputsUsed = r.uses` (RL1 finding 3). Doc status: `rejectedProven` → REJECTED_PROVEN,
  `screenedOut` → UNDECIDED ("would fire, ranking only"), `undecided` → UNDECIDED, `pass` → PASS.
* `decideCore`: gates, then the first rule with floor ∧ provOK; else `screenedOut` if any floor
  holds, else `pass`. `decideT = applyPolicy (blocked x) (decideCore x)`; a removal with a
  blocking tag becomes UNDECIDED. `unknown` tags block (fail-closed, RL1 finding 7b);
  `informational` tags do not.

Theorems (all generic in `K`):

* **(i) `decide_rejectedProven`**: `decideT = rejectedProven r` ⇒ `wf ∧ gvalid ∧ ¬blocked ∧
  floorR r ∧ ∃ p, p.tier = fullDFPT ∧ ∀ n ∈ r.uses, provOf x n = some p` (every used input
  fullDFPT, same run). `rejectedProven_assumptions`, `iso_records_mIso`, `nmG_records_certNM` (NM-G records [ME-MODEL, CERT-H, ENV-SOUND, M-ISO, CERT-NM]).
* **(ii)** `decide_of_not_wf`, `decide_of_not_gvalid` (⇒ UNDECIDED);
  `decide_not_proven_of_input` (a used input absent / not fullDFPT ⇒ not that rule);
  `decide_not_proven_of_lowTier` (no fullDFPT input ⇒ never REJECTED_PROVEN);
  `decide_not_proven_of_mixed_runs`; `decide_of_blocked` (⇒ UNDECIDED or PASS).
* **(iii)** `decide_downgrade`, `decide_downgrade_rank` (more tags: same or UNDECIDED);
  `applyPolicy_eq`, `policy_rank_mono`, `decideCore_tags`, `blocked_mono`.
* Envelope direction: `Ahi_ge`, `lamHi_ge`, `mHi_ge`, `lamHi_pos` (from `wf`).
* Core-rank lemmas: `core_rank_two`, `core_rank_one`, `core_of_rank_two`, `core_of_rank_one`.

## E. `TargetLogic.Monotone` — worse data never upgrades (RL1 findings 2, 6, 7c)

`Worse x x'`: same validity data, tags, nominal values and provenance; δ_A, δ_λ, δ_m, δ_𝒜, δ_w,
ε_sys, ε_ω not smaller; q-spread arbitrary; matrix / DOS data may be dropped.
**`decide_worse`**: if `x` passes both gates and `cfg.gUp` is antitone on `[0,∞)`, then
`rank(decideT x') ≤ rank(decideT x)`. Supporting: `env_mono`, `floor_mono` (each rule),
`provOf_worse`, `provOK_mono`, `kResolved_mono`.

## F. `TargetLogic.Sound` — semantic soundness per rule (over ℝ, constants pinned)

Every theorem assumes `CertH`, `cfg.a0 = a0Cert`, `0 < cfg.Tstar`,
`decideT cfg x = rejectedProven r`, and ONLY the ENV-SOUND / model hypotheses of that rule:

| rule | model (assumption) | ENV-SOUND hypothesis | conclusion |
|---|---|---|---|
| `kappa_sound` | `MSpec md.B α`, the data's band/k index (M-DISCRETE), `N_i > 0` | `Amat ≤ calAhi` entrywise | `MStable` ∀T ≥ T*, every PSD Coulomb matrix and cutoff; also ∃ `0 < T' < T*` stable on `[T', ∞)` |
| `rowmax_sound` | same | same | same |
| `dos_sound` | `MSpec dd.n α`, sheets (M-SHEET-ISO, SHEETS-CORRECT) | `w_i(1−δ_w) ≤ wfrac_i`, `A_iso ≤ A_hi` | same |
| `iso_sound` | any `MSpec B α` | every row sum `≤ A_hi` (M-ISO + ENV) | same |
| `nmCert_sound` | single band `TcR1.Spectrum` (M-ISO) | `m ≤ m_hi`, `λ ≤ λ_hi` | `PipelineLogic.Stable` ∀T ≥ T*, every μ* ≥ 0, cutoff |
| `nmG_sound` | single band (M-ISO) | `m ≤ m_hi`, `λ ≤ λ_hi`; plus `CertNMSet` and `cfg` = the certified table (`Config.ofCert`) | same (via `gUp_cert`: the lookup lands on a certificate with `λ_hi ≤ λ_k` and `gUp = m0·T*²`) |

The doc claims Tc ≤ T* only; the strict form (∃ T' < T*) is an extra for the multiband
rules (RL1 finding 6). Helpers: `stable_of_cw`, `exists_T'`.
Non-vacuity (RL1 findings 1, 8): `kappa_env_satisfiable`, `iso_env_satisfiable`,
`nm_env_satisfiable`, `dos_env_satisfiable` (the sheet model `N_i = w_i` with the zero spectrum; uses the new
gate Σ w_i(1−δ_w) ≤ 1) — on every well-formed input the zero spectrum (B ≥ 1) satisfies the
ENV-SOUND hypotheses, so the gate never admits an input whose rejection is only vacuously sound.

## G. `TargetLogic.Regression` — RL1 counterexamples, now negated

Config: `a0Cert`, T* = 300. `ex_proven` (positive control, N1-ISO fires on a doc-conformant input);
`ce_negA_fixed`, `ce_signImag_fixed`, `ce_imag_fixed`, `ce_signShift_fixed`, `ce_shift_fixed`,
`ce_eps0_fixed`, `ce_epsOmega0_fixed` (⇒ UNDECIDED); `ce_mix_fixed` (κ input from ML: no rule gives
REJECTED_PROVEN), `ce_mixRun_fixed` (matrix from another run); `unknown_tag_blocks`.

## H. Export for the Python mirror (`export/decide_cases.json`)

`ExportCases.lean` runs the same `decideT` at `K = ℚ` (exact). 563 cases, config = `Config.ofCert 4.33148145 300` (certified table): a 480-case product
grid (A × m × n_sheets × matrix option × DOS option × tier of A; λ chosen as the smallest of 1, 2.5,
8, 20 that makes the nominal moments consistent), 25 λ/m cases around the NM table, and 58
one-factor cases (each gate item, domain violations, provenance, every tag kind, margin
clause). Each case has the full inputs, `expected` (Lean outcome, doc status, rule,
assumptions, inputs used with provenance), the gates, per-rule `floor` / `prov_ok`, and the
derived numbers (A_hi, λ_hi, m_hi, floor, NM-CERT floor, gUp, κ_hi, rowmax_hi, w_min,lo) as exact
rationals `"p/q"`. `moments_consistent_nominal` flags whether `A² ≤ λm/4` holds exactly for the nominal values
(WELL-FORMED allows a 1 % tolerance); the one-factor cases include violations.

## Modelling choices / deviations (read before relying on `decideT`)

1. Provenance is declared data: Lean proves "a rule whose inputs are not all fullDFPT from one
   run never gives REJECTED_PROVEN"; it cannot check that the declared tier and run are true.
2. No `allow_unverified_rejection` override (stricter than the core `certify_policy`).
3. Tier `converged` of revision 1 is gone; the doc's four tiers are used.
4. ε_sys is checked in both gates (WELL-FORMED table and G-VALID YAML list it); harmless.
5. The 25 % margin clause and the symmetry / Σw / moment tolerances are policy; no proof uses them
   (the Σ w_i(1−δ_w) ≤ 1 gate is what makes N1-DOS non-vacuous, `dos_env_satisfiable`).
7. `NMG.lean` still contains the superseded `NMHyp` / `gLookup` / `nmg_*` lemmas (pre-certificate
   route); `decideT` and `nmG_sound` no longer use them.
6. RANK-AIM is not modelled. The σ-exclusion test (RL1 finding 7a) is still a precomputed Bool.

## Not proved / open

* CERT-H, CERT-NM and the ENV-SOUND / model assumptions (hypotheses by design).
* Perron–Frobenius existence (not needed: κ_hi is a Collatz–Wielandt bound).
* General-measure (non-atomic) versions.
* That the Python mirror computes `decideT` (test it against `export/decide_cases.json`).
