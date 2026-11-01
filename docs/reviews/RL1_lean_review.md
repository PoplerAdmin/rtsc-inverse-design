# RL1 — Adversarial review of `lean/` (library `TargetLogic`)

Reviewer: agent R-L1-review-lean, 2026-10-28.
Target: `lean/TargetLogic/{Matrix,Moments,NMG,Decide}.lean`, `lean/CheckTarget.lean`, `lean/README.md`,
measured against `docs/theory/01_eliashberg_targets.md` §0.2–§0.4, §5.1 and
`docs/reviews/R12_architecture_review.md` findings 1–4, 13.
Scratch Lean (all compile against the project build, no `sorry`, standard axioms only):
`docs/reviews/rl1_scratch/RL1Decide.lean`, `RL1NMG.lean`, `RL1HProp.lean`; compiler/axiom output in
`docs/reviews/rl1_scratch/axioms_scratch.txt`.

Severity: **critical** = the specification lets a proven status through on inputs the theory document
excludes, or a soundness claim is wrong; **major** = a stated guarantee is weaker than the README
says or a doc rule is not faithfully encoded; **minor** = wording, missing side-lemma, hygiene.

---

## A. Summary verdict

* **The proofs are not worthless.** Build and axiom audit reproduce exactly (2203 jobs; 49/49
  TargetLogic theorems on `[propext, Classical.choice, Quot.sound]`; no `sorry`/`axiom`/
  `native_decide`). The main semantic theorems are **non-vacuous**: I proved in Lean that the
  numerical hypotheses `H` and `PropN` are jointly satisfiable (`RL1H.H_and_PropN_consistent`,
  witness tp = 1, a0 = 1/64, η = 1/2), and discharged *every* hypothesis of
  `decide_rejectedProven_sound` and of `n1_dos_stable` on concrete inputs and models
  (`RL1H.fully_instantiated`, `RL1H.n1_dos_inst`). The conclusion `MStable` is a genuine
  restriction (`RL1.ex_unstable`: an unstable one-band model exists for every `tp`). `NMHyp` is
  non-vacuous too — a weak version of it is even a *theorem* from `H` (`RL1N.nmhyp_of_H`), and
  `nmg_rule_stable` is instantiated with a firing rule (`RL1N.nmg_rule_inst`).
* **The multiband definitions are faithful** (𝒜 = `Amat`, row sums, `A_iso`, `w_i`, the
  Perron/Collatz–Wielandt input; §B.1 below).
* **The decision layer `decideT` is not a faithful encoding of §0.4.** It returns
  `rejectedProven` on inputs §0.2 excludes (negative A, sign-flipped imaginary modes or
  anharmonic shifts, ε_sys = ε_ω = 0, mixed provenance). Each is proved in Lean
  (`RL1.ce_*`). For these inputs `decide_rejectedProven_sound` gives no guarantee, because its
  hypothesis `SoundMB` is unsatisfiable or assumes what the gate should have checked. The
  combinatorial theorem `decide_rejectedProven` still holds for them, so it is the weaker of
  the two "soundness" statements and should not be presented as one.
* The semantic soundness theorem moves **all physics** (envelope correctness, M-ISO, sheet
  granularity, n_sheets correctness) into one hypothesis, `SoundMB`. The outcome carries none of
  these assumptions, although the YAML requires `assumptions_recorded`.

---

## B. Findings

### 1. [critical] G-VALID / ENVELOPE have no domain or sign checks; `rejectedProven` is reachable on impossible or sign-flipped data

Evidence (all proved, `rl1_scratch/RL1Decide.lean`, config a0 = 4.33148145, T\* = 300):

* `ce_negA`: `A_nom = −1` (a signed, non-admitted α²F, R12 finding 2), everything else
  doc-conformant ⇒ `decideT = rejectedProven` (N1-ISO, since A_hi < 0 < floor).
  `ce_negA_no_model`: for this input **no** one-band model satisfies `SoundMB`, so
  `decide_rejectedProven_sound` says nothing about it. `decide_rejectedProven` still
  "certifies" it: proof grade ∧ GValid ∧ no blocking tag ∧ FiresCertified.
* `ce_sign`: `maxImagModeK = −500` passes `maxImagModeK ≤ 30`. QE and phonopy print imaginary
  modes as **negative** frequencies, so this is the default representation. Also
  `anharmonicShift = −0.30` (30 % softening) with H present, an anharmonic check done and α²F
  **not** recomputed passes `0.10 < shift → recomputed`. The doc tests `max|imag mode|` and
  "shifts > 10 %" (magnitudes; ε_ω uses `|anharmonic shift|`). Softening raises λ ∝ ω⁻² and
  A ∝ ω⁻¹, so ignoring it is anti-conservative. Result: `rejectedProven`.
* Not separately encoded, same mechanism: `lamNom < 0` or `omegaMax < 0` make the adiabatic
  test `λ_hi ω_max ≤ 0.1 E_F` hold trivially. `kappaHi`, `rowMaxHi` or `wMin`-derived ratios
  < 0 fire N1-KAPPA / N1-ROWMAX / N1-DOS. `mNom < 0` fires NM-G.

Required change. Add to `GValid` (or to a separate `WellFormed` gate that comes first):
`0 ≤ Anom, 0 < lamNom, 0 ≤ mNom, 0 < omegaMax`; `0 ≤ k` for `kappaHi = some k`, `0 ≤ r` for
`rowMaxHi`; state the validity fields as magnitudes (`|maxImagModeK| ≤ 30`,
`|anharmonicShift| > 0.10 → recomputed`), or rename them `…Abs` and add `0 ≤` checks. Then prove a
lemma that `GValid x` makes `SoundMB x sp dos` satisfiable for some model with `B ≥ 1`, i.e. the
gate never admits an input whose rejection can only be vacuously sound. Add the four
counterexamples above as regression `example`s with expected outcome `undecided`.

### 2. [major] The ENVELOPE convention is not enforced: ε_sys = 0.10 and ε_ω ≥ max(0.02, |shift|, q-spread) are free parameters

Evidence. `Envelope.WF` only asks `0 ≤ epsSys`, `0 ≤ epsOmega < 1`. §0.2(b) fixes
`ε_sys = 0.10` and `ε_ω = max(q-spread of ω_log, |anharmonic shift|, 0.02)`.
`RL1.ce_eps0`: A_nom = 1250, δ = 0.03, ε_sys = ε_ω = 0, k-resolved check done ⇒ A_hi = 1287.5 <
1299.44 ⇒ `rejectedProven`. `RL1.ce_epsDoc_not_fired`: with the doc values (0.10, 0.02) the same
data give A_hi ≈ 1445 and no certified floor fires. The Lean spec therefore allows a Python
mirror with a smaller envelope than the theory document requires. Nothing ties ε_ω to
`anharmonicShift` either.

Required change. Add `epsSys = 0.10` (or `0.10 ≤ epsSys`), `0.02 ≤ epsOmega` and
`|anharmonicShift| ≤ epsOmega` (when a check exists) to `GValid`, and (optionally) a field for the
q-spread of ω_log with `qSpreadOmegaLog ≤ epsOmega`. Prove the monotonicity lemmas
"increasing any of δ, ε_sys, ε_ω never turns a non-rejection into a rejection". They are true for
the current formulas but not stated (see finding 7).

### 3. [major] `SoundMB` bundles every physical assumption; the outcome records none of them

Evidence. `decide_rejectedProven_sound` is proved from `H`, `PropN` and `SoundMB x sp dos`.
`SoundMB` contains, for one model `(sp, dos)`:
(a) the DFPT numbers are upper envelopes (κ_hi certificate, row-sum bound, `A_iso ≤ A_hi`);
(b) `wMin ≤` every true DOS fraction;
(c) **M-ISO**: `nSheets = 1 → every true row sum ≤ A_hi`;
(d) implicitly, that the model's band index *is* the right granularity (M-DISCRETE /
M-SHEET-ISO) and that the reported `nSheets` is correct.
(c) is the assumption §5.1 says can fail ("intra-sheet anisotropy … can still make κ > A_iso").
The 25 % margin clause, which the doc uses as a hedge against M-ISO, plays **no** role in the
proof: `FiresISO`'s margin conjunct is discarded (`⟨h1, hlt, -⟩` in `decide_rejectedProven_sound`).
`Outcome` has no assumption list, so a `rejectedProven` from N1-ISO (resting on M-ISO) cannot be
told apart from one from N1-KAPPA (resting only on envelope soundness and M-DISCRETE). The YAML
asks for `assumptions_recorded`, and the core's `annotate_model` adds
`M-ISO-SINGLEBAND`/`DISCRETISATION` (R12 finding 1). The README sentence "REJECTED_PROVEN is a
stability theorem" is true only relative to (a)–(d).

Required change. (i) Return the fired rule and its recorded assumptions with the outcome (e.g.
`rejectedProven (rule : CertRule)` with `CertRule.assumptions`). Prove one soundness theorem per
rule with only the `SoundMB` fields that rule uses: N1-KAPPA needs `kappa`, N1-ISO needs `iso`, and
so on. (ii) Split `SoundMB` into `EnvelopeSound` (numerical) and `ModelAssumption`
(M-ISO / M-SHEET-ISO / M-DISCRETE). (iii) In the README, state that N1-ISO rejections are proven
*assuming M-ISO*, and that the margin clause is policy, not part of any proof.

### 4. [major] One `tier` for all numbers: mixed provenance gives `rejectedProven`

Evidence. `Inputs.tier` is a single field covering `env` **and** `mb`. `RL1.ce_mix`: a full-DFPT
isotropic record (`nSheets = 2`, band-resolved flag true) plus a `kappaHi = 1000` of any origin
(ML regression, coarse run, hand entry) ⇒ `rejectedProven` via N1-KAPPA.
`not_rejected_of_not_proofGrade` is true but only speaks about the single declared label.
Rejections from coarse or ML data are therefore reachable whenever the multiband inputs come from
a different run than the isotropic α²F. That is the natural workflow: the README itself notes that
DOS fractions come "cheaply" from a band-structure run. R12 finding 4 asks for exactly this binding.

Required change. Give each rule input its own provenance (`kappaHi : Option (ℝ × Tier)`, etc., or a
`Prov` record per field). Make `FiresCertified` count a rule only when *all* its inputs are proof
grade, and re-prove `not_rejected_of_not_proofGrade` as "if any input used by the fired rule is not
proof grade, then no rejection".

### 5. [major] N1-DOS and N1-KAPPA/ROWMAX need envelopes the theory document never defines; the Lean exposes the needed direction but nothing supplies it

Evidence. `SoundMB.dosw` needs `wMin ≤ wfrac_i` for every *true* sheet fraction, i.e. a **lower**
envelope on w_min. §0.4 lists `w_i = N_i/N` with no envelope, and ENVELOPE (§0.2b) covers only A,
λ, m. A smearing-dependent DOS at E_F can err in the anti-conservative direction. The same holds
for κ_hi and A_i,hi: the YAML says "upper envelope" but ENVELOPE gives no formula (A-type quantities
should scale like A_hi: (1+δ)(1+ε_sys)/(1−ε_ω)).

Required change (doc 01 + Lean). Add ENVELOPE formulas for κ_hi, A_i,hi and w_min,lo (e.g.
w_min,lo = min_i w_i·(1−δ_w) from the smearing ladder of the DOS). Add `dW` to `MBInputs` with
`FiresDOS` using `w·(1−dW)`, and put the corresponding `δ ≤ 0.5`, `δ ≥ 0` checks into G-VALID.

### 6. [minor] "Stable at every T ≥ T\*" gives Tc ≤ T\*, not Tc < T\*

Evidence. README: "Stable at every T ≥ T\*" means `Tc < T*` with `Tc := sup{T : unstable}`. That
is false as stated: stability on `[T*, ∞)` gives `sup ≤ T*`. Unstable temperatures can accumulate
at T\* from below. For the certified branches a strict statement *is* available
(`n1_rowmax_stable` also returns `bound/a0 < T*`, and `tc_le_rowmax` gives `T ≤ bound/a0`).
`decide_rejectedProven_sound` drops it, though. For NM-G only `≤` is available. Practically
irrelevant (Tc = 300 K exactly), but the README overstates.

Required change. Either write "Tc ≤ T\*", or add to `decide_rejectedProven_sound` the conclusion
`∃ T' < T*, ∀ T ≥ T', MStable …` (take T' = bound/a0, available in every certified branch).

### 7. [minor] `decide_downgrade` covers list inclusion on `tags` only; the tag semantics and other "downgrade" directions are not formalised

Evidence. `decide_downgrade` is correct and covers every tag *addition* (`ts ⊆ ts'` as lists,
including duplicates and reordering; `decideCore` is literally independent of `tags`). But:
(a) `sigmaExcluded r c` carries a pre-computed Bool. The core's `_sigma_exclusion_directions` logic
("excluded smearing's λ or A above the envelope-widened range", R12 finding 1) that decides `r`
is not modelled, so the Lean guarantees nothing about the most delicate tag (the Pb case).
(b) The core's `blocking_tags` matches strings by prefix and regex. Unknown spellings map to `other`, which
**does not block** (fail-open). A typo such as `A-FILE-INCOMPLET` silently loses its veto.
(c) There is no monotonicity theorem for non-tag degradations: larger δ/ε, smaller w_min, larger
κ_hi, `kResolvedCheckDone := false`, fewer smearings or q-meshes. By inspection all of these are
direction-safe in the current formulas, but nothing proves it, and finding 2's fix will touch them.

Required change. (a) Model the σ-exclusion test as a function of the excluded and admitted
(A, λ) values and the envelope, and prove the core's intent. (b) Make unknown tags block rejection
(fail-closed), or have the Python mapper reject unknown tags with an error. (c) Prove
`decideT` rank-monotone under each of the degradations listed in (c).

### 8. [minor] The empty model satisfies `SoundMB` trivially; nothing ties `B` to the material

Evidence. `RL1.ex_B0_sound`, `RL1.ex_B0_stable`: for B = 0, `SoundMB xNeg sp dos` holds (all
`∀ i` are vacuous), and every B = 0 model is `MStable`. Harmless, since the theorem is
universally quantified over models. But a reader may take "the true material" to be pinned down,
and it is not: no hypothesis says `B ≥ 1`, `B = nSheets` (sheet model) or that the model is the
DFPT discretisation.

Required change. Add `0 < B` (or `Nonempty (Fin B)`) to `SoundMB`, and document which model `(sp,
dos)` is meant for each rule (sheet model for N1-DOS, k/band discretisation for N1-KAPPA/ROWMAX,
one band for N1-ISO).

### 9. [minor] `tp` and `a0` are not pinned; the H witness shows the constants matter

Evidence. Every theorem holds for any real `tp` (MultiBand: "tp stands for 2π"). `Config` fixes
neither `tp` nor `a0 = 4.33148145`. `H 1 (1/64) (1/2)` is provable (`RL1H.H_witness`), so a
mis-instantiation of `H` (wrong `tp`, wrong units) would yield "proven" stability for a different
kernel. The core certificate is for one specific `(tp, a0, η)`.

Required change. Define `tpPhys := 2 * Real.pi` and `a0Cert := 4.33148145` in Lean. Make
`decide_rejectedProven_sound` take `H tpPhys cfg.a0 η` with `cfg.a0 = a0Cert` (or `cfg.a0 ≥ a0Cert`,
which is direction-safe because a larger a0 lowers no floor). Record in README which η the
certificate gives (core: 1e-9).

### 10. [minor] NM-G: `cfg.gUp` is an arbitrary function; the table, the s2 scaling and `gLookup` are not connected to `decideT`

Evidence. `gLookup_le` is proved, but `decideT` uses `cfg.gUp`, and the soundness theorem simply
assumes `∀ lh l, 0 ≤ l ≤ lh → gUp lh ≤ g l`. `λ_hi` "unknown → g(∞)" cannot be represented
(`lamHi` is always computed). The scaling `g·s2` with s2 = (T\*/300)² is not tied to `Tstar`.

Required change. Set `cfg.gUp := gLookup g up 100 gInf` scaled by `(Tstar/300)^2` in a concrete
config, and derive `hgUp` from `gLookup_le` inside the theorem. The hypotheses should then be the
numerical ones only (`NMHyp`, antitone, `g∞ ≤ g`).

### 11. [minor] Smaller faithfulness points

* `FiresKappa`/`FiresRowMax` do not check `bandResolvedAvailable`. `FiresDOS` checks
  `dosFractionsKnown` but no rule checks that the supplied `Option` and the flag agree. Make the
  flags derived (`bandResolvedAvailable := kappaHi.isSome ∨ rowMaxHi.isSome`), not independent.
* E_F "or half the bandwidth if smaller" and "ladder on the finest k-mesh" are inputs, not
  modelled. That is acceptable, but the README should list them under "not checked".
* SCREENED_OUT carries no criterion, e_c or source id (R12 finding 3a).

### 12. [info] NM: a weak version of NMHyp is a theorem; the numerical part matters only at λ ≳ 1

`RL1N.nmhyp_of_H`: from `H` (a0, η > 0) and T\* > 0,
`NMHyp κ gT T*` holds with `gT(λ) = 4(a0T*)²/max(λ, ε)` (TcR1 `stable_of_ge` + L-MOM2
`A² ≤ λm/4`; the zero spectrum is handled by `RL1N.stable_of_zero`). At T\* = 300 K,
4(a0T\*)² = 6.754×10⁶ K². Compared with the §0.3 table: λ = 0.8: 8.44 vs 9.453; λ = 1: 6.75 vs
6.851; λ = 2: 3.38 vs 4.015; λ = 3: 2.25 vs 3.469 (×10⁶ K²). So around λ ≈ 1 the [NUM] content of
NM-G adds under 2 %. This is the certified "λ → ∞ end" of §3.2 in another form. Suggest: a check
`m_hi < 4(a0T*)²/λ_hi` could be a fifth certified floor (`rejectedProven`, assumptions M-ISO only).
Caveat: with exact values it is implied by N1-ISO (A² ≤ λm/4). It adds something only when the A
envelope is looser than the λ and m envelopes.

---

## C. Faithfulness checks that passed

* **Multiband definitions.** `Amat i j = (Σ_s W_sij)/N_i = ∫α²F_ij` with `W_ij = N_i α²F_ij`
  symmetric. This is 𝒜 of §5.1 (S_ij = N_i 𝒜_ij symmetric). `rowA` = YAML `A_i = Σ_j calA_ij`.
  `wfrac` = `N_i/N`. `Aiso` = A of `α²F_iso = Σ_ij (N_i/N) α²F_ij` (`Aiso_eq`). The YAML's
  `PerronRoot(sqrt(N_i/N_j) calA_ij)` is a diagonal similarity of 𝒜 and has the same spectrum. A
  Collatz–Wielandt certificate `𝒜v ≤ κ_hi v, v > 0` bounds it from above, so `SoundMB.kappa` is
  the right (and checkable) form.
* **N1-DOS proof.** The row-sum route `w_i A_i ≤ Σ_k w_k A_k` is correct and needs no symmetry.
  The doc's S-route is also correct; same constant.
* **Direction of every floor comparison.** All rules compare *upper* envelopes (`A_hi`, `κ_hi`,
  `rowMax_hi`, `m_hi`, `λ_hi` inside a decreasing g) with `<`. `Ahi_ge`/`lamHi_ge`/`mHi_ge` are
  correct for well-formed envelopes and nonnegative nominals (see finding 1 for the missing
  nonnegativity).
* **G-VALID thresholds.** These match the YAML literally (30 K, 0.1, 0.10, 3, 2, 0.5). The
  `n_sheets > 1` clause matches; N1-ISO/NM-G require `nSheets = 1`.
* **Provenance and policy.** `decide_rejectedProven`, `decide_of_not_gvalid`,
  `decide_of_not_proofGrade`, `decide_of_blocked`, `decide_downgrade(_rank)` say what the README
  says. The policy is strictly safer than the core's (no `allow_unverified_rejection`).
* **L-MOM1..3** carry the explicit hypotheses R12 findings 2/13 asked for (`Ω_j ≤ ω_max`,
  `λ > 0`) and are not used by any rule.

## D. Non-vacuity evidence (Lean, `rl1_scratch/`)

| theorem under review | instance | file / name |
|---|---|---|
| `H ∧ PropN` (all of A, B of README) | tp = 1, a0 = 1/64, η = 1/2 (uniform bound \|F⁺(ge)\| ≤ 10π\|x\|² via arctan telescoping + MVT + MultiBand.Odd row sums) | `RL1HProp.lean` `H_witness`, `PropN_witness`, `H_and_PropN_consistent` |
| `decide_rejectedProven_sound` | **all hypotheses discharged**: cfg (a0 = 1/64, T\* = 64000), doc-conformant input, one-band model W = 500 | `RL1HProp.lean` `fully_instantiated` |
| same, with the doc constants a0 = 4.33148145, T\* = 300 | all but `H`, `PropN` discharged | `RL1Decide.lean` `ex_sound_inst` |
| `n1_dos_stable` | all discharged, two bands, w_min = 1/2, A_iso = 200 | `RL1HProp.lean` `n1_dos_inst` |
| `nmg_rule_stable` | all discharged except `H` (`NMHyp` derived from `H`), g = 4(a0T\*)²/max(λ,1), rule fires | `RL1NMG.lean` `nmg_rule_inst` |
| conclusion `MStable` non-trivial | unstable B = 1 model at T = 300 for every tp | `RL1Decide.lean` `ex_unstable` |

## E. Build and audit reproduction

Re-ran in WSL (`nice 19`, 2 CPUs, existing build dir `~/rtsc_targeted_work/L1-lean/build`):
`lean/scripts/build.sh` → `Build completed successfully (2203 jobs).` (only warnings: Mathlib
module-system import notes in `Moments.lean`). `lean/scripts/check.sh` →
Check 34, CheckPipeline 75 + 7 axiom-free, CheckMultiBand 67, CheckTarget 49, no disallowed
axiom, no error/sorry lines, source scan clean. All 49 theorems declared in
`TargetLogic/*.lean` appear in `CheckTarget.lean` (checked by name diff). **Claims confirmed.**

## F. What I could not verify

* That `H tpPhys 4.33148145 η` itself holds. My witness proves only consistency of the
  hypothesis set, at a much smaller a0. The certified value rests on the core's interval-arithmetic
  certificate, which I did not re-run.
* `PropN` at the certified a0 (wMB §4.2), same.
* Python mirror of `decideT`: none exists yet to test against.
