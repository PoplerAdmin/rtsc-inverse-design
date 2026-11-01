# R30 — Adversarial review of Campaign 1 (`docs/design/30_campaign1.md`, `runs/campaign1/`)

Reviewer R30, 2026-10-30. Read-only on compute. The only computation was small offline Python
(nice 19, 2 threads) over the repo export (`runs/campaign1/export/*.jsonl`, `tables/*`) and the
cached Alex-MP reference queries (`~/rtt_scratch/front_cache/alexmp_ref/ref_<chemsys>.json`). The
analysis scripts are not in the repo; every number below can be re-derived from the export with
the joins named in each finding. Labels: **(meas.)** computed by me from the export;
**(cited)** with an id; **UNVERIFIED** not checked.

Scope note: another agent is fixing the coarse-screen cutoff/σ issue (Mg₂IrH₆ coarse
A_mid 341 K vs ≈ 990 K in doc 15). I do not repeat that diagnosis. Finding F2 lists the
**other** campaign decisions that depend on it.

---

## 0. Verdict in five lines

1. **Every removal in the campaign is numerically safe.** No candidate was lost that plausibly had
   A ≳ 300 K. Two removal rules, though, were applied **outside their calibrated or stated domain**:
   - the S3 H7-ml window (34 G4 non-hydride removals);
   - the S1 `n_atoms ≤ 20` rule (19 lane-D removals).

   The one S5 removal (LaH₂) rests on a σ ladder that admission truncated to σ ≥ 0.06 Ry. That is
   exactly the regime where the Mg₂IrH₆ control shows a ×3 under-estimate (F1, F2, F5).
2. **"Meaningful MatterGen output" is shown only for composition targeting.** It is **not** shown
   for the target property:
   - The two best "generated" DFT results (NaIrH₃, BaPtH₃) are **re-generated Alexandria entries**.
   - S2 dedupe silently folded 25 lane-D known compounds *into* G1 candidates.
   - G1 does not beat lane D on metallicity (64 % vs 67 % at S4) or on A_coarse (F3, F4).
3. **The novelty claims are half right.**
   - BaRhH₆ (and BaIrH₆) are new compositions vs Alex-MP. No literature hit was found.
   - **LaRhH₂ is a new polymorph of a composition already in Alexandria** (4-atom R-3m entry
     `msft-2023063021282812988979740863-GGA`, E_hull 0.114 eV/atom). The generated cell is P4/mmm
     and does not match it.
   - BaRhH₆ is mislabelled "C2/m-like (sg 8)": it is R-3m at 0.1 Å, and sg 8 is Cm (F6).
4. **NaIrH₃ was the wrong S6 choice for discovery.** λ ≈ 0.15–0.25, S5 status ERROR, and S6
   admission is at risk at small σ. It is defensible only as an end-to-end chain test, which then
   ends in `undecided` by construction. The S6 slot should go to a **Mg₂IrH₆ S6-lite**, the
   positive control that calibrates every hydride number (F7).
5. **No generation strategy on offer plausibly reaches A ≥ 1300 K at ≤ 1 GPa.** The only family
   within reach of the floor is the A₂MH₆ (K₂PtCl₆-type) [MH₆] complex family: Mg₂IrH₆ is at
   ≈ 76 % at small σ. Campaign 2 should map that family's **electron-count / A-site-alloy**
   neighbourhood systematically, with a deterministic prototype lane plus quaternary-conditioned
   MatterGen. It should report honestly that the expected outcome is ≤ ~1000 K (§3).

---

## 1. Findings, ranked by severity

Severity: **S1** changes a conclusion or violates the contract; **S2** biases the reported
comparison or wastes ≥ 1 day of compute; **S3** documentation or accuracy.

### F1 (S1) — Two heuristic removal rules were applied outside their calibration domain (contract R2/R7)

**(a) S3 H7-ml window on non-hydrides.**
- The threshold `ehull_self_max_eV = 0.60` (`src/rtt/configs/tiers.toml:10`) is built from three
  pieces, all of them hydride-specific:
  - MatterSim-vs-DFT hull errors on **160 H-entries of H-Ir-Mg / H-Li-Mg / H-Ca-Pd**;
  - a "Mg offset";
  - an **H ZPE** allowance.
- 34 of the 38 S3 removals are G4 **non-hydride** B–C–Li / B–C–Mg / B–C–N cells, plus one lane-D
  B–C cell (meas., `s3_table.jsonl`, `S3 = SCREENED_OUT`).
- No calibration of MatterSim hull error exists for borocarbides. In one case the stage's own
  second opinion disagrees:
  - **c-0b81a0e54dbb** (B–C–N, 6 atoms): e_hull_ml_self 0.632 > 0.60, but e_hull_ml_dftref 0.524
    < 0.60.
- Three Li–B–C cells sit within 0.02 eV/atom of the cut:
  - c-32bf9da2e451: 0.603;
  - c-45edc6dce0e0: 0.607;
  - c-5b81a26468d4: 0.616.
- Physics impact is small, because B–C cells > 0.6 eV/atom above hull are poor candidates. The
  contract issue is real, though: the outcome's `source_id` claims a calibration that does not
  cover these candidates.
- **Fix:** for compositions outside the calibration set, emit `DEFERRED(uncalibrated_window)`
  instead of `SCREENED_OUT`, or calibrate separately (≈ 50 B–C–N–Li–Mg Alex-MP entries,
  ~10 GPU-min).
- Doc 12 line 197 still states the window as "0.35 + 0.07". It should cite 0.60 / doc 21.

**(b) S1 `n_atoms_max = 20` on lane D.**
- The rule's stated source is "doc10#1 (MatterGen training ≤ 20 atoms)" (`tiers.toml:5`). That is
  a generator-domain argument, and it has no force for database entries.
- 19 lane-D entries were `SCREENED_OUT(n_atoms_primitive)` for this reason, among them:
  - known Na₃IrH₆ (mp-707836, 40 atoms);
  - Y₁₀H₂₁, LaH₂₀ (21 atoms);
  - several AlH₃/MgH₂/YH₃ polymorphs;
  - LiB₁₃C₂, B₄C, Mg(BC)₂.
- Most are closed-shell insulators, so impact ≈ 0.
- **Fix:** lane D should get `DEFERRED(cost_class)` with a cost source, not a MatterGen-domain
  screen.

**(c) Not a finding — checked and sound: the H1 removals.**
- Eight S4 `EXCLUDED_ASSUMPTION(H1)` outcomes, all with PBE gaps 1.68–3.69 eV:
  - SrPtH₆, BaPtH₂, Ba₂Pt₂H₈, NaAuH₂, BaPtH₆, and AlH₃ ×3.
- After the B2 fix, `decide_s4` (`src/rtt/dft/electronic.py:611–638`) removes only when:
  - the tetrahedron N(E_F) is at or below `nef_min`;
  - there are no sheets;
  - the gap exceeds `gap_max + marginal band`.

  Anything in the marginal band is `DEFERRED`. The Gaussian/tetrahedron definition question does
  not reach these cases, because a 1.7 eV PBE gap is far outside any smearing artefact.
- Independent check (meas.): the valence-electron parity of each cell separates the S4 outcomes
  perfectly.
  - All 4 G1 H1-exclusions have **even** valence-electron counts: 18, 14, 32, 14.
  - Every passing odd-count cell is a necessarily-metallic NM band structure: NaIrH₃ 13, BaPtH₃ 15,
    BaRhH₆ 17, La₂Rh₂H₇ 31, Ba₂IrH₆ 19, Mg₂IrH₆ 19, BaRh₂H₇ 27.
- The removals are sound. The S2/S3 stages could have predicted them for free (see F4, plan §3.2).

**(d) S2 DUPLICATE merges** are sound as equivalences: tight matcher, `scale=False`, `stol=0.1`,
plus a bonding guard (`src/rtt/structure/dedupe.py:30`). The choice of *representative* is wrong,
though; see F3.

### F2 (S1) — Campaign decisions that depend on the coarse-σ/cutoff issue (beyond Mg₂IrH₆ itself)

The S5 profile uses ecut 60 Ry and a σ ladder 0.01–0.10 Ry, and builds A_mid over σ ≥ 0.03
(`src/rtt/dft/profiles.py:178–182`). Doc 15 §1 says SG15 H needs ≥ 80 Ry.

Two meas. comparisons from `s5_table.jsonl` against doc 15 locate the problem:
- At matched σ = 0.03 Ry, campaign q 4³ Mg₂IrH₆ gives **A = 518 K**, against doc 15's 520 K
  (q 3³, 80 Ry). At equal σ the two agree.
- The campaign's "341 K" is the median over σ ∈ [0.03, 0.10], where A falls (N(E_F) halves
  between σ 0.005 and 0.05, doc 15 §3).

So the dominant effect is the **σ statistic**, not the cutoff. The other agent owns that
diagnosis. Its consequences elsewhere:

1. **LaH₂ `SCREENED_OUT(coarse-A envelope)`** (c-37743845859b), the only S5 removal.
   - Mode-positive admission failed for σ ≤ 0.05 on **both** meshes. Only σ 0.06–0.10 were
     admitted, so δ_σ = 0.096 is blind to the small-σ regime.
   - That is exactly the Mg₂IrH₆ failure mode, and doc 13 §0 warns about truncated ladders.
   - Doc 13's calibration set (Al, Pb, Nb, H3S; doc 13 §6) contains no 0-GPa hydride.
   - **Numerically the removal still stands.** Even the Mg₂IrH₆ ratio A(σ 0.005)/A_mid ≈ 2.9
     applied to (1+e_c)·A_hi = 65 K gives ≈ 190 K ≪ 1299 K.
   - Re-tag it `PROVISIONAL_MARGIN + TRUNCATED_LADDER`. Do not cite it as calibrated.
2. **ScH₂ "rule would remove (132 K)"**: same situation (σ ≥ 0.07 admitted); same conclusion.
3. **The headline "no candidate carries more than ~110 K"** (doc 30 lines 28–31, 370–375).
   - The bound is a σ-averaged statistic. A small-σ re-evaluation could raise A up to ≈ 3×, which
     would put the best non-control at ≲ 0.25 of the floor.
   - The conclusion "far below 1299 K" survives. The quoted numbers and the "8 %" do not; they
     should be restated at the corrected σ statistic.
   - NaIrH₃ is the exception that goes the other way: its A *rises* with σ (101 → 122 K), so its
     small-σ value is lower.
4. **The §6 ranking and the S6 choice** use A_mid. NaIrH₃ (108 K), LaRhH₂ (96 K q 3³), BaPtH₃
   (85 K) and ScH₂ (74 K) lie within the σ-statistic uncertainty of one another, so the order
   among them is not meaningful.
5. **Doc 30 misstates the Mg₂IrH₆ control margin** (lines 237–238 and 285–286: "PASS … only
   through the e_c margin", "close to removing").
   - `s5_table.jsonl` gives A_nom 351.4 K > 346.5 K, so the control passes **without** e_c.
   - The removal test value is (1+e_c)·A_hi = **1639 K** against 346.5 K, i.e. 4.7× clear.
   - The screen did not nearly remove the control. The real problem is that the control's *A*
     is wrong by ×2.9, which the screen's own δ_σ = 1.575 happened to cover.
   - R8 recall is therefore not at risk from this run. The control's value as a calibration of
     *A* is.
6. **Not dependent:** S1, S2, S3, H1/H2 (gaps ≥ 1.7 eV; moments 0) and H7-heur. S4 also uses 60 Ry,
   but doc 15's cutoff test (8.7 meV/atom, −2.7 kbar at 60 Ry) is far inside the 0.56 eV/atom
   H7-heur window.

### F3 (S1 for the MatterGen claim) — The funnel credits the generator with rediscovered training data

Measured from `duplicates.jsonl` × `outcomes.jsonl` (S2 `value.novel`):

- **25 G1 candidates are canonical duplicates of lane-D Alex-MP entries.** S2 processed G1 first,
  so the D copy became `DUPLICATE` and the known compound continued under lane G1.
- They include every G1 cell that reached S5 apart from LaRhH₂:
  - **NaIrH₃** (c-0ae3f133a308 ≙ D c-cc63c9c69fce, Alexandria `msft-2023070218340169597120730953`);
  - **BaPtH₃** (c-d985155f945a ≙ D c-58dcb188f839).
- They also include G1 S4 entries: SrPtH₆, BaPtH₂, NaAuH₂, SrPtH₄ (Sr₂Pt₂H₈), LaRh₂H₇ and LaRhH₄.
- **8 of the 15 G1 S3 PASSes are not novel.**
- Of the 7 novel ones:
  - 3 are S4 metals: LaRhH₂, BaRhH₆, La₂Rh₂H₇;
  - 1 is an insulator: BaPtH₄;
  - 1 is an ERROR: BaIrH₆;
  - 2 were not run: Ba(RhH₄)₂, Ba₂IrH₇.

Consequences:
- Lane D lost 25 of its best entries to G1. The D-vs-G comparison of doc 30 §3–§4 is confounded,
  and the D lane is understated.
- The S2 representative was chosen by insertion order, not by cost or symmetry. Example:
  NaPtH₃ kept the 20-atom sg 35 cell (c-20a77d23a970, C3) and dropped the 5-atom Pm-3m G1 and D
  cells (C1).
- **Fix:**
  - Run D before G (or prefer the D/DFT-relaxed copy as canonical).
  - Record generated duplicates of known entries as `rediscovered = true` on the D row.
  - Choose the canonical copy by (lane D first, then fewest atoms, then highest symmetry).
  - Re-cut the funnel.

### F4 (S2) — "Meaningful MatterGen output" per funnel stage (meas.)

| Stage | G1 (chem.sys + E_hull) | G5 (base) | G4 (non-H) | D (Alex-MP) |
|---|---|---|---|---|
| S2 valid | 384/384 (0 INVALID) | 41/41 after S1 | 80/80 | 487/487 |
| unique (not DUP) | 369/384 (96 %) | 41/41 | 80/80 | 456/487 (25 lost to G1, F3) |
| novel structure vs Alex-MP | 336/369 (91 %) | 31/41 | 78/80 | 0 |
| novel **composition** vs Alex-MP | 247/369 (67 %) | – | – | – |
| E_hull,ML ≤ 0.1 (self / DFT-ref hull) | 185 / 122 of 369 | 28 / 24 of 41 | 5 / 1 of 80 | 367 / 320 of 456 |
| **SUN** (novel ∧ E_hull ≤ 0.1), ML-self / DFT-ref | **155 (42 %) / 101 (27 %)** | 21 (51 %) / 16 (39 %) | 5 / 1 | – |
| SUN with H | 154 | 0 | 0 | – |
| SUN with H_f ≥ 0.5 | 42 | 0 | 0 | (known: 261 with H_f ≥ 0.5) |
| SUN with odd valence-e count (NM metal guaranteed) | 44 of 155 | n/a | n/a | n/a |
| S4 metallic / decided | 7/11 (64 %); novel 3/4 | – | – | 8/12 (67 %) |
| best A_coarse | 108 K (NaIrH₃, **not novel**) | – | – | 74 K (ScH₂) |

Reading:
- **Demonstrated:** composition control.
  - 98 % of G1 cells are in the requested system.
  - G1 produced 154 novel, near-hull, H-containing structures in 24 systems, 110 of them in new
    compositions. G5 produced none.
- **Not demonstrated:** enrichment for the target property.
  - G1's S4 metal fraction equals lane D's. n is tiny, and no comparison of a proxy for A is
    possible: 1 novel S5 result, and that one informal.
  - The G1 E_hull,ML distribution is *worse* than unconditioned G5 (SUN 42 % vs 51 %). This is
    expected: G1 asks for E_hull U(0, 0.3).
- **SUN counts are upper bounds.**
  - Five named G1 SUN cells sit 60–90 meV/atom *below* the DFT-reference hull: Ba₂IrH₇ −0.086,
    Ba(RhH₄)₂ −0.090, Na₃IrH₆ (P-3m1) −0.086, BaIrH₆ −0.069, BaRhH₆ −0.062.
  - A generated structure that far below a DFT hull of 30–100 entries is more likely MatterSim
    over-binding (the window calibration's p95 is 0.092 eV/atom) than a discovery.
  - S3's score rewards low E_hull and therefore preferentially promotes these ML-energy errors.
    Mark them `ML_BELOW_DFT_HULL` and do not use "stable" in any statement until S7a.
- **G4** is non-informative: 5/80 within 0.1 eV/atom, and the window issue of F1(a) applies.
  **G5** loses 36 % of its cells to `no_pseudo`, so it is no longer an unbiased control.
- **All 992 S2 entries were valid.** For raw MatterGen output 0 INVALID is plausible but
  unusually clean. Check that the validity test (min distance, density) actually ran on lane G
  (outcome `origin = rtt.structure.stage`; I did not audit the thresholds).

### F5 (S2) — The S6 launch on NaIrH₃ bypassed the stage gate and is low-value for discovery

- **Gate bypass.** `S6.requires = ("S5",)` (`src/rtt/engine/stages.py:321`), but `rtt step`
  calls `run_job` directly and does not check `requires` (`src/rtt/cli.py:51–68`). NaIrH₃'s
  current S5 outcome is **ERROR**. Its S6 ran anyway (tier_runs 1087/1088), and doc 30 does not
  list this as a recorded deviation (§3 lists only `CAMPAIGN_GENERATED_QUOTA`).
  - **Fix:** `cmd_step` should refuse unless `st.ready(rec)`, or require `--override <reason>`,
    which writes a tag.
- **Physics.**
  - λ = 0.153–0.247 on the admitted q 4³ σ (σ 0.02–0.10), ω_log 770–935 K, A 101–122 K rising
    with σ. That is a weak-coupling metal; Allen–Dynes would give a few K at most (est.).
  - Doc 30 §5.6 already expects `undecided`: the harmonic-only H compound fails G-VALID 2. The
    ≈ 124 core-h (36 h on 4 ranks) buys one chain test whose decision layer is pre-determined.
- **Admission risk.**
  - The S6 ladder starts at σ = 0.005 Ry (`profiles.py:185`).
  - On q 3³, NaIrH₃'s negative per-mode α²F fraction grew monotonically toward small σ
    (3·10⁻⁶ at 0.10 → 2.7·10⁻³ at 0.01 Ry), and every σ failed.
  - If q 6³ behaves like q 3³, the finest mesh may admit < 3 σ. The run would then end
    `ERROR`/`undecided` on G-VALID 4, not on physics. This is not verified; it is a plausible
    36-hour failure.
- **Recommendation:**
  - Let the running job finish (sunk cost, and it is the first S6 → certify → decide chain on
    real QE output).
  - Do not launch further S6 on λ < 0.5 candidates.
  - The next S6-class slot goes to the Mg₂IrH₆ control (F7).

### F6 (S2) — Named-candidate statements: corrections

| Doc 30 claim (line) | Problem | Corrected statement |
|---|---|---|
| LaRhH₂ "novel" (176, 316) | composition LaRhH₂ exists in Alexandria (4-atom R-3m, `msft-2023063021282812988979740863-GGA`, E_hull 0.114); StructureMatcher (default ltol/stol, scale = True) does not match it to either the raw or the relaxed G1 cell | **new P4/mmm polymorph of a known hypothetical composition.** E_hull,ML 0.029, DFT-ref 0.080. Its λ ≈ 1.2 comes from a soft mode at ω_log 41 K (near-unstable). Not a high-A structure |
| BaRhH₆ "C2/m-like (sg 8)" (177, 317) | sg 8 is **Cm**; spglib at 0.1 Å gives **R-3m** (raw and MatterSim-relaxed) | R-3m BaRhH₆ with isolated RhH₆ octahedra (min Rh–H 1.63 Å, min H–H 1.89 Å); new composition vs Alex-MP; no literature hit (web search 2026-10-30; ICSD/OQMD/GNoME **not checked**) |
| BaRhH₆ / BaIrH₆ H2 "non-magnetic" | formally [M^IV H₆]²⁻, low-spin d⁵ with 17 valence e per cell. The S4 spin check starts only from a ferromagnetic moment of 0.5 with scalar-relativistic SG15 (no SOC, no U). Ir⁴⁺/Rh⁴⁺ d⁵ octahedra are the textbook J_eff = ½ / Mott case in oxides; whether that carries over to hydrides is **UNVERIFIED** | H2 "checked" is too strong. Flag `H2_FM_START_ONLY`; add AFM and +SOC scf before any S5 on these two |
| NaIrH₃ "#1 named candidate" (313) | ranked first by "furthest stage", which reflects a manual queue choice (F5); not novel; λ ≈ 0.2 | rank by A_coarse at a common σ statistic, with novelty shown. Then Ba₂IrH₆ (D, Fm-3m, 19 e, isoelectronic with Mg₂IrH₆, 1 sheet, H_DOS 0.60) is the most informative open item, not NaIrH₃ |
| "Mg₂IrH₆ passes S5 by < 2 %" (237, 285–286) | wrong; see F2 item 5 | passes with A_nom > floor and 4.7× test margin; the problem is A accuracy, not recall |
| "insulators passing H1 … 5" (B2, 384) and §5.1 | §5.1 lists 8 H1 exclusions; B2 says 5 passed before the fix | consistent after reading, but state which 3 were caught by the tetrahedron path from the start |
| "Generated 'novel' 336/369" (137) | structure-level novelty with a tight, unrelaxed matcher (`scale=False`, `stol=0.1`) | also report composition-level novelty (247/369) and the 25 rediscoveries hidden by S2 order |

Every named compound's evidence statement does mark novelty "checked". For LaRhH₂ that check
does not support the word "novel" in the doc text.

### F7 (S2) — What to run next with ≤ 48 h per job (12 usable threads)

Cost anchors (meas., doc 30 §5.3 and doc 15 §4):
- C1 5d hydride S5: 3–6.5 core-h.
- C2 (9-atom Fm-3m) S5: ≈ 30 core-h at 60 Ry. At 80 Ry with a small-σ ladder, ≈ 1.5–2× (est.).
- Mg₂IrH₆ q 4³ / k 12³ / k_fit 48³ at 80 Ry: ≈ 34 h on 4 ranks (doc 15 §4, extrapolated).

Priority, each ≤ 48 h:

1. **Mg₂IrH₆ S5 re-run with the corrected profile** (owned by the other agent). This gates every
   hydride number.
2. **Mg₂IrH₆ "S6-lite"**: q 3³ (≈ 4 h) and q 4³ (≈ 34 h), each a separate job on 4 ranks, at
   80 Ry with k_fit a multiple of k_scf.
   - This gives the two q meshes that G-VALID asks for, within the 48 h cap.
   - It is the first A5→6-type calibration point on a 0-GPa hydride (doc 13 §6 asks for exactly
     this).
   - The full q 6³ profile (≈ 9 days) is out of reach. The S6 tiers need a cap-aware mesh pair
     for C2 cells, recorded as a deviation.
3. **Ba₂IrH₆ S5** (c-539acaa4a6e7, corrected profile, ≈ 10–15 h on 4 ranks). It is isoelectronic
   with Mg₂IrH₆ and is a direct test of the A-cation effect.
4. **Finish the running LaRhH₂ q 4³.** Do not promote it further: its soft mode makes it an
   instability candidate, not a high-A one.
5. **BaRhH₆ S4 re-check** (AFM + SOC scf, < 3 core-h) **before** its S5.

### F8 (S3) — Smaller items

- **B11 (missing `targets/chemsys_list.txt`) is still unfixed.** Campaign 2 must not run G1 with a
  silent fallback. Make the loader fail hard.
- **`score_s3` should gain two free features** (ranking only, R2-compliant):
  - valence-electron parity / distance from an 18-electron count per MHₙ unit (doc 02 S7). It
    would have pre-sorted every S4 H1 outcome of this run correctly;
  - an `ML_BELOW_DFT_HULL` penalty (F4).
- **B3 symmetrization idealises to the 0.01 Å label.** For BaRhH₆ that is Cm, while the cell is
  R-3m at 0.1 Å. That costs k/q symmetry in QE. Prefer the highest group within the S2/S3
  tolerance, after confirming that MatterSim forces stay small.
- **Doc 30 §4 mixes denominators:** E_hull "≤ 0.1: 185/353" against "of 369". 16 G1 cells have
  no ML hull. Say why.
- **The export snapshot predates two running jobs.** The S6 tier_runs (1087/1088) and LaRhH₂
  (1080/1081) are `started`. Regenerate before any commit that cites them.

---

## 2. Is there a path to A ≥ 1300 K at ≤ 1 GPa? (strategic)

### 2.1 Where the evidence stands

- **Floor:** a0·T\* = 1299 K (doc 01 §0.4).
- **Best 0-GPa hydride with a validated calculation on this machine:** Mg₂IrH₆, A ≈ 990 K at
  σ = 0.005 Ry, q 3³ (doc 15 §3), i.e. ≈ 76 % of the floor. This is consistent with Dolui et al.
  (λ ≈ 2.3–2.5, Tc 160–175 K): K. Dolui et al., PRL 132, 166001 (2024), arXiv:2310.07562, cited
  in doc 15.
- **Other A₂MH₆ predictions:**
  - Sanna et al. (npj Comput. Mater. 10, 44 (2024), arXiv:2310.06804) give Mg₂XH₆
    (X = Rh, Ir, Pd, Pt) at 45–80 K, and > 100 K for electron-doped Mg₂PtH₆.
  - Zheng et al. (Mater. Today Phys. 42, 101374 (2024), arXiv:2401.17024; abstract read
    2026-10-30) screened X₂MH₆ cubic hydrides: ~26 dynamically stable with strong coupling;
    Mg₂RhH₆, Mg₂IrH₆, Al₂MnH₆ and Li₂CuH₆ highlighted at Tc > 50 K.
- **The ceiling literature** (Gao et al., Nat. Commun. 16, 8253 (2025), arXiv:2502.18281, doc 02
  §3.3) puts designable ambient Tc at 100–120 K and calls 300 K "extremely unlikely".
- **Synthesis:** Mg₂IrH₆ itself is not made. Only Mg₂IrH₅ (< 28 GPa) and Mg₂IrH₇ (> 40 GPa) form:
  Hansen et al. arXiv:2406.09538; Sinha et al. arXiv:2602.23675; doc 02 §3.2.
- **Kinetic instability** of the lightest-cation variants (Li₂AuH₆, Li₂AgH₆; doc 02 §3.2,
  Ding/Chen/Shi) removes the most obvious "lighter A raises ω" route.
- Dihydrides such as LaH₂ reach ~11 K at most even in the most favourable metastable H
  configuration: Kim et al., arXiv:2311.17328 (preprint). That is consistent with this campaign's
  LaH₂ λ ≈ 0.17.

### 2.2 What would have to happen to close the remaining ~30 %

A = λ⟨ω⟩/2. From Mg₂IrH₆ (λ ≈ 2.4, A ≈ 990 K), A ≥ 1300 K needs one of:
- **(i)** about +30 % λ at fixed ⟨ω⟩. This is limited by H4 (λ ≤ 4.7, flag at 3.7) and by the
  T1g-libration instability that already carries much of λ (doc 15 §0.3).
- **(ii)** about +30 % ⟨ω⟩ at fixed λ. That means a stiffer [MHₙ] unit, a lighter M, or more of
  the coupling in the H stretch band (λ_>600 cm⁻¹ ≈ 1.19 of 2.42).
- **(iii)** both, partly.

Mechanisms consistent with the H1–H7 profile (doc 02 §5, S4/S7):

1. **E_F tuning onto the DOS peak** (doc 02 S7, "doped neighbour of an insulator"). Mg₂IrH₆'s
   N(E_F) is strongly σ-dependent, i.e. E_F sits near a sharp feature. Ordered A-site alloys
   in doubled Fm-3m cells move E_F by ±½ e per MH₆:
   - Mg/Li, Mg/Na (hole side);
   - Mg/Al, Mg/Sc (electron side).

   A rigid-band scan from the S4 DOS costs nothing. The 18-atom alloy cells are C3 (S5 ≈ 50–150
   core-h, est.).
2. **M-site 4d/5d substitution and mixed-M** at fixed electron count, e.g. Mg₂Ir₁₋ₓRhₓH₆ or
   Mg₂Ir₀.₅Os₀.₅H₆+δ. The stiffness of the M–H σ bond sets the stretch-band frequency.
3. **More H per complex at the same metallic count.** Candidates are A₂MH₇- and A₂MH₈-type cells
   with [MH₆] + interstitial H, or higher-coordination complexes ([ReH₉]²⁻-type; known salts are
   closed-shell insulators, **UNVERIFIED** for metallic analogues). This raises H_f and H_DOS.
   G1 produced 42 SUN cells with H_f ≥ 0.5, but none in this motif class was checked.

**Honest expectation.** The literature (§2.1) and this campaign both place the reachable region
at A ≈ 0.5–0.8 × floor, i.e. Tc ≲ 150–200 K. Campaign 2 should be designed to:
- map the A₂MH₆ neighbourhood with calibrated A;
- find whether *any* member exceeds Mg₂IrH₆;
- deliver a better-founded negative statement if none does.

Claiming a route to 300 K is not supported by anything in hand.

---

## 3. Campaign-2 plan

Budget envelope:
- GPU ≈ 6 h (shared).
- CPU ≈ 6 days × ~250 core-h ≈ 1500 core-h.
- Every QE job ≤ 48 h, ≤ 4–6 ranks, machine total ≤ 12.

### 3.1 Preconditions (no campaign-2 S5 removal before these hold)

| # | Item | Cost |
|---|---|---|
| P1 | Corrected S5 profile (σ statistic, cutoff ≥ 80 Ry for H, k_fit % k_scf = 0), re-run Mg₂IrH₆ and require A_coarse within ±25 % of doc 15's 990 K (owned by the other agent) | ≈ 40 core-h |
| P2 | Mg₂IrH₆ S6-lite (q 3³ + q 4³ at full settings, two ≤ 48 h jobs) → first hydride entry in `calibration`; re-derive e_c for hydrides | ≈ 150 core-h |
| P3 | Fix: CLI gate bypass (F5); dedupe order/representative (F3); fail-hard chemsys list (B11); `DEFERRED` instead of `SCREENED_OUT` outside calibration domains (F1) | code |
| P4 | Re-cut the campaign-1 funnel with D-first dedupe and composition-level novelty | minutes |

### 3.2 Lanes

| Lane | What | Systems | Size | GPU / CPU |
|---|---|---|---|---|
| **P (new, deterministic)** | K₂PtCl₆ Fm-3m prototype substitution: A₂MH₆ with A ∈ {Li, Na, Mg, Ca, Sr, Ba, Al, Sc, Y, La} and M ∈ {Mn, Fe, Co, Ru, Rh, Os, Ir, Ni, Pd, Pt, Cu, Ag, Au}; plus ordered A/A′ alloys in doubled cells, at valence counts 18 ± ½, ± 1 per MH₆ | the Mg₂IrH₆ neighbourhood | ≈ 130 + ~300 alloys | S3 ≈ 1 GPU-h |
| G1′ | MatterGen `chemical_system_energy_above_hull`, **quaternary** A–A′–M–H systems around the P-lane hits (Mg–Li–Ir–H, Mg–Na–Ir–H, Mg–Al–Ir–H, Mg–Sc–Ir–H, Mg–Ir–Rh–H, Mg–Ir–Os–H, Ca–Mg–Ir–H …); E_hull conditioning U(0, 0.15) | 12–16 systems | ~500 | ≈ 3 GPU-h |
| G2 / G3 | `space_group` ∈ {225, 221, 229}; `dft_band_gap = 0` (needs the user to approve the checkpoint download, doc 30 §10.4) | G1 ternaries | ~300 | ≈ 1.5 GPU-h |
| D | Alex-MP, all X₂MH₆ / XMH₃ / X₂MH₇ entries in the systems above, any E_hull ≤ 0.2; **run S2 first** | same | ~200 | – |
| C | Mg₂IrH₆ (+ Ba₂IrH₆ as a second anchor once run), Al, Nb | – | – | – |
| G4 | **paused** until a non-hydride hull calibration exists (F1a) | – | 0 | – |

### 3.3 Stage rules for campaign 2

- **S3 ranking (no removal):**
  - add valence-e parity / 18-e distance;
  - add `ML_BELOW_DFT_HULL`;
  - per-lane quotas: P 40 %, G1′ 30 %, D 20 %, G2/G3 10 %;
  - rank cut 40 for S4.
- **S4:** ≈ 40 candidates × ≤ 3 core-h ≈ 120 core-h.
  - Add the AFM + SOC scf for d⁵ M⁴⁺ complexes (H2).
  - Record a rigid-band ΔN(E_F) for ±0.5 e from the S4 DOS (ranking only).
- **S5:** corrected profile, ≤ 15 runs.
  - C1/C2 cost ≈ 30–60 core-h each (est.); 18-atom alloy cells only for the top 3 by rigid-band
    N(E_F).
  - ≈ 700 core-h total.
  - **Removal only inside the re-derived hydride e_c (P2).** Until then every S5 outcome is
    `DEFERRED(uncalibrated)` or `PASS`.
- **S6-lite:** only for candidates with A_coarse ≥ 0.75 × A_coarse(Mg₂IrH₆) at matched settings.
  At most 2 runs, ≈ 300 core-h, each mesh ≤ 48 h.
- **S7a** (own-DFT hull + ZPE): for every S6-lite candidate and for any S5 result above Mg₂IrH₆.
  ≈ 20 core-h per system.

### 3.4 Success criteria (pre-register before running)

- **Primary:** does any non-control compound have A_coarse (corrected) > A_coarse(Mg₂IrH₆) at
  matched settings, with E_hull(DFT, ZPE) ≤ 0.10 eV/atom?
  - If yes → S6-lite and S7a.
  - If no → publish the mapped A(composition, electron count) surface of the A₂MH₆ family as the
    campaign's negative result.
- **Secondary, for MatterGen:** G1′ vs P vs D on (S4 metal fraction, A_coarse distribution),
  Mann–Whitney as in doc 12 §3.3, reported and not gated. Report composition-level novelty and
  rediscovery counts.

---

## 4. Summary

- **Removals.** None of the 38 + 19 + 8 + 1 removals lost a plausible high-A candidate. Two rules
  were used outside their stated calibration:
  - S3 H7-ml window on 34 G4 non-hydrides; one of them is under the window on the DFT-reference
    hull;
  - S1 n_atoms on lane D (19 entries).

  The sole S5 removal (LaH₂) sits on a σ-truncated ladder of the kind the Mg₂IrH₆ control shows to
  under-estimate A by ×2.9. It survives only because its margin is ~20×. The H1 removals are sound
  (gaps ≥ 1.7 eV, and electron parity agrees in every case).
- **Dependence on the σ/cutoff issue.** It touches the LaH₂ and ScH₂ envelopes, every quoted
  A_coarse (the "8 %" and "~110 K"), the §6 order and the S6 choice. It does not touch S1–S4. Doc
  30's "Mg₂IrH₆ passes S5 by < 2 %" is wrong: the margin is 4.7×. The problem is A accuracy, not
  recall.
- **MatterGen.** Composition conditioning works (98 % in-system; 155 SUN-ML of 369; 110 with new
  compositions). Property enrichment is not shown:
  - the two best "generated" DFT results are Alexandria rediscoveries;
  - S2 order hid 25 such cases from the funnel;
  - G1's metal fraction (64 %) equals lane D's (67 %).
  - Several "stable" SUN cells sit 60–90 meV/atom below the DFT hull, which suggests MatterSim
    error.
- **Novelty.** BaRhH₆ and BaIrH₆ are new compositions (no literature hit; ICSD/OQMD/GNoME not
  checked). BaRhH₆ is R-3m, not "sg 8 C2/m-like", and as a d⁵ Rh⁴⁺ complex it needs an AFM/SOC
  check. LaRhH₂ is a new polymorph of a composition already in Alexandria.
- **S6.** NaIrH₃ (λ ≈ 0.2, not novel, S5 ERROR, gate bypassed) is a chain test, not a candidate.
  The next ≤ 48 h jobs:
  - corrected Mg₂IrH₆ S5;
  - Mg₂IrH₆ S6-lite (q 3³ + q 4³ as separate jobs);
  - Ba₂IrH₆ S5;
  - BaRhH₆ magnetic re-check.
- **Strategy.** Only the A₂MH₆ family is within ~25 % of the floor (Mg₂IrH₆ ≈ 990 K). Campaign 2
  should map it by deterministic prototype substitution plus quaternary-conditioned MatterGen,
  with E_F tuning via A-site alloys, under calibrated hydride margins. The expected result is
  A ≲ 0.8 × floor. A route to 300 K at ≤ 1 GPa is not supported by any evidence in hand.
- **Not verified:**
  - ICSD/OQMD/GNoME novelty;
  - the J_eff = ½ analogy for hydrides;
  - the S6 admission risk for NaIrH₃;
  - all cost extrapolations marked (est.);
  - the validity-test thresholds of S2;
  - claims of a 2025 ML paper (arXiv:2511.03865) about doped Mg₂IrH₆ that surfaced in a search
    summary but are **not** in its abstract. Not used.
