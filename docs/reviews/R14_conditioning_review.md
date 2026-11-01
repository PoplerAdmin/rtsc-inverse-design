# R14 — Adversarial review of doc 14 (conditioning MatterGen on A = ∫α²F dω)

Reviewer agent R14, 2026-10-29. Read-only on compute: no GPU, no QE, CPU work at nice 19 with
≤ 2 threads. Nothing was committed.

**Reviewed:**
- `docs/design/14_property_conditioning.md` (cited as "D14:line")
- `data/LABELS_SOURCES.md`
- `data/labels/{labels_all.tsv.gz, manifest_log10A.tsv.gz, summary_log10A.json, adherence_summary.json}`
- `src/ml/labels/*`
- the D1 scratch outputs `~/rtsc_targeted_work/D1-labels/{labels,gen,raw}`, read only

**Reviewer scripts:** `docs/reviews/r14_scripts/`, run with `~/.rtt/envs/py/bin/python -I`.
- `smear.py` recomputes A, λ and m = λ⟨ω²⟩ at all 10 smearings for all 8,249 entries of S1 (the
  Materials Cloud 3kbt5-r3n56 set). It took 10 min on 1 thread. Output:
  `~/rtsc_targeted_work/R14-review/smear.jsonl`.
- `a2.py`: leakage and adherence statistics.
- `a4.py`: smearing tail.
- `a5.py`: z-scores and variance decomposition.
- `a6.py`: cross-source agreement and novel-only kNN.
- `a7.py`: correlations with DOS(E_F).

**Verdict in one line.** The main negative conclusion is correct and robust: no open 0–1 GPa
structure+α²F label comes near A = 1300 K, and the adapter cannot reach that target. Four claims
are overstated:
- the headline ceiling of 346 K (it is a smearing convention);
- the z-scores of the targets (arithmetic error);
- "every proxy is monotone up to 300 K" (the best proxy plateaus at 100 K);
- the statistical support.

The doc also skips the acceptance test that doc 12 requires before the adapter can be used.
The strategy section (§7 here) recommends dropping A as the *driving* condition.

---

## 1. Findings, ranked by severity

### F1 (High): the label ceiling of 346 K depends on the smearing convention, and the convention contradicts the pipeline's own decision rule

**What the doc does.** It labels every S1 row at smearing #6 (σ = 0.030 Ry), because that
reproduces the authors' `McMillan.dat` (D14:77–79; `compute_labels.py:96–101`). The "label
uncertainty" it reports covers only σ #4–#8, i.e. 0.020–0.040 Ry (D14:80–81;
`compute_labels.py:100`). It never reports σ = 0.005–0.015 Ry.

**Why that conflicts with the pipeline.** The pipeline's own decision rule is the upper ENVELOPE
over admitted smearings (`docs/theory/01_eliashberg_targets.md` §0.4; doc 13 §0). Its positive
control shows exactly this sensitivity: Mg2IrH6 gives A ≈ 990 K at σ = 0.005 Ry, and λ falls
from 2.42 to 0.91 between 0.005 and 0.05 Ry (doc 15 §0, item 4).

**Recomputed tail of the 8,058 S1 training rows, from the archive's own a2F.dos1–10 (`a4.py`):**

| σ (Ry) | max A (K) | n(A ≥ 300) | n(A ≥ 500) | n(A ≥ 1000) | max m (K²) |
|---|---|---|---|---|---|
| 0.005 | **859** (Li2ScCuO4) | 19 | **6** | 0 | 2.16×10⁶ |
| 0.010 | 671 | 13 | 2 | 0 | 1.48×10⁶ |
| 0.020 | 416 | 10 | 0 | 0 | 1.07×10⁶ |
| 0.030 (used) | 346 | 5 | 0 | 0 | 9.35×10⁵ |
| 0.050 | 345 | 2 | 0 | 0 | 8.57×10⁵ |

- A(0.005)/A(0.030): 5/50/95/99 % = 0.72 / 1.09 / 1.48 / 1.85. The median is 1.14 for hydrides
  and 1.09 for the rest.
- Rows whose A is far above 300 K at 0.005 Ry but not at 0.030 Ry: NaH3Pd 744 vs 271, NaH3Au
  645 vs 297, NaC6 565 vs 331, LiH3Pt 532 vs 287, K(HPd)2 500 vs 124.

**Reading.**
- At σ = 0.005 Ry and S1's k-grids, the values are probably under-converged. Small-σ values are
  an upper end of the label uncertainty, not truth.
- Even so, the doc's phrases "largest among training rows: 346 K", "zero entries have A ≥ 500 K"
  and "a factor of 3.8" (D14:22–23, 109) hold only under one convention.
- Under the envelope convention the ceiling is ~860 K. The gap to 1299 K becomes a factor of
  **1.5**, not 3.8.
- The conclusion "no label ≥ 1000 K" survives under every convention.
- **Fix:** store A at all 10 σ, apply mode-positive admission (as in the core's Pb fix / doc 13),
  and report both the σ = 0.030 label and the admitted upper envelope. Retraining on the envelope
  is option 7.2-A below.

### F2 (High): "every proxy is monotone 30 → 100 → 300 K" is false for the best proxy; saturation starts at about 100 K, not 300 K

D14:31–37 and 204–206 claim monotone adherence. The per-structure data in
`gen/proxies_comp.json` say otherwise (`a2.py`, `a6.py`).

**Composition kNN Â, the best proxy** (validation Spearman 0.69):

| Target | mean kNN Â (K) | bootstrap 95 % CI | novel formulas only |
|---|---|---|---|
| A100 | 82.7 | 76–89 | 84.0 |
| A300 | 76.5 | 69–85 | 77.4 |

- A300 vs A100: Mann–Whitney **p = 0.40**.
- A1500 vs A300: p = 0.63.
- The doc shows medians (79 → 81 K), which hide that the means *decrease* from A100 to A300.

**x_H** (fraction of H atoms): A300 vs A100 p = 0.12, so it does not separate either.

**What still moves from 100 to 300 K:**
- x_light (fraction of H, Li, Be, B, C, N, O, F atoms): p = 1.4×10⁻⁶;
- mean atomic mass;
- the MatterSim ω_max.

These are composition variables whose Spearman correlation with A is only 0.22 (x_light) in the
training set (`a5.py`). From 100 to 300 K the sampler is therefore drifting in *chemistry*, not
in coupling.

**Metal-free samples.** Under a plain "contains no metal" rule, 16 of the 64 A300 samples are
metal-free (the doc reports 11/64; the doc does not define its rule).

**Fix:** state that conditioning has an effect only up to roughly the 95th percentile
(~100 K), and that above that it only lightens the chemistry.

### F3 (High): the doc mis-states the z-scores of its own targets

D14:109 gives z = +2.1 (300 K), +3.8 (800 K), +6.2 (1500 K). The training statistics are mean
1.6392 and sd 0.2487 (`summary_log10A.json:833–836`), which give:

| Target | z |
|---|---|
| 300 K | **+3.38** |
| 800 K | **+5.10** |
| 1500 K | +6.20 (correct) |
| 1299 K | +5.94 |
| 2060 K | +6.74 |

**Consequences:**
- "A = 300 K" is not "inside the label range". Only 5–7 of 8,675 rows sit at or above it. The
  "in-range vs out-of-range" framing of §4.2 rests on the wrong z.
- **Mechanism the doc misses: NoiseLevelEncoding aliases.**
  - The property embedding (`rtt_log10A.yaml`; `mattergen/diffusion/model_utils.py:63–73`)
    encodes the *standardised* value z as sin/cos(z·k) with k = 10000^(−2i/d).
  - The fastest channel has k = 1, so its period in z is 2π ≈ 6.28.
  - At z = 6.20 (1500 K) that channel equals its value at z = −0.08, which is the training mean.
  - At z = 5.10 (800 K) it equals z = −1.18.
  - Training covers z ∈ [−3, +3.6] only. Far targets are therefore partly encoded as *ordinary*
    labels.
  - This is the likely mechanism of the plateau and fall-back at 800/1500 K in D14:211–215. It
    also means no γ setting can make this encoder extrapolate to 1300 K. Asking for 1300 K is
    ill-posed for this encoder, not merely uncertain.

### F4 (High): the adapter skipped its acceptance test, and the baseline used is the wrong one

Doc 12 §3.3 (lines 265–274) sets the acceptance test:
- target q_hi = 90–95th percentile, "never 1300 K";
- **500 samples**, compared against lane G1 (chemical-system conditioning) at equal GPU cost;
- outcome measured on S5 coarse-DFPT A, ≤ 10 runs per lane.

Doc 14 instead:
- uses 64 samples;
- compares against its own fine-tuned unconditional model;
- runs no DFPT.

The two strongest effects (x_light, H fraction) are exactly what G1 gets for free by
conditioning on X–Y–H systems. Doc 14 does not show the adapter beats G1 on anything. Under
doc 12 the adapter is **not accepted** and must not replace any G1 share yet. A small naming
mismatch: doc 12 calls the property `rtt_A_K`; the implementation uses `rtt_log10A`.

### F5 (Medium-High): full fine-tuning damaged the prior, and outputs are partly retrieval

Fraction of generated formulas that also appear in the training set (`a2.py`):

| Model | formulas in train | kNN distance |
|---|---|---|
| Base model, 64 samples (`gen_uncond64`) | **1.6 %** | 0.88 |
| Fine-tuned, unconditional | **18.8 %** | 0.42 |
| Fine-tuned, conditioned (by target) | 14–22 % | — |

- The doc notes the shift (D14:159–164) but not the 12-fold rise in exact training-formula
  reproduction.
- With lr 5×10⁻⁶ and about 12k optimiser steps on 7.8k rows, the network has partly memorised
  S1. S1 is itself an ML-preselected, near-hull Alexandria subset (D14:53).
- Novelty is what the pipeline needs, and full fine-tuning removes it.
- The adapter-only ablation (D14:256–258) should come before any further use.

### F6 (Medium): the statistics overstate support

- **Multiple comparisons.**
  - The adherence run makes 7 composition tests × 6 sets plus 3 phonon tests × 5 sets ≈ 57
    Mann–Whitney tests (`gen_proxies.py:179–181`; `gen_phonon_proxy.py:251–253`). Bonferroni
    α = 0.05/57 ≈ 9×10⁻⁴.
  - Reported "effects" that fail this threshold: A100 x_H (p 0.02), mass (0.04), ω_max (0.03)
    and phonon proxy (0.005); A1500 x_H (0.002); A800 x_H (0.035).
  - D14:259–261 says "several in-range effects … survive any reasonable correction". That is true
    for the kNN effect at A30/A100 and the x_light effect at A100/A300, and false for most of the
    table.
- **Proxies are not independent.** x_light, x_H, mass, ω_max and the phonon proxy are all driven
  by composition. The phonon proxy is literally linear in x_light and x_H (`calib_proxy.py:32`).
  "Every proxy moves" amounts to about two signals, not seven.
- **Survivorship in the phonon columns.** Medians are taken over MatterSim-stable subsets whose
  size changes from 48/64 (uncond) to 25/64 (A300). The sets being compared are different
  sub-populations. Doc 15 shows the MatterSim stability filter removes the soft, strongly
  coupled hydrides preferentially.
- **Seed and atom counts.** One seed gives identical atom counts in every set (MWU p = 1.000 for
  nsites in `adherence_summary.json`). That is not a confound *between* sets. It does mean all
  results rest on a single draw of atom counts and a single noise realisation. Seed-to-seed
  variance is unknown, so the bootstrap CIs above are lower bounds on the true uncertainty.
- **Calibration reuse.** The phonon-proxy coefficients used to score generated sets are a full
  fit on the same 370 test rows that report the CV Spearman (`calib_proxy.py:76`). This is minor,
  but no held-out set is left.

### F7 (Medium): the kNN proxy is circular and cannot show what the doc needs

- kNN Â is the median *training* label of the 5 nearest training compositions
  (`gen_proxies.py:128–131`). By construction it cannot exceed 346 K.
- It measures whether the generator reproduces the composition-to-label association it was
  trained on. That is the training objective, not physical A. Agreement is expected even if the
  generated structures are not metals or not stable. Metal-free cells still get a kNN value from
  their nearest metallic neighbours.
- **Label noise between sources is as large as the proxy's error.** For the 275 (formula, nsites)
  pairs present in both S1 and JARVIS, A_JARVIS/A_S1 has a 5–95 % range of **0.79–1.56**
  (Spearman 0.89; `a6.py`). MgB2, for example, is 197 K in S1 and 302 K in JARVIS. That is the
  same size as the kNN MAE of 0.13 dex.
- kNN adherence can only rank compositions inside the known label cloud. The doc says this
  (D14:182–185), but then still counts kNN as adherence evidence in §0 item 4.

### F8 (Low-Medium): leakage between splits is real in structure but does not inflate the metrics

- The split is by reduced formula (`build_dataset.py:297–305`). Polymorphs of one formula share a
  split, and so do cross-source duplicates (285 formulas). That part is sound.
- Near-duplicates are not excluded: **98.8 %** of test rows have a single-element-substitution
  analog in train (same anonymised formula, all elements but one shared). 47 % share a chemical
  system with a train row.
- But a chemical-system-disjoint resplit gives kNN Spearman **0.70** (MAE 0.127 dex vs baseline
  0.179). That is the same as the formula split (0.69). Performance also falls only gently with
  distance to the nearest training row: Spearman 0.78 / 0.68 / 0.61 for L1 < 0.25, 0.25–0.5 and
  > 0.5 (`a2.py`).
- So prototype leakage does not change the reported numbers. It would matter for a learned
  surrogate. Use chemical-system or prototype grouping for any future surrogate validation.
- Cross-source duplicates are double-weighted in training (MgB2 appears 3 times). Prefer one
  source per material, or average them.

### F9 (Low): label arithmetic and gates — mostly correct, with caveats

**Verified correct:**
- **Units.** The QE a2F.dos header says "frequencies in Rydberg" (checked in the raw archive).
  1 Ry = 157,887.5 K, matching `compute_labels.py:34–35,91`. S1 A matches the authors' integral.
- **JARVIS frequency unit.** Here λ is the wrong check, because λ does not depend on the
  frequency unit. ω_log is the right check: recomputed/published ω_log has median 1.003 and
  5–95 % range 0.94–1.16 over stable entries. So meV is right.
- **Integration.** Trapezoid integration on the native grid, ω > 0, negative α²F clipped. The
  negative-frequency weight is negligible (A_neg/A median 0, 99th percentile 0.55 %).
- **Definitions.** A, λ and m = 2∫ωα²F dω follow theory doc 01 §0 (line 32).

**Caveats:**
- **The ≤ 20-atom gate is vacuous.** No row in either source exceeds 20 atoms.
- **The stability gates are inconsistent.** S1 uses the minimum frequency along the dispersion
  path, ≥ −30 cm⁻¹. 409 kept rows lie between −30 and −5 cm⁻¹. JARVIS uses the authors' own flag.
  For soft-mode rows λ is unreliable: AgH3Pt has recomputed λ 3.0 vs published 1.2 and
  ω_log 1.3 K. It is excluded and λ is not a label, but the same soft modes are where A is most
  σ-sensitive.
- **Pressure.** All rows are 0 GPa, but in PBEsol. PBEsol cells are typically a few per cent
  smaller in volume than PBE cells; the effective PBE pressure is **UNVERIFIED**. Downstream
  DFPT here is PBE (SG15). This is a small systematic difference, not a bug.

### F10 (Low): the data search missed a few things, none of which changes the conclusion

| Source | Licence | Content | Effect on the ceiling |
|---|---|---|---|
| Antiperovskite EPC set, Materials Cloud `2hg1p-63427` (doi:10.24435/materialscloud:6g-az) | CC BY 4.0 | 384 DFPT calculations; reported max Tc ≈ 18 K | None. It sits in the D1 scratch record dump but is not listed in `LABELS_SOURCES.md` |
| High-pressure hydride set, figshare article 25270012 (arXiv:2312.12694; doi:10.1088/2752-5724/ad4a94) | CC BY 4.0 | 8.25 GB `raw-dft-hydride.tar.gz`; ~900 hydrides at 0–500 GPa; whether the archive holds α²F files is **UNVERIFIED** | The only open source likely to contain A ≥ 1000 K labels with structures, but all at high pressure (see §7.2-H) |
| BEE-NET `database.json` (the public BEE-NET GitHub repository) | — | 5,241 entries, all `agm…` ids, i.e. a **subset of S1** | Resolves the "open release not located" item in `LABELS_SOURCES.md`: no new labels |
| MC3D-BCS 2025, Materials Cloud `0kqen-w8h90` (doi:10.1103/sb28-fjc9) | CC BY 4.0 | 250 EPW α²F inside `EPW-screening.aiida` (8.4 GB); experimentally known compounds | Better label quality, but no tail. D14's "15 GB" is the QE-screening archive; the EPW archive alone is 8.4 GB |
| GNoME ambient hydride search, arXiv:2508.19781 (Commun. Phys. 9, 94 (2026)) | — | max Tc 17 K; no data release located | None |
| Hydride Tc compilation, arXiv:2512.20228 | — | 2,059 records; Tc only, no α²F | Not usable as labels |

**Conclusion of F10.** Item 3 of the brief ("cannot reach A ≥ 1300 K because no labels exist
there") is **correct**: no open 0–1 GPa α²F with a structure has A ≥ 1000 K. At ambient
pressure, the highest literature values are single points such as Mg2IrH6 (PRL 132, 166001
(2024); in-house coarse A ≈ 990 K at σ = 0.005 Ry, doc 15). A deeper reason: an open,
structure-resolved 0 GPa α²F with A ≥ 1300 K would itself be close to the answer this project
is looking for. Its absence is expected, not a gap in the search.

---

## 2. Physics finding that drives the strategy

Three facts from the training data (`a5.py`, `a7.py`):

1. **A is controlled by coupling, not by phonon scale.**
   - Var(log A) = 0.062, Var(log λ) = 0.069, Var(log⟨ω⟩) = 0.024, with a negative covariance
     (2·cov = −0.031).
   - Spearman correlation with A: λ 0.81, ⟨ω⟩ 0.19, x_H 0.23, x_light 0.22.
   - This matches doc 02 §1 (λ = η/(M⟨ω²⟩)): the variable that moves A is the Hopfield factor
     η = N(E_F)⟨I²⟩. The generator does not see it, and MatterSim cannot compute it.
2. **λ⟨ω²⟩ does not help as a conditioning label.**
   - Training z of the rigorous floor: m ≥ 2.70×10⁶ K² (NM rule, μ* = 0) sits at z = +5.96.
     The A floor sits at +5.94.
   - The A bound is the binding one in 100 % of training rows: A/4.33 < 0.1827·√m everywhere.
   - The bound-Tc ceiling of the training set is 80 K.
   - Conditioning on m is the same extrapolation problem.
3. **An electronic descriptor carries more of A than any phonon proxy.**
   - The JARVIS `alex_supercon` DOS(E_F), joined one-to-one to S1: Spearman 0.43 with A, and 0.57
     between λ and DOS(E_F) per atom.
   - The MatterSim phonon proxy gets 0.29; DOS(E_F)·(1 + x_light) gets 0.47.
   - DOS(E_F) needs one static DFT, minutes per cell, not DFPT.

What the generator *can* control well is chemistry: x_light goes from 0.16 to 0.63. What
decides A is electronic coupling, which no current label supplies at volume.

---

## 3. Compute costs on this machine (basis)

Measured unless marked:

| Step | Cost |
|---|---|
| Full fine-tune, 7.8k rows × 200 epochs (D14:136) | 1.7 GPU-h |
| Conditional generation (doc 12 §3.3; D14:157, on a shared GPU) | ≈ 10–19 s per structure (64 in 633–1,220 s) |
| MatterSim relax + phonons (doc 10 §3) | 2–6 s per structure |
| Coarse DFPT, Mg2IrH6, 9-atom cubic cell, q 3³ (doc 15 §1) | **3.0 h on 4 ranks** |
| Coarse DFPT, pipeline budget (doc 12 §1) | **2–6 per day** |
| Coarse DFPT, 6–10-atom low-symmetry generated cells | **UNVERIFIED**; expect several times more irreducible q and k → 10–30 h on 4–6 ranks |

The CPU is the binding resource. On this machine DFPT, not the GPU, sets the pace of any loop
that makes labels.

## 4. Strategy options under a ~350 K (σ = 0.030) / ~860 K (envelope) label ceiling

| # | Option | GPU | CPU | Expected value | Verdict |
|---|---|---|---|---|---|
| A | Re-label with the admitted-σ upper envelope plus σ-spread weights, then retrain | 1.7 h | 10 min (done here) plus admission | Tail grows from 5 to ~19 rows ≥ 300 K and the ceiling to ~860 K. Still < 1300, and the small-σ values are noisy | **Do**, cheap. Keep both labels |
| B | Condition on λ⟨ω²⟩ instead of A | 1.7 h | 0 | None (z +5.96 vs +5.94; A always binding) | **Skip** |
| C | Decomposed chemistry targets (x_H, x_light, H-network motif) | 0 (the shipped `chemical_system` checkpoint already does it) | 0 | Hits the necessary-condition chemistry directly with no label ceiling | **Do** (it is G1/G4 in doc 12) |
| C′ | ω_max as a learned condition (MatterSim phonons on ~50k Alex-MP-20 cells, then fine-tune) | ~2–3 GPU-days of labels + ~3 h fine-tune | 0 | Low: ω_max explains little of A (Spearman 0.25), and MatterSim softens H modes 15–20 % and gets libration signs wrong (doc 15) | **Filter only**, do not condition |
| D | Adapter-only fine-tune (`adapter.full_finetuning=false`), same labels | ~1.5 h | 0 | Restores the base prior (formulas in train 1.6 % vs 18.8 %). Needed before any use of the A-adapter | **Do** |
| E | Guidance and seed sweep on D: γ ∈ {1, 2, 4} × targets {100, 150, 250 K} (z ≤ +3) × 2 seeds × 128 | ≈ 6–7 h + ≈ 2–4 h MatterSim | 0 | Gives a real dose–response with seed variance. Do not test targets above z ≈ 3.5 (F3 aliasing) | **Do**, as one batch |
| F | Doc-12 acceptance test: adapter (D, best γ) vs G1 at equal GPU, 500 samples each, top 10 per lane to coarse DFPT | ≈ 3 h | 20 coarse runs ≈ **4–10 days** | The only *decisive* test of whether A-conditioning adds anything over chemistry conditioning | **Do**, gate for option G |
| G | Active learning, generate → coarse DFPT → add labels → retrain MatterGen | ~6 h per round | 2–6 labels/day → **60–180 labels/month** | Tens of labels barely move an 8.7k-row diffusion model, even with up-weighting. Breaking the tail would take ≳ 100 labels at ≥ 300 K, i.e. months, if such compounds are found at all | **Not for MatterGen retraining yet.** Feed the labels to the doc-12 G1 system-reweighting and to a ranking surrogate instead |
| H | Add high-pressure hydride α²F (figshare 25270012) with a pressure covariate, sample at P = 0 | ~3 h (after the 8 GB download and parsing) | parsing | High-A rows exist only at high P. Classifier-free guidance cannot separate "high A" from "compressed cell", and the base model's density prior is ambient | **Skip** for generation. At most use the set to train or validate a surrogate |
| I | Cheap electronic label at volume: DOS(E_F) and H-projected DOS(E_F) from one static DFT, as a ranking feature or second condition | fine-tune ~2 h once ≥ 5k labels exist | S4 ≈ 25/day (doc 12); the 8,253 `alex_supercon` DOS(E_F) values exist already | Targets the variable that dominates A (§2: 0.43–0.57 Spearman). Its high tail is not capped by EPC data | **Do as a ranking feature now**; conditioning later |

## 5. Recommendation

1. **Drop A as the driving condition. Make chemistry conditioning plus filters the backbone.**
   - Use the shipped checkpoints: G1 chemical-system X–Y–H and G4 light-element non-hydrides, with
     energy_above_hull and band_gap = 0 where available.
   - Then MatterSim relax and hull filter (never MatterSim stability alone for hydrides, doc 15),
     an S4 electronic screen (DOS(E_F), H-DOS share), and S5 coarse DFPT.
   - Cost: generation only, ≈ 8.6k structures per GPU-day. The bottleneck is S5 at 2–6 runs per day.
2. **Keep the A-adapter only as a minor lane, and only after D + E + F.**
   - Retrain adapter-only (≈ 1.5 GPU-h).
   - Sweep γ and seeds at targets ≤ 250 K (≈ 10 GPU-h total including MatterSim).
   - Run the doc-12 acceptance test against G1 (≈ 3 GPU-h + 4–10 CPU-days).
   - Use it only if it beats G1 on A_coarse. Never ask it for A > ~300 K: the encoder aliases far
     targets back toward the mean (F3).
3. **Relabel with the σ-envelope (option A)** for honesty about the ceiling. Retrain only if
   option F says the adapter is worth keeping.
4. **Spend CPU on labels that change the decision, not on retraining.** Every coarse-DFPT label
   goes into three places:
   - (i) the doc-12 G1 system-weight update;
   - (ii) a held-out validation set for any surrogate;
   - (iii) a small ranking surrogate on composition + DOS(E_F) + MatterSim ω features, validated
     with a chemical-system-grouped split (F8).
   - Revisit MatterGen retraining when ≥ ~100 in-house labels at A ≥ 300 K exist. At 2–6 labels
     per day, and assuming such compounds are found, that is ~1–3 months.
5. **Fix the doc** before it is cited:
   - the σ-convention ceiling (F1);
   - z-scores (F3);
   - the monotonicity claim (F2);
   - the multiple-comparison statement (F6);
   - the missing-source list (F10);
   - an explicit statement that the adapter has not passed the doc-12 acceptance test (F4).

## 6. What this review could not verify

- Whether the small-σ (0.005 Ry) S1 values are converged. The S1 k-grids differ per entry, and no
  per-mode admission was applied here.
- Whether figshare 25270012 contains α²F files. The archive was not downloaded (8.25 GB).
- Coarse-DFPT cost for 6–10-atom low-symmetry generated cells (estimate only).
- The units and provenance of `alex_supercon` `dosef`. It was used only as a rank correlate.
- Seed-to-seed variance of the adherence results. No GPU was used.
