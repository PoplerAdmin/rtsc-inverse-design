# RC1 — Adversarial code review: front end (S1–S3) and back end (S4, S5, report)

Reviewer: agent RC1-review-code, 2026-10-29.
Scope: `src/rtt/source/`, `src/rtt/structure/`, `src/rtt/mlff/` (B2, doc 21);
`src/rtt/dft/`, `src/rtt/score/s4.py`, `src/rtt/report/` (B3, doc 22).
Core modules (`records`, `ledger`, `decide`, `certify`, `engine`) were read for context and
**not edited**. Where a finding comes from how the engine calls B2/B3 code, the finding says so.

Spec used: `docs/design/12_architecture.md` rev 2 §3.1 and §4 (R2, R3, R5, R7, R8);
`docs/theory/01_eliashberg_targets.md` §0.4 and §5.1; `docs/theory/02_materials_profile.md` §4.1
and §5(a)/(b); `docs/design/10_ml_capabilities.md` §3.3–3.4.

Question asked of every path: can a candidate that might be a real high-Tc superconductor be
**removed** wrongly (SCREENED_OUT, EXCLUDED_ASSUMPTION, DUPLICATE, or a DEFERRED/ERROR that
nothing re-runs)? Or can it **pass** with a stronger claim than the data support?

**Tests.** `tests/review/test_rc1_findings.py` holds 17 tests (19 cases); each is named
`test_F<n>_...` after its finding. All 19 **fail** on the code as reviewed. They run on the CPU
in about 5 s:
`OMP_NUM_THREADS=2 nice -n 19 ~/.rtt/envs/py/bin/python -m pytest -q tests/review`.
The existing `tests/back` and `tests/front` suites still pass (79 passed, 4 skipped), so the new
tests do not disturb them. Each test encodes the behaviour the spec requires. It should pass once
the finding is fixed.

Severity:
- **Critical**: a hard filter is silently inactive, or the evidence statement asserts a check
  that was never made.
- **High**: a plausible route to wrongly removing a target-class candidate.
- **Medium**: wrong in a narrower set of cases, or an accounting error that hides removals.
- **Low**: robustness or provenance issues.

---

## Summary table

| # | Sev. | Where | One line | Test |
|---|---|---|---|---|
| F1 | **Critical** | `dft/electronic.py::_th`, `decide_s4`; engine `S4Electronic.decide` | H1 and H2 are never applied in the wired pipeline, and the evidence statement still says "checked" | `test_F1_...` |
| F2 | High | `dft/electronic.py::analyse_bands`, `decide_s4` | H1 uses "a band crosses E_Fermi on the mesh" instead of doc 02's N(E_F) > 0, so doped insulators with E_F at a band edge are EXCLUDED | `test_F2_...` |
| F3 | Medium | `decide_s4` | gap_eV_max = 0.0 with strict `>`: a 10 meV mesh-limited gap is EXCLUDED rather than deferred to the hybrid check doc 02 asks for | `test_F3_...` |
| F4 | High | `decide_s4` (H2) | a metastable ferromagnetic SCF solution excludes the candidate even when the NM energy is lower; E_spin is recorded but never compared | `test_F4_...` |
| F5 | Medium | `decide_s4` (H2), targets `moment_muB_tol` | the tolerance is per **cell** (0.05 μB, unsourced, provisional), so it tightens with cell size | `test_F5_...` |
| F6 | Medium | `decide_s4` (H2) | a ferrimagnet (total ≈ 0, absolute 4 μB) is a plain PASS; the evidence reads "H2 checked" | `test_F6_...` |
| F7 | Low | `decide_s4` | a PASS decided with the provisional H2 tolerance carries no provisional mark, so S9 reports "none used" | `test_F7_...` |
| F8 | Low | `decide_s4` (H7-heur) | removes on runs flagged `RELAX_PRESSURE_MISMATCH` (cutoff not converged) | `test_F8_...` |
| F9 | High | `dft/coarse.py::decide_coarse` | imaginary modes drive coarse removal indirectly: the weight at ω ≤ 0 is dropped (A-DROP-NONPOS) and the A margin still removes | `test_F9_...` ×3 |
| F10 | High | `decide_coarse` | multiband (n_sheets > 1) candidates are screened on the isotropic A, which doc 01 §5.1 shows is invalid without bound | `test_F10_...` |
| F11 | Low | `decide_coarse` | a missing ε_sys silently becomes 0, which is less conservative | `test_F11_...` |
| F12 | High | `mlff/screen.py::run_s3` (+ engine scheduler) | DEFERRED (ml_failed, rank_cut, relax_not_converged) is never revisited, so it is a silent permanent removal | `test_F12_...` |
| F13 | High | engine `scheduler.next_stage` (reference only) | ERROR at S4/S5 is terminal; matdyn NaN from imaginary modes ends as ERROR, which is the F9 route again | `test_F13_...` |
| F14 | Medium | `report/summary.py::control_recall` | DEFERRED and ERROR controls count as "passed", so R8 auto-disable cannot fire for them | `test_F14_...` |
| F15 | High | tiers `S3.ehull_self_max_eV` = 0.42 | the 0.07 eV/atom ML margin is smaller than one measured MAE (0.074 for H-containing entries) and ignores ZPE | `test_F15_...` |
| F16 | Medium | `structure/dedupe.py` | StructureMatcher on **unrelaxed** cells merges molecular H₂ (0.74 Å) with dissociated H (1.60 Å); the second becomes DUPLICATE | `test_F16_...` |
| F17 | Low | `dft/collect.py::collect_run` | the spectrum hash is taken after parsing, from the file on disk, not from the bytes that were parsed | `test_F17_...` |
| N1–N6 | Low–Med | various | findings without tests (see the end of this document) | – |

---

## Findings

### F1 (Critical) — H1 and H2 are inactive in the wired pipeline, and the evidence says "checked"

**Evidence.** `engine/stages.py` `S4Electronic.decide` calls

```python
decide_s4(res.payload, gap_max=hard["H1_metallic"], moment_tol=hard["H2_nonmagnetic"], ...)
```

where `hard = tgt.get("hard")` is the raw TOML. So `gap_max = {"gap_eV_max": 0.0, "source": ...}`
and `moment_tol = {"moment_muB_tol": 0.05, ...}`.

`electronic._th` handles a dict with `x.get("value")`, which returns **None**. In `decide_s4`:
- the H1 check `g is not None and ...` is skipped;
- the H2 check is skipped;
- the tags `H1_NOT_CHECKED` and `H2_NOT_CHECKED` are added only when `res.gap_eV is None` or the
  filter is disabled, **not** when the threshold is missing.

So an insulating ferromagnet (gap 2 eV, 3 μB) gets `PASS` with tags `['H7_HEUR_NOT_CHECKED']`.
`report/evidence.py` then prints "H1 metallic (S4): checked (gap_eV = 2; …)" and the same for
H2. Doc 22 §8.7 warned about these tables, but the engine was wired to them unchanged.

**Effect.** A passed check is claimed where none was made (contract R3, §3.6 "an absent
computation can never read as a passed check"). Insulators and magnets travel on to S5/S6, which
also wastes their core-hours.

**Required fix.**
1. `_th` must raise when given a dict without `value` (or accept the named keys
   `gap_eV_max` / `moment_muB_tol`).
2. `decide_s4` must add `H1_NOT_CHECKED` / `H2_NOT_CHECKED` whenever the threshold is None.
3. The engine (B1) must pass `Threshold(hard.H1_metallic.gap_eV_max, source)` etc.
4. The cleanest repair is to make targets.toml use `value =` for these two entries so that
   `targets.threshold("hard.H1_metallic")` works.

### F2 (High) — H1 is implemented as "a band crosses E_Fermi on the mesh", not "N(E_F) > 0"

**Evidence.**
- Doc 02 H1 reads "N(E_F) > 0 at the PBE level". It has **no** clause on distance to a gap
  edge, and it states that E_F near a gap is only flagged, Mg₂IrH₆ being one such case.
  Doc 02 S7 expects every ambient candidate to be a doped insulator.
- `analyse_bands` sets `gap = 0` only if some band has `min_k ε < E_Fermi < max_k ε`.
  Otherwise `gap = CBM − VBM`.
- With the profile smearing (Marzari–Vanderbilt 0.02 Ry ≈ 0.27 eV), a lightly electron-doped
  insulator has its E_Fermi a few tens of meV **below** the sampled CB minimum. The test puts
  E_Fermi 0.05 eV below the CB minimum. N(E_F) comes out at 1.37 states/eV/atom, yet
  `gap_eV` = 0.55 eV, so the outcome is `EXCLUDED_ASSUMPTION(H1)`.
- A coarse nscf mesh that misses a small pocket gives the same result.

**Effect.** Priority family 1 (metallic A₂MH₆ / AMH₃, doc 02 §5(b)) is the family most exposed.
EXCLUDED_ASSUMPTION is terminal.

**Required fix.**
- Decide H1 on electron counting: N_el is odd per spin-degenerate band, or the occupation of
  bands above the gap is ≥ 1e-3 e/cell. Alternatively use N(E_F) with a sourced lower
  threshold.
- When E_Fermi lies within a few smearing widths of a band edge, give `DEFERRED(H1-marginal)`.
  Doc 02 asks for a hybrid check there.
- Record `nelec` and the occupations.

### F3 (Medium) — marginal gaps are EXCLUDED, not deferred

**Evidence.**
- `gap_eV_max = 0.0` and the comparison is `res.gap_eV > g`. A 10 meV gap gives EXCLUDED, even
  though the code itself adds `GAP_SMALL_MESH_LIMITED` for gaps under 0.05 eV.
- Doc 02 H1: "check with a hybrid functional if the gap is marginal".
- PBE also underestimates gaps, so a small PBE gap could move either way under a hybrid, and a
  small gap is also where mesh errors live.

**Required fix.** A sourced `marginal_gap_eV` band (for example the mesh's energy resolution plus
the smearing width) inside which the outcome is `DEFERRED(H1-marginal)`. EXCLUDED only above it.
The threshold should be compared with its uncertainty, not with an exact 0.

### F4 (High) — H2 excludes on a metastable ferromagnetic solution

**Evidence.**
- Doc 02 H2 says "non-magnetic **ground state**".
- The spin step is a single SCF from `starting_magnetization = 0.5` on every species. Such a
  start can converge to a metastable moment-bearing solution.
- `ElectronicResult.e_total_spin_Ry` is parsed but never used. `decide_s4` excludes whenever
  |M| > tol, even when the FM energy is above the NM energy (the test puts it 0.136 eV/cell
  higher).
- Contrary to §3.1 ("spin-polarised **relax**"), the code runs a spin-polarised scf on the
  NM-relaxed cell. That is a deviation in its own right.

**Required fix.**
- Exclude only if E_spin < E_NM − ε_E (sourced, e.g. 1 meV/atom) **and** |M| > tol.
- If the moment is above tol but E_spin ≥ E_NM, record the result as NM ground state with a
  flag.
- Note the scf-vs-relax deviation in doc 22, or implement the relax.

### F5 (Medium) — the H2 tolerance is per cell, unsourced and provisional

**Evidence.**
- `moment_muB` is the |total magnetization| per cell.
- `moment_muB_tol = 0.05` comes from `doc02#5a:H2`, but doc 02 gives no number ("moment ≈ 0").
  The value is marked provisional.
- Applied per cell, the tolerance shrinks per atom as the cell grows. A 20-atom cell with
  0.003 μB/atom (smearing-level residual) is excluded.

**Required fix.**
- Define the tolerance per atom or per magnetic site, with a written source. Doc 02 or doc 11
  should give the number and its rationale (for example the smearing-induced residual measured
  on known NM controls).
- Test it on Pd-containing controls: Pd is a near-ferromagnetic paramagnet in the S10 palette.
- Below about 2× tol, give `DEFERRED(H2-marginal)`.

### F6 (Medium) — ferrimagnets and antiferromagnets pass H2 with "checked"

**Evidence.**
- `ANTIFERRO_OR_FERRI_SUSPECTED` is raised when |M_tot| < 0.05 < |M|_abs, but `decide_s4` only
  copies the flag into the tags.
- The outcome is PASS, and `evidence.py` prints "H2 non-magnetic: checked (moment_muB = 0)".
- The isotropic singlet Eliashberg machinery (and so every certified bound downstream) assumes a
  non-magnetic normal state.

**Effect.** This is a stronger claim than the data support.

**Required fix.** If |M|_abs exceeds the per-site tolerance, then:
- give `EXCLUDED_ASSUMPTION(H2)`, or at least `DEFERRED(H2-afm)`;
- the evidence statement must not print "checked" while that flag is present.

### F7 (Low) — the provisional H2 tolerance is not recorded on PASS

**Evidence.** When `moment_tol.provisional` is True, a PASS carries neither `provisional = True`
nor a `PROVISIONAL_*` tag. As a result, `evidence.py` "provisional margins used upstream" and
`summary.provisional_outcomes` both omit it.

**Required fix.** Tag every S4 outcome decided with a provisional threshold, PASS included. S5
already does this for M_c.

### F8 (Low; dormant) — H7-heur removes on runs whose cutoff is not converged

**Evidence.**
- The lower-bound argument (formation energy vs elements ≤ E_hull) is correct. But it assumes
  converged total energies at the same profile.
- `RELAX_PRESSURE_MISMATCH` flags exactly the opposite case (doc 22 §2: SG15 semicore files at
  60 Ry). `decide_s4` still removes there.
- Doc 02 H7 is also ZPE-corrected, while this check is not.
- Today the engine never passes `e_form_eV`, so the path is dormant.

**Required fix.**
- With `RELAX_PRESSURE_MISMATCH`, give `DEFERRED(H7-heur-unconverged)`.
- Add a sourced ZPE allowance to `ehull_heur_max_eV` before the path is wired.

### F9 (High) — imaginary coarse modes drive SCREENED_OUT(coarse-A) indirectly

**Evidence.**
- R2 says "coarse imaginary modes never remove a candidate".
- The core `atoms_to_spectrum` **drops** every α²F atom at ω ≤ 0 and records
  `A-DROP-NONPOS (... lowers A, anti-conservative for REJECTED)`.
- `collect_qe_eph` itself warns that "matdyn's a2F is unreliable/NaN/zero for dynamically
  unstable structures".
- `decide_coarse` reads none of this: no `NEGATIVE_FREQ_A2F_WEIGHT` or `IMAGINARY_MODES` flag,
  and no assumptions. It applies the margin and only adds `NEEDS_PHONON_CHECK`.
- A soft mode near instability is exactly where λ is largest (doc 02 §4.2). Its per-mode A
  contribution ∝ η/(Mω), which is large when the anharmonic frequency is small but real.
- Parametrised test: min ω −500 cm⁻¹, or the `IMAGINARY_MODES` flag, or the
  `NEGATIVE_FREQ_A2F_WEIGHT` flag, each with A = 400 K. Every case gives `SCREENED_OUT`.

**Required fix.** If `needs_phonon_check`, any `IMAGINARY_MODES` / `NEGATIVE_FREQ_A2F_WEIGHT`
flag, or any spectrum assumption `A-DROP-NONPOS` / `A-NEG-CLIP` is present, `decide_coarse` must
return `DEFERRED(imaginary_modes)`, never `SCREENED_OUT`.

### F10 (High) — multiband candidates are screened on the isotropic A

**Evidence.**
- Doc 01 §5.1 shows that for n_sheets > 1 the band-averaged A_iso falls by the DOS fraction w
  **without limit**. Its decoupled two-band example has A_iso = 806 K for a material with
  Tc = 300 K.
- Only κ, max_i A_i or A_iso/w_min are valid N1 inputs there. Doc 01 G-VALID item 5 makes such
  materials UNDECIDED without N1-MB data.
- M_c = 2 covers w ≥ 0.5 at most.
- `decide_coarse` has no multiband input; S4 `n_sheets` is never consulted. The doc 02 5(c)
  channel (MgB₂-like σ-hole borides) is multiband by design, and its control (MgB₂) is blocked,
  so control recall cannot catch the problem.

**Required fix.** Pass S4 `n_sheets` and the sheet DOS fractions into `decide_coarse`.
- For n_sheets > 1, screen on A_iso / w_min, using the minimum DOS fraction over crossing bands
  from the S4 nscf.
- If w is unknown, give `DEFERRED(multiband)`.
- Doc 12 §3.1 S5 should state this; it is a spec gap as well as a code gap.

### F11 (Low) — a missing ε_sys is treated as 0

**Evidence.** `one_eps = 1 + Fraction(str(es if es is not None else 0))`. In the test, A = 600 K
is SCREENED_OUT with ε_sys None but PASS with 0.10.

**Required fix.** Raise, as is already done for M_c, or return ERROR.

### F12 (High) — DEFERRED at S3 is a silent permanent removal

**Evidence.**
- `run_s3` processes only candidates with **no** S3 outcome.
- The scheduler (`engine/scheduler.next_stage`) advances only PASS-like outcomes. No code
  anywhere selects DEFERRED rows for re-processing.
- Therefore all of these never return:
  - `DEFERRED(ml_failed)`, including a whole-batch CUDA OOM, which defers the entire batch;
  - `DEFERRED(rank_cut)`, which is per `run_s3` call and never re-ranked against later days;
  - `DEFERRED(relax_not_converged)`.
- The rank cut orders by `score_s3`, which includes the MatterSim-phonon ω₂ term (weight 0.25,
  computed from a DOS with imaginary modes). Through rank_cut, MatterSim phonons therefore *do*
  remove candidates permanently, against R2 and doc 10 §3.4.
- Doc 12 §3.5's audit A3→4 samples some deferred candidates, but it is a measurement, not a
  re-entry path.

**Required fix.**
- A re-entry mechanism: `run_s3` and `next_stage` should treat DEFERRED as "not reached" when
  its reason is retryable (`ml_failed`, `relax_not_converged` with more steps, `budget`).
- A global rank queue for `rank_cut` (re-ranked daily with the new scores) instead of a per-call
  cut.
- Bounded retries counted in `Outcome.retry`.
- A summary line counting DEFERRED by reason and age.

### F13 (High) — ERROR is terminal (engine; reference only)

**Evidence.**
- `next_stage` returns None once a stage has an ERROR outcome, and `Outcome.retry` is never
  incremented anywhere.
- QE failures are therefore permanent, for example a spin scf that does not converge for a
  frustrated system, or `run_profile` diagnostics.
- At S5, `spectra_from_run` and `collect_run` errors (NaN α²F from matdyn on unstable
  structures, or a failed smearing admission) give `col.ok = False` → ERROR → terminal. This is
  the F9 imaginary-mode route ending as a removal again.

**Required fix (B1).**
- Retry policy: ERROR is re-queued with `retry + 1` up to a sourced limit, then becomes
  `DEFERRED(error_limit)` and is listed in S9.
- S5 collection failures caused by imaginary modes should be `DEFERRED(imaginary_modes)`.

### F14 (Medium) — control recall counts DEFERRED and ERROR controls as passed

**Evidence.** `control_recall` defines REMOVING without `DEFERRED` or `ERROR`. In the test, all
six controls are DEFERRED(ml_failed) at S3, and `control_recall()["S3"]["passed"] == 6`.

**Effect.**
- The S9 recall line reads 6/6 when none passed.
- R8 (auto-disable below recall 1) can never fire for controls that a screen stalls.
- Given F12/F13, a stalled control *is* removed.

**Required fix.** Count only PASS-like outcomes as passed. Report DEFERRED and ERROR separately
and treat them as failures of recall.

### F15 (High, threshold arithmetic) — the S3 H7-ml margin is smaller than the measured ML error

**Evidence.**
- `tiers.toml S3.ehull_self_max_eV = 0.42`, sourced "0.35 + 0.07". Units are consistent:
  eV/atom on both sides, and `ehull_self` returns pymatgen eV/atom.
- The 0.07 comes from doc 10 §3.3's "treat |E_hull| < ~0.07 as indistinguishable from 0". That
  is a resolution near zero, **not** an error bound.
- The same section measured per-entry MAEs of 0.074 eV/atom for H-containing H-Ir-Mg entries,
  0.064 overall, and a systematic +0.10 eV/atom offset on MP-sourced Mg entries, with MgH₂
  scattering up to +0.17. A hull distance is a difference of two such errors (≈ √2 σ).
- Doc 02 H7 is ZPE-corrected, while the ML hull is not. ZPE differences are "tens of meV/atom"
  (doc 02 §4.1).
- So a candidate with true ZPE-corrected E_hull 0.30 eV/atom (Li₂AgH₆ sits at 0.319) can
  plausibly read above 0.42 on the ML hull.
- The test requires the margin to be ≥ 2√2 × 0.074 ≈ 0.21 eV/atom, which puts the threshold at
  ≈ 0.56. It would fail at any margin below one MAE.
- One mitigation exists: the "MatterSim-only hull" partly cancels the systematic MP2020-mixing
  errors. But no tail statistic was ever measured for it.

**Required fix.**
- Derive the margin from a measured high quantile (e.g. 95th percentile) of
  |E_hull^ML − E_hull^DFT| on H-containing structures. The A3→4 audit pairs are the data
  source.
- Add a sourced ZPE allowance.
- Until then, use ≥ 0.35 + 0.2 and mark it provisional. The current entry is not tagged
  provisional, although it should be.

### F16 (Medium) — dedupe on unrelaxed cells merges physically different polymorphs

**Evidence.**
- S2 dedupes with `StructureMatcher(ltol 0.2, stol 0.3, angle 5°, scale=True)` on the
  **unrelaxed** generator output, before S3 relaxes anything.
- With `scale=True`, stol is relative to (V/N)^(1/3) (≈ 2.4 Å in the test). Mg + H₂ with
  H–H = 0.74 Å (molecular, penalised by doc 02 S4) and H–H = 1.60 Å (dissociated, a different
  bonding class) match.
- The later one becomes `DUPLICATE`, a terminal outcome, and its fate follows the first one.
- Doc 10 §3.4 also reports two distinct Fm-3m Mg₂IrH₆ minima (Ir–H 1.64 vs 1.73 Å). These would
  be merged too.
- A duplicate is never re-examined if its representative is later SCREENED_OUT, ERROR or
  DEFERRED (F12).

**Controls.** Controls are never DUPLICATE (good). But a generated cell matching a control
becomes DUPLICATE of that control, which is only evaluated at its surrogate T\*_s (80 K for
Mg₂IrH₆). That structure therefore never gets a 300 K decision.

**Required fix.**
- Dedupe only after S3 relaxation, on relaxed structures, with tighter tolerances (stol ≈ 0.1,
  ltol ≈ 0.1, `scale=False` within the same P).
- Or keep the S2 dedupe as a *tentative* link and re-check after relaxation.
- Re-activate a duplicate when its representative is removed by anything other than a decision
  that applies equally to both.
- A duplicate of a control must be evaluated at the real T\*.

### F17 (Low) — the spectrum hash is not bound to the parsed bytes

**Evidence.**
- `collect_run` first calls `wf.collect()`, which parses `a2F.dos*`. Only then does it hash the
  files on disk (`run_artifact_paths`, `register`).
- A rewrite in between attaches a `source_sha256` that is not that of the parsed spectrum. Two
  ways this can happen: a second matdyn in the same workdir, or a re-run that is adopted.
- R4 relies on this hash.
- The source lookup is also by file *name* only (`by_name[src.name]`). This is safe today only
  because every run has its own workdir.

**Required fix.**
- Hash before parsing and re-hash after. Mismatch → ERROR.
- Or parse from an in-memory copy whose hash is recorded.
- Match by resolved path, not by basename.

---

## Further findings (no test; lower confidence or outside B2/B3)

- **N1 (Medium).** The S4 decision has no hashed artifacts. `S4Electronic.run` registers no
  artifact: not `spin.out`, not the nscf XML, not `electronic_result.json`. So
  EXCLUDED_ASSUMPTION(H1/H2), a removal, cannot be re-verified from the ledger. R3/R4 should
  extend to S4 removals: register `scf.out`, `nscf.out`, `spin.out`, `data-file-schema.xml` and
  `electronic_result.json`.
- **N2 (Low).** Evidence ranking: `gather()` reads `S6.metrics.A_mid_K` and `A_lo_K`. The engine's
  S6 metrics hold only `ladder`, so A_mid is always None and the third ranking key (§3.6)
  silently does nothing.
- **N3 (Low).** The ZPE from `phonon.dos` drops imaginary modes and renormalises the positive
  part to 3N modes. For a reference phase with soft modes this biases its ZPE, and with it the
  S7 H7 removal, in an unknown direction. The unstable-phase case should be flagged.
- **N4 (Low).** S1 `outside_mattergen_domain` is applied to lane D (database) and control
  structures too, where MatterGen's domain is irrelevant. It is harmless for current lanes. The
  rule should be gated to lanes G\*.
- **N5 (Low, design).** The coarse plateau (σ ≥ 0.03 Ry ≈ 0.41 eV) over-smooths sharp DOS
  features at E_F. Those features are exactly the van Hove / high-symmetry cases doc 02 S6
  favours, and there A at large σ can be biased low. The all-σ rule and M_c = 2 partly cover
  this, but M_c is calibrated only on 1–4-atom elemental and H3S cases (doc 12 §3.4). The A5→6
  audit must include a vHs-at-E_F case.
- **N6 (Low).** The control lane passes `T*_s` through `targets_for`. The S5 `evidence` step
  uses `ctx.targets_for(rec)` correctly. `S5Coarse.decide` receives `tgt` from the engine, which
  was not checked here (engine code); confirm that the surrogate floor reaches `decide_coarse`
  for controls.

## What was checked and found sound

- **S3 decision.** `decide_s3` reads only `e_hull_ml_self_eV`, `error` and `relax.converged`.
  Nothing else in S3 can remove (except via the rank cut, F12): not MatterSim phonons, not
  ALIGNN, not the DFT-reference hull, not features or novelty. An unconverged relax never
  removes; a missing hull or P > 1 GPa gives no removal.
- **S3 hull direction.** Skipping reference phases (size cap, unconverged, failed) can only lower
  E_hull, the lenient direction. Units are eV/atom throughout.
- **S5 margin rule.**
  - `decide_coarse` uses exact `Fraction` comparisons and the all-σ rule.
  - It adds `PROVISIONAL_MARGIN` and `PROVISIONAL_EPS_SYS`.
  - NO_PLATEAU gives PASS, and M_c is never passed to the core engine.
  - Units are consistent: A in K (core `to_kelvin_factor("Ry")`), σ in Ry from the `lambda`
    file, floor = a0·T\* in K.
  - ZPE = 1.5⟨ħω⟩ per atom is dimensionally correct (meV/atom).
- **S4 H7-heur bound.** Formation energy vs elemental references is a valid lower bound on
  E_hull; a non-ground-state elemental reference only makes it lower. `refhull` with missing
  phases is lenient, as stated.
- **`prune.py`.**
  - Every deletion is preceded by a sha256 manifest line, fsync'ed to every manifest.
  - `*.save`, dyn files, `elph_dir` and `a2Fsave` are never touched by per-q or wfc pruning.
  - A q is pruned only after the *next* dyn file exists. In ph.x the `elphon` call of q N
    finishes before q N+1's dynmat is written, so the wavefunctions are no longer needed.
  - No evidence-loss path was found.
- **Controls** cannot be DUPLICATE themselves (but see F16).

## Files

- Review: `docs/reviews/RC1_code_review.md` (this file).
- Tests: `tests/review/conftest.py`, `tests/review/test_rc1_findings.py`. 19 test cases are
  expected to fail until the fixes land.
