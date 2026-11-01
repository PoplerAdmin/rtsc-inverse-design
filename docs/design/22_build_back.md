# 22. Build notes: back end (`rtt.dft`, `rtt.score.s4`, `rtt.report`)

Status: built and tested, 2026-10-29. Builder: agent B3-back. Spec: `12_architecture.md` rev 2
(§3.1 S4–S9, §3.2, §3.4, §3.6, §5.1, §6, §7.2, §8). Core interfaces (`rtt.records`,
`rtt.config`, `rtt.ledger`) are B1-core's (`20_build_core.md`); this code imports them and the
the core DFPT code from `rtt.core.dfpt` (no second copy).

Labels: **(meas.)** = measured here, with the script; **UNVERIFIED** = not checked.

---

## 1. What exists

| Module | Purpose | Main API |
|---|---|---|
| `dft/profiles.py` | the four §7.2 profiles; SG15-NC-only pseudopotential table; tiers.toml overrides | `PROFILES`, `Profile.eph_settings(q, pressure_GPa, do_relax)`, `Profile.with_overrides(**kw)`, `profiles_from_tiers(tiers)`, `sg15_pseudos(elements)`, `available_elements()` |
| `dft/workflow.py` | profile runner around the **unchanged** core `EPhWorkflow`; generic `Chain` runner for pw.x/pp.x chains | `run_profile(struct, profile, workdir, launcher, *, material_id, q, pressure_GPa, do_relax, policy) -> DFPTRun`, `run_full(...) -> list[DFPTRun]`, `default_launcher(ranks)`, `RunPolicy`, `wait_for_ranks` |
| `dft/prune.py` | hash-then-delete pruning, manifests, verification, disk guard | `prune_q_scratch(workdir, q_done=None, with_dvscf=False) -> int`, `prune_final(workdir, mode)`, `verify_manifest(manifest, base)`, `inventory(root)`, `disk_guard(path, min_free_gb)`, `ScratchWatcher` |
| `dft/electronic.py` | S4 chain and derived quantities; S4 decision | `run_electronic(struct, workdir, *, profile, launcher, pressure_GPa, policy) -> ElectronicResult`, `analyse_bands`, `decide_s4(res, gap_max=, moment_tol=, ehull_heur_max=, e_form_eV=) -> Outcome`, `formation_energy_per_atom_eV` |
| `dft/belli.py` | Belli φ, H_DOS, Tc_φ (port of `src/scoring/phi_elf.py`) | `phi_and_hdos(workdir) -> BelliResult`, `networking_value`, `h_fraction_voronoi`, `tc_phi` |
| `dft/refhull.py` | own-DFT reference energies per profile (cached), hull with optional ZPE | `ref_hull(chemsys, profile, phases, workroot, *, cache_dir) -> RefHull`, `RefHull.e_above_hull(comp, e_per_atom_eV, zpe_meV_atom, with_zpe)`, `phase_energy`, `h2_molecule`, `zpe` |
| `dft/collect.py` | Spectrum objects through core parsers + `admit_smearings`; sha256 `ArtifactRef`s (ledger-registered when a ledger is given) | `spectra_from_run(run, *, ledger, tier_run_id) -> (list[Spectrum], list[ArtifactRef])`, `collect_run(...) -> Collected`, `spectrum_digest(sp)` |
| `dft/coarse.py` | S5 metrics and the coarse margin rule | `coarse_metrics(collected, run) -> CoarseMetrics`, `decide_coarse(metrics, *, floor_K, margin_M_c, eps_sys) -> Outcome`, `zpe_meV_per_atom` |
| `score/s4.py` | S4 ranking score (ranking only) | `score_s4(electronic, phi=None, *, T_star_K) -> float`, `components(...)` |
| `report/evidence.py` | §3.6 evidence statement | `evidence_statement(cid, ledger) -> str`, `gather(ledger, cid)`, `render(EvidenceInputs)`, `items(...)` |
| `report/summary.py` | daily summary (R5, R8) | `daily_summary(ledger) -> str`, `control_recall(ledger)`, `provisional_outcomes(ledger)` |

`dft/coarse.py` is an addition to the §5.1 list. The S5 margin rule had to live somewhere I
own; `rtt.decide` (B1) stays the only producer of decide outcomes.

Tests: `tests/back/` (42 unit tests, about 3 s, plus one opt-in real QE smoke test).

## 2. Profiles (§7.2)

All profiles use SG15 ONCV PBE-1.2 files only (`<El>_ONCV_PBE-1.2.upf`). `sg15_pseudos` refuses
any file whose UPF header is not `pseudo_type="NC"` and any fully relativistic file (W4: QE 7.6
e-ph segfaults with US/PAW). ecutrho = 4·ecutwfc.

| Profile | kind | ecut (Ry) | k (1/Å incl. 2π) | q | e-ph σ ladder | prune | Source |
|---|---|---|---|---|---|---|---|
| `s4` | electronic | 60 | scf 0.25, nscf 0.15 | – | – | wfc at end | doc 11 §3.3 (Belli chain at 60 Ry) |
| `coarse` | dfpt, tier `coarse` | 60 | k_scf 0.25, k_fit 0.125 | tiers `S5.q` (4³) | 0.01 Ry × 10, mode-positive, plateau σ ≥ 0.03 | per-q wfc + dvscf; wfc at end | doc 11 §3.1 (H3S: gives exactly k 12³ / k_fit 24³, tested) |
| `full` | dfpt, tier `full` | 90 | 0.17, fit factor 2 | tiers `S6.q_meshes` | 0.005 Ry × 10, mode-positive | per-q wfc only (dvscf kept for EPW) | core defaults; H3S reference |
| `ref_hull` | hull | 90 | 0.20 | – | – | wfc at end | §3.1 S7a |

Override any field from `tiers.toml` with a table `[profiles.<name>]`. Overrides are recorded
in `Profile.overrides` and change `Profile.sha256()`, which keys the reference-energy cache.

**Finding (meas.):** SG15 `Al_ONCV_PBE-1.2` has `z_valence = 11` (2s2p semicore). At 25 Ry the
pressure at the experimental lattice is −6799 kbar. It is −374 kbar at 40 Ry and −8 kbar at
60 Ry. A vc-relax at 25 Ry shrank the cell to 9.4 Å³/atom, compared with 16.5 Å³/atom
experimentally. The 60 Ry coarse/s4 cutoff is adequate for Al, but other SG15 files with
semicore states may need more. Two guards now catch this:
- S4 runs a fresh-basis scf on the relaxed cell and flags
  `RELAX_PRESSURE_MISMATCH(P=…)` when |P − P_target| > max(10 kbar, 5 %).
- `coarse_metrics` does the same check on the DFPT scf.

## 3. Runner behaviour (`run_profile`)

- **Steps.** The core `EPhWorkflow` is driven one step at a time with
  `wf.run(until=step, collect=False)`. The core token chain (`dfpt_state.json`) and its
  adopt-if-identical-input rule stay in force. Changing a step's input reruns it and every step
  after it.
- **Whole-run adoption.** After success, `profile_done.json` stores the matdyn token and the
  sha256 of every output. A later call with the same inputs adopts the run when every hash still
  matches. This works even after `tmp/` has been pruned, which plain core adoption cannot do
  because it needs `tmp/pwscf.a2Fsave`. Tested on real QE: the second call re-runs nothing.
- **Timing.** One row per step goes to `timing.jsonl`, with wall time (monotonic clock),
  children user+sys CPU, core-h, ranks, time spent waiting for ranks, time paused, and bytes
  pruned. The core `RunRecord.elapsed_s` uses `time.time()` and came out **negative**
  (−0.10 s) under WSL clock adjustment (meas.). The wrapper therefore measures the step itself.
- **Rank budget.** Before each parallel step, `wait_for_ranks` polls `pgrep -c -x` over the QE
  executables until running + needed ≤ `max_machine_ranks` (default 12), because other agents
  share the machine.
- **Disk.** `disk_guard` runs before each step (`machine.toml [disk] min_free_GB`). During ph.x,
  `ScratchWatcher` does two things every `watcher_interval_s`:
  - It prunes finished q points. A q point counts as finished when its dyn file is complete and
    the next dyn file exists. Selection: the wavefunctions, plus dvscf when
    `prune_dvscf_per_q` is set.
  - Below the free-space floor it SIGSTOPs, and above the floor + 10 GB it SIGCONTs, **only
    descendants of this Python process**. Other agents' QE processes are never signalled.
- **Manifests.** Every pruned file is hashed and its line fsync'ed to `<run>/PRUNED_MANIFEST.tsv`
  (plus a central copy if `RTT_PRUNE_MANIFEST_DIR` is set) before deletion. Paths are
  `~/…`-relative and the machine label is neutral. `verify_manifest` re-hashes restored or
  regenerated files (MATCH / MISMATCH / MISSING): this is the outside-verification path.
- **Failure.** `run_profile` never raises on QE failure. `DFPTRun.status = "failed"` and the
  core diagnosis is in `error`; the stage maps this to `ERROR`.
- **S6.** `run_full` runs `full` on each `q_meshes` entry in `q<n>x<n>x<n>/`. The first mesh
  re-relaxes at the full cutoff and k (R12 #14), and the later meshes reuse that cell. The
  force/stress-check alternative in §3.1 is **not** implemented.

## 4. S4 electronic (`run_electronic`)

Chain: vc-relax → scf → ELF (`pp.x` 8) → nscf → ILDOS (`pp.x` 10, only with H) → spin-polarised
scf (nspin 2, starting_magnetization 0.5 on every species, own prefix).

- **n_sheets:** the number of bands (per spin, summed) with min_k ε < E_Fermi < max_k ε on the
  nscf mesh. A band with several pockets counts once, so this is a lower bound on the number of
  sheets.
- **E_F (G-VALID item 1, literal):** over the crossing bands, min(E_Fermi − bottom, width/2).
  `crossing_bands` records the depth, height and width of each crossing band, and
  `E_F_hole_eV` adds the band-top distance for information.
- **gap:** 0 if any band crosses. Otherwise CBM − VBM on the mesh, flagged
  `GAP_SMALL_MESH_LIMITED` when below 0.05 eV.
- **N(E_F):** Gaussian with σ = 0.136 eV, both spins, per eV per atom. It is also evaluated at
  2σ and flagged `NEF_MESH_SENSITIVE` when the two differ by more than 25 %.
- **moment:** |total magnetization| of the spin scf. Flagged
  `ANTIFERRO_OR_FERRI_SUSPECTED` when the total is small but the absolute magnetization is not.

`decide_s4` returns an `rtt.records.Outcome`. The order is ERROR → H1
(`EXCLUDED_ASSUMPTION`) → H2 (`EXCLUDED_ASSUMPTION`) → H7-heur (`SCREENED_OUT`) → PASS. Each
filter can be switched off (R3). Any check not made adds a tag to the PASS outcome:
`H1_NOT_CHECKED`, `H2_NOT_CHECKED` or `H7_HEUR_NOT_CHECKED`.

The heuristic hull is the formation energy against elemental references at the *same* profile
(`refhull.ref_hull(..., profile=s4 profile)`). This is a lower bound on E_hull, because the
hull lies at or below the elemental tie line, so removing on it is conservative. No X-HULL is
applied at S4 (R12 #8).

Belli φ is our own percolation code (6-connectivity, bisection), not Critic2. H_DOS is the
share of the ILDOS in H Voronoi cells, because SG15 has no PSWFC for projwfc. Both deviations
from the paper are UNVERIFIED against its table (as stated in doc 11 §3.3). For a structure
without H, H_f = 0 and Tc_φ = None.

`score_s4` is for ranking only. Weights for H compounds: 0.6 · Tc_φ/T\* + 0.25 · N/(N + 0.3) +
0.15 · H_DOS window score. Without H, the N(E_F) term alone. These weights are conventions to
re-fit against the A4→5 audit.

## 5. S5 coarse (`collect.py`, `coarse.py`)

- `collect_run` rebuilds the core `EPhWorkflow` from the `DFPTRun` and calls its `collect()`
  (`collect_qe_eph` with `smearing_admission="mode_positive"`, `admit_smearings`,
  `qe_a2f_dos_to_spectrum`), so the spectra and their provenance are exactly the core ones.
- Each spectrum gains `source_sha256` (its `a2F.dos` file), `source_artifact_id` (with a ledger),
  `profile`, `spectrum_tier`, `q` and `spectrum_digest`.
- Every input, output, dyn, matdyn, state, timing and manifest file becomes an `ArtifactRef`.
- `spectra_from_run` raises `CollectError` when admission fails or any requested smearing is
  unusable (the core rule: no verdict from a subset).
- `coarse_metrics` reports, per smearing: σ, admitted, A, λ, ω_log, ω₂. Excluded smearings are
  reported too, for information. It also reports:
  - the plateau (admitted σ ≥ `plateau_sigma_min_Ry`);
  - the minimum frequency on the matdyn mesh and per irreducible q;
  - the ZPE from `phonon.dos`;
  - the scf pressure.
- `decide_coarse` emits **only** `SCREENED_OUT(coarse-A)`, `PASS` or `ERROR`. It removes iff
  A·(1 + ε_sys) < floor/M_c on **every** plateau σ, compared exactly in `Fraction`.
  - Tags: `PROVISIONAL_MARGIN` on every outcome made with a provisional M_c;
    `PROVISIONAL_EPS_SYS` likewise for ε_sys.
  - An empty plateau never removes (`NO_PLATEAU`).
  - Imaginary modes never remove; they add `NEEDS_PHONON_CHECK` (|ω| > 30 K).
  - M_c is never passed to the core engine.

## 6. Evidence statement (`report/evidence.py`)

All 22 items of the §3.6 template are always printed, each as `checked (value; record ids)` or
`NOT CHECKED (what is missing)`. The items are:
- structure, lane, novelty;
- H1, H2;
- S6 full-grid stability and the Γ diagnostic;
- G-VALID 1–5, each with its fields; a missing field reads "G-VALID fails";
- ENVELOPE inputs;
- decide outcome and rule;
- the core verdict with recheck;
- ZPE hull and synthesis route;
- anharmonic check;
- MD/PIMD;
- μ\* at 0.10/0.13/0.15; each value without a verdict is listed as NOT CHECKED;
- multiband status;
- provisional thresholds used upstream;
- audits and control recall.

The ranking keys follow §3.6: decide outcome, then T_low of the latest full/converged verdict at
μ\*_hi = 0.13 (0 if missing), then A_mid. The S7 metric names it reads are listed in the module
docstring; whoever builds S7 must write those keys. `daily_summary` gives:
- the current outcomes per stage;
- core REJECTED_PROVEN by spectrum tier;
- compute per stage/profile;
- the provisional outcomes;
- control recall per stage, flagging any recall below 1.

## 7. Real QE smoke test (meas.)

`RTT_QE_SMOKE=1 pytest tests/back/test_qe_smoke.py`. Setup: QE 7.6, 2 MPI ranks, nice 19, run
while other agents' ph.x used 4–10 ranks.

- **System:** fcc Al from a = 4.04 Å.
- **s4 profile** (ecut 60, k spacing 0.5, nscf 0.15 → 18³):
  - wall 15 s, 0.005 core-h for relax + scf + ELF + nscf + spin;
  - gap 0, moment 0, n_sheets 2 (bands 5 and 6 of the semicore-including set);
  - N(E_F) 0.52 /eV/atom (2σ: 0.44; free-electron value 0.39);
  - φ = 0.578, scf P = −1.0 kbar;
  - S4 outcome PASS.
- **coarse profile** (ecut 60, k 6³ / k_fit 12³, q 2³, matdyn 8³) on the S4-relaxed cell:
  - wall 39 s, 0.016 core-h; ph.x 23 s;
  - per-q pruning removed 10 MB during ph.x, with 7 manifest lines;
  - λ = 1.08 / 0.68 / 0.46 / 0.37 / 0.35 at σ = 0.01–0.05 Ry, ω_ln ≈ 300–318 K. The
    small-σ blow-up is the coarse-k artefact doc 11 describes. At σ ≥ 0.03, λ = 0.35–0.46,
    close to Al's ≈ 0.4.
  - The artifacts were registered in an in-memory ledger.
  - `decide_coarse` at the real target gave `SCREENED_OUT(coarse-A)` with
    `PROVISIONAL_MARGIN`; at a T\* = 1 K surrogate it gave PASS.
  - A second `run_profile` / `run_electronic` call adopted everything and re-ran nothing.

Two earlier runs failed in useful ways:
- At 25 Ry the cell collapsed (§2).
- With nscf at 9³, N(E_F) came out at 0.72 /eV/atom and E_F(G-VALID) at 0.075 eV, from
  under-sampled band-6 pockets. That run motivated the 2σ check and the `crossing_bands` record.

## 8. Open problems and what I could not verify

1. **G-VALID item 1 E_F definition (theory owner).** "Bottom of the shallowest crossing band"
   makes E_F the depth of the smallest pocket. For Al this depends on the mesh: 0.075 eV at 9³
   and 0.31 eV at 18³. With λ ≈ 0.4 and ω_max ≈ 0.04 eV, λω_max/E_F is 0.2 at 9³, which fails
   the 0.1 gate, and 0.05 at 18³, which passes. A simple metal can therefore fail the
   non-adiabatic gate through a shallow pocket. Doc 01 should state whether E_F is per sheet
   (with N1-MB) or should exclude pockets carrying a small DOS share, and which mesh it needs.
   I implemented the literal definition and record every band's depth.
2. **N(E_F)** uses Gaussian sampling of nscf eigenvalues. Tetrahedra (`dos.x`) would be more
   accurate. The 2σ flag is a guard, not a fix.
3. **S6 force/stress-check alternative** to re-relaxing is not implemented (`run_full` always
   re-relaxes on the first mesh).
4. **ZPE from `phonon.dos`** is normalised to 3N modes from matdyn's DOS. It was not compared
   with a direct mode sum (UNVERIFIED accuracy at coarse q).
5. **`refhull`** is tested with hand-made energies and the cache key only. No real multi-phase
   QE hull has been run, and the H₂-in-a-box reference (10 Å box, Γ only, fixed cell) is
   untested on QE. The choice of reference phases is left to the caller (database lane D).
6. **Belli φ and H_DOS:** not validated against the paper's table (UNVERIFIED, as in doc 11).
7. **No `Stage` classes.** Building S4/S5 `Stage` objects with `run`/`decide` is the engine's
   job (B1). They should call:
   - S4: `run_electronic` → `decide_s4` → `score_s4`;
   - S5: `run_profile(coarse)` → `collect_run(ledger, tier_run_id)` → `coarse_metrics` →
     `decide_coarse`, plus `rtt.decide.decide_t(tier=coarseDFPT)` for the decide record.

   Thresholds come from `targets.toml` `hard.H1_metallic.gap_eV_max` and
   `hard.H2_nonmagnetic.moment_muB_tol` (inline tables without `value`, so pass
   `Threshold(value, source)` built from them), and from tiers `S4.ehull_heur_max_eV`,
   `S5.margin_M_c` and targets `envelope.eps_sys`.
8. **Dispatch:** `run_profile` dispatches `s4` to `run_electronic`, but `ref_hull` profiles must
   go through `refhull.ref_hull`.
9. The core `RunRecord.elapsed_s` can be negative under WSL (see §3). The core state
   files therefore carry wrong elapsed times; the wrapper's `timing.jsonl` is the record to use.

---

## 9. Response to RC1 (`docs/reviews/RC1_code_review.md`) and the B1 follow-ups

Status after the fixes: `tests/back` (54 passed + the QE smoke), `tests/review` and `tests/front`
all pass: 106 passed, 4 skipped (opt-in tests), CPU only. The real QE smoke was re-run with the
new chain (Al, 2 ranks): it passes in 66 s.

| Finding | Fix | Where |
|---|---|---|
| F1 (back half) | `th_parts` reads `value` **or** exactly one named key (`gap_eV_max`, `moment_muB_tol`, `moment_muB_per_atom_tol`, `window_eV`, …) from inline TOML tables. A table with neither raises `ValueError`, so a threshold can never become None silently. A threshold that is None, or a filter that is disabled, adds `H1_NOT_CHECKED` / `H2_NOT_CHECKED` / `H7_HEUR_NOT_CHECKED`. `evidence.py` prints NOT CHECKED for those tags and for any `DEFERRED(H1-…/H2-…)`. | `dft/profiles.py::th_parts`, `decide_s4`, `report/evidence.py` |
| F2 | H1 is decided on states at E_F (doc 02 "N(E_F) > 0"), not on band crossing alone. The S4 chain gained `nscf_tet` (`occupations='tetrahedra_opt'`) and `dos.x`. `N_EF_per_eV_atom` is now the **tetrahedron** DOS at E_F. The Gaussian σ and 2σ values are kept as `N_EF_gauss_*`. A candidate counts as metallic if any estimate exceeds `nef_min` (0.01 states/eV/atom, provisional), or if a band crosses. `nelec` is recorded. | `electronic.py` |
| F3 | A gap is EXCLUDED only above `gap_eV_max + marginal_gap` (0.3 eV ≈ the MV smearing 0.02 Ry, provisional). Inside that band the outcome is `DEFERRED(H1-marginal)`: a hybrid functional or denser mesh check is required. | `decide_s4` |
| F4 | H2 excludes only when the spin solution is lower than the NM one by more than `mag_energy_tol` (1 meV/atom, provisional). Same cell, k, cutoff and smearing; `magnetic_energy_meV_atom` is recorded. A higher-energy magnetic solution passes with `H2_METASTABLE_MAGNETIC_SOLUTION(dE, m)`. Missing energies give `DEFERRED(H2-energy-unknown)`. Not implemented: a spin-polarised *relax*; the spin scf runs on the NM-relaxed cell, as before. | `decide_s4` |
| F5 | The moment tolerance is **per atom**. The engine still passes targets `moment_muB_tol` (0.05, provisional). `decide_s4` reads it as μB/atom. Proposal to B1/theory: rename it `moment_muB_per_atom_tol` and give it a source (e.g. the residual measured on NM Pd-containing controls). Moments within 2× the tolerance give `DEFERRED(H2-marginal)`. | `decide_s4`, `ElectronicResult.moment_per_atom_muB` |
| F6 | The moment per atom is max(\|M_tot\|, \|M\|_abs)/nat, so ferrimagnetic and antiferromagnetic solutions go through the same energy test and are never a plain PASS. | `decide_s4` |
| F7 | Every S4 outcome decided with a provisional threshold, PASS included, gets `provisional=True` and a `PROVISIONAL_THRESHOLD(H1|H2|H7)` tag. | `decide_s4` |
| F8 | H7-heur gives `DEFERRED(H7-heur-unconverged)` instead of SCREENED_OUT when the run is flagged `RELAX_PRESSURE_MISMATCH`. An optional `zpe_allowance_eV` argument exists, but see (c) below: the allowance now sits in the tiers value. | `decide_s4` |
| F9 | `coarse_metrics` turns the core `A-DROP-NONPOS` assumption into the flag `A_DROP_NONPOS`. `A-NEG-CLIP` is recorded as `A_NEG_CLIP`; it raises A and so is not a removal hazard. If a removal would fire while any of `needs_phonon_check`, `IMAGINARY_MODES…`, `NEGATIVE_FREQ_A2F_WEIGHT` or `A_DROP_NONPOS` is present, the outcome is `DEFERRED(imaginary_modes)` with `NEEDS_PHONON_CHECK`. | `dft/coarse.py` |
| F10 | `decide_coarse(..., n_sheets, w_min)`: <br>• `n_sheets == 1` screens on A; <br>• `n_sheets > 1` with a known `w_min` screens on A/w_min (doc 01 §5.1); <br>• otherwise, including n_sheets unknown, a removal becomes `DEFERRED(multiband)`. <br>S4 now provides `sheet_w_min`: the smallest per-band DOS fraction at E_F, min over the σ and 2σ Gaussians. **Engine (B1):** `S5Coarse.decide` must pass S4 `n_sheets` and `sheet_w_min`. Until it does, every coarse removal is DEFERRED, which is conservative. | `coarse.py`, `electronic.py` |
| F11 | A missing ε_sys raises `ValueError`. | `decide_coarse` |
| F14 | `control_recall` had already been changed (PASS-like only). The summary table now also shows deferred and error counts and flags recall < 1 whenever passed < reached. | `report/summary.py` |
| F17 | `collect_run` hashes every file the parser reads (`a2F.dos*`, `lambda`, `freq.dos`) **before** parsing and again after. The spectra carry the pre-parse hash, and any change is a collection error. Sources are matched by resolved path, not basename. `source_artifact_id` is set only when the registered hash equals the parsed one. | `dft/collect.py` |
| N1 | `electronic_artifacts(res, ledger=, tier_run_id=)` registers the files an S4 decision rests on: scf/nscf/nscf_tet/spin in/out, `pwscf.dos`, `data-file-schema.xml`, the relaxed structure, state, timing and `electronic_result.json`. **Engine (B1):** call it in `S4Electronic.run` before `finish_tier_run`. | `electronic.py` |
| N2 | `coarse_metrics` now writes `A_lo_K` (min over admitted σ) and `omega_max_K` next to `A_mid_K`. For a full run, the §3.6 ranking keys are available once S6 stores these metrics. | `coarse.py` |

**(b) S6 full-grid stability (`min_freq_K`).** `coarse_metrics` works for any DFPT run, coarse
or full. It now writes:
- `min_freq_K`: signed, in K; negative means imaginary, the QE sign convention. It is the
  minimum over every irreducible q of the run's own DFPT mesh (the `dyn*` dynamical matrices of
  the run's final relaxed structure) and the matdyn interpolation mesh (uniform
  k_matdyn = `matdyn_factor`·q, read directly from `freq.dos`, so it survives a failed α²F
  collection).
- `max_imag_mode_K` = max(0, −min_freq_K).
- `min_freq_sources`: the number of irreducible-q dyn files used, and whether the matdyn mesh
  was used.

The mesh is uniform, not a high-symmetry path. For S6, B1 should pass the finest mesh's
`metrics["min_freq_K"]` into `decide.inputs` (which already expects `min_freq_K`).

**(c) `tiers.toml` `S4.ehull_heur_max_eV`: 0.42 → 0.56 (provisional).** I changed this one line
only; B1's other edits to tiers.toml are untouched. The value is sourced in the file as:
- 0.35 (H7 window);
- + 0.10 ZPE allowance (H7 is ZPE-corrected, the heuristic is not);
- + 0.10 allowance for the formation-energy error of the s4 profile (60 Ry, coarse k).

Both allowances are UNVERIFIED and should be measured against `ref_hull`. The F15 argument
about ML error does not apply here, because S4 uses DFT, and E_form vs elements is already a
lower bound on E_hull. The remaining risks are profile convergence and ZPE.

**Measured with the new chain (Al smoke, nscf 18³):**
- Tetrahedron N(E_F) = 0.397 /eV/atom, against 0.39 for a free-electron metal. The earlier Gaussian
  value was 0.517.
- In a run over the full energy range, the integrated DOS at E_F was 11.0 = nelec, so `dos.x`
  reports the spin-summed DOS.
- That full-range `dos.x` took **470 s** for Al because the SG15 semicore states sit near
  −92 eV. `dos.x` now runs on E_F ± 1 eV and takes 8.7 s.
- sheet_w_min = 0.23.
- E_spin − E_NM = +0.01 meV/atom, so NM is the ground state.
- With the S4 sheet data, the coarse rule for Al gives SCREENED_OUT or PASS. Without the S4
  sheet data, it gives DEFERRED(multiband).

**Still open:**
- The spin-polarised relax (F4) is not implemented.
- The provisional numbers have no measured basis yet: `nef_min`, `marginal_gap`,
  `mag_energy_tol`, the per-atom moment tolerance and the two H7 allowances.
- The per-band DOS fractions use Gaussian sampling of eigenvalues, not per-band tetrahedra.
- Engine wiring (B1): pass `n_sheets` / `sheet_w_min` into `decide_coarse`; call
  `electronic_artifacts`; feed `min_freq_K` into S6.

---

## 10. S5 rule replaced by doc 13 T2-SCREEN

Doc 13 §5 found that the factor-2 rule on max_σ A_coarse is not safe enough. S5 now uses T2-SCREEN.

**Rule** (`dft/coarse.py::decide_coarse`, `coarse_envelope`; a port of
`src/scoring/coarse_envelope.py`):
- **Two meshes** (`metrics` = q 4³, `metrics_coarser` = q 3³):
  SCREENED_OUT(`coarse-A envelope`) iff (1 + e_c)·A_hi,coarse < A_target.
  - A_target = `Targets.floor` (a0·T\*; the surrogate floor for controls).
  - A_hi,coarse = A_nom(1 + δ_σ + δ_q)(1 + ε_sys)/(1 − ε_ω), with:
    - X_nom = the value on q 4³ at the smallest σ admitted on **both** meshes;
    - δ_σ taken over the σ admitted on q 4³;
    - δ_q = |A₄ − A₃|/A_nom at that σ;
    - ε_ω from the q-spread of ω_log at that σ, floor 0.02, cap 0.9.
  - e_c = tiers `S5.e_c` = 0.5.
- **One mesh** (fallback): SCREENED_OUT(`coarse-A max`) iff M_c·(1 + ε_sys)·max over admitted σ
  of A < A_target, with M_c = 2.
  - Only at q ≥ 4³; otherwise `DEFERRED(single_mesh_q_below_4)`.
  - The (1 + ε_sys) factor is **stricter** than doc 13's 2·max A. It only makes removal harder.
    It is kept so that ε_sys stays a required input (RC1 F11).
- A, λ and ω_log are always integrated from the α²F atoms (core `Spectrum` methods), never
  read from QE's `lambda` footer, which is wrong when there are imaginary modes.

**Preconditions.** If any of these fails, a removal that would fire becomes
`DEFERRED(<first blocker>)`, and every blocker is tagged `BLOCKED:<name>`:
- no mesh below q 3³ (q 2³ is never used);
- mode-positive admission succeeded on every mesh (`CoarseMetrics.admission_ok`);
- at least 3 admitted σ on q 4³;
- a non-empty σ plateau;
- S4 `n_sheets == 1` (`w_min` is recorded only; the A/w_min variant of §9 F10 was dropped,
  following the coordinator's instruction, because doc 13 did not calibrate it);
- no imaginary modes, `A_DROP_NONPOS`, `IMAGINARY_MODES` or `NEGATIVE_FREQ_A2F_WEIGHT` on any
  mesh (also tagged `NEEDS_PHONON_CHECK`).

If the envelope cannot be formed (no common admitted σ), the outcome is PASS with
`ENVELOPE_INCOMPLETE`: no removal.

The rule never uses λ, ω_log, m or imaginary modes as a removal criterion. `PROVISIONAL_MARGIN`
stays on every outcome made with a provisional e_c or M_c. In tiers.toml I set `e_c`
`provisional = true`, because doc 13 calls it a judgment on n = 4 materials, and re-sourced
`margin_M_c` as the doc-13 fallback. B1's other `[S5]` edits are untouched (`q_meshes =
[[3,3,3],[4,4,4]]`, `e_c`).

**Runner.** The coarse profile now has `q_meshes = ((3,3,3),(4,4,4))`, and `profiles_from_tiers`
reads `S5.q_meshes`. `workflow.run_coarse` runs q 3³ first, with the relax, then q 4³ on the
same relaxed cell, so δ_q and ε_ω compare one structure.

**Engine (B1).** `S5Coarse` must:
- call `run_coarse` instead of `run_profile`;
- collect and build `coarse_metrics` per mesh;
- call `decide_coarse(m_q4, metrics_coarser=m_q3, e_c=th.threshold("S5.e_c"),
  eps_sys=tgt.threshold("envelope.eps_sys"), n_sheets=<S4 n_sheets>, floor_K=tgt.floor)`.

With today's single-mesh wiring, the fallback applies at q 4³ and every removal without S4
`n_sheets` is DEFERRED, so the current wiring is conservative.

Tests: 11 new unit tests in `tests/back/test_collect_coarse.py` cover:
- the envelope formula;
- the two-mesh screen;
- each precondition blocker;
- the fallback restricted to q ≥ 4³;
- the required inputs.

The real QE smoke test now runs Al at q 3³ + q 4³.

**Measured on QE (Al smoke, 2 ranks, nice 19, 268 s in total).**
- ph.x took 35 s at q 3³ and 164 s at q 4³. The q 4³ run reused the q 3³-relaxed cell and did
  not relax again.
- T2-SCREEN at the real target, treating Al as single-sheet:
  - A_nom = 64.4 K, δ_σ = 0.030, δ_q = 0.020, ε_ω = 0.068;
  - (1 + e_c)·A_hi = 119.8 K < 1299.4 K, giving `SCREENED_OUT(coarse-A envelope)` with
    `PROVISIONAL_MARGIN`.
- With the real S4 `n_sheets = 2`, the outcome is `DEFERRED(multiband)`. At a T\* = 1 K
  surrogate it is PASS.
- `tests/back`, `tests/review` and `tests/front`: 114 passed, 4 skipped (opt-in tests).
