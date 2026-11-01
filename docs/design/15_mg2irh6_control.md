# 15 — Mg2IrH6 positive control: PBE DFPT vs MatterSim

Agent P1-control, 2026-10-28/29. Everything here was **computed** on this machine with QE 7.6
(`~/micromamba/envs/qe/bin`, 4 MPI ranks, nice 19, OMP_NUM_THREADS=1) unless it is marked
*cited*. Scripts, inputs, small outputs and a sha256 manifest of the large outputs are in
`runs/controls/mg2irh6/`. Raw scratch lives in `~/rtsc_targeted_work/P1-control/` and is not in the repo.
Every e-ph number is **tier-2, informational**: a coarse spectrum, not a verdict.

*Cited* reference: K. Dolui, L. J. Conway, C. Heil, T. A. Strobel, R. Prasankumar, C. J. Pickard,
"Feasible route to high-temperature ambient-pressure hydride superconductivity",
*PRL* 132, 166001 (2024), arXiv:2310.07562. The numbers below were read from the arXiv HTML version
on 2026-10-29: main-text a = 6.66 Å. Coarse DFPT (QE, PBE, USPP, 50/500 Ry, q 2³, k 4³, e-ph k 8³)
gives Allen–Dynes Tc 134.1 K at 0.005 Ry Gaussian broadening and 28.9 K at 0.05 Ry
(their Table 4; μ* not checked, **UNVERIFIED**). The converged run (QE, 80/640 Ry, q 6³, k 18³,
MP smearing 0.04 Ry) gives λ ≈ 2.5 at 0.005 Ry. EPW anisotropic gives λ = 2.3 and Tc = 160–175 K
(μ* = 0.1–0.16). The brief's a = 6.655 Å is consistent with "6.66".

## 0. Verdict

1. **PBE harmonic DFPT: Fm-3m Mg2IrH6 is dynamically stable at 0 GPa.** Every directly computed
   frequency is real at Γ, X, L (q 2³) and at the four irreducible points of q 3³. The q-grid
   Fourier interpolations (24³ and 32³ meshes) also have no imaginary branch. The lowest optical
   frequencies are 152.8 cm⁻¹ at Γ, 93.7 cm⁻¹ at X and 77.0 cm⁻¹ at L.
2. **MatterSim's instability is an ML artefact.** At the *same* PBE geometry, MatterSim-5M puts
   the Γ **T1g** mode at **−216 cm⁻¹** (−245 cm⁻¹ after it relaxes x_H itself). This is the
   rigid rotation of the IrH6 octahedra. DFPT puts the same T1g mode at **+242 cm⁻¹**. Identity
   of the mode: phonopy irreps on the MatterSim force constants; QE symmetry analysis on DFPT.
   The error is not small noise: the mode is off by ~460 cm⁻¹ and has the wrong sign of ω².
   MatterSim also softens every H mode by 15–20 % (A1g stretch 1456 vs 1834 cm⁻¹; T1u stretch
   1184 vs 1414 cm⁻¹).
3. Worse for the pipeline: in DFPT this T1g libration is the **most strongly coupled mode at Γ**
   (λ_Γ,ν ≈ 2.55 per mode at σ = 0.005 Ry). The harmonic ML force field fails on exactly the
   mode that carries the superconductivity. A MatterSim stability filter would therefore remove
   precisely the soft, strongly coupled hydrides we want.
4. Coarse e-ph (q 3³; two valid runs, k_scf/k_fit = 8³/24³ and 6³/36³): λ = 2.42–2.53 at
   σ = 0.005 Ry, falling to 0.91 at 0.05 Ry; ω_log ≈ 576–620 K; A = ∫α²F dω ≈ 990 K at
   σ = 0.005 Ry. Allen–Dynes Tc (μ* = 0.10) is 128–129 K at σ = 0.005 Ry and 101–109 K at
   0.010 Ry. This matches Dolui's coarse
   Allen–Dynes value (134 K at 0.005 Ry) and their λ ≈ 2.5. Allen–Dynes underestimates
   Eliashberg Tc at λ > 2, which is consistent with their 160–175 K. The control reproduces the
   literature at coarse level. **The pipeline's QE e-ph chain (coarse_dfpt.py) is validated on an
   ambient-pressure hydride.**

**Rule (reaffirms 10_ml_capabilities §3.4, now with DFT evidence):** MatterSim harmonic phonons
must never reject a hydride with MH6-type (or other) rotational/librational H modes. An
imaginary MatterSim mode means "run a DFPT Γ check". That check costs ~11 min for a 9-atom cell
on 4 ranks (§2). Use Mg2IrH6 as the positive control of every ML stability filter: a filter
that rejects it is miscalibrated.

## 1. Settings and convergence

Structure: Fm-3m (K2PtCl6 type), ibrav = 2, 9-atom primitive cell. Ir sits on 4a (0,0,0),
Mg on 8c ±(¼,¼,¼), H on 24e (x,0,0) and permutations. The start is the Alex-MP DFT entry
`msft-2023062809350933252138247157-GGA` (a = 6.633 Å, x_H = 0.2606).
Pseudopotentials are SG15 ONCV PBE-1.2, norm-conserving, scalar-relativistic, with no NLCC:
Mg (10 e⁻), Ir (17 e⁻), H (1 e⁻), 43 valence electrons, 26 bands.
Other settings are the core conventions (`rtt/core/dfpt/qe_inputs.py`): non-spin-polarised,
Marzari–Vanderbilt smearing 0.02 Ry, conv_thr 1e-10, mixing_beta 0.5, ecutrho = 4·ecutwfc, and
tr2_ph 1e-16.

**Cutoff** (scf, k 6³, Alex-MP geometry; `ecut/ecut.json`). ΔE is relative to 120 Ry:

| ecutwfc (Ry) | ΔE (meV/atom) | P (kbar) | F_x(H) (Ry/bohr) | wall (s) |
|---|---|---|---|---|
| 40 | 44.9 | −37.6 | −0.00341 | 4 |
| 50 | 15.8 | −4.4 | −0.00387 | 5 |
| 60 | 8.7 | −2.7 | −0.00401 | 5 |
| 70 | 4.1 | −1.8 | −0.00410 | 9 |
| **80** | **1.6** | **−0.4** | **−0.00411** | 12 |
| 90 | 0.5 | +0.6 | −0.00412 | 10 |
| 100 | 0.1 | +1.1 | −0.00412 | 14 |
| 120 | 0 | +1.3 | −0.00412 | 18 |

Production uses **80 Ry / 320 Ry**: 1.6 meV/atom, 1.8 kbar and 0.2 % in the H force from the
120 Ry values. The relaxed lattice at 100 Ry differs by 0.005 Å (below). SG15 H needs ≥ 80 Ry;
50–60 Ry, typical for USPP, would be under-converged here. That is the same wavefunction
cutoff as Dolui's converged run (80 Ry with USPP).

**k grid** (scf, 80 Ry; `kconv/kconv.json`): ΔE relative to 16³ is 0.60/0.33/−0.05/−0.02 meV/atom
at 6³/8³/10³/12³, and P stays within 0.25 kbar. 8³ (29 irreducible k) is used for DFPT (but see the Γ phonon k test in §2: the T1g mode needs ≥ 12³) and 10³
for the relaxation.

**vc-relax** (0 GPa, cell_dofree = 'ibrav', BFGS, press_conv_thr 0.2 kbar, forc 1e-4 Ry/bohr;
`relax_e80`, `relax_e100`):

| Run | a (Å) | x_H | d(Ir–H) (Å) | final P (kbar) | wall |
|---|---|---|---|---|---|
| 80 Ry, k 10³ (production) | **6.6521** | 0.25862 | 1.720 | 0.08 | 138 s |
| 100 Ry, k 10³ | 6.6567 | 0.25846 | 1.720 | 0.08 | 150 s |
| 80 Ry, started from the MatterSim "second minimum" (a 6.528, x 0.2512, Ir–H 1.64 Å) | 6.6510 | 0.25868 | 1.720 | 0.37 | 278 s (2 ranks) |
| Dolui et al. (*cited*, PBE USPP) | 6.655 / "6.66" | — | — | — | — |
| Alex-MP entry (PBE PAW, *cited* via reference set) | 6.633 | 0.2606 | 1.729 | — | — |
| MatterSim-5M (10_ml_capabilities §3.4) | 6.713 | — | — | — | — |

The PBE/SG15 lattice agrees with Dolui to **−0.003 Å (−0.05 %)** at 80 Ry and +0.002 Å at
100 Ry. The relax started from MatterSim's second Fm-3m minimum returns to the same structure:
ΔE < 0.01 meV/atom and Δa = 0.001 Å. PBE has **one** minimum along this path, so MatterSim's
double well in the IrH6 complex is also an artefact. MatterSim overestimates a by +0.9 %.

## 2. Harmonic DFPT phonons (stability)

`phon_q2_e80_k8`: scf 80 Ry, k 8³, at a = 6.6521 Å and x = 0.25862, residual P = 0.3 kbar;
ph.x with ldisp and q 2³, phonons only. Frequencies in cm⁻¹ are before the acoustic sum rule
(ASR); the Γ acoustic modes are at 5.3 cm⁻¹ before ASR.

| q | lowest modes | top | wall (4 ranks) |
|---|---|---|---|
| Γ | 5.3 (×3, acoustic) · **152.8 T1u** · 203.6 T2g · **242.4 T1g** · 711.9 T2g · 786.5 T2u · 853.3 T1u | 1413.5 T1u, 1520.8 Eg, 1834.4 A1g | 10.7 min |
| L (½,½,½) | **77.0** (×2, Eu), 131.1 (A2u), 175.6, 182.1 | 1628.4 | 26.5 min |
| X (0,1,0) | **93.7** (×2), 128.7, 142.6, 144.9 | 1578.2 | 20.1 min |

Total: 57 min. After q2r (zasr = crystal) and matdyn with ASR on a 24³ mesh, the minimum over
13 824 q is 0.0 cm⁻¹ (acoustic at Γ) and no frequency is below −1 cm⁻¹. The q 3³ run (§3) adds
four more irreducible points, all real: Γ; (−⅓,⅓,−⅓)·2π/a, i.e. ⅔ of Γ–L, at 73.9 cm⁻¹;
(0,⅔,0)·2π/a, ⅔ of Γ–X, at 83.0 cm⁻¹; (⅔,0,⅔)·2π/a, on Γ–K, at 98.3 cm⁻¹. Its 32³ matdyn mesh
minimum is also 0.0 cm⁻¹. All q-points are in units of 2π/a with a the conventional cubic
lattice constant; ph.x prints them in units of 2π/a_prim.

**Γ k-convergence** (`phon_G_e80_k{12,16}`, 2 ranks, 51 and 103 min; k 6 comes from the
e-ph run). Only the T1g libration depends noticeably on k:

| k_scf | T1u | T2g | **T1g** | T2g | T2u | T1u | T1u | Eg | A1g |
|---|---|---|---|---|---|---|---|---|---|
| 6³ | 153.2 | 203.8 | **268.7** | 724.2 | 785.9 | 851.9 | 1416.9 | 1530.0 | 1834.7 |
| 8³ (production) | 152.8 | 203.6 | **242.4** | 711.9 | 786.5 | 853.3 | 1413.5 | 1520.8 | 1834.4 |
| 12³ | 152.9 | 203.7 | **261.3** | 715.3 | 788.1 | 854.0 | 1414.7 | 1532.3 | 1834.5 |
| 16³ | 152.9 | 203.7 | **259.6** | 714.9 | 788.0 | 854.0 | 1414.6 | 1531.4 | 1834.5 |

The converged T1g frequency is ≈ 260 cm⁻¹. The k 8³ value is 7 % soft. This is the
Fermi-surface-sensitive (strongly coupled) mode. Every other mode is within 1 %. Stability is
not in question at any k. The k 8³ e-ph run overweights this mode's λ by about
(260/242)² ≈ 1.15.

MatterSim-5M at the same geometry (`scripts/mattersim_gamma_irreps.py`, CPU, Γ exact in the 1×1×1 cell,
finite displacement 0.01 Å):

| Γ irrep | DFPT (cm⁻¹) | MatterSim at QE geometry | MatterSim, x_H relaxed (0.25567) |
|---|---|---|---|
| T1g (IrH6 rotation) | **+242.4** | **−216.1** | **−244.7** |
| T1u | 152.8 | 115.7 | 114.5 |
| T2g | 203.6 | 155.7 | 155.6 |
| T2g | 711.9 | 530.4 | 507.5 |
| T2u | 786.5 | 644.7 | 651.5 |
| T1u | 853.3 | 739.1 | 732.1 |
| T1u | 1413.5 | 1183.7 | 1242.4 |
| Eg | 1520.8 | 1286.7 | 1380.7 |
| A1g | 1834.4 | 1456.0 | 1493.7 |

At the PBE geometry the MatterSim residual forces are 0.16 eV/Å and its pressure is 0.56 GPa.
Its Γ T1g value agrees with 10_ml_capabilities §3.4 (−245 cm⁻¹ at a = 6.655 Å, from 3×3×3
supercells). That section's X and L instabilities were not re-run here. The DFPT X and L
values above are real.

## 3. Coarse electron–phonon run

Chain: `src/scoring/coarse_dfpt.py`, the same QE sequence as `rtt.core.dfpt`. Steps:
scf_fit (la2F, k_fit) → scf (k 8³) → ph.x (electron_phonon = 'interpolated', ldisp, q 3³,
el_ph_sigma 0.005 Ry × 1..10) → q2r (la2F) → matdyn (la2F, dos, 32³, ndos 200). Template:
`scripts/make_template.py 6.652064 0.2586234 80`.

**Timings** (`eph_q3_k8_kf24/coarse_run.json`, 4 ranks, shared with another agent's 6 ranks):

| step | wall |
|---|---|
| scf_fit, k 24³ (413 irreducible k) | 209 s |
| scf, k 8³ | 20 s |
| ph.x: Γ / ⅔Γ–L / ⅔Γ–X / Γ–K point (240, 200, 328 k+q listed) | 15.4 / 52 / 37 / 65 min = **2.82 h** (elphsum, 4 calls, 16.4 min of it) |
| q2r + matdyn (32³, 10 σ) | 2 s + 332 s |
| **total** | **2.98 h** |

**Two k_fit = 36³ attempts are INVALID and were discarded.** Both used k_scf = 8³:
(a) `scripts/chain_kf36.sh`, a re-sum that reused the q 3³ dvscf with trans = .false.;
(b) a fresh coarse_dfpt.py chain. Both gave unphysical couplings: Γ acoustic γ of 195–9244 GHz,
a triply degenerate T1u mode split into λ = 0.22/0.82/0.46, and λ = 43.9 for the 152.8 cm⁻¹ mode
in (b). The cause is QE elphsum's integer ratio nti = nk1fit/nk1 = 36/8. QE's elphsum needs
k_fit to be a multiple of k_scf; the core `dfpt/grids.py` enforces this, and
`src/scoring/coarse_dfpt.py` does not (§3.1). The valid k_fit-convergence point is the
commensurate run **k_scf 6³ / k_fit 36³** (`eph_q3_k6_kf36`) below. The k_fit 24³ run satisfies
both QE hard constraints (24 % 8 = 0, 24 % 3 = 0). Only the "usual practice" k_scf % nq = 0
is not met, and QE does not require it. Its internal checks pass: Γ acoustic γ ≤ 0.02 GHz, and
degenerate modes agree to < 1 %.

**α²F moments** (`scripts/a2f_table.py`). σ is the el_ph_sigma Gaussian double-delta width.
A = ∫α²F dω. "A/4.33" is the float value of the universal bound k_B·Tc ≤ A/4.33148.
Allen–Dynes uses μ* = 0.10 and 0.13. The negative α²F weight fraction is ≤ 1.2·10⁻⁴ in every
smearing, so all ten pass the "mode_positive" admission rule.

k_fit 24³:

| σ (Ry) | λ | ω_log (K) | ω₂ (K) | A (K) | A/4.33 (K) | Tc_AD μ*=0.10 | Tc_AD μ*=0.13 |
|---|---|---|---|---|---|---|---|
| 0.005 | **2.423** | 603 | 1011 | 990 | 229 | **129.0** | 117.6 |
| 0.010 | 1.995 | 615 | 1025 | 830 | 192 | 109.2 | 99.2 |
| 0.015 | 1.711 | 619 | 1031 | 716 | 165 | 93.9 | 84.8 |
| 0.020 | 1.513 | 620 | 1035 | 636 | 147 | 82.2 | 73.6 |
| 0.025 | 1.359 | 621 | 1038 | 572 | 132 | 72.3 | 64.2 |
| 0.030 | 1.233 | 621 | 1040 | 520 | 120 | 63.8 | 56.0 |
| 0.035 | 1.130 | 621 | 1042 | 477 | 110 | 56.4 | 48.9 |
| 0.040 | 1.044 | 621 | 1044 | 441 | 102 | 50.0 | 42.8 |
| 0.045 | 0.972 | 622 | 1045 | 411 | 95 | 44.5 | 37.6 |
| 0.050 | 0.913 | 623 | 1047 | 387 | 89 | 39.9 | 33.3 |

k_scf 6³ / k_fit 36³ (`eph_q3_k6_kf36`; commensurate; 1240 irreducible fit k). Walls:
scf_fit 611 s, ph.x 3.10 h (elphsum 4602 s of it), matdyn 408 s; total 3.3 h. The lowest
q 3³ frequencies are 6.0 (Γ acoustic), 67.5, 76.6 and 91.3 cm⁻¹, all real. The Γ T1g mode is at
268.7 cm⁻¹ here.

| σ (Ry) | λ | ω_log (K) | ω₂ (K) | A (K) | A/4.33 (K) | Tc_AD μ*=0.10 | Tc_AD μ*=0.13 |
|---|---|---|---|---|---|---|---|
| 0.005 | **2.525** | 576 | 973 | 994 | 230 | **128.1** | 116.9 |
| 0.010 | 1.949 | 583 | 982 | 775 | 179 | 101.2 | 91.9 |
| 0.015 | 1.679 | 584 | 986 | 670 | 155 | 87.0 | 78.5 |
| 0.020 | 1.490 | 585 | 990 | 596 | 138 | 76.1 | 68.1 |
| 0.030 | 1.220 | 584 | 993 | 489 | 113 | 59.1 | 51.9 |
| 0.050 | 0.908 | 586 | 999 | 365 | 84 | 37.2 | 30.9 |


The λ values match QE's own footer in each `a2F.dosN` (2.4225 at σ1). λ depends strongly on σ
because N(E_F) itself changes with σ: QE's `lambda` file gives N(E_F) = 18.96 → 8.51 states/Ry
from σ = 0.005 to 0.05. This is the known sharp DOS structure at E_F in this compound. The
small-σ end is the physical limit, but it is the least converged in k_fit. The two valid runs
change k_scf and k_fit together (8³/24³ → 6³/36³). They agree to within 0–4 % in λ at every
σ ≥ 0.005 Ry: 2.42 vs 2.53 at 0.005, 2.00 vs 1.95 at 0.010, 0.91 vs 0.91 at 0.05. Tc_AD at
σ = 0.005 Ry is 129 vs 128 K. A ≈ 990 K in both. ω_log differs by 4 % (603 vs 576 K). The
remaining σ dependence is therefore mostly physical: the sharp N(E_F) structure, not k_fit
noise. The coarse-q (3³) and k_scf errors are not covered by this test. The matdyn `lambda`
file gives λ 2–4 % above the binned α²F integral (2.47 vs 2.42 and 2.64 vs 2.53). That is the
200-bin DOS discretisation.
ω_log is almost σ-independent (600–620 K). The spectral weight splits about evenly between
ω < 600 cm⁻¹ (Ir/Mg modes and H librations: λ ≈ 1.23 of 2.42 at σ1) and the H bend/stretch
bands above 600 cm⁻¹ (λ ≈ 1.19).

**Comparison with Dolui et al.** Their coarse q 2³ λ is ≈ 2.5 at 0.005 Ry; ours is 2.42 with q 3³.
Their coarse Allen–Dynes Tc is 134 K at 0.005 Ry; ours is 129 K at μ* = 0.10. At 0.05 Ry they
get 28.9 K and we get 39.9 K. Their converged EPW λ is 2.3. We agree to within the
coarse-sampling spread. Our ω_log is ≈ 600 K; Dolui's ω_log was not found in the text read,
**UNVERIFIED**. For the targeted pipeline, **A ≈ 1000 K** at the physical (small-σ) end, so
k_B·Tc ≤ A/4.33 ≈ 230 K. Mg2IrH6 lies below the A ≳ 1300 K that 300 K would require, which is
expected: it is a 160 K material.

### 3.1 Two bugs in shared scoring code (not fixed here: owned by other agents)

**(i) `src/scoring/coarse_dfpt.py` grid check.** It checks only `kfit % nq == 0`. It must also
require `kfit % kscf == 0`, or QE silently produces garbage λ (above). The core
`rtt/core/dfpt/grids.py` checks this together with `kscf % nq == 0`. Suggested one-line fix next
to the existing `ap.error`:
`if a.kfit % a.kscf: ap.error("k_fit must be a multiple of k_scf (QE elphsum nti = nk1fit/nk1)")`.

**(ii) `src/scoring/a2f_summary.py` integration.**

`summarise()` first drops every point with α²F ≤ 0 and then integrates with
`np.gradient(w)` on the remaining **non-uniform** frequencies. On a gapped spectrum (the q 3³
phonon DOS of Mg2IrH6 is zero over 562–645 and 857–1300 cm⁻¹), the points at the gap edges get a dw of about half
the gap width. On this run it reports λ = 3.68 instead of 2.42 (QE: 2.4225), ω_log = 877 K
instead of 603 K, A = 2138 K instead of 990 K and Tc_AD = 267 K instead of 129 K. It inflates
λ, A and Tc for every material with phonon gaps (most hydrides). **Fix:** integrate on the full
uniform grid with `dw = w[1] - w[0]` and clip negative α²F to 0 instead of deleting points.
That is what `runs/controls/mg2irh6/scripts/a2f_table.py` does; it imports `allen_dynes` and
`RY_TO_K` from a2f_summary. Any number already produced with a2f_summary on gapped spectra
should be recomputed.

## 4. Cost of a full (converged) run

Measured scaling: ph.x wall ≈ 0.2 min per listed (k, k+q) point per irreducible q at 80 Ry,
4 ranks. Γ costs ≈ 0.37 min per k. The scf_fit and elphsum costs scale with k_fit³.

| Run | irreducible q | k_scf | k_fit | estimated wall, 4 ranks |
|---|---|---|---|---|
| this run | 4 (q 3³) | 8³ / 6³ | 24³ / 36³ | 3.0 h / 3.3 h (measured) |
| intermediate | 8 (q 4³) | 12³ | 48³ | ≈ 8 × 0.2 min × ~1100 k ≈ 30 h + ~1 h scf_fit + ~2 h elphsum ≈ **1.4 days** |
| Dolui-like | 16 (q 6³) | 18³ | 36³–72³ | ≈ 16 × 0.2 min × ~3900 k ≈ 210 h ≈ **9 days** (+ elphsum); RAM per rank stays < 1 GB |
| EPW (Wannier, Dolui's anisotropic Tc) | — | — | 36³ fine | needs the q 6³ (or 4³) DFPT above plus Wannierization; not attempted |

The estimates assume linear scaling with the number of k-points and no I/O bottleneck. They are
extrapolated, not measured. The practical step is a q 4³ run with k_fit 48³, which fits in
~1.5 days on 4 ranks. A Γ-only stability check costs ~11 min and the q 2³ phonon check ~1 h
for a 9-atom hydride; that is the right "DFPT veto" for MatterSim imaginary modes.

## 5. Files (`runs/controls/mg2irh6/`)

- `scripts/qe_mg2irh6.py` — ecut / kconv / relax / phonon (q grid) driver.
- `scripts/make_template.py` — ibrav = 0 template for coarse_dfpt.py.
- `scripts/chain_eph.sh` — the valid k 8³/k_fit 24³ chain; `scripts/chain_kf36.sh` — the INVALID non-commensurate re-sum (kept as a record). The valid k 6³/k_fit 36³ run is a direct `coarse_dfpt.py --kscf 6 --kfit 36 --nq 3` call.
- `scripts/a2f_table.py` — corrected α²F moments.
- `scripts/mattersim_gamma_irreps.py` — MatterSim Γ phonons + irreps on CPU.
- `scripts/collect.sh` — copies small outputs here and writes `sha256_manifest.txt`.
- Per run directory: `*.in`, `*.json`, `dyn*`, `a2F.dos*`, `lambda`, `phonon.dos`, `freq.mesh.gp`,
  and `*.out.summary` (energies, pressures, frequencies, λ_qν, walls grepped from the QE
  stdout). Large files (wavefunctions, dvscf, a2Fsave, full stdout) stay in scratch and are
  listed with sha256 in `sha256_manifest.txt`.
- `mattersim_gamma.json` — the MatterSim Γ result.

## 6. Not verified / open

- **Anharmonicity and zero-point motion** were not computed. Dolui (*cited*) report that SSCHA
  keeps the structure stable and hardens the H T1g modes; this was not checked here.
- Only q 3³ e-ph. k_fit was tested only jointly with k_scf (8³/24³ vs 6³/36³; λ agrees to 4 %), and k_scf affects the T1g frequency by up to 10 % (Γ table). q convergence was not tested. Tc here
  is Allen–Dynes, not Eliashberg; the certified Eliashberg bounds (`rtt.core.certify`) were not run on these spectra.
- MatterSim's X and L imaginary modes were not re-identified by irrep (only Γ was).
- The SG15 Ir pseudo has no semicore 4f in valence (17 e⁻: 5s5p5d6s). No PAW/all-electron
  cross-check was made beyond the lattice-constant agreement with Dolui and Alex-MP.
- Dolui's μ* for Table 4 and their ω_log were not read; marked **UNVERIFIED**.
- The machine was shared with another agent's QE run (6 ranks) during these timings, so
  walls on an idle machine may be ~10–30 % shorter (**UNVERIFIED**).
