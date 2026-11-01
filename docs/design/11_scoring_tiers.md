# 11. Scoring tiers between MatterSim screening and certified DFPT

Status: design note with measurements. Author: agent T4-surrogates. Date: 2026-10-28.
Scripts: `src/scoring/` (all CPU-only; see §6). Raw outputs live in the agent scratch dir
`~/rtsc_targeted_work/T4-surrogates/` (not in the repo).

Labels used below: **[measured]** = run on this machine in this work, script given;
**[cited]** = from the named paper or repository page, read during this work;
**[est.]** = my own estimate, with the logic shown; **UNVERIFIED** = recalled or not checked.

---

## 0. Answer in one table

| Tier | What runs | Wall time per structure (≤ 10 atoms) on this machine | Expected error in the target regime (λ ≳ 2, ω_log ≳ 1000 K) | May decide | May NOT decide |
|---|---|---|---|---|---|
| T0 | MatterSim relax + hull + finite-displacement phonons (exists, `rtt.mlff`) | seconds to minutes (GPU or CPU) | n/a for e-ph | DEFER (ML imaginary modes, far above hull) | anything about Tc; exclusion |
| T1a | Composition/structure heuristics (H fraction, light-element fraction, √(m_H/M), N(E_F)) | < 1 s | ordering only; no calibrated error | reorder | DEFER on its own; any verdict |
| T1b | ML Tc/λ/ω_log predictors (ALIGNN 2.0 checkpoints, BETE-NET) | < 1 s on CPU **[measured, ALIGNN]** | **general models: Tc 0.2–3 K and λ ≈ 0.3–0.5 for H3S / CaH6 / LaH10 (DFT: 200–280 K, λ_H3S = 2.4); hydride model: 2× low on H3S (which is in its training set), ~8× low on LaH10 [measured]**; useless as a magnitude estimate here | weak reorder *within* low-Tc chemistries | DEFER of hydride / high-ω candidates; any verdict |
| T1c | Belli networking value φ (needs SCF + ELF + H share of N(E_F)) | ~11 s (H3S, 4 atoms, 2 ranks) **[measured, §3.3]**; ~5–30 min for 10 low-symmetry atoms **[est.]** | ±60 K on Tc in its own training set **[cited]**; H3S: 273 K vs ~200–230 K **[measured]**; no phonons, so blind to dynamical (in)stability and to soft-mode λ | reorder hydrides | DEFER; any verdict |
| T2 | Coarse DFPT α²F (same QE chain as certification, smaller cutoff, k, k_fit and q grids, 2 ranks); **score on A, not λ** | 4 atoms, O_h: 36 s (q 2³), 4.5 min (q 3³), 11 min (q 4³), 35 min (q 6³) **[measured]**; 10 atoms, low symmetry: ~1–15 h **[est.]** | q 2³: λ −55…−60 %, A −50 %; q 3³ / 4³ / 6³ (σ = 0.03–0.05 Ry): λ +40 % / +10…15 % / +90 %, ω_log −25 % / −17 % / −71 %, **A +11…13 % / −1…+2 % / +13…14 %**; one spurious imaginary mode at q 6³ **[measured, one material]** | reorder; DEFER below a factor-2 margin on A (§4); choose which candidates get T3 | any verdict (not converged; no error bars that the certifier accepts); rejection on coarse imaginary modes |
| T3 | Converged DFPT α²F + certified bounds (`rtt.core.certify`, Arb/Lean-backed) | H3S (4 atoms): ph.x 1 h 58 min on 8 ranks **[measured, earlier benchmark run]**; 10 low-symmetry atoms: days **[est.]** | certified bracket (given α²F and stated assumptions) | REJECTED_PROVEN / CERTIFIED_CANDIDATE / UNDECIDED | — |

The one-line recommendation: **use T2 (coarse DFPT at q 3³–4³, judged on A = ∫α²F dω) as the
real filter, T1c and T1a as cheap orderers, and do not use any public ML Tc model to defer
high-Tc candidates.** The public models fail exactly where this project needs them (§2.3).
Coarse λ and ω_log are unreliable (spurious soft modes); coarse A was within −1 % … +14 % of
the converged value in every q ≥ 3³ run at σ ≥ 0.03 Ry (one material; at σ = 0.01 Ry, up to +56 %).

The governing principle: ML and cheap estimates only reorder or DEFER;
only certified computations (T3) issue REJECTED_PROVEN / CERTIFIED_CANDIDATE (contract S2:
nothing in the scoring code can produce a Status; scores only reorder candidates, and the budget
planner can only mark candidates as not computed, DEFERRED).

---

## 1. Why a middle tier is needed (cost anchors on this machine)

[measured, earlier benchmark runs of this pipeline]:

| Material | atoms | settings | ph.x wall (8 MPI ranks) | end-to-end |
|---|---|---|---|---|
| Al fcc | 1 | 29 irr. q (8³), k 16³ | 3271 s | 55 min |
| Pb fcc | 1 | 29 irr. q (8³), k 16³ | 8189 s | 2 h 18 min |
| H3S Im-3m, 200 GPa | 4 | 16 irr. q (6³), k 18³, k_fit 36³, ecut 90 Ry | ~7100 s (1 h 58 min) | ~2 h |

pw.x steps are < 1 min in all cases; ph.x is ≥ 98 % of the cost. Screening (MatterSim) is
seconds per structure. A generator like MatterGen produces thousands of structures per run, so
there is a gap of about four orders of magnitude between T0 and T3.

Reference spectrum used below (H3S from an earlier benchmark run, σ = 0.01 Ry, integrated by
`src/scoring/a2f_summary.py`, which reproduces matdyn's own λ = 2.3946 to 4 digits):
λ = 2.395, ω_log = 1250 K, ω₂ = 1478 K, **A = ∫α²F dω = 1648 K**, A/4.33148 = 380 K,
Allen–Dynes Tc = 229 K (μ* = 0.13). Over σ = 0.01…0.10 Ry, λ runs 2.40 → 1.51 and A runs
1648 → 1084 K: the smearing spread alone is ±20 % in A.

---

## 2. ML predictors of α²F / λ / ω_log / Tc

### 2.1 Survey

| Model | Reference | Training data (size, domain) | Reported accuracy | Weights / code / license | Runs on CPU here? |
|---|---|---|---|---|---|
| BETE-NET (CSO, CPD, FPD variants) | J. B. Gibson, A. C. Hire, P. M. Dee, O. Barrera, B. Geisler, P. J. Hirschfeld, R. G. Hennig, *npj Comput. Mater.* **11**, 7 (2025), doi:10.1038/s41524-024-01475-4; arXiv:2401.16611 **[cited]** | 818 dynamically stable materials, α²F from QE/PBEsol (recomputed Choudhary set + small MP metals ≤ 8 atoms); 27 unary, 456 binary, 333 ternary, 2 quaternary **[cited]**. Hydride content not stated. | Test MAE: CSO (structure only) λ 0.21, ω_log 45 K, Tc 2.5 K; FPD λ 0.18, ω_log 29–30 K, Tc 2.1 K; linear baseline Tc 2.84 K **[cited]**. Authors note under-predicted low-ω peaks → under-predicted λ for high-λ materials; weaker for C, N, O **[cited]**. No ensemble uncertainty reported. | github.com/henniggroup/BETE-NET, **Apache-2.0**, 100 checkpoints × 19.7 MB per variant (~2 GB each) in the repo **[checked via GitHub API]**. Pinned stack: torch 1.10, e3nn 0.4.2, torch-geometric 2.0.2 (CUDA 11.3 install line) | Probably (pure PyTorch), but needs a separate legacy env; **not tested** (UNVERIFIED) |
| BEE-NET (successor) | J. B. Gibson, A. C. Hire, P. Prakash, P. M. Dee, B. Geisler, J. S. Kim, Z. Li, J. J. Hamlin, G. R. Stewart, P. J. Hirschfeld, R. G. Hennig, *npj Comput. Mater.* **12**, 95 (2026), doi:10.1038/s41524-026-01964-8; arXiv:2503.20005 **[cited]** | 5241 α²F from the Cerqueira et al. set (subset of ~7000) **[cited, repo README]** | CPD: Tc MAE 0.87 K, λ MAE 0.073, ω_log MAE 17.4 K, true-negative rate 0.991–0.994 **[cited]**. Pipeline: 1.3 M structures → 741 DFT-confirmed Tc > 5 K; 2 synthesized and confirmed **[cited]** | github.com/henniggroup/BEE-NET: **no license file**, `models/` "not tracked in git", no download link found **[checked]** | No (weights not obtainable; no license) |
| ALIGNN (JARVIS e-ph set) | K. Choudhary, K. Garrity, *npj Comput. Mater.* **8**, 244 (2022), doi:10.1038/s41524-022-00933-1; arXiv:2205.00060 **[cited]** | 1058 e-ph calculations after BCS pre-screen of 1736 JARVIS materials; 105 stable with Tc > 5 K **[cited]**; q grid criticised by Cerqueira et al. as "2x2x2, severely underconverged" **[cited]** | Predicting α²F first beats direct Tc **[cited]**; MAE numbers not read (UNVERIFIED) | ALIGNN: github.com/usnistgov/alignn (→ atomgptlab/alignn), NIST public-domain-style terms **[checked]**. ALIGNN 2.0 checkpoint `tc_supercon_radius` (test MAE 1.64) on Figshare **[checked, registry]** | **Yes — measured** |
| ALIGNN 2.0 "Alexandria superconductor" models (`alex_supercon_{Tc,la,wlog,debye,dosef}`) | ALIGNN 2.0, arXiv:2609.19487 (2026) **[cited, abstract only]**; checkpoint registry in `alignn/pretrained.py` | 8241 Alexandria (`agm…`) entries, split 6592/824/825 **[measured from ids file]**. Presumably the Cerqueira–Sanna–Marques e-ph data (UNVERIFIED) | Registry test MAE: Tc 0.88 K, λ 0.071, ω_log 20.3 K **[cited, registry]** | Figshare, pure-torch (no DGL) **[checked]** | **Yes — measured** |
| ALIGNN 2.0 hydride Tc (`tc_supercon_hydride_radius`) | Data: K. Wines, K. Choudhary, *Mater. Futures* **3**, 025602 (2024), doi:10.1088/2752-5724/ad4a94; arXiv:2312.12694 **[cited]** | 953 hydride structures (763 train), 0–500 GPa, JARVIS DFT **[measured from ids file]**; >900 hydrides, 122 stable with Tc > 39 K **[cited]**. Contains H3S and CaH6 at 0–300 GPa in the **training** split **[measured]** | Registry test MAE 10.1 K; leaderboard (hydride + bulk) MAE 7.96 **[cited]** | Figshare **[checked]** | **Yes — measured** |
| Cerqueira–Sanna–Marques ML | T. F. T. Cerqueira, A. Sanna, M. A. L. Marques, *Adv. Mater.* **36**, e2307085 (2024); arXiv:2307.10728 **[cited]** | ~7000 e-ph calculations (QE, PBEsol, PseudoDojo NC; q = half the k grid, double-δ MP smearing 0.05 Ry); 200 k metals screened; 541–545 with Tc > 10 K (abstract versions differ) **[cited]** | modnet-type ensemble (λ, ω_log, Tc jointly); validation Tc MAE 2.74 K at mean Tc 8.07 K (predict-the-mean: 5.56 K) **[cited]** | Data: Materials Cloud, doi:10.24435/materialscloud:qv-bq, **CC BY 4.0** **[cited]**. Code "will be released on github" — release not found (UNVERIFIED) | Model: no. Data: yes (useful as a calibration set) |
| Tran & Vu | H. Tran, T. N. Vu, *Phys. Rev. Mater.* **7**, 054805 (2023); arXiv:2211.03265 **[cited]** | 584 structures (λ, ω_log) | not read (UNVERIFIED) | not checked | — |
| Hutcheon, Shipley, Needs | *Phys. Rev. B* **101**, 144505 (2020); arXiv:2001.09852 **[cited]** | literature binary hydrides; NN used to steer structure search | not read | not checked | — |
| Stanev et al. | V. Stanev et al., *npj Comput. Mater.* **4**, 29 (2018); arXiv:1709.02727 **[cited]** | ~12 000 SuperCon entries (experimental), composition only | ~92 % accuracy for Tc ≷ 10 K classifier **[cited]** | — | — |
| Hamidieh | K. Hamidieh, *Comput. Mater. Sci.* (2018); UCI dataset 464, 21 263 entries, CC BY 4.0 **[cited]** | SuperCon (experimental), 81 composition features | — | — | — |
| ML e-ph matrix elements (Hamiltonian + gradient from an equivariant GNN) | Y. Zhong, S. Liu, B. Zhang, Z. Tao, Y. Sun, W. Chu, X.-G. Gong, J.-H. Yang, H. Xiang, *Nat. Comput. Sci.* (2024), doi:10.1038/s43588-024-00680-x; arXiv:2302.00439 **[cited]** | per-system training (tests: H₂O, MoS₂, …) | "close agreement" with DFPT, orders-of-magnitude speed-up **[cited]** | not checked | Not a universal model: must be trained per chemical system |

**Universal e-ph ML models (2025–2026).** I did not find a released, universal (any-chemistry)
model of electron–phonon matrix elements or α²F beyond the BETE-NET/BEE-NET and ALIGNN families.
Two searches found related but not universal work: Zhong et al. 2024 (above, per-system),
Novoa, Belli, Bianco, Contreras-García, Errea, arXiv:2507.06749 (2025; real-space e-ph
descriptors for hydrides, no fitted Tc model in the abstract) **[cited]**, and a 2026 npj
boron/carbon study (doi:10.1038/s41524-026-01962-w; 417 compounds, chemistry-specific)
**[cited, search snippet only]**. Absence of evidence: a model released after my search index,
or under a name I did not search, may exist (UNVERIFIED).

**Why SuperCon-trained models (Stanev; Hamidieh; and similar composition models) are less
suitable.** (i) Labels are experimental Tc of mostly unconventional superconductors (cuprates,
Fe-based, heavy fermions) mixed with conventional ones; the mechanism is not part of the label,
so the model learns "cuprate-like composition → high Tc", which is irrelevant to an e-ph target.
(ii) Composition-only inputs cannot tell polymorphs apart (e.g. H3S Im-3m vs R3m, or a hydride
at 0 vs 200 GPa) — the pipeline's candidates differ precisely in structure. (iii) No high-pressure
hydrides with Tc > 200 K exist at 1 atm in SuperCon, so the target is outside the label range.
(iv) The output is a number with no link to α²F, so it cannot feed the A-based bounds used by
certification.

### 2.2 Measurement: public ALIGNN checkpoints on reference structures [measured]

Script: `src/scoring/alignn_screen.py` (ALIGNN 2026.8.11 pure-torch path, CPU, 2 threads,
`CUDA_VISIBLE_DEVICES=` empty). 5 models × 7 structures in 4.6 s wall including model loading,
peak RSS 1.0 GB; inference 0.03–0.34 s per structure for all 5 models. Primitive and conventional
cells agree to within ~1 K except LaH10 (26 vs 32 K).

| Structure | Reference (DFT / exp.) | `alex_supercon_la` (λ) | `alex_supercon_wlog` (K) | `alex_supercon_Tc` (K) | `tc_supercon_radius` (K) | `tc_supercon_hydride` (K) |
|---|---|---|---|---|---|---|
| Al fcc | exp 1.2 K | 0.44 | 277 | 2.2 | 1.1 | 4.1 |
| Pb fcc | exp 7.2 K | 1.54 | 42 | 7.3 | 2.0 | 5.3 |
| Nb bcc | exp 9.2 K | 1.01 | 161 | 13.3 | 3.4 | 6.7 |
| MgB2 | exp 39 K; isotropic ME ~20–25 K (UNVERIFIED) | 0.26 | 580 | 1.0 | 1.1 | 6.0 |
| H3S Im-3m (200 GPa cell) | bench DFPT: λ 2.40, ω_log 1250 K, AD Tc 229 K | **0.52** | 545 | **2.8** | 0.8 | 100 |
| CaH6 Im-3m (a = 3.50 Å, ~150 GPa) | DFT Tc ~215–250 K (Shipley et al. 2021 give 216–253 K at ~100 GPa) | 0.27 | 579 | 0.2 | 1.2 | 174 |
| LaH10 Fm-3m (a = 5.10 Å, ~200 GPa) | DFT Tc ~250–280 K (UNVERIFIED) | 0.43 | 292 | 1.4 | 1.3 | **32** |

(Lattice constants of CaH6 and LaH10 and the LaH10 H(32f) x = 0.12 are approximate,
UNVERIFIED; minimum H–H distances 1.24 Å and 1.06 Å, physically plausible.)

### 2.3 Failure modes for out-of-distribution high-Tc candidates

1. **Range compression.** The general models (Alexandria, JARVIS) are trained on sets whose Tc
   distribution is concentrated below ~20 K (Cerqueira et al. validation mean Tc 8.07 K). For
   H3S they predict λ = 0.52 instead of 2.40 and Tc ≈ 1–3 K instead of ~200 K. A model that
   has never seen λ > 2 regresses to the mean. Their good MAE (0.9 K) and high true-negative
   rate are statements about the *training distribution*, not about the tail we search.
2. **Misses even in-sample.** The hydride model has H3S and CaH6 (at 0–300 GPa) in its training
   split, and still returns 100 K for H3S (DFT ~200 K). For LaH10 (only LaH7 in training)
   it returns 26–32 K. A candidate that is "new" in the sense the pipeline wants will be at
   least as far out of distribution as LaH10.
3. **No pressure variable.** The hydride model was trained on the same structures relaxed at
   0/100/200/300/500 GPa; it sees pressure only through bond lengths. For a 1-atm target with
   compressed-looking H sublattices the model has no way to know whether the structure is
   dynamically stable.
4. **No stability awareness.** None of these models know about imaginary modes; T0 (MatterSim)
   and T2/T3 must do that.
5. **Light elements are the documented weak spot** (BETE-NET authors: C, N, O; low-ω peak
   magnitudes under-predicted → λ of high-λ materials under-predicted) **[cited]**.
6. **Cerqueira et al.'s own warning** that the coarse first pass yields false positives
   (caught later) and may miss false negatives (never caught) **[cited]** — exactly the
   asymmetric risk that forbids ML deferral of high-Tc-looking candidates.

Consequence: **a public ML model may not be used to DEFER any candidate that T1a/T1c rank as
hydrogen-rich or light-element-rich.** At most it can reorder within chemistries that are
inside its training distribution (intermetallics, nitrides, borides with Tc < 40 K), which is
not where a Tc ≥ 300 K candidate is expected to be.

What would change this: fine-tuning an α²F-predicting model (BETE-NET's architecture is
Apache-2.0 and its weights are available) on our own T2/T3 spectra of hydrides. That is a
later, GPU task; until it is validated on held-out T3 spectra it stays in the "reorder only"
class.

---

## 3. Cheap physics estimators

### 3.1 λ from coarse DFPT [measured]

Script: `src/scoring/coarse_dfpt_probe.sh` (same QE 7.6 chain as `rtt.core.dfpt`:
`scf_fit (la2F)` → `scf` → `ph.x` (`electron_phonon='interpolated'`) → `q2r` → `matdyn`),
H3S Im-3m in the bench's relaxed 200 GPa cell, SG15 ONCV PBE pseudopotentials,
**2 MPI ranks, 1 thread each, nice 19**. Spectra integrated with `src/scoring/a2f_summary.py`.

| Run | ecut (Ry) | k_scf / k_fit | q grid (irr. q) | ph.x wall | λ (σ=0.01 / 0.03 / 0.05 Ry) | ω_log (K, σ=0.03) | A (K, σ=0.01 / 0.03 / 0.05) | AD Tc μ*=0.13 (σ=0.03) |
|---|---|---|---|---|---|---|---|---|
| reference (bench, 8 ranks) | 90 | 18³ / 36³ | 6³ (16) | 7100 s | 2.40 / 2.10 / 1.86 | 1287 | 1648 / 1480 / 1316 | 210 K |
| c1 | 60 | 8³ / 16³ | 2³ (3) | **36 s** | 0.98 / 0.94 / 0.97 | 1490 | 827 / 763 / 786 | 84 K |
| c2 | 60 | 12³ / 24³ | 3³ (4) | **271 s** | 4.64 / 2.95 / 2.59 | 972 | 2576 / 1671 / 1467 | 217 K |
| c3 | 60 | 12³ / 24³ | 4³ (8) | **635 s** | 2.82 / 2.32 / 2.14 | 1065 | 1778 / 1462 / 1347 | 194 K |
| c4 | 60 | 12³ / 24³ | 6³ (16) | **2099 s** | 5.32 / 4.06 / 3.57 | 370 | 2179 / 1674 / 1496 | 119 K |

pw.x steps: 2–8 s; q2r + matdyn: ~22 s, in every run. (λ, ω_log, A from
`a2f_summary.py` on matdyn's `a2F.dosN`; QE's own `lambda` file, which sums λ_qν on the q grid,
differs by a few %, e.g. c3 σ=0.01: 2.95 vs 2.82.)

> **Erratum (doc 15 finding).** `a2f_summary.py` used to drop α²F ≤ 0 points and then take
> `np.gradient` on the uneven remainder, inflating integrals across phonon gaps. Fixed; the
> coarse rows recomputed with the fixed script (σ = 0.01 / 0.03 / 0.05 Ry):
>
> | Run | λ | ω_log (K, σ=0.03) | A (K) | AD Tc μ*=0.13 (σ=0.03) |
> |---|---|---|---|---|
> | c1 | 0.87 / 0.97 / 0.95 | 1497 | 712 / 786 / 773 | 89 K |
> | c2 | 3.29 / 2.45 / 2.07 | 962 | 1861 / 1389 / 1186 | 184 K |
> | c3 | 2.46 / 2.06 / 1.82 | 1059 | 1556 / 1298 / 1147 | 173 K |
> | c4 | 4.44 / 3.38 / 2.87 | 385 | 1816 / 1427 / 1235 | 104 K |
>
> The reference row was not produced with this script and is unchanged. The doc 13 calibration
> (e_c) does not use `a2f_summary.py` and is unaffected; production `rtt` integrates via the
> core `rtt.core.spectrum`.

Lowest phonon frequencies on the q grid (cm⁻¹): reference 10 (Γ acoustic), then 180 (q2),
226 (q10); c3: 243, 270, 304; **c4: −190 (an imaginary mode at q13), 119 (q2), 199 (q10)**.
The c4 run uses the reference's own q grid, but the lower cutoff and k grid make H3S —
which sits close to a soft-mode instability at this pressure — spuriously unstable at one q and
softer elsewhere. That moves weight to low ω: λ is +93 % and ω_log −71 % at σ = 0.03, yet
**A is only +13 %**, because A = ∫α²F dω carries no 1/ω weight.

Observations (one material only — treat as an anecdote, not a calibration):

- **q 2³ is not usable.** It samples only Γ, H and one more point; H3S's coupling is carried
  by modes away from those points. λ is low by 60 %, A by 50 %, ω_log high by 16–20 %.
  Any method that samples only Γ or a few high-symmetry q (including frozen-phonon
  deformation potentials, §3.2) inherits this error.
- **Small σ with a coarse k_fit is unstable.** With k_fit 24³, σ = 0.01 Ry gives λ = 4.6
  (q 3³) and 2.8 (q 4³), i.e. noise from too few Fermi-surface points. At σ = 0.03–0.05 Ry
  λ is high by +40 % (q 3³) and +10…+15 % (q 4³) against the reference at the same σ.
  Coarse runs must report a σ-plateau value, never the smallest σ.
- **ω_log is biased low by 15–25 %** at q 3³–4³ and by 71 % in c4 (soft modes);
  λ and ω_log errors partly cancel in A (A = λ⟨ω⟩/2 with ⟨ω⟩ = ∫α²F dω / ∫α²F/ω dω), which at
  σ = 0.03 is −1 % (q 4³), +13 % (q 3³) and +13 % (q 6³, c4) against the reference.
  **A is the robust coarse quantity; λ and ω_log are not.** This is convenient: A is also
  the quantity in the universal bound k_B Tc ≤ A/4.33148 used by certification.
- **Refining one coarse parameter can make things worse.** c4 (finer q than c3) is further from
  the reference in λ than c3, because its denser q grid hits the q point where the
  under-converged cutoff/k grid produce an instability. Coarse settings must be
  calibrated as a set, not assumed to converge monotonically.
- **Coarse DFPT can report spurious instabilities.** An imaginary coarse mode is therefore no
  reason to reject or defer; it is a reason to run a converged phonon check.
- Cost: per core-second, q 3³ is ~105× cheaper than the reference
  (271 s × 2 ranks vs ~7100 s × 8 ranks), q 4³ ~45×, and c4 (q 6³ at low cutoff/k) only ~14×.
  Most of the reference's cost is in the q count, so a cheap tier has to cut q.

**Cost for ≤ 10-atom cells [est.].** ph.x time scales roughly as
N_irr-q × (3 N_at) × N_k(q) × N_PW log N_PW × N_bands, with N_k(q) the k points kept by the
small group of q. Going from H3S (4 atoms, O_h, 12 modes) to a 10-atom, low-symmetry hydride:
modes ×2.5, N_PW ×~2.5 (volume), bands ×~2.5, and symmetry loss multiplies both the number
of irreducible q (×~3–6 at q 4³) and k per q (×~3–10). Combined: ×~150–2000. With the c3
anchor (635 s) that gives **~1–15 h on 2 ranks** for a q 3³–4³ coarse run on a 10-atom,
low-symmetry cell, and **minutes to ~1 h** for small high-symmetry cells (≤ 5 atoms, cubic or
hexagonal). The 15 GB RAM limit also bites for large k_fit with 10 atoms (UNVERIFIED; not
measured). So T2 is "cheap" only relative to T3; it must be applied to tens, not thousands,
of candidates.

### 3.2 Hopfield / deformation-potential estimates

- **Gaspari–Gyorffy / rigid muffin-tin (RMTA).** λ = Σ_atoms η_i / (M_i ⟨ω²⟩_i), with the
  McMillan–Hopfield parameter η_i from the electronic structure only (phase shifts and
  partial DOS at E_F). Used for hydride screening by Shipley, Hutcheon, Needs, Pickard,
  *Phys. Rev. B* **104**, 054501 (2021), arXiv:2105.02296, as a DOS-based Tc model to pick
  121 candidates for converged e-ph calculations **[cited]**; and for H3X hydrides
  (arXiv:1903.03255) **[cited, snippet]**. A plane-wave/pseudopotential RMTA implementation
  appeared in 2025 (Radevych, Shishidou, Weinert, Margine, Kolmogorov, Mazin,
  arXiv:2508.21281), validated against LAPW η, not against DFPT λ **[cited]**; no public code
  found (UNVERIFIED). QE has no RMTA. Cost if implemented: one SCF + DOS (minutes) plus
  ⟨ω²⟩_H from MatterSim phonons. Weakness: RMTA assumes a rigidly moving spherical potential,
  which is least justified for covalent H networks (H3S, clathrates) — the systems we care
  about. No large-scale benchmark of RMTA λ vs DFPT λ for hydrides was found (UNVERIFIED).
  **Not recommended for implementation now**: no code, no calibration in our regime.
- **Frozen-phonon deformation potentials at Γ or a few q.** For a phonon at commensurate q,
  displace ±u along the eigenvector in a supercell, compute the Fermi-surface average of
  |∂ε_nk/∂u|² (intraband at q = Γ) and assemble λ_qν. At q = Γ the intraband coupling is not
  the quantity that enters λ for a metal with a 3D Fermi surface (it misses all
  inter-sheet / finite-q scattering); at a few zone-boundary q it needs supercells (2× to
  8× the atoms) and a dense k-mesh nscf per displacement. Cost per q is therefore comparable
  to or above one DFPT q point at the same k grid, and the sampling error is the same as
  the q 2³ row of §3.1 (λ −60 % for H3S). **Not recommended**: coarse DFPT dominates it
  in both cost and accuracy on this machine.

### 3.3 Belli networking value φ [cited + measured cost]

F. Belli, T. Novoa, J. Contreras-García, I. Errea, *Nat. Commun.* **12**, 5381 (2021),
doi:10.1038/s41467-021-25687-0 **[cited]**: φ is the highest ELF value whose isosurface
percolates through the crystal in all three directions; with hydrogen fraction H_f and
H fraction of N(E_F) H_DOS, Tc ≈ (750 φ H_f H_DOS^{1/3} − 85) K, scatter ~±60 K, fitted on
~135 hydrogen-based compounds whose Tc came from the literature (mixed methods and μ*);
DOS-based descriptors give false positives for molecular hydrides **[cited]**. No phonons are
needed. Inputs: one SCF, ELF (QE `pp.x plot_num=8`, norm-conserving pseudopotentials — SG15
qualifies), and a projected DOS (`projwfc.x`).

Implementation here: `src/scoring/phi_elf.py` (my own percolation code: 6-connectivity on the
FFT grid, periodic wrap-vector rank, bisection on the ELF threshold; not validated against the
paper's table — UNVERIFIED agreement); `src/scoring/belli_probe.sh` times the chain.

**Measured on H3S (200 GPa cell, ecut 60 Ry, 2 ranks):** scf (k 12³) 2.8 s, ELF (`pp.x`)
1.8 s, nscf (k 16³) 3.9 s, ILDOS around E_F (`pp.x plot_num=10`) 1.8 s, φ 0.3 s: **~11 s in
total**. Result: φ = 0.729, H_f = 0.75, H_DOS = 0.664, Tc_φ = 273 K, against the DFPT-level
~200 K (bench Allen–Dynes 210–229 K at σ = 0.01–0.03) — within the paper's ±60 K band only at
its upper edge. Two deviations from the paper: (i) SG15 pseudopotentials have no atomic
wavefunctions, so `projwfc.x` fails ("Cannot project on zero atomic wavefunctions"); H_DOS
here is the fraction of the E_F ± 0.136 eV integrated LDOS inside H Voronoi cells, not a
Löwdin projection (UNVERIFIED equivalence). Using pseudopotentials with PSWFC (e.g.
PseudoDojo) would allow the paper's projection. (ii) The ELF is on the 25³ FFT grid of a
60 Ry cutoff; φ's dependence on grid density was not tested. Cost for a 10-atom
low-symmetry cell **[est.]**: SCF/nscf scale roughly with N_k × N_PW × N_bands, so ×~50–200,
i.e. **~5–30 min on 2 ranks**.

Use: ordering of hydrogen-containing candidates only. It is blind to dynamical stability
and to the phonon side of λ (H3S vs H3Se is the paper's own counter-example raised by
Novoa et al. 2025 **[cited]**), and the fit is anchored on high-pressure hydrides.

### 3.4 Allen–Dynes with an ML ω_log

ω_log can be estimated cheaply: either from an ML model (Alexandria ω_log MAE 20 K
in-distribution **[cited]**, but 545 K vs 1250 K for H3S **[measured]**) or from the MatterSim
phonon DOS F(ω) as a stand-in for α²F(ω)/λ (exact only if coupling is mode-independent;
for H3S the reference α²F ω_log, 1250 K, sits on the H-dominated part of F(ω)). The hard
part of Tc is λ, not ω_log: at fixed ω_log = 1250 K, the Allen–Dynes Tc for λ = 1.5 / 2.0 /
2.5 is 144 / 195 / 237 K (μ* = 0.13, ω₂ = 1478 K, computed with `a2f_summary.allen_dynes`). An ML ω_log therefore adds little unless λ is also known. Allen–Dynes itself is
only a fit to Eliashberg solutions; it is used here for ordering, while T3 uses the full
isotropic Eliashberg bounds.

---

## 4. Recommended stack and decision rules

Pipeline order (cheapest first). Every tier writes its numbers to the candidate's ledger
with its tier label; no tier overwrites another's numbers.

1. **T0 — MatterSim** (exists): relax, E_hull, finite-displacement phonons. DEFER if
   E_hull > threshold or ML phonons show imaginary modes beyond tolerance (`rtt.mlff`
   policy). Not a rejection.
2. **T1 — ordering scores** (< 1 min each):
   - T1a: composition heuristics (H fraction, light fraction,
     √(m_H/M_eff), N(E_F) if available).
   - T1c: Belli Tc_φ for hydrogen-containing candidates (needs one small SCF + ELF + projwfc).
   - T1b: ALIGNN/BETE-NET numbers **recorded but weight 0 for hydrides / light-element
     candidates**; nonzero weight only inside their training domain.
   - Output: a rank. **Allowed action: reorder; DEFER only by the budget planner (i.e. "not
     reached yet"), never by a threshold on an ML Tc.**
3. **T2 — coarse DFPT** (q 3³–4³, k_fit ≥ 24³-equivalent density, ecut ~2/3 of the converged
   value, σ-plateau value reported, 2 ranks per job so several run side by side):
   - Output: α²F_coarse, λ, ω_log, A at several σ, plus the min frequency of the coarse
     phonons (a DFPT-level stability hint).
   - **Allowed DEFER rule:** defer from T3 only if A_coarse(σ) < A_target / 2 for **every**
     σ in the plateau (factor-2 margin, where A_target ≥ T* · 4.33148 for target T* and
     A_target is whatever the theory documents set; T* = 300 K gives 1299 K). Justification: in
     the one measured case, coarse A ranged from −50 % (q 2³ only, which is excluded) to +56 %
     (q 3³ at σ = 0.01 Ry). For q ≥ 3³ and σ ≥ 0.03 Ry it stayed within −1 % … +14 %. Deferral
     is only hurt by *under*-estimates of A, and a factor-2 margin covers the worst one seen
     (q 2³) as well. This margin must be re-checked as more T3
     references (Nb, MgB2, Al, Pb, H3S) are compared with their own coarse runs.
   - Candidates with imaginary coarse-DFPT modes are not deferred on that basis alone
     (coarse grids create spurious soft modes); they are queued for a T3 phonon-only
     check first.
4. **T3 — converged DFPT + certification (`rtt.core.certify`)**: the only tier that issues REJECTED_PROVEN,
   CERTIFIED_CANDIDATE or UNDECIDED.

Never allowed below T3: REJECTED_PROVEN, EXCLUDED_*, CERTIFIED_CANDIDATE, or any statement
"X is not a superconductor above T". A DEFERRED candidate is reported as "not computed",
with the tier and rule that deferred it, so a reviewer can re-queue it.

---

## 5. Open problems and what I could not verify

- The T2 error statement rests on **one** material (H3S). It must be repeated on the other
  bench materials (Al, Pb, Nb, MgB2) and on at least one ambient-pressure hydride before the
  factor-2 DEFER margin is trusted. A one-command repeat is `coarse_dfpt_probe.sh` with a
  different cell (needs small edits: species block and positions are hard-coded for H3S).
- BETE-NET was not run (legacy torch 1.10 / e3nn 0.4.2 stack would need its own env; ~2 GB
  of weights per variant). Its training set (818 materials, no stated hydrides) makes the same
  OOD failure likely, but that is an inference, not a measurement.
- BEE-NET weights and license are not available; the Cerqueira et al. ML code release was not
  found.
- ALIGNN JARVIS-paper MAE numbers and the identity of the Alexandria superconductivity
  training set were not read from the papers (UNVERIFIED).
- `phi_elf.py` is untested against Belli et al.'s published φ values; H_DOS uses a Voronoi
  ILDOS split instead of the paper's projected DOS; ELF grid dependence untested.
- The cause of the spurious c4 instability (60 Ry cutoff vs k 12³ / k_fit 24³ vs the
  0.02 Ry scf smearing) was not isolated. One more run (ecut 90, k 12³, q 6³) would separate
  them; it was not done because of the CPU budget.
- The DEFER margin is tied to the target A. Its value (≈ 1300 K for T* = 300 K) must come
  from the theory documents, not from this note.
- No universal e-ph ML model was found; the search may have missed 2026 preprints.
- The 10-atom cost estimate in §3.1 is a scaling argument, not a measurement.

---

## 6. Scripts (all in `src/scoring/`)

| File | What it does | How it was run |
|---|---|---|
| `alignn_screen.py` | ALIGNN 2.0 pure-torch checkpoints (Tc, λ, ω_log, hydride Tc) on reference or user POSCARs; CPU only | alignn + jarvis-tools installed with `pip --target` into the scratch dir (main venv untouched), `PYTHONPATH` set, 2 threads, nice 19 |
| `coarse_dfpt_probe.sh` | Coarse QE DFPT chain with per-step timings (`<workdir> <ecut> <kscf> <kfit> <nq> [ph_timeout]`) | 2 MPI ranks, nice 19, `setsid nohup` |
| `a2f_summary.py` | λ, ω_log, ω₂, A, A/4.33148, Allen–Dynes Tc from matdyn `a2F.dosN` | float, informational |
| `phi_elf.py` | Belli networking value φ from an ELF cube (+ Tc_φ if H_DOS given) | float, informational |
| `belli_probe.sh` | scf → ELF → nscf → ILDOS at E_F → φ and Tc_φ for H3S, with timings | 2 MPI ranks, nice 19 |

Sources (all accessed 2026-10-28): arXiv 2401.16611, 2503.20005, 2205.00060, 2307.10728,
2312.12694, 2105.02296, 2508.21281, 2507.06749, 2302.00439, 2001.09852, 2211.03265,
1709.02727, 2609.19487; github.com/henniggroup/BETE-NET; github.com/henniggroup/BEE-NET;
github.com/atomgptlab/alignn (`alignn/pretrained.py`); pages.nist.gov/jarvis_leaderboard;
PMC8446067 (Belli et al. full text); Crossref records for doi:10.1038/s41524-024-01475-4
and doi:10.1038/s41524-022-00933-1; b2find record of doi:10.24435/materialscloud:qv-bq.
