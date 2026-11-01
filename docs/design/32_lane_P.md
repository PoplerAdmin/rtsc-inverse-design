# 32 — Lane P: deterministic prototype substitution in the A₂MH₆ family (campaign 2)

Agent P2-lane, 2026-10-30. Code: `src/rtt/source/prototype.py` (lane module),
`src/rtt/configs/lane_p.toml` (lane config), `tests/front/test_lane_p.py` (tests), driver
`scripts/campaign2_lane_p.py` (S1–S3), `scripts/campaign2_lane_p_report.py` (tables, S4
selection, S5 queue), `scripts/lane_p_band_count.py` and `scripts/lane_p_rigid_band.py` (§1
checks), `scripts/lane_p_barhh6_magnetism.py` (§6). Tables: `runs/campaign2_laneP/`. Ledger
(not in the repo): `~/rtsc_targeted_work/P2-lane/ledger.sqlite` (run id `campaign2-laneP`); the
campaign-1 ledger was only read.

Labels: **(meas.)** read from a run of this agent; **(cited)** with an id; **(est.)** my
extrapolation with the assumption stated; **UNVERIFIED** recalled, not checked.

## 0. Summary

1. **Rule.** Per [MH₆] complex, N = g(M) + 6 + Σox(A)/n(M); N = 18 is the closed shell.
   - Mg₂IrH₆ has **N = 19**: one electron per complex in the e_g* band, i.e. *electron*-doped
     relative to the 18-electron Mg₂FeH₆-type insulators. It is "hole-doped" only relative to
     Mg₂PtH₆ (N = 20) or Mg₂IrH₅.
   - The PBE bands confirm this: a 3.3 eV gap above the 18-electron manifold, with E_F 1.07 eV
     into the next band.
   - E_F already sits on the e_g* DOS peak (rigid band: ±½ e per complex lowers N(E_F) 2–3×).
2. **Lane P.** `src/rtt/source/prototype.py`, `configs/lane_p.toml`, 7 test functions (10 cases), lane `P` in the
   ledger. It covers four prototypes (A2MH6, rock-salt AA′MH6, A3A′(MH6)2, A4MM′H12). From 689
   enumerated cells, **347 were planned**: 62 ternaries, 165 AA′ at N = 19, and 120 half-electron
   alloys.
3. **S1–S3.**
   - 347 PASS at S1 and at S2. S3 made no removal: 40 PASS at the rank cut, 307 in the rank
     queue.
   - 311 novel structures and 266 novel compositions vs Alex-MP.
   - 253 at E_hull,ML ≤ 0.10 eV/atom.
   - 2.6 GPU-h.
4. **S3 honours the doc-15 rule**: MatterSim phonons cannot remove a hydride. In this run
   MatterSim distorted the Mg₂IrH₆ control to P1 (the fake T1g mode) and the control still passed.
5. **S4 top 10 by S3 score**: all Li-containing rock-salt AA′MH₆ at N = 19. The S3 score
   saturates; see §2.2.
   - **10/10 PASS**: metallic, non-magnetic, 1 sheet, E_form −0.37 to −0.49 eV/atom.
   - **N(E_F) 0.14–0.20 states/eV/atom, against 0.386 for Mg₂IrH₆.** No cell beats the
     control's N(E_F).
   - Best: LiScIrH₆, LaLiCoH₆, LiScCoH₆. LiScIrH₆ has a rigid-band DOS peak of 0.43 at
     N = 19.5.
6. **S5** was not run. The queue is ranked by N(E_F) with cost estimates: 26–52 core-h at 60 Ry
   and 40–80 core-h at 80 Ry for the leader; all 10 together 450–900 core-h at 80 Ry. Run at
   most LiScIrH₆ and LiScCoH₆. 18-atom alloys are out of the 48 h cap.
7. **BaRhH₆**: no FM or AFM moment survives (PBE, two AFM bipartitions, MV 0.02 and
   0.005 Ry); SOC-AFM decays to ≤ 0.04 µB/Rh (not converged). An SOC restart is running in the
   background.
8. **Honest expectation unchanged** (R30 §2.2): nothing in lane P is above the Mg₂IrH₆ anchor at
   the S4 level.

---

## 1. The substitution rule and its electron count

### 1.1 Counting

The parent is the K₂PtCl₆ (Fm-3m, sg 225) structure of Mg₂FeH₆ and Mg₂IrH₆: isolated [MH₆]
octahedra on an fcc lattice (4a), A cations on the simple-cubic 8c sublattice, H on 24e (x,0,0).
Count each H as hydride H⁻ (a two-electron σ donor), the A cations at their group valence, and
M by its group number g(M) (dⁿs⁰ count of the neutral atom). Per [MH₆] complex:

    N = g(M) + 6 + q,        q = Σ ox(A) / n(M)        (complex charge [MH₆]^q−)

This is the textbook 18-electron count of an octahedral complex. M has oxidation state 6 − q,
so its d count is n_d = g(M) − (6 − q), and N = 12 (six H⁻ pairs) + n_d = g(M) + 6 + q. For a
low-spin octahedron the six
M–H σ-bonding orbitals (a₁g + t₁u + e_g) and the non-bonding t₂g hold exactly 18 electrons:

| N | occupation | example | expectation |
|---|---|---|---|
| 18 | σ⁶ … t₂g⁶, e_g*⁰ (closed shell) | Mg₂FeH₆, K₂PtH₆; campaign-1 SrPtH₆, BaPtH₆ | insulator |
| 17 | one t₂g hole (d⁵ low spin) | Na₂IrH₆, BaRhH₆-like [M^IV H₆]²⁻ | metal; moment possible (§5) |
| 19 | one electron in the M–H σ* (e_g) band, ¼ filled | **Mg₂IrH₆**, Ba₂IrH₆, Ca₂RhH₆ | metal |
| 20 | e_g*² (½ filled) | Mg₂PtH₆, Ba₂PdH₆ | metal or high-spin moment |
| 18.5 / 19.5 | half an electron less / more than Mg₂IrH₆ | Mg₃Li(IrH₆)₂, Mg₂(Ir½Pt½)H₆ | metal |

The total valence count of a cell is n(M)·N (the same counting; H = 1). An odd count forces a
non-magnetic metal (Kramers-degenerate bands cannot be filled by an odd number of electrons).
This is the parity feature of R30 F1(c), which separated every campaign-1 S4 outcome correctly.

### 1.2 Mg₂IrH₆ is electron-doped relative to the 18-electron insulators, not hole-doped

The brief describes Mg₂IrH₆ as "hole-doped relative to the 18-electron Mg₂FeH₆-type
complexes". By the count above it is the opposite: N = 9 + 6 + 4 = 19, i.e. **one electron per
complex beyond the closed shell**, in the antibonding M–H e_g* band. "Hole-doped" is correct only
relative to N = 20 (Mg₂PtH₆, whose e_g* band holds two electrons) or relative to the
five-coordinate Mg₂IrH₅ parent, where adding an H⁻ ligand creates a new bonding level that takes
two electrons while the H atom brings one (net: one hole in the Ir d manifold, d⁸ → d⁷).
Both statements describe the same state; the band check below decides which bands hold E_F.

### 1.3 Band-level check (meas.)

`scripts/lane_p_band_count.py` on the campaign-1 S4 nscf of the Mg₂IrH₆ control (PBE, SG15,
60 Ry, 56 irreducible k, 43 valence electrons in the pseudopotentials = 24 semicore Mg 2s2p / Ir
5s5p + 19). The 18-electron closed shell is 42 electrons = 21 bands:

| band | range (eV) |
|---|---|
| 19 | 5.14 – 6.20 |
| 20 | 5.99 – 6.59 |
| 21 (top of the 18-e manifold) | 5.99 – **6.67** |
| 22 (bottom of e_g*) | **9.94** – 13.23 |
| 23 | 12.01 – 13.97 |
| E_F | **11.00** |

There is a 3.27 eV gap above band 21, and E_F lies 1.07 eV inside band 22. Mg₂IrH₆ is a
closed-shell 18-electron framework plus one electron per complex in the bottom of the e_g*-derived
conduction band. The rule's N is therefore the band-filling coordinate. The gap also explains the
campaign-1 H1 exclusions of the N = 18 complexes SrPtH₆ and BaPtH₆ (R30 F1(c): the eight H1
exclusions have PBE gaps of 1.7–3.7 eV).

### 1.4 Where E_F sits on the DOS: rigid-band scan (meas.)

`scripts/lane_p_rigid_band.py` on the same run's tetrahedron DOS (k 11³, dos.x window E_F ± 1 eV,
0.1 eV box average), with N(E) = 19 + ΔIDOS per complex:

| N per complex | E − E_F (eV) | DOS (states/eV/atom) |
|---|---|---|
| 18.50 | −0.28 | 0.127 |
| 18.75 | −0.11 | 0.203 |
| 18.90 | −0.03 | 0.304 |
| **19.00** | 0.00 | **0.386** |
| 19.15 | +0.05 | 0.407 (max) |
| 19.25 | +0.08 | 0.363 |
| 19.50 | +0.17 | 0.195 |
| 19.75 | +0.46 | 0.047 |

Mg₂IrH₆ already sits on the e_g* DOS peak: the maximum is +0.15 e per complex away and only 5 %
higher. In a rigid band, ±½ electron per complex **lowers** N(E_F) by a factor 2–3. Ordered
half-electron alloys (N = 18.5, 19.5) are therefore not expected to raise λ through N(E_F). Their
only route is a changed band shape or stiffer/softer modes, which a rigid band cannot show.
Isoelectronic (N = 19) substitutions keep E_F on the peak while changing the cation mass, size
and the M–H bond. This motivates the priority order below. (Caveat: one k mesh, MV-smeared
scf; the σ-dependence of N(E_F) in doc 15 §3 is the same feature.)

### 1.5 Enumeration (deterministic; `configs/lane_p.toml`)

- **A** ∈ {Li, Na, K, Rb, Cs, Mg, Ca, Sr, Ba, Sc, Y, La}.
- **M** ∈ {Mn, Re, Fe, Ru, Os, Co, Rh, Ir, Ni, Pd, Pt}. Tc is excluded (radioactive).

Prototypes (`rtt.source.prototype`):

| key | cell | atoms | ordering |
|---|---|---|---|
| A2MH6 | Fm-3m primitive | 9 | parent |
| AA'MH6 | F-43m primitive | 9 | the A simple-cubic sublattice split into its two fcc sublattices (rock-salt order of A/A′, 8c → 4c + 4d) |
| A3A'M2H12 | (a/√2, a/√2, a) two-complex cell | 18 | one A of four is A′ (1 : 3); only heterovalent pairs, half-integer N |
| A2MM'H12 | same cell | 18 | M/M′ alternate on (001) planes of the fcc M lattice (L1₀-type); only heterovalent pairs |

Starting geometry: a from the mean A–H contact at x_H = ¼ (Shannon CN8 radius + 1.46 Å
effective H⁻ radius, fitted so that Mg₂IrH₆ gets a = 6.65 Å; radii recalled, **UNVERIFIED** to
0.01 Å, used only as a relax start) and x_H = d(M–H)/a with d = 1.58 Å (3d) or 1.70 Å (4d/5d).
MatterSim relaxes cell and positions at S3.

Targets: ternaries at N ∈ {17, 19, 20}, alloys at N ∈ {18.5, 19, 19.5}. N = 18 is not enumerated:
those cells are closed-shell insulators by §1.3. Priority:
1. |N − 19|;
2. |N_parent − 19|, where the parent is the majority-A ternary, so a small doping of a
   19-electron analogue comes first;
3. distance from Mg₂IrH₆ (non-Mg A fraction + non-Ir M fraction);
4. cell size.

Caps per prototype: A2MH6 200, AA'MH6 200, A3A'M2H12 60, A2MM'H12 60 (GPU budget; judgment).
The full enumeration has 689 cells (62 / 165 / 312 / 150). The plan keeps **347**: all 62
ternaries, all 165 rock-salt AA′ (all N = 19), and the first 60 of each 18-atom alloy family.

Provenance per candidate (`Origin.conditioning`):
- lane P, generator `rtt.source.prototype/1`;
- prototype, parent, substitution {A, A′, M, M′}, formula;
- N, parent N, filling class, Δ vs Mg₂IrH₆;
- cell valence count and parity.

Flags: `PROTOTYPE_SUBSTITUTION`, `PROTOTYPE=…`, `N_PER_COMPLEX=…`, `ODD_VALENCE_PARITY`.
The lane is a *proposal* rule. It removes nothing from any other lane (R2). Ledger: `P` was
added to `records.LANES` and to the `candidates.lane` CHECK in `ledger/schema.sql` (new ledgers
only; existing ledgers keep their table definition).

---

## 2. S1–S3 on the campaign-2 lane-P ledger (meas.)

Driver: `scripts/campaign2_lane_p.py --phase all --rank-cut 40`. The stages are the standard
ones: `ingest_raw`/`prefilter` (S1), `run_s2` (validity, dedupe, Alex-MP novelty) and
`run_s3` + `rerank_s3` (S3). Two run options differ from campaign 1, and neither changes a
decision rule:
- **`phonons=False`.** MatterSim phonons are dropped from the score for every H compound
  (`rank_metrics`), and `decide_s3` never reads them, so they are not computed.
- **An entry-level hull cache (`hull_fn`).** Each Alex-MP reference entry is relaxed with the
  same MatterSim model and settings exactly once and reused across chemical systems. The shared
  per-system cache is read but never written. The numbers are the same as `build_self_hull`'s:
  every reference is relaxed independently.

The Mg₂IrH₆ control (lane C) ran through every stage alongside the lane.

### 2.1 Funnel

| stage | lane P | notes |
|---|---|---|
| S1 | 347 PASS / 347 | every element has an SG15 file; cells ≤ 18 atoms |
| S2 | 347 PASS, 0 INVALID, 0 DUPLICATE | |
| novel structure vs Alex-MP | 311 / 347 | the 36 non-novel are all ternary A2MH6 (Fm-3m entries already in Alexandria; one `mp-` entry, Na₂IrH₆) |
| novel composition vs Alex-MP | 266 / 347 | by prototype: A2MH6 13/62, AA'MH6 133/165, A3A'M2H12 60/60, A2MM'H12 60/60. Known AA′ compositions are hexagonal Alexandria polymorphs, e.g. LaLiIrH₆ P-6m2 at E_hull,DFT 0.036, CaLiPtH₆ P312 at 0.098 eV/atom |
| S3 | 0 SCREENED_OUT; 40 PASS (rank cut 40) + 307 DEFERRED(rank_cut) | control: PASS (exempt) |
| E_hull,ML-self ≤ 0.05 / ≤ 0.10 eV/atom | 208 / 253 of 303 with a hull | |
| no ML hull | 44 | every one contains Ru or Re (next item) |
| DFT-reference hull ≤ 0.10 | 246 | |
| `ML_BELOW_DFT_HULL` (E_hull,DFT-ref < −0.05) | 44 | R30 F4: treat as likely MatterSim over-binding |
| MatterSim relax keeps the prototype space group (0.1 Å) | 337 / 347 | 7 AA′ cells and 3 alloys lose symmetry |

By prototype (E_hull,ML-self: median, number ≤ 0.10 eV/atom):

| prototype | n | with hull | median (eV/atom) | ≤ 0.10 |
|---|---|---|---|---|
| A2MH6 | 62 | 52 | 0.006 | 48 |
| AA'MH6 | 165 | 150 | 0.033 | 122 |
| A3A'M2H12 | 60 | 57 | 0.008 | 49 |
| A2MM'H12 | 60 | 44 | 0.019 | 34 |

GPU: 2.57 GPU-h of recorded S3 compute (348 cells, including about 8 000 newly relaxed hull
references). This is 2.5× R30's ≈ 1 GPU-h estimate. The hull references dominate, especially
the 30–64-atom reference cells that do not converge in 1000 FIRE steps. Peak GPU memory was kept
below 15 GB with a `set_per_process_memory_fraction(0.6)` cap. An earlier two-worker attempt
reached 23.9 GB and was killed; its 73 `started` S3 tier_runs remain as orphan rows in the
ledger.

### 2.2 S3 findings that matter for campaign 2

1. **The S3 score saturates.** The hull term is clipped at E_hull ≤ 0, and the H, mass and
   metallicity terms are equal across the family. 67 lane-P cells tie at 0.8211 and 9 at 0.830.
   The top of the ranking is decided by the *light-element fraction* (Li counts as light), so
   every top-10 cell is a Li-containing rock-salt AA′ cell. Past the ties, the rank cut 40 falls
   by candidate id. Fix (owner of `score/s3.py`): do not clip below 0, or add a tie-breaker such
   as the DFT-reference hull or the R30 18-electron distance. Ranking only.
2. **Ru and Re have no ML hull.** Their elemental Alex-MP references do not converge under the
   MatterSim relax, so `ehull_self` returns None ("incomplete hull") for all 44 Ru/Re cells.
   They PASS with `NO_ML_HULL` and a neutral hull term. That is lenient and correct under R2,
   but they are ranked blind. Fix: accept the elemental reference at its last step when
   ΔE < 1 meV/atom, or take the elemental energy from a fixed-cell relax.
3. **A known Fm-3m ternary reproduces its own reference.** Its Alexandria entry is in the
   reference set, re-relaxed by the same model, so its E_hull,ML is that entry's hull distance
   (20 of the 30 known cells with a hull sit at |E_hull| < 0.005 eV/atom). This is not
   independent evidence of stability.

## 3. Does S3 honour "MatterSim phonons never reject a hydride"? (task 2)

**Yes, in code and in this run.**
- **Code.** `decide_s3` (`src/rtt/mlff/screen.py`) reads only `e_hull_ml_self_eV`, `error` and
  the relax-converged flag. Its docstring and code contain no phonon access. The only removal is
  `SCREENED_OUT(H7-ml)`, and outside the hydride calibration domain that becomes
  `DEFERRED(uncalibrated_window)`. `rank_metrics` deletes the `phonons` block for every
  composition with H/Li/Be before scoring, so MatterSim phonons cannot reach a hydride through
  the rank cut either (RC1 F12). `tiers.toml [S3] ml_phonons_remove = false`. The only other
  `ml_phonons` caller is `mlff/calibrate.py` (ZPE for the window calibration, not a decision).
- **This run.** Phonons were not computed. There were 0 SCREENED_OUT. The control passed S3.
- **The MatterSim artefact is visible anyway.**
  - The **control** started at the Alex-MP geometry (x_H = 0.2606). MatterSim relaxed it to
    **P1** (sg 1 at 0.1 Å) with an unconverged relax, at E_hull,ML −0.016 eV/atom.
  - The **lane-P copy** started from the idealised prototype (x_H = 0.2556). It stayed **Fm-3m**
    at E_hull,ML +0.001 eV/atom.
  - The ≈ 17 meV/atom difference is MatterSim rotating the IrH₆ octahedra along its spurious
    T1g instability (doc 15 §2). The rotation *lowers* the ML energy. It cannot remove the
    control; it only makes its ML hull distance look better.
- **Downstream protection.** S4 does not inherit the distortion. `dft_structure` uses the
  candidate's own cell symmetrised to its label, not the MatterSim-relaxed cell, and QE's
  vc-relax keeps the symmetry.

One residual risk, not a violation: a hydride whose *prototype* cell MatterSim distorts during
relaxation gets an S3 E_hull for the distorted cell. That is lenient (lower E), but the
recorded `relaxed_structure` is not the cell S4 computes. 10 lane-P cells changed space group
in S3 (§2.1).

## 4. S4 (DFT electronic, H1/H2) on the top 10 by S3 score (meas.)

**Selection.** The top 10 lane-P cells by stored S3 score: the 9 cells at 0.830 and the single
cell at 0.8277. This set is unambiguous, because the next score is 0.8261. All ten are S3 PASS,
all are metallic by count (N = 19, odd cell parity), and all are rock-salt-ordered **AA′MH₆
(F-43m, 9 atoms) with A = Li**. That is the score saturation of §2.2(1): the light-element term
picks Li.

**Run.** `rtt step S4 <cid> --ranks 2` on the lane-P ledger; profile `s4` (60 Ry, vc-relax,
scf, nscf, tetrahedron DOS, spin check, ELF, H7-heur E_form). Wall time was 04:29–07:15 on 2
ranks, about 5.5 core-h including the new elemental references (Li, La, Y, Sc, Co, Rh, Ir,
Ba, Pt, Ca). With the references cached, a run takes 3.5–7 min.

Columns `gap` and `E_F − e_g* bottom` come from `scripts/lane_p_band_count.py`; `DOS max ±½ e`
from `scripts/lane_p_rigid_band.py`. Raw values are in `runs/campaign2_laneP/s4_bandcheck.json`.

| cid | formula | S4 | N(E_F) tet (st/eV/atom) | H_DOS | moment | sheets | E_form (eV/atom) | a_conv DFT (Å) | gap above 18-e manifold (eV) | E_F − e_g* bottom (eV) | DOS max within N ± ½ (at N) | E_hull ML self / DFT-ref | novel comp. |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| c-43d572803f9d | LiScIrH₆ | PASS (`NEF_MESH_SENSITIVE`) | 0.201 | 0.40 | 0 | 1 | −0.44 | 6.652 | 3.13 | 0.95 | **0.425 (19.48)** | −0.012 / −0.068 | yes |
| c-5ff054e14807 | LaLiCoH₆ | PASS | 0.200 | 0.39 | 0 | 1 | −0.45 | 6.966 | 1.48 | 0.70 | 0.228 (18.51) | −0.025 / +0.018 | no (P-6m2 known) |
| c-7c5823b24da2 | LiScCoH₆ | PASS | 0.197 | 0.33 | 0 | 1 | −0.43 | 6.336 | 1.55 | 0.92 | 0.237 (18.71) | −0.044 / +0.084 | yes |
| c-b91240dd4cd0 | LaLiIrH₆ | PASS | 0.195 | 0.49 | 0 | 1 | −0.48 | 7.239 | 2.99 | 0.84 | 0.212 (18.59) | −0.063 / −0.141 | no (P-6m2 known) |
| c-b16887414468 | LaLiRhH₆ | PASS | 0.178 | 0.48 | 0 | 1 | −0.46 | 7.174 | 2.58 | 0.90 | 0.200 (18.59) | −0.018 / −0.098 | no (P-6m2 known) |
| c-7a746dec73c4 | LiYCoH₆ | PASS | 0.176 | 0.39 | 0 | 1 | −0.47 | 6.633 | 1.54 | 0.87 | 0.212 (18.62) | −0.032 / −0.016 | yes |
| c-33f78121fcea | BaLiPtH₆ | PASS | 0.163 | 0.63 | 0 | 1 | −0.37 | 7.709 | 3.35 | 0.94 | 0.183 (18.51) | −0.000 / +0.094 | yes |
| c-3ba28b0824f3 | LiYRhH₆ | PASS | 0.156 | 0.48 | 0 | 1 | −0.47 | 6.852 | 2.69 | 1.11 | 0.191 (18.69) | −0.012 / +0.052 | yes |
| c-9bb28830217f | LiYIrH₆ | PASS | 0.154 | 0.50 | 0 | 1 | −0.49 | 6.913 | 3.12 | 1.07 | 0.194 (18.65) | −0.046 / +0.021 | yes |
| c-6c2d8da02768 | CaLiPtH₆ | PASS | 0.144 | 0.59 | 0 | 1 | −0.41 | 7.030 | 2.93 | 1.32 | 0.164 (18.74) | +0.005 / +0.076 | no (P312 known) |
| *reference* c-9a8fa100a45b | Mg₂IrH₆ (campaign-1 S4, same profile) | PASS (`NEF_MESH_SENSITIVE`) | **0.386** | 0.59 | 0 | 1 | −0.24 | – | 3.27 | 1.07 | 0.407 (19.15) | – | known |
| *reference* c-539acaa4a6e7 | Ba₂IrH₆ (campaign-1 D) | PASS | 0.148 | 0.60 | 0 | 1 | – | – | – | – | – | – | known |

What S4 shows:

1. **10/10 metallic and non-magnetic. The rule's band picture holds in every cell.** Each has a
   1.5–3.4 eV gap above the 18-electron manifold, with E_F 0.7–1.3 eV inside the e_g*-derived
   band. Prediction and outcome agree 10/10, as the parity rule predicted. No H1 or H2
   exclusion.
2. **No top-10 cell reaches the Mg₂IrH₆ N(E_F).** They sit at 0.14–0.20 states/eV/atom against
   0.386, i.e. 37–52 % of the control, with Ba₂IrH₆ at 0.148 in the same range.
   - The e_g* band is wider and E_F is not on its peak. In 9 of 10 the rigid-band DOS maximum
     lies on the *hole* side (N ≈ 18.5–18.7) and is only 5–20 % higher.
   - Mg₂IrH₆ is special in this family because its E_F lies on the e_g* peak (§1.4).
3. **One rigid-band lead: LiScIrH₆.** Its DOS peak is at N ≈ 19.5, with 0.425 states/eV/atom,
   above Mg₂IrH₆'s 0.39.
   - Reaching it needs +½ e per complex (q = 4.5), e.g. (Li½Mg½)ScIrH₆ = LiMgSc₂(IrH₆)₂ in a
     doubled cell. That cell has three A species, which none of the four prototypes enumerates.
   - The S4 tag `NEF_MESH_SENSITIVE` on this cell says the peak is sharp. Rigid band only; not a
     result.
4. **All E_form are −0.37 to −0.49 eV/atom.** That is deeper than Mg₂IrH₆ (−0.24), but E_form is
   against the elements, not the hull. Four of the ten have Alexandria hexagonal polymorphs at
   E_hull,DFT 0.04–0.10 eV/atom (§2.1), and the F-43m cell is not shown to be lower. S7a would
   decide.
5. **What to read from this.** Isoelectronic Li-based rock-salt alloys are good metals but do
   not beat the control on the one quantity S4 measures that λ scales with. The campaign-2 S3
   score should not promote them by Li content alone.

## 5. Ranked S5 queue (not run; the coarse screen is being recalibrated, doc 16)

**Cost model (est.).** The anchor is Mg₂IrH₆ S5 ≈ 30 core-h at 60 Ry (campaign 1, measured;
9 atoms, Fm-3m, 43 SG15 electrons, V₀ = 73.6 Å³). It is scaled by:
- (ecut/60)^1.5;
- × V/V₀ (plane waves);
- × nat/9 (perturbations);
- × n_el/43 (bands);
- × 48/|G| (irreducible k+q).

The F-43m cells have |G| = 24 and no inversion. The *low* end credits time reversal (effective
48); the *high* end does not. The true factor is probably in between: per-q k+q lists grow
without inversion. V is the S4 DFT-relaxed cell. Model and data: `scripts/campaign2_lane_p_report.py`,
`runs/campaign2_laneP/s5_queue.json`.

**Ranking (ranking only).**
1. N(E_F) from the S4 tetrahedron DOS, descending. λ ∝ N(E_F)⟨I²⟩/(M⟨ω²⟩), and N(E_F) is the
   only factor S4 measures.
2. H_DOS.

| rank | cid | formula | N(E_F) | H_DOS | V (Å³) | core-h at 60 Ry | core-h at 80 Ry | wall on 4 ranks at 80 Ry (h) | note |
|---|---|---|---|---|---|---|---|---|---|
| 1 | c-43d572803f9d | LiScIrH₆ | 0.201 | 0.40 | 73.6 | 26–52 | 40–80 | 10–20 | `NEF_MESH_SENSITIVE`; rigid-band peak at N = 19.5 (§4 item 3); new composition |
| 2 | c-5ff054e14807 | LaLiCoH₆ | 0.200 | 0.39 | 84.5 | 30–59 | 46–91 | 11–23 | 3d M (lighter, stiffer M–H?); hexagonal polymorph known |
| 3 | c-7c5823b24da2 | LiScCoH₆ | 0.197 | 0.33 | 63.6 | 22–45 | 34–69 | 9–17 | smallest cell, cheapest; DFT-ref hull +0.084 |
| 4 | c-b91240dd4cd0 | LaLiIrH₆ | 0.195 | 0.49 | 94.8 | 33–67 | 51–102 | 13–26 | `ML_BELOW_DFT_HULL` (−0.141) |
| 5 | c-b16887414468 | LaLiRhH₆ | 0.178 | 0.48 | 92.3 | 32–65 | 50–100 | 12–25 | `ML_BELOW_DFT_HULL` (−0.098) |
| 6 | c-7a746dec73c4 | LiYCoH₆ | 0.176 | 0.39 | 73.0 | 26–51 | 39–79 | 10–20 | |
| 7 | c-33f78121fcea | BaLiPtH₆ | 0.163 | 0.63 | 114.5 | 40–80 | 62–124 | 15–31 | highest H_DOS |
| 8 | c-3ba28b0824f3 | LiYRhH₆ | 0.156 | 0.48 | 80.4 | 28–56 | 43–87 | 11–22 | |
| 9 | c-9bb28830217f | LiYIrH₆ | 0.154 | 0.50 | 82.6 | 29–58 | 45–89 | 11–22 | |
| 10 | c-6c2d8da02768 | CaLiPtH₆ | 0.144 | 0.59 | 86.9 | 31–61 | 47–94 | 12–23 | MatterSim relax distorted it to Cm (S4 used F-43m) |

All 10 together cost 300–600 core-h at 60 Ry and 450–900 core-h at 80 Ry. The campaign-2 S5
envelope is ≈ 700 core-h for ≤ 15 runs (R30 §3.3).

**Recommendation.** Run at most 2–3 of these after the doc-16 profile lands:
- **#1 LiScIrH₆**, the best N(E_F) and a 5d M like the control;
- **#3 LiScCoH₆**, the cheapest, and a 3d-M test of the "stiffer M–H bond" route of R30 §2.2(2);
- optionally **#2 LaLiCoH₆**.

The prior expectation is A below Mg₂IrH₆, because N(E_F) is ≈ ½ of the control's. A result
above the control would therefore be informative, and one below it would confirm that the
N(E_F) peak, not the A-site chemistry, sets the Mg₂IrH₆ coupling. The R30 anchors outside lane P
(Mg₂IrH₆ corrected-profile re-run, Ba₂IrH₆ c-539acaa4a6e7) stay ahead of this queue.

**18-atom alloy cells** (A3A'M2H12, A2MM'H12; none reached S4 here). By the same model (P4/mmm
|G| = 16 or P-4m2 |G| = 8, nat 18, V ≈ 2 V₀, n_el ≈ 79), a cell like LiMg₃Ir₂H₁₂ costs
≈ 650–1350 core-h at 60 Ry (est.). That is beyond the 48 h × 4-rank cap (192 core-h). Half-electron
alloys therefore need either a reduced S5 (q 2³, or Γ + X only) or a virtual-crystal /
rigid-band proxy. §1.4 also predicts that a rigid band at N = 18.5 or 19.5 *lowers* N(E_F) for
Mg-based parents, which argues against spending S5 on them.

## 6. BaRhH₆ magnetism re-check (R30 F6; task 3)

Cell: the campaign-1 S4 vc-relaxed BaRhH₆ (c-336a1451e874; 8 atoms; R-3m at 0.1 Å). Rh⁴⁺ is
formally low-spin d⁵, N = 9 + 6 + 2 = **17** (one t₂g hole per complex), and the cell has an odd
count of 17 electrons. Campaign 1 started its spin check only from a uniform +0.5 moment.

Script: `scripts/lane_p_barhh6_magnetism.py`. All runs use the S4 settings: SG15 PBE, 60/240 Ry,
MV 0.02 Ry, k 7×7×4 or 4×7×7. They run in 16-atom supercells doubled along a3 or a1, with Rh
split into two species Rh1 and Rh2, at 2 ranks.
- **Doubling along a3** makes the two shortest Rh–Rh contacts, a1+a3 and a2+a3 (4.63 Å),
  antiparallel.
- **Doubling along a1** makes a1+a3 and a1−a2 (5.04 Å) antiparallel.

The Rh net is near-triangular, so no collinear AFM is unfrustrated. These are the two simplest
bipartitions.

| run | start | final total / absolute moment (µB/cell) | E − E_NM (meV/f.u.) | wall |
|---|---|---|---|---|
| nm_a3 | nspin 1 | – | 0 | 3.6 min |
| fm_a3 | +0.8 / +0.8 | 0.00 / 0.00 | +0.03 | 9.5 min |
| afm_a3 | +0.8 / −0.8 | 0.00 / 0.00 | +0.03 | 11.3 min |
| nm_a1 | nspin 1 | – | 0 | 6.8 min |
| afm_a1 | +0.8 / −0.8 | 0.00 / 0.00 | +0.03 | 18.1 min |
| soc_nm_a3 | noncollinear + SOC, no moment | – | (SOC reference) | 17.9 min |
| soc_afm_a3 | SOC, AFM ±0.8 along z | −0.03 / 0.08 after 57 iterations, **not converged** (scf accuracy stalled at 4·10⁻⁷ Ry) | −0.25 vs soc_nm (unconverged) | 2.1 h, killed |

SOC uses FR pseudopotentials for Rh and H (SG15 FR-1.0). Ba has no SG15 FR file and is treated
scalar-relativistically (QE: "Not all pseudopotentials have spin-orbit data"). The +0.03 meV/f.u.
of the collinear spin-polarised runs is the nspin = 1 vs 2 smearing offset: both collapse to
the same zero-moment solution.

**Result.**
- At the PBE level (no U), with MV 0.02 Ry smearing, BaRhH₆ has **no stable FM or AFM moment**
  from either collinear bipartition. Every start collapses to 0 µB.
- With SOC, an AFM start decays to ≤ 0.04 µB per Rh. Its energy is within 0.25 meV/f.u. of
  SOC-NM, but that run did not reach conv_thr.
- The H2 outcome stays "non-magnetic". The tag changes from `H2_FM_START_ONLY` to
  `H2_FM_AFM_SOC_CHECKED(PBE, no U)`.

**Caveats.**
- 0.02 Ry smearing can suppress a weak itinerant moment.
- The J_eff = ½ Mott picture of R30 F6 needs a U. No DFT+U was run, so the oxide analogy is not
  tested.

Follow-ups: the a3 NM/FM/AFM triple at **MV 0.005 Ry** with ±1.0 starts (done, below) and the
SOC-AFM restart with `mixing_beta` 0.1 and `local-TF` mixing (running, §8).

**Small-smearing check (meas., finished 07:35).** In the a3 cell at MV 0.005 Ry, with ±1.0 µB
starts:
- NM: −551.42733938 Ry;
- FM: −551.42733528 Ry, 0.00 / 0.00 µB;
- AFM: −551.42733528 Ry, 0.00 / 0.00 µB.

All three converged, with walls of 3.5, 7.4 and 8.9 min. Both magnetic starts again collapse to
the zero-moment state, so the smearing caveat does not change the result.

**H2 for BaRhH₆ is "non-magnetic (PBE, collinear FM + 2×AFM at 0.02 and 0.005 Ry; SOC-AFM
not converged)".** The SOC-AFM restart is still running (§8).

## 7. Open problems and what was not verified

1. **The S3 score cannot rank this family.** See §2.2(1): saturation at E_hull ≤ 0 plus a
   light-element bonus made the top 10 all Li-AA′. The ranking needs a non-saturating hull term
   and/or a band-filling feature. The R30 18-electron distance is constant (N = 19) across the
   lane, so it does not help *within* lane P; a rigid-band N(E_F) proxy would, but it needs DFT.
   Owner: `score/s3.py`.
2. **Ru/Re cells have no ML hull** (§2.2(2)). 44 cells are ranked without the hull term.
3. **Composition novelty against ICSD, OQMD and GNoME was not checked.** Only Alex-MP was. The
   AA′ F-43m cells are new *structures*. 4 of the 10 S4 compositions are known in Alexandria as
   hexagonal polymorphs at E_hull,DFT 0.04–0.10 eV/atom, and whether the F-43m cell lies below
   them was not computed (S7a).
4. **The rule's DOS statements rest on one k mesh** (S4 11³ tetrahedron, 60 Ry). Mg₂IrH₆ and
   LiScIrH₆ both carry `NEF_MESH_SENSITIVE`, so the peak heights (0.39–0.43) are uncertain at the
   ~20 % level (UNVERIFIED).
5. **Starting radii** (Shannon CN8) were recalled, not checked. They only seed the relax.
6. **Not enumerated:**
   - three-A cells (e.g. LiMgSc₂(IrH₆)₂, §4 item 3);
   - 1 : 3 cells with integer N;
   - (111)-ordered (L1₁) M/M′;
   - A₂MH₇ / A₂MH₈ (R30 §2.2(3)).
7. **No DFT+U and no hybrid functional on BaRhH₆** (§6). The SOC-AFM run did not converge in its
   first attempt; a restart is queued (§8).
8. **The S4 set is homogeneous** (all Li-AA′, N = 19) because of item 1. No half-electron alloy,
   no N = 17 or 20 ternary and no Mg-based isoelectronic cell (e.g. NaMgPtH₆, MgScRuH₆) reached
   S4. The next S4 slots should be taken by hand from these classes:
   - **MgYOsH₆** c-279f7a83a52f (E_hull,ML −0.023);
   - **CaMgCoH₆** c-0d00ea0c4fe8 (+0.002);
   - **LiMg₃Ir₂H₁₂** c-66a04423bb33 (18.5 e, −0.018);
   - **K₂IrH₆** c-1b53417679e6 (17 e, −0.006);
   - **Sr₂NiH₆** c-036f8de605cc (20 e, −0.001).
   All are S3 PASS or rank-queue entries with E_hull,ML ≤ 0.

## 8. Running background jobs (at hand-off)

| job | what | log |
|---|---|---|
| `~/rtsc_targeted_work/P2-lane/queue3.sh` | BaRhH₆ SOC-AFM restart (`soc_afm_a3b.in`, ≤ 150 iterations, ~1 min each); 2 ranks, nice 19; then exits. The MV 0.005 Ry runs before it are done | `~/rtsc_targeted_work/P2-lane/queue3.log`; outputs `~/rtsc_targeted_work/P2-lane/barhh6_mag/{nm,fm,afm}_a3_s005.out`, `soc_afm_a3b.out` |

No GPU job is running. The lane-P ledger is idle (S1–S4 complete). Re-generate the tables with
`scripts/campaign2_lane_p_report.py --ledger ~/rtsc_targeted_work/P2-lane/ledger.sqlite --out runs/campaign2_laneP`.
