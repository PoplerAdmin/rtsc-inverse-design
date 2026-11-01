# R02 — Adversarial review of `docs/theory/02_materials_profile.md`

Reviewer: R2-review-materials. Method: primary sources fetched on-line (arXiv abstracts/HTML,
Europe PMC full text of Nat. Commun. papers), independent recomputation of every number
(POSIX awk, Appendix R-A), and a cross-check of the profile against the compounds the document
itself cites. No DFT was run. I did not rewrite the target document.

Severity: **critical** = the error would make the pipeline exclude valid candidates or would
invalidate a headline conclusion; **major** = wrong or unsupported claim that changes a
threshold or a citation used as a premise; **minor** = cosmetic, wording, or a small number.

## Summary verdict

The literature synthesis is mostly accurate. The Belli formula, the Gao et al. 2025 numbers, the
Sun 2016 statistics, the Hansen/Sinha Mg₂IrH₅/Mg₂IrH₇ results, the Ding et al. PIMD numbers and
all derived arithmetic (Hopfield, Belli thresholds, Allen–Dynes grid) check out. The document
fails in **Section 5**. Several rows meant as necessary conditions are not necessary. Four
rows (8, 15, 11, 1) would **reject Mg₂IrH₆ and Li₂AuH₆**, the document's own best existing
candidates and its priority class 1. Rows 3 and 4 and the closing "non-hydrides cannot reach
300 K" rest on a misreading of the Semenok bound, which the document's own numbers contradict.
One citation is misattributed (Li₂MgH₁₆) and one is misdescribed (arXiv:1307.1323).

---

## A. Citation checks

| # | claim in doc | verdict | evidence |
|---|---|---|---|
| A1 | Belli et al., Nat. Commun. 12, 5381 (2021): Tc = (750·Φ_DOS − 85) K, Φ_DOS = φ·H_f·H_DOS^{1/3}, "within 60 K", 178 collected / 43 discarded, no RMSE; φ(PdH)=0.19, φ(YH₄)=0.43, φ(H₃S)=0.68; YH₁₀ 326 K, H 300 K, H₃S 200 K; charge transfer 0–0.25 e/H; isolated H ≤ 40 K; H–H > 1.55 Å → ≤ 50 K; DOS alone not a descriptor | **Verified** (Europe PMC full text PMC8446067) | Note: Belli state the fit is "valid to estimate the superconducting critical temperature in hydrogen-based superconductors". Their Discussion also uses 1.55 Å for the *H–host* distance, so the paper is internally inconsistent; the doc uses only the H–H reading. Also "electride ≤ 50 K" has the exception Si₂H₆ 153 K, and "ionic highest 45 K (BeH₂)" conflicts with BeH₂ 97 K elsewhere in the same paper. |
| A2 | Pickard, Errea & Eremets, Annu. Rev. CMP 11, 57 (2020): S = Tc/√(Tc,MgB₂² + P²) | **Plausible / partially verified.** Vol. 11, pp. 57–76 is confirmed. Errea 2022 (arXiv:2201.09500) gives the same formula and attributes it to its Ref. 16; I could not render either PDF to confirm Eq. 8. The quote "more stable structures … lower densities of electronic states" was **not verified**. | |
| A3 | Gao et al., Nat. Commun. 16, 8253 (2025), arXiv:2502.18281 | **Verified**: >20 000 metals; ω_log "rarely exceeds 1800 K"; RT "extremely unlikely"; λ=2, ω_log=1800 K → 260 K (McMillan); "100–120 K" for designed hypotheticals; Li₂AgH₆ 108.8/83.0 K, Li₂AuH₆ 116.1/91.0 K; hull 0.319/0.172 eV/atom; ClB₂C₈ λ≈3.3, ω_log 425 K, 55 K; Al₂OsH₇ ω_log ≈ 300 K, 38 K; LiMoN₂ defects; isotropic MgB₂ ~2× low; no E_hull cutoff stated | "Practical estimate ~100–120 K" is fair, but in the paper it is the Tc of hypothetical compounds one "can design", not a derived limit (see D1). |
| A4 | Semenok, Altshuler & Yuzbashyan, arXiv:2407.12922 | **Verified with caveats**: λ* = 3.69 (Einstein, the *minimum* over shapes), 4.72 (Debye); λ_c = λ*/2; Tc < 0.20 ω_max (stable), < 0.32 ω_max (metastable); Tc < 600 K/√A (stable), < 950 K/√A (metastable) with ω₀ = 3000 K; "stable room-temperature … only … hydrogen (A=1) and deuterium (A=2)"; the MgB₂ ω_max = 1624 K in Table S1 is theirs | The abstract states only "λ ≲ 4". arXiv v2 (Aug 2025) lists no journal ref, so "Adv. Mater. 2025" is **unverified**. The √A form is a heuristic (ionic plasma frequency ∝ M^-1/2 with ω₀ fixed by metallic H), not a theorem. See C1. |
| A5 | Dolui et al., PRL 132, 166001 (2024) | **Mostly verified**: 160–175 K at ambient, λ ≈ 2.52, vHs from flat bands along L–W, ~4 eV crystal-field splitting, significant λ from H bands near 0.1 and 0.2 eV, quantum-anharmonic hardening of 30–80 meV and ~15 meV softening above 150 meV, non-zero-spin structures removed | **Not found:** "keep structures within 50 meV/atom" (Sec. 4.1). The arXiv text says only "close to the convex hull", with ternary hulls at 1 and 10 GPa. A different paper (arXiv:2401.17024) uses < 80 meV/atom. The paper does not say the 0.1/0.2 eV modes contribute "most of λ". |
| A6 | Mg₂IrH₆ ~60 meV/atom at 0 GPa marked UNVERIFIED | **Verified** from the Carnegie-group synthesis paper (Hansen et al., arXiv:2406.09538) and Sinha et al. 2026 ("∼60 meV/atom above the convex hull at ambient pressure"). | Drop the UNVERIFIED tag and cite Hansen et al. |
| A7 | Ding, Chen & Shi, arXiv:2604.12367 | **Verified**: Li₂AgH₆ lattice collapse; Li₂AuH₆ fluorite Li–Au sublattice intact, H partially dimerises and diffuses, Tc = 22 K | Now carries a journal DOI 10.1103/w7rx-sknt, so it is no longer only a preprint. It reports kinetic instability (a structural change), not an "anharmonic correction" of the same structure, so the doc's Sec. 4.2 item 3 miscategorises it. |
| A8 | Li₂MgH₁₆ "Sun, Lv, Xie, Liu & Ma, PRL 123, 097001 (2019); arXiv:2305.04875", 473 → 277 K "with anharmonicity" | **MISATTRIBUTED.** arXiv:2305.04875 is **Chen & Shi, "Coexistence of superconductivity and superionicity in Li₂MgH₁₆", PRB 109, L140505 (2024)**: 277 K at 260 GPa from a stochastic path-integral treatment of the *superionic* phase, versus 473 K harmonic. Sun et al. 2019 is the original 473 K prediction. The doc merges the two and calls the effect "anharmonicity". | |
| A9 | Gao, Lu & Xiang: Li₂B₃C, Li₃B₄C₂ "Tc > 50 K; arXiv:1307.1323" | **MISDESCRIBED.** arXiv:1307.1323 is Gao, Lu & Xiang, PRB 91, 045132 (2015), on Li₂B₃C only. Its abstract claims only a Tc "might be even higher than" MgB₂ and gives no 50 K figure. Li₃B₄C₂ comes from another paper (not located). | |
| A10 | Sanna et al., npj CM 10, 44 (2024) | **Verified**: 45–80 K, >100 K for e-doped Pt, "thermodynamically stable" | Exec. summary item 3 lists Sanna among "weakly metastable ≤ 70–80 K". Sanna et al. claim *stable*, so the correct wording is "on-hull, ≤ 80 K". |
| A11 | Cerqueira et al., Adv. Funct. Mater. (2024), 150k hydrides, "max 70" | **Slight misquote**: the abstract reads "≈50 systems > 20 K, some even reaching above 70 K", and most systems are slightly unstable. | |
| A12 | Cerqueira, Sanna & Marques, Adv. Mater. 36 (2024) | **Verified** (e2307085): ~7000 e-ph calculations, ~200 000 metals on/near hull, 541 with Tc > 10 K, LiMoN₂ > 38 K | |
| A13 | Hansen et al. arXiv:2406.09538; Sinha et al. arXiv:2602.23675 | **Verified** (Mg₂IrH₅ up to ~28 GPa; Mg₂IrH₇ > ~40 GPa, persists to ~20 GPa, reverts to Mg₂IrH₅) | |
| A14 | Ouyang et al., PRB 111, L140501 (2025) | **Verified** (~140 K at ambient) | |
| A15 | Sun et al., Sci. Adv. 2, e1600225 (2016) | **Verified**: 29 902 phases, median 15 ± 0.5, 90th percentile 67 ± 2 meV/atom; nitrides 63/195 | |
| A16 | Trachenko et al., arXiv:2406.08129 | **Verified** ("10²–10³ K"). No journal reference appears on arXiv; the cited JPCM DOI was not checked. | |
| A17 | Kong et al. cond-mat/0102499 | **Verified** (λ = 0.87, σ 42 % of N(0), ~590 cm⁻¹, π contributes 0.25) | |
| A18 | Choudhary & Garrity npj CM 8, 244 (2022); Gibson et al. (BETE-NET) | **Verified** (1736 → 1058 → 105; 818 materials, Tc MAE 2.5/2.1 K, ~5× precision) | |
| A19 | SrB₃C₃ (Zhu et al.) 27–43 K, onset ~20 K at 40 GPa | **Verified** from arXiv:1708.03483 | Rb₀.₅Sr₀.₅B₃C₃ "78 K" remains **unverified** (not in that abstract). |
| A20 | Moussa & Cohen PRB 74, 094520 | Content verified from arXiv (MgB₂ at its bound; B-diamond far from its bounds); journal ref not confirmed. Related: the same authors' 2008 APS abstract puts the maximum B-diamond Tc at ~55 K (20–30 % B, μ* = 0.12). | |

No fabricated papers were found. The two misattributions (A8, A9) and the unsupported
"Dolui 50 meV" (A5) are the citation defects.

---

## B. Recomputation of derivations and numbers (all reproduced; see R-A)

* **Hopfield (Sec. 1.3).** k = M_H ω² with ħω = k_B·1500 K = 0.1293 eV gives 64.5 N/m =
  **4.03 eV/Å²**. η_H = 2.5–3 × 4.03 = **10.1–12.1 eV/Å²** ✓. At 1000 K: 4.5–5.4 ✓.
  N_H = 10/25 = 0.4 ✓. I also checked whether the η ≳ 10 floor depends on the chosen point.
  Along the whole Tc = 300 K Allen–Dynes curve, η = λ·M·ω₂² is ≈ 10.8 eV/Å² at
  (λ=4, ω_log≈1020 K) and ≈ 16 eV/Å² at (λ=2, ω_log≈1760 K). The floor of ~10 eV/Å² per
  coupled H is therefore *robust* in the single-band isotropic model, which is good.
* **A = (λ/2)⟨ω⟩** ✓ (Allen–Dynes ω₁ definition). λ⟨ω⟩ ≥ 2600 K, 2600/3 = 867 K ✓.
  MgB₂ A ≈ 0.435 × 700–750 K ≈ 300–330 K ✓ (Semenok's table gives ω_log = 680 K).
* **Belli thresholds.** (300+85)/750 = 0.513 ✓; 325/750 = 0.433 ✓; with φ = 0.7: 0.733 ✓;
  H_f = 0.85 → H_DOS ≥ 0.642 ✓; H_f = 0.75 → H_DOS ≥ 0.935 (doc 0.92, minor);
  Mg₂XH₆ at φ = 0.7, H_DOS = 1 → 265 K ✓.
* **Allen–Dynes grid** reproduced to the kelvin. ω_log needed for 300 K at μ* = 0.10,
  ω₂/ω_log = 1.2: 1757 / 1451 / 1258 / 1123 / 1022 K for λ = 2/2.5/3/3.5/4 ✓.
  Sensitivity the doc omits: at ω₂/ω_log = 1.5 → 1699/1382/1179/1036/929 K; at μ* = 0.15 →
  2038/1661/1431/1272/1155 K. The ±15 % band should be stated.
* **S values.** Target 300/√(39²+1²) = **7.69**, not 7.6 (minor). H₃S (203 K, 155 GPa) = 1.27 ✓.
  LaH₁₀ (250 K, 170 GPa) = 1.43, not ≈1.3 (minor; the value depends on which Tc/P is used).

The arithmetic is sound. The problems are in how the numbers are *used* (Sections C and D).

---

## C. Findings on the logic of the ranked profile

### C1 (critical) — "H required" / "non-hydrides cannot reach 300 K" misreads Semenok and contradicts the doc's own numbers
*Doc:* Sec. 3.3 says the bound "excludes non-hydrides (MgB₂ ω_max ≈ 1624 K → Tc < 325 K
stable; Li/B/C/N-only frameworks have no route)". Row 3 says "for 300 K: H required (Semenok
bound)". Sec. 5 ends with "Non-hydride frameworks cannot reach 300 K by the Semenok bound and
the A requirement".
*Evidence:*
(i) 0.20 × 1624 K = 325 K **> 300 K**, so the doc's own example does *not* exclude an
MgB₂-frequency material, even in a stable phase.
(ii) The ω_max form allows 300 K whenever ω_max ≥ 1500 K (stable) or ≥ 940 K
(metastable). Diamond's optical phonon is ≈1918 K (1332 cm⁻¹), which gives 384 K stable and
614 K metastable. c-BN and sp²/sp³ C–C/B–C stretches are similar.
(iii) The √A form is a heuristic with ω₀ = 3000 K fixed by metallic hydrogen. It gives
ω_max(C) = 866 K, which strongly underestimates real stiff covalent lattices. Even taking it at
face value, the **metastable** bound is Li 359 K, Be 317 K, B 286 K, C 274 K. Li- and Be-based
metastable compounds are therefore *not* excluded, and the doc's own profile admits metastable
phases up to 200 meV/atom.
(iv) Semenok et al. themselves restrict their room-temperature claim to *stable* phases.
(v) The A requirement λ⟨ω⟩ ≥ 2600 K is satisfiable by C/B frameworks if λ ≈ 2 with
⟨ω⟩ ≈ 1300 K. That is unprecedented (the best is p-doped graphane at λ ≈ 1.3 and ~90 K;
Savini, Ferrari & Giustino, PRL 105, 037002 (2010)), but it is a coupling-strength question, not
a bound.
*Required change:* Remove "H required" as a necessary condition. Re-state: "No proven bound
excludes Li/Be/B/C/N-carried 300 K superconductivity in metastable phases. H is a strong
*prior* because only H reaches A ≥ 1300 K with λ ≤ 3." Keep the non-hydride channel as a
minority exploration slot (e.g. 10–20 % of the generation budget: metallised sp³ B–C–N
frameworks, Li/Be-rich intermetallics and borocarbides, hole-doped graphane-like sheets), not
merely as "calibration/controls". Correct "Li/B/C/N-only frameworks have no route".

### C2 (critical) — Row 8 (H bonding motif) rejects the document's own priority class 1
*Doc:* Row 8 requires "extended network: H–H 1.0–1.55 Å … or M–H–M bridging; ELF networking
φ ≥ 0.6; … reject isolated H⁻".
*Evidence:* In K₂PtCl₆-type Mg₂IrH₆ / Li₂AuH₆, H sits in *isolated* [MH₆] octahedra. Each H
binds one metal, and the shortest H–H is ≈ √2·d(M–H) ≈ 2.3–2.4 Å, so there is no H–H
< 1.55 Å and no M–H–M bridge. The doc itself says their φ "is surely lower" than 0.7. Belli's
own "isolated" class tops out at 40 K, yet Dolui et al. compute 160–175 K for Mg₂IrH₆. **Belli's
class rule is falsified at ambient pressure by the very compound the doc ranks first.** Row 8
would discard Mg₂IrH₆, Li₂AuH₆, Mg₂PtH₆ (Sanna) and every A₂MH₆/AMH₃ candidate in "priority
class 1".
*Required change:* Downgrade row 8 to a ranking feature, not a filter. Explicitly allow
"M–H σ-bonded [MH_n] complexes with H-1s/M-e_g states at E_F". Note that Belli's class
statistics were fitted on compressed hydrides and are contradicted by Mg₂IrH₆.

### C3 (critical) — Row 15 (≲ 5–6 Å³/H) is a high-pressure number that excludes every ambient hydride
*Evidence (reviewer computation from textbook lattice constants):* H₃S Im-3̄m at 155 GPa
(a ≈ 3.09 Å, 6 H/cell) gives ≈ 4.9 Å³/H, so the threshold is the *compressed* value.
LiH (a = 4.08 Å) gives 17 Å³/H. Mg₂FeH₆ (a = 6.44 Å, 24 H per conventional cell) gives
≈ 11 Å³/H, and Mg₂IrH₆ is similar or larger. At 0–1 GPa essentially no hydride reaches 6 Å³/H.
*Required change:* Delete row 15, or replace it with a calibrated, ambient-referenced quantity.
For example, compute the H-volume distribution for the known ambient candidates and set the
threshold at their 10th percentile. Never use it as a filter.

### C4 (critical) — E_hull ≤ 50 meV/atom as a MatterGen condition guarantees the pipeline reproduces the known ceiling
*Doc:* Row 5 sets E_hull ≤ 50 meV/atom as the "synthesizable band", conditioned directly in
MatterGen, with ≤ 150–200 as a flagged exploratory window.
*Evidence:* The doc's own Sec. 3.3/4.1 (and Gao et al. Fig. 4) say Tc rises along a Pareto
front with rapidly increasing E_hull. Every 100–170 K candidate sits at 60 (Mg₂IrH₆), 172
(Li₂AuH₆) or 319 (Li₂AgH₆) meV/atom, so **all of them fail row 5**. Li₂AgH₆ fails even the
exploratory window. Other problems:
(a) Sun 2016 percentiles describe phases made at ambient or near-ambient conditions. Nitrides
already have a 90th percentile of 195 meV/atom, and high-pressure-quenched phases (SrB₃C₃,
diamond-anvil hydrides) are not represented.
(b) DFT hulls for hydrides omit ZPE, which differs between competing H phases by tens of
meV/atom.
(c) The "Dolui 50 meV" support is not in the paper (A5).
(d) Conditioning MatterGen on low E_hull steers sampling into the near-hull region that
Choudhary, Cerqueira ×2, Sanna and Gao have already exhausted (40–80 K). A 300 K material, if
one exists, is expected off that region.
*Required change:* Make E_hull a ranking/Pareto variable, not a necessary condition. For
generation, sample the conditioning value across 0–0.3 eV/atom (or do not condition on it).
Post-filter at ≤ 0.35 eV/atom (ZPE-corrected for H), with a stated synthesis route required
above ~0.1. Remove "Dolui 50 meV" or mark it unverified.

### C5 (major) — Row 11 hard-conditions on five space groups and excludes known high-Tc structure types
*Evidence:* The list {221, 223, 225, 229, 191} omits P6₃/mmc (194: CeH₉, YH₉, (La,Ce)H₉),
Fd-3̄m (227: Li₂MgH₁₆, which the doc cites), I4/mmm (139: YH₄, many ternaries), R-3̄m (166:
LiMoN₂, the doc's own on-hull candidate), Fm-3̄m-derived tetragonal distortions, and all
lower-symmetry ternaries. The doc concedes "symmetric phases are also the ones prone to
distortion". The "every known high-Tc SC is high-symmetry" argument is survivorship: only
high-symmetry structures were prioritised for expensive e-ph calculations. In addition, the
MatterGen `space_group` checkpoint conditions on *one* space group at a time.
*Required change:* Use space group as a soft prior. Spread the budget over the ≥ 10 groups
above, plus an unconditioned share. Do not post-filter on symmetry.

### C6 (major) — Rows 1 and 13 contradict each other
Row 1 rejects E_F "within ~0.2 eV of a gap edge". Row 13 asks for "metallic composition
adjacent to a closed-shell insulating parent" (doped Mg₂IrH₅, σ-hole MgB₂, doped LiBC). Those
materials have E_F near a gap edge by construction. Mg₂IrH₆, with a vHs next to the
Ir-t₂g/e_g gap, would fail row 1.
*Required change:* Drop the 0.2 eV clause from row 1, or turn it into a flag
("doping-sensitive"). Keep row 13.

### C7 (major) — Rows 4/6 and exec. summary item 2: Belli regression used as a necessary condition
*Evidence:*
(i) A fit with ±60 K scatter, calibrated on ~135 compressed hydrides, gives no *bound*. Only
two data points (YH₁₀ 326 K, metallic H 300 K) lie at ≥ 300 K, so 300 K is at the extrapolation
edge, where the prediction interval is wider than ±60 K.
(ii) Mg₂IrH₆ is a direct ambient counterexample to the class rule (C2), so transferability is
already doubtful.
(iii) "Forces H_f ≳ 0.75" holds only for φ ≤ 0.7. At φ = 0.8 the floor is 0.64, and at φ = 1
it is 0.51.
(iv) **Internal inconsistency:** the doc's priority class 1 (A₂MH₆, H_f = 0.667; AMH₃,
H_f = 0.60) fails its own "≥ 0.75 target". By the doc's own Belli arithmetic, class 1 cannot
reach 300 K unless φ ≳ 0.77 (A₂MH₆) or ≳ 0.86 (AMH₃).
*Required change:* Present H_f and H_DOS as ranking scores. Use H_f ≥ 0.5 at most as a soft
floor, and only for the H-channel. State that the Belli fit is untested at 1 atm and is
contradicted by Mg₂IrH₆. Reconcile the class-1 priority with rows 4/6. Either class 1 is a
calibration family for ≤ 200 K, or the H_f target goes.

### C8 (major) — The A ≥ 1300 K "necessary condition" holds only in the single-band isotropic model; multiband/anisotropic enhancement is absent from the doc
*Evidence:* The multiband notes of the certified core (`rtt.core`, wMB)
states that the DOS-averaged A can underestimate the Tc bound, with a counterexample exceeding
it by 9 %. The multiband bound uses κ, the Perron root of the symmetrised band-coupling matrix.
Gao et al. note that isotropic calculations underestimate MgB₂ by ~2×. The Gao/Cerqueira
"ceiling" datasets are isotropic, so they are biased low for exactly the σ-band/multiband
materials the doc dismisses as ≤ 80 K.
*Required change:* Add a multiband caveat to Secs. 1.4/3.3. Use κ (band-resolved) in place
of A wherever a candidate has distinct coupled sheets. Add a "multiband σ-hole" class
(MgB₂-like with stiffer bonds, Li/Be/B/C) to the exploration share.

### C9 (major) — Row 10 caps λ at 3.7 and contradicts its own data
λ* = 3.69 is Semenok's *minimum*, for Einstein spectra. For Debye-like shapes λ* = 4.72, and
the abstract says only "λ ≲ 4". The doc cites Li₂AgH₆ at λ ≈ 4 (and Li₂MgH₁₆ at λ = 3.7).
Its own Allen–Dynes grid shows that (λ = 4, ω_log ≈ 1020 K) reaches 300 K.
*Required change:* Use λ ∈ [2, 4.7]. Treat λ > 3.7 as a "near-instability" flag rather than a
rejection.

### C10 (major) — MatterGen "C" column overstates what can be conditioned jointly
*Evidence:* The installed checkpoints are `chemical_system`,
`chemical_system_energy_above_hull`, `dft_band_gap`, `dft_mag_density`,
`dft_mag_density_hhi_score`, `ml_bulk_modulus`, `space_group` (plus the base models). There is
no joint band_gap + mag_density + E_hull + space_group model, and `energy_above_hull` exists
only paired with `chemical_system`. Rows 1, 2, 5 and 11 marked "C" cannot be applied together
without fine-tuning. Training data for these models (Alex-MP) contains few H-rich metals, so
H_f ≥ 0.6 is also out of distribution.
*Required change:* Mark which rows are jointly conditionable with the available checkpoints.
Specify the generation plan (e.g. `chemical_system_energy_above_hull` over enumerated X–Y–H
and Li/Be–B–C systems, and everything else as post-filter). Flag fine-tuning as a prerequisite
for anything else.

### C11 (minor) — Row 9 "H bands ≥ 1000–1500 K carry ≥ 50 % of λ" is stricter than A requires
A weights α²F without 1/ω, so modes below 1000 K can carry most of λ while A ≥ 1300 K still
holds (e.g. λ = 4 at ω_log ≈ 1000 K). Required change: state the phonon criterion directly as
A ≥ 1300 K (or κ, see C8), and ω_max ≥ 940 K for metastable (≥ 1500 K for stable) as a flag.

### C12 (minor) — Sec. 3.4 "the required corner (λ ≥ 2.5 and ω_log ≥ 1250 K)"
This is not the 300 K region: (2.5, 1250 K) gives ~260 K, and (2.0, 1760 K) gives 300 K. Use the
iso-Tc curve and include the μ*/ω₂ sensitivity from Section B. "λ ≳ 3 occurs with
ω_log ≲ 500 K" is supported by two examples only; label it as such.

### C13 (minor) — Mislabelled corrections in Sec. 4.2 item 3
Li₂MgH₁₆ 473 → 277 K comes from superionic path-integral treatment (Chen & Shi 2024), and
Li₂AuH₆ 140 → 22 K comes from kinetic H₂ formation (a different structure). Neither is a
"20–40 % anharmonic correction". Both downward results come from one group and one method
(the Shi group's stochastic path integral). Required change: separate "anharmonic
renormalisation (SSCHA, typically −10 to −40 %)" from "kinetic instability (structure
changes)", and note the single-method provenance.

### C14 (minor) — Small items
S target 7.69 (not 7.6). LaH₁₀ S ≈ 1.4. H_DOS ≥ 0.935 (not 0.92). Ding et al. is now
published (DOI 10.1103/w7rx-sknt). Mg₂IrH₆ 60 meV is verified (A6). Cerqueira AFM "above 70 K",
not "max 70". Sanna is on-hull (A10). Semenok's journal ref is unverified. Row 14's
"avoid 3d magnetic" would exclude Mg₂FeH₆-type and NaNiH₃ parents (Gao et al.'s high-ω_log
example). Use row 2 (non-magnetic) to handle magnetism and do not exclude elements.

---

## D. Is "realistic expectation ≤ ~150 K" fair?

### D1 (major) — Fair as a prior, overstated as evidence of a ceiling
*For the doc:* Every computed ambient dataset tops out at 40 K (on hull), ~80 K (near hull)
and 100–175 K (metastable, unsynthesised, partly destroyed by PIMD). Gao et al. call RT
"extremely unlikely". As an expectation for *this pipeline's* output, ≤ 150 K is reasonable.
*Against generalising it to a ceiling:*
(a) **Non-independence.** Cerqueira ×2, Sanna and Gao share the Marques/Sanna/Errea group,
workflow and Alexandria-derived structure space. The "systematic searches" are largely one
sampling of one prior.
(b) **Sampling bias.** The searches enumerate known prototypes, small cells and near-hull
compositions. "ω_log rarely exceeds 1800 K" describes that sample, not physics. Gao's own
statement of 100–120 K is about compounds one can design from it.
(c) **Isotropic, harmonic, adiabatic e-ph** throughout. Isotropic calculations underestimate
MgB₂ ~2×, and anharmonicity cuts both ways (it stabilises H₃S and LaH₁₀).
(d) The **ambient record has moved 39 K → predicted 160 K in a few years**, so the distribution's
upper edge is not converged.
(e) The pipeline as specified (C4, C5, C2) would *enforce* the ceiling by construction.
*Required change:* Re-word as "Prior expectation for outputs: ≤ ~150 K. No fundamental bound
excludes 300 K for metastable H-based (and, by the metastable Semenok bound, Li/Be-based)
compounds. The evidence for a practical ceiling comes from correlated, isotropic searches of
near-hull prototype space." Keep the "≤ 150 K" expectation, but separate it from the claim of
impossibility.

### D2 (minor) — Non-hydride exclusion leans on folklore-level evidence
Hole-doped LiBC failing is marked UNVERIFIED. Li₂B₃C has never been tested. B-doped diamond
reached 4–11 K against a ~55 K Moussa–Cohen estimate. Graphane (~90 K predicted) has never been
doped. These are *synthesis* failures, not physics ceilings, and they argue for keeping a
modest non-hydride channel (C1).

---

## E. Priority list of required changes

1. **C2, C3, C4** — remove or soften rows 8, 15 and 5 as filters, because they reject Mg₂IrH₆,
   Li₂AuH₆ and every ambient hydride.
2. **C1** — remove "H required" and "non-hydrides cannot reach 300 K", and fix the 325 K
   self-contradiction.
3. **C7, C9, C5, C6** — convert the Belli H_f/H_DOS, λ ≤ 3.7, five-space-group and gap-edge
   rules into ranking features or flags. Reconcile priority class 1 with rows 4/6.
4. **C8** — add the multiband (κ) caveat.
5. **A8, A9, A5** — fix the Li₂MgH₁₆ attribution, the Li₂B₃C citation and the
   "Dolui 50 meV" claim.
6. **C10** — make the MatterGen conditioning column match the available checkpoints.
7. **D1** — re-word the bottom line so the expectation is separate from impossibility.

## R-A. Recomputation script (POSIX awk)

```sh
awk 'function ad(l,wl,w2,mu){L1=2.46*(1+3.8*mu);L2=1.82*(1+6.3*mu)*(w2/wl);
 f1=(1+(l/L1)^1.5)^(1/3);f2=1+(w2/wl-1)*l*l/(l*l+L2*L2);
 return f1*f2*wl/1.2*exp(-1.04*(1+l)/(l-mu*(1+0.62*l)))}
function solve(l,r,mu, lo,hi,m,i){lo=100;hi=5000;for(i=0;i<60;i++){m=(lo+hi)/2;
 if(ad(l,m,r*m,mu)<300) lo=m; else hi=m}; return m}
BEGIN{for(l=2;l<=4.01;l+=0.5) printf "lam=%.1f wlog300(mu.1,r1.2)=%.0f (r1.5)=%.0f (mu.15)=%.0f\n",
  l,solve(l,1.2,0.1),solve(l,1.5,0.1),solve(l,1.2,0.15);
 kB=8.617333e-5;hbar=6.582120e-16;M=1.6735e-27;w=1500*kB/hbar;k=M*w*w;
 printf "k=%.2f N/m=%.3f eV/A2\n",k,k/16.0218;
 for(A=1;A<=12;A++) if(A==1||A==7||A==9||A==11||A==12)
  printf "A=%d stable<%.0f metastable<%.0f\n",A,600/sqrt(A),950/sqrt(A)}'
```
Output: ω_log(300 K) = 1757/1451/1258/1123/1022 K (μ* 0.1, r 1.2); 1699/1382/1179/1036/929
(r 1.5); 2038/1661/1431/1272/1155 (μ* 0.15). k = 64.54 N/m = 4.03 eV/Å². Semenok
stable/metastable: H 600/950, Li 227/359, Be 200/317, B 181/286, C 173/274 K.

## Not verified by this reviewer
Pickard–Errea–Eremets Eq. 8 and the "more stable → lower DOS" quote (PDF not renderable here).
Semenok's journal ref. Trachenko's JPCM DOI. Li₃B₄C₂ source. Rb₀.₅Sr₀.₅B₃C₃ 78 K.
Hole-doped LiBC negative experiments. Gao et al.'s statement that Li₂AgH₆ has λ ≈ 4. Absolute
η_H values of real hydrides (H₃S, Mg₂IrH₆), which should be computed to calibrate row 7.
Mg₂IrH₆ volume per H: estimated from the isostructural Mg₂FeH₆, not computed.
