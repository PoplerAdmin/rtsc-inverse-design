# 02 — From α²F targets to a materials target profile (ambient pressure)

Agent: T2-materials. **Revision 2**, written after adversarial review R02
(`docs/reviews/R02_materials_profile_review.md`). The response to each finding is in the last
section.

Status: literature synthesis plus small derivations. No DFT was run. Every number is one of:
derived here (shown), computed (script in Appendix A), or cited.
- **UNVERIFIED** marks citations whose bibliographic details I recalled and did not check
  on-line.
- "(per R02)" marks items the reviewer verified on-line that I adopted without re-fetching.

Sibling document (α²F level, T1) sets the target. It requires A = ∫α²F(ω)dω ≥ 4.33148·k_B·300 K
≈ **1300 K (≈112 meV)**. This comes from the certified single-band isotropic bound
k_B·Tc ≤ A/4.33148 in `rtt.core.certify`. T1 also gives λ ≈ 2–4 and ω_log ≳ 1000 K.

This document asks two questions:
1. What must a crystal look like for that α²F to exist at 0–1 GPa?
2. What does the literature say about whether such crystals exist?

---

## 0. Executive summary

1. **λ = η/(M⟨ω²⟩), with η = N(E_F)⟨I²⟩** (McMillan–Hopfield). Light atoms raise the frequency
   prefactor, but the coupling itself must come from a large Hopfield parameter η on those
   atoms. Derived (Sec. 1.3): if H carries the coupling, 300 K needs η_H ≳ 10 eV/Å² per H atom.
   This floor holds along the whole 300 K iso-Tc curve in the single-band model.
2. **Necessary conditions, within stated model assumptions (Sec. 5a).** Hydrogen is not on
   this list. The conditions force light atoms or very stiff bonds carrying strong coupling.
   - Metallic and non-magnetic.
   - Dynamically and kinetically stable on the quantum, zero-point-inclusive energy surface.
   - A ≥ 1300 K. For multiband materials, a band-resolved analogue κ replaces A.
   - λ ≤ λ* ≈ 3.7–4.7. Together with A ≥ 1300 K this forces ⟨ω⟩ ≥ ~550 K.
   - ω_max ≥ 940 K for metastable phases, ≥ 1500 K for stable ones (Semenok et al.).
3. **Soft priors (Sec. 5b).** Everything else is used for ranking or conditioning, not
   rejection: H fraction, Belli networking value, H share of the DOS at E_F, symmetry, volume
   per H, bonding motif, and E_hull.
   - E_hull is a Pareto axis. A zero-point-aware post-filter window extends up to ~0.35 eV/atom.
   - Several of these priors fail at ambient pressure for the best existing candidate,
     Mg₂IrH₆. It has isolated [IrH₆] units, H_f = 0.67, E_hull ≈ 60 meV/atom, and a predicted Tc
     of 160–175 K.
4. **Minority non-hydride channel (Sec. 5c), 10–20 % of the generation budget.** This covers
   Li, Be, B, C and N compounds and multiband σ-hole materials. No proven bound excludes them at
   300 K in metastable phases. The obstacle is coupling strength.
5. **Ambient-pressure state of the art (Sec. 3).**
   - Measured record: MgB₂, 39 K.
   - Predicted on or near the hull: about 40 K (LiMoN₂), and ≤ 80 K for Mg₂XH₆ (Sanna et al.,
     claimed stable).
   - Predicted, slightly unstable hydrides: some "above 70 K" (Cerqueira et al., Adv. Funct.
     Mater. 2024).
   - Predicted, strongly metastable hydrides: 100–175 K (Mg₂IrH₆, Li₂AuH₆, Li₂AgH₆). None has
     been synthesised.
     - Mg₂IrH₆ attempts produced Mg₂IrH₅ and Mg₂IrH₇ instead.
     - Path-integral molecular dynamics (PIMD) finds Li₂AgH₆ collapses and Li₂AuH₆ forms H₂,
       which drops its Tc to 22 K.
6. **Expectation vs. impossibility (Sec. 3.5).**
   - *Prior expectation for this pipeline's outputs: ≤ ~150 K.*
   - *Separately, no fundamental bound excludes 300 K for metastable H-, Li- or Be-based
     compounds.* The metastable Semenok bound gives H 950 K, Li 359 K and Be 317 K. The ω_max
     form of the same bound allows any phase with ω_max ≥ 940 K.
   - The evidence for a practical ~100–150 K ceiling comes from searches that are correlated
     (largely one group), isotropic, harmonic, and confined to near-hull prototype space.

---

## 1. Microscopic decomposition of λ

### 1.1 McMillan–Hopfield form

For a monatomic metal (McMillan 1968, Phys. Rev. 167, 331; Hopfield 1969, Phys. Rev. 186, 443;
page numbers **UNVERIFIED**):

  λ = N(E_F)⟨I²⟩ / (M⟨ω²⟩) ≡ η / (M⟨ω²⟩).

- N(E_F) is the DOS per spin per atom at E_F.
- ⟨I²⟩ is the Fermi-surface average of the squared electron–ion matrix element.

For compounds whose species have well-separated modes, λ ≈ Σ_j η_j/(M_j⟨ω²⟩_j). This is the
Gaspari–Gyorffy-type estimate used as a hydride Tc pre-screen by Shipley, Hutcheon, Needs &
Pickard (PRB 104, 054501 (2021), arXiv:2105.02296).

| factor | wants | physical origin | price |
|---|---|---|---|
| N(E_F) on the coupled atoms | large | flat bands or van Hove singularities (vHs) of bonding or antibonding character at E_F | high N(E_F) favours Stoner, charge-density-wave or distortion instabilities; stable phases tend to open a gap |
| ⟨I²⟩ | large | E_F states whose energy changes quickly with bond length: σ bonds, H 1s–H 1s, H 1s–M d | partly emptied strong bonds cost energy |
| M | small | H ≫ Li, Be, B, C, N | zero-point motion, anharmonicity, diffusion (Sec. 4) |
| M⟨ω²⟩ = k | small for λ | — | λ = η/k does not depend on M at fixed k; M enters Tc only through ω ∝ √(k/M) |

So a light element is needed for high frequencies, but it does not by itself supply coupling.

### 1.2 Mechanisms that make η large

**Metallized covalent bonds (σ-band metallicity)**
- MgB₂ has partly emptied B sp² σ bands. The σ holes couple to the in-plane B–B stretch near
  ~590 cm⁻¹ (Kong, Dolgov, Jepsen & Andersen, PRB 64, 020501(R) (2001),
  arXiv:cond-mat/0102499).
  - Total λ = 0.87.
  - The σ bands carry 42 % of N(E_F) but dominate λ. The π electrons contribute about 0.25.
- An & Pickett (PRL 86, 4366 (2001), arXiv:cond-mat/0102391) call this "covalent bonds driven
  metallic".
- Gao, Lu & Xiang propose σ-metallization as a general route (Physics (Wuli) 44, 421 (2015),
  doi:10.7693/wl20150701; author list **UNVERIFIED**).
  - Their stoichiometric Li₂B₃C paper (PRB 91, 045132 (2015), arXiv:1307.1323) claims only
    that Tc "might be even higher than" MgB₂.
  - The Li₃B₄C₂ prediction comes from a paper I did not locate.
- Hole-doped LiBC was predicted far above MgB₂ (Rosner, Kitaigorodsky & Pickett, PRL 88, 127001
  (2002), arXiv:cond-mat/0111592). Experiments reportedly found no superconductivity
  (**UNVERIFIED**). These are failures of synthesis and doping, not physics limits.

**H-derived states at E_F**
- In H₃S, LaH₁₀ and Mg₂IrH₆, H 1s states contribute strongly at E_F.
- For Mg₂IrH₆, Dolui et al. 2024 report:
  - E_F sits near a vHs-like feature from flat bands along L–W.
  - H-1s/Ir-e_g hybridisation splits the crystal field by about 4 eV.
  - λ ≈ 2.3–2.5, with significant contributions from H bands near 0.1 and 0.2 eV.
- In Mg₂IrH₆ the H is **not** networked. It sits in isolated σ-bonded [IrH₆] octahedra, and
  the shortest H–H distance is about 2.3–2.4 Å (derived from the geometry).

**Networks of weakened covalent bonds (compressed hydrides)**
- Belli et al. 2021 find the highest Tc in weak-covalent H–H networks, with charge transfer of
  0–0.25 e per H.
- Tc falls below ~50 K when H–H > ~1.55 Å.
- Their "isolated" class stays ≤ 40 K.
- These statistics come from compressed hydrides. Mg₂IrH₆ beats the isolated-class cap by
  about 4× at ambient pressure. **So Belli's class rule does not transfer to ambient pressure,
  and I use it only as a ranking feature.**

### 1.3 Derived requirement on the Hopfield parameter (this work)

Take λ_H ≈ 2.5–3, with the coupled H modes at ⟨ω²⟩^½ ≈ 1500 K (ħω = 0.129 eV):

  M_H⟨ω²⟩ = 1.674×10⁻²⁷ kg × (1.96×10¹⁴ s⁻¹)² = 64.5 N/m ≈ **4.0 eV/Å² per H**.

That gives η_H ≈ **10–12 eV/Å² per H**.

R02 checked this along the whole Allen–Dynes 300 K curve (μ* = 0.1, ω₂ = 1.2·ω_log):
- η = λ·M·ω₂² ≈ 10.8 eV/Å² at (λ = 4, ω_log ≈ 1020 K).
- η ≈ 16 eV/Å² at (λ = 2, ω_log ≈ 1760 K).

**So the ~10 eV/Å² floor is robust whenever H carries the coupling, in the single-band model.**

Assume ⟨I_H²⟩^½ ≈ 5 eV/Å (order of magnitude, **UNVERIFIED**). Then N_H(E_F) ≈ 0.4
states/(eV·spin·H). Absolute η_H for H₃S and Mg₂IrH₆ should be computed to calibrate this
(open item).

### 1.4 What A = ∫α²F dω means for materials, and the multiband caveat

The first-moment identity gives A = (λ/2)⟨ω⟩.
- A ≥ 1300 K therefore requires λ⟨ω⟩ ≥ 2600 K.
- With λ ≤ λ* = 3.69–4.72 (Sec. 3.3), this forces **⟨ω⟩ ≥ 550–700 K**.
- The coupled modes must therefore lie mostly above ~550 K (47 meV). This is a necessary
  condition, derived here.

For comparison, MgB₂ has A ≈ 0.435 × 700–750 K ≈ 300–330 K (derived, approximate).

**Multiband caveat (R02-C8).** A ≥ 1300 K is necessary only in the single-band isotropic
model.
- The multiband notes of `rtt.core` (wMB) say the DOS-averaged A can underestimate the
  multiband Tc bound. Its counterexample exceeds the A-bound by 9 %.
- The multiband bound uses κ, the Perron root of the symmetrised band-coupling matrix.
- Isotropic calculations underestimate MgB₂ by about 2× (Gao et al. 2025).

For candidates with distinct coupled Fermi sheets, the hard condition is the κ-based analogue,
not A. The isotropic "ceiling" datasets (Sec. 3) are biased low for exactly these σ-band and
multiband materials.

---

## 2. Quantitative descriptor models from the literature

### 2.1 Belli, Novoa, Contreras-García & Errea, Nat. Commun. 12, 5381 (2021), arXiv:2103.07320

**Definitions**
- φ (networking value): the highest ELF value at which an isosurface percolates through the
  crystal in all three directions.
- H_f: the H atom fraction.
- H_DOS: the H share of N(E_F).
- Φ_DOS = φ·H_f·H_DOS^{1/3}.

**Fit:** **Tc = (750·Φ_DOS − 85) K**, stated accuracy "within 60 K".
- No RMSE or correlation coefficient is reported.
- Data: 178 compounds collected and 43 discarded, so about 135 used. All are compressed
  hydrides.
- Example φ values: PdH 0.19, YH₄ 0.43, H₃S 0.68.

**Derived thresholds (this work)**
- Tc = 300 K requires Φ_DOS ≥ **0.513**, or ≥ 0.433 at the −60 K error edge.
- The H_f floor (with H_DOS = 1) depends on φ: 0.51 at φ = 1, 0.64 at φ = 0.8, 0.73 at φ = 0.7.
- At φ = 0.7: H_f = 0.85 needs H_DOS ≥ 0.64, and H_f = 0.75 needs H_DOS ≥ **0.935**.
- Even with H_DOS = 1, A₂MH₆ (H_f = 0.667) needs φ ≳ 0.77, and AMH₃ (H_f = 0.60) needs
  φ ≳ 0.86.

**Status: a ranking feature, not a bound.**
1. Only two training points (YH₁₀ and metallic H) are at ≥ 300 K, so 300 K is at the
   extrapolation edge.
2. The fit is calibrated at high pressure.
3. Mg₂IrH₆ is an ambient-pressure counterexample to the class rule (Sec. 1.2).

### 2.2 Figure of merit S (Pickard, Errea & Eremets, Annu. Rev. Condens. Matter Phys. 11, 57–76 (2020), arXiv:1910.00385)

  S = Tc / √(Tc,MgB₂² + P²), with Tc in K, P in GPa, and Tc,MgB₂ = 39 K.

I read this as Eq. 8 of the arXiv PDF, extracted to text with pdftotext.
- The review's own values:
  - S(MgB₂) = 1.
  - A 390 K superconductor at ambient pressure would have S = 10.
  - S(H₃S) = S(LaH₁₀) = 1.3, and S(Hg-1223) = 3.5.
- Recomputed: H₃S at 203 K and 155 GPa gives 1.27, and LaH₁₀ at 250 K and 170 GPa gives 1.43.
- Our target (300 K at 1 GPa) has **S = 7.69**. Mg₂IrH₆ (160 K at 0 GPa) would have S = 4.1.

The same extracted text contains this sentence: "more stable structures or compositions will
typically have lower densities of electronic states at the Fermi level". Errea reuses S in
J. Phys.: Condens. Matter 34 (2022), doi:10.1088/1361-648X/ac5b46, arXiv:2201.09500.

### 2.3 Other descriptor and ML models, with their reported errors

**Screening workflows and ML models**
- **Shipley et al. 2021:** pre-screen on DOS plus Gaspari–Gyorffy, then 121 converged
  electron–phonon calculations.
- **Choudhary & Garrity, npj Comput. Mater. 8, 244 (2022):** pre-screen on Debye temperature
  and N(E_F). 1736 candidates → 1058 electron–phonon calculations → 105 with Tc ≥ 5 K.
- **Gibson et al., BETE-NET** (npj Comput. Mater. 2025, arXiv:2401.16611): 818 materials, Tc
  MAE 2.5 K (2.1 K with phonon-DOS features). Trained on low-Tc data, so use it for triage only.
- **Cerqueira, Sanna & Marques, Adv. Mater. 36, e2307085 (2024):** about 7000
  electron–phonon calculations and about 200 000 metals screened. 541 have Tc > 10 K; LiMoN₂
  reaches > 38 K.

**Qualitative reviews.** None of these gives a closed-form Tc model.
- Hilleke & Zurek, J. Appl. Phys. 131, 070901 (2022), arXiv:2111.03697.
- Boeri et al., "The 2021 room-temperature superconductivity roadmap", J. Phys.: Condens.
  Matter 34, 183002 (2022).
- Duan et al., Natl. Sci. Rev. 4, 121 (2017) (details **UNVERIFIED**).

**Reported and observed errors**
- Belli: ±60 K.
- BETE-NET: 2–2.5 K MAE, only in the low-Tc regime.
- Spread between methods for a single compound: Mg₂IrH₆ ranges from 65 to 175 K.
- Allen–Dynes: the ω_log needed for a given Tc shifts by ±15 % across reasonable μ* and
  spectral-shape choices (Sec. 3.4).

---

## 3. Ambient-pressure state of the art

### 3.1 Measured (1 atm, conventional)

- **MgB₂, 39 K** (Nagamatsu et al., Nature 410, 63 (2001); details **UNVERIFIED**).
- **B-doped diamond, about 4–11 K** (Ekimov et al., Nature 428, 542 (2004), **UNVERIFIED**).
  - Moussa & Cohen (arXiv:cond-mat/0607832; PRB 74, 094520 (2006), journal ref **UNVERIFIED**)
    put MgB₂ at their doping bound and B-diamond far below its bound.
  - A 2008 estimate by the same authors gives a B-diamond maximum of about 55 K (per R02).
- **SrB₃C₃** (Zhu et al., Phys. Rev. Research 5, 013012 (2023); arXiv:1708.03483):
  - Predicted Tc 27–43 K.
  - Synthesised at about 50 GPa and quenchable.
  - Measured onset about 20 K at 40 GPa.
  - The Rb₀.₅Sr₀.₅B₃C₃ value of 78 K is **UNVERIFIED**.
- **p-doped graphane:** predicted λ ≈ 1.3 and Tc ~90 K (Savini, Ferrari & Giustino, PRL 105,
  037002 (2010); per R02). Never made.
- **Fullerides:** excluded, because their pairing mechanism is debated.

### 3.2 Predicted ambient-pressure conventional superconductors

| compound | predicted Tc (K) | E_hull / stability | source / status |
|---|---|---|---|
| LiMoN₂ | ~38–40 | near hull | Cerqueira et al. Adv. Mater. 2024; Gao 2025: intrinsic defects destroy superconductivity |
| Mg₂XH₆ (X = Rh, Ir, Pd, Pt), Sanna | 45–80; >100 with electron doping (Pt) | claimed thermodynamically stable | Sanna et al., npj Comput. Mater. 10, 44 (2024), arXiv:2310.06804 |
| ML-searched hydrides (>150k compounds) | ~50 systems > 20 K, "some even reaching above 70 K" | mostly slightly above hull | Cerqueira, Fang, Errea, Sanna & Marques, Adv. Funct. Mater. (2024), doi:10.1002/adfm.202404043 |
| ClB₂C₈ | 55 (λ ≈ 3.3, ω_log 425 K) | metastable | Gao et al. 2025 |
| Al₂OsH₇ | 38 (λ ≈ 3.3, ω_log ≈ 300 K) | metastable | Gao et al. 2025 |
| **Mg₂IrH₆** | **160–175** (anisotropic Migdal–Eliashberg, μ* 0.10–0.16); λ ≈ 2.3–2.5 | **≈ 60 meV/atom** at 0 GPa | Dolui et al., PRL 132, 166001 (2024); E_hull from Hansen et al. arXiv:2406.09538 and Sinha et al. arXiv:2602.23675 (per R02) |
| Li₂AuH₆ | ~140 (Ouyang et al., PRB 111, L140501 (2025)); 116 Eliashberg / 91 SCDFT (Gao 2025) | 172 meV/atom | Ding, Chen & Shi (arXiv:2604.12367; published, DOI 10.1103/w7rx-sknt, per R02): kinetic instability, H partly dimerises and diffuses, Tc → **22 K** |
| Li₂AgH₆ | 109 Eliashberg / 83 SCDFT; λ ≈ 4 | 319 meV/atom | Ding et al.: **lattice collapse** |
| Li₂MgH₁₆ | 473 harmonic (Sun, Lv, Xie, Liu & Ma, PRL 123, 097001 (2019)); 277 at 260 GPa, superionic, stochastic path integral (Chen & Shi, PRB 109, L140505 (2024), arXiv:2305.04875) | needs ~250 GPa | not an ambient material; listed only for the λ and anharmonicity discussion |
| LaBH₈, BaSiH₈ | ~120 (LaBH₈) | metastable down to ~40 GPa / ~3 GPa | cited in Errea 2022 |
| "LiMgH" (named in the brief) | — | — | not found; no claim made |

**Mg₂IrH₆ synthesis attempts**
- Up to ~28 GPa, only the charge-balanced insulator Mg₂IrH₅ forms (Hansen et al.).
- Above ~40 GPa, cubic insulating Mg₂IrH₇ forms. It persists down to ~20 GPa on
  decompression, then reverts to Mg₂IrH₅ (Sinha et al. 2026).

**Pattern.** Hydrides prefer charge-compensated, semiconducting compositions. The metallic
high-Tc stoichiometry is a doped neighbour of a stable insulator.

### 3.3 Bounds and ceiling papers

**Gao, Cerqueira, Sanna, Fang, Dangić, Errea, Wang, Botti & Marques, Nat. Commun. 16, 8253
(2025), arXiv:2502.18281.** Electron–phonon calculations for more than 20 000 metals.
- ω_max reaches ~5400 K in hydrides, but ω_log "rarely exceeds 1800 K".
- λ and ω_log trade off against each other (they are two moments of the same α²F).
- Tc against E_hull forms a steep Pareto front.
- They judge the "optimal" α²F (a narrow high-frequency peak with λ ≈ 2) to be unphysical.
  McMillan at λ = 2, ω_log = 1800 K, μ* = 0.1 gives 260 K; I reproduce 259 K.
- 100–120 K is the Tc of hypothetical compounds one "can design" from this space. It is not a
  derived limit.
- They call room temperature "extremely unlikely".
- The calculations are isotropic and harmonic. No E_hull cutoff is stated.

**Semenok, Altshuler & Yuzbashyan, arXiv:2407.12922** (journal ref **UNVERIFIED**).
- λ < λ*: 3.69 for Einstein spectra (the minimum over shapes) and 4.72 for Debye spectra. The
  abstract says only "λ ≲ 4".
- Tc < 0.20·ω_max for stable phases and < 0.32·ω_max for metastable ones.
- A heuristic √A form (with ω₀ = 3000 K, from metallic hydrogen) turns this into a bound by
  the mass number A of the lightest element:

| lightest element | H | Li | Be | B | C |
|---|---|---|---|---|---|
| stable bound, 600 K/√A | 600 | 227 | 200 | 181 | 173 |
| metastable bound, 950 K/√A | 950 | 359 | 317 | 286 | 274 |

- **Corrected consequences:**
  - The √A form excludes 300 K only for stable non-H phases and for metastable B- or
    C-carried phases. It is a heuristic and underestimates ω_max for stiff covalent lattices.
  - The ω_max form allows 300 K whenever ω_max ≥ 1500 K (stable) or ≥ 940 K (metastable).
    MgB₂ (ω_max ≈ 1624 K) gives 325 K, and diamond (≈ 1918 K) gives 384 K stable and 614 K
    metastable.
  - **So the bound does not exclude non-hydrides; it excludes soft lattices.** Semenok et al.
    restrict their "only H/D reach room temperature" statement to stable phases.

**Trachenko, Monserrat, Hutcheon & Pickard, arXiv:2406.08129** (JPCM DOI
10.1088/1361-648X/adbc39, **UNVERIFIED**). Fundamental constants bound Tc at roughly
10²–10³ K. This does not restrict a 300 K target.

**Allen & Cohen, PRL 29, 1593 (1972).** Strong coupling softens phonons. This is the origin of
the old ~30 K folklore, which MgB₂ and the hydrides disproved.

### 3.4 Where the T1 targets sit: the 300 K iso-Tc curve (computed, Appendix A)

ω_log (K) needed for Tc = 300 K by Allen–Dynes:

| λ | 2.0 | 2.5 | 3.0 | 3.5 | 4.0 |
|---|---|---|---|---|---|
| μ* = 0.10, ω₂/ω_log = 1.2 | 1757 | 1451 | 1258 | 1123 | 1022 |
| μ* = 0.10, ω₂/ω_log = 1.5 | 1699 | 1382 | 1179 | 1036 | 929 |
| μ* = 0.15, ω₂/ω_log = 1.2 | 2038 | 1661 | 1431 | 1272 | 1155 |

The band is about ±15 %. For example, (λ = 2.5, ω_log = 1250 K) gives only about 260 K.

No published ambient-pressure data point lies on or above this curve:
- The highest-ω_log compounds (~1800 K, e.g. NaNiH₃) have small λ.
- The two λ ≈ 3.3 examples (ClB₂C₈, Al₂OsH₇) have ω_log ≤ 425 K. That is two examples, not a
  law.
- Inverting Allen–Dynes, Mg₂IrH₆ (λ ≈ 2.4, Tc 160–175 K) implies ω_log ≈ 800–900 K (derived,
  approximate).

### 3.5 Expectation vs. impossibility (separated, per R02-D1)

**Prior expectation for this pipeline's outputs: ≤ ~150 K.**
- Every computed ambient dataset tops out at about 40 K on the hull, ≤ 80 K near it, and
  100–175 K for strongly metastable phases. The last group is unsynthesised, and PIMD has
  partly destroyed it.
- Gao et al. call room temperature "extremely unlikely".
- Our generator, MatterGen (trained on Alex-MP), samples a similar prior.

**No fundamental bound excludes 300 K**:
- for metastable H-based compounds;
- for metastable Li- and Be-based compounds, by the metastable Semenok forms;
- for any stiff light-element framework with ω_max ≥ 940 K, by the ω_max form.

**The practical-ceiling evidence is weaker than it looks.**
1. Cerqueira (×2), Sanna and Gao share one group, one workflow and one Alexandria-derived
   structure space.
2. The searches sample prototypes, small cells and near-hull compositions.
3. The calculations are isotropic, harmonic and adiabatic. Isotropic MgB₂ comes out about 2×
   too low, and anharmonicity can raise or lower Tc.
4. The ambient frontier moved from 39 K (measured) to 160 K (predicted) in a few years, so the
   upper edge has not converged.

---

## 4. The stability problem

### 4.1 Thermodynamic windows

**What known compounds look like.** Sun et al., Sci. Adv. 2, e1600225 (2016), over 29 902 ICSD
phases:
- Median E_hull 15 meV/atom; 90th percentile 67 meV/atom.
- Nitrides: median 63 meV/atom; 90th percentile 195 meV/atom.
- These are phases made at or near ambient conditions. Pressure-quenched phases are absent.

**Thresholds in use**
- Dolui et al. 2024 keep candidates "close to the convex hull" at 1 and 10 GPa. My earlier
  "50 meV/atom" attribution came from an automated summary and is not in their text (R02-A5).
  **It is withdrawn.**
- Another hydride search (arXiv:2401.17024) uses < 80 meV/atom (per R02).
- MatterGen's S.U.N. metric conventionally uses 0.1 eV/atom (sibling agent to confirm).

**Zero-point energy.** DFT hulls for hydrides leave out zero-point energy (ZPE). ZPE differs
between competing H phases by tens of meV/atom (order of magnitude; **UNVERIFIED** for specific
systems).

**Where the high-Tc candidates sit.** The 100–175 K candidates are at 60, 172 and 319 meV/atom,
and Tc rises with E_hull along a Pareto front (Gao 2025).

**Conclusion**
- E_hull is a **Pareto axis**, not a necessary condition.
- The outer post-filter is ≤ 0.35 eV/atom, ZPE-corrected for H-containing phases.
- Above about 0.1 eV/atom, the candidate needs a plausible synthesis route, for example:
  - a doped neighbour of a stable insulator;
  - a pressure-quench parent;
  - topotactic H insertion.

### 4.2 Dynamical stability and the coupling–instability tension

Schematically (Allen & Cohen 1972; Moussa & Cohen 2006; Semenok et al.):

  Mω² = k_bare − c·η, so λ = η/(k_bare − c·η).

λ therefore grows fastest next to a soft-mode instability. In a harmonic calculation that
instability shows up as imaginary modes. In reality it shows up as a gap-opening distortion or
as decomposition.

**Practical rules**
1. Require no imaginary modes on the quantum energy surface, computed with the stochastic
   self-consistent harmonic approximation (SSCHA) or PIMD. Treat harmonic imaginary H modes as
   "undecided", not rejected.
   - SSCHA stabilises H₃S (Im-3̄m) and LaH₁₀ (Fm-3̄m) about 100 GPa below the classical
     instability (Errea 2022).
   - In Mg₂IrH₆, quantum anharmonicity hardens the 30–80 meV modes and softens modes above
     150 meV by about 15 meV (Dolui et al.).
2. Flag a λ dominated by modes below ~550 K. Those modes add little to A and signal an
   approaching instability (Sec. 1.4).
3. Keep two different downward effects apart:
   - **Anharmonic renormalisation of the same structure** (SSCHA-type). Typically −10 to −40 %
     on harmonic Tc for H-rich phases. This range is **UNVERIFIED** as a general rule; the
     effect is compound-dependent and can change sign.
   - **Kinetic instability, where the structure itself changes.** Examples:
     - Li₂AuH₆ forms H₂ (Tc 140 → 22 K).
     - Li₂AgH₆ collapses.
     - Superionic Li₂MgH₁₆ drops from 473 to 277 K at 260 GPa.

     All three results come from one group using one family of methods (the Shi group's
     stochastic path-integral and PIMD). Independent confirmation is pending.
4. Test kinetic persistence at the operating temperature separately from dynamic stability at
   0 K. Run MD or PIMD at ≥ 300 K and look for H diffusion and H₂ formation.

### 4.3 Kinetic protection and chemical precompression

- **Covalent cages** (B–C clathrates such as SrB₃C₃) survive quenching from 50 GPa because
  sp³ bonds give large barriers. They reach high ω, but so far λ < 1.
- **[MH₆]ⁿ⁻ complexes** (K₂PtCl₆-type) give rigid H environments at 1 atm. They are insulating
  at the stable electron count, so the metallic versions are doped and metastable.
- **Chemical precompression** (Hilleke & Zurek 2022; arXiv:2205.11569) has lowered the
  stability pressure of clathrate hydrides to tens of GPa, not yet to 0–1 GPa.

---

## 5. Materials target profile

**Use labels**
- **HARD:** a necessary condition within the stated model assumptions. Applied as a
  post-filter at the stage where it can be evaluated. Never applied implicitly through
  MatterGen.
- **SOFT:** a ranking feature and/or a MatterGen conditioning prior. Never a reason to reject.
- **MG (MatterGen conditioning).** Only the installed checkpoints, each used on its own:
  - single property: `chemical_system`, `dft_band_gap`, `dft_mag_density`, `ml_bulk_modulus`,
    `space_group`;
  - pairs: `chemical_system_energy_above_hull` (the only way to condition on E_hull) and
    `dft_mag_density_hhi_score`.

  No checkpoint combines band gap, magnetism, E_hull and space group. Any other combination
  needs fine-tuning first. H-rich metals are sparse in the Alex-MP training data, so
  conditioning toward high H fraction is out of distribution.
- **PF (post-filter):** MatterSim pre-screen, then a DFT, DFPT or EPW verdict.
- **RF:** ranking feature.

### 5(a) HARD necessary conditions (within the stated assumptions)

| # | condition | threshold | assumption / justification | use |
|---|---|---|---|---|
| H1 | Metallic | N(E_F) > 0 at the PBE level; check with a hybrid functional if the gap is marginal. **No** clause on distance to a gap edge; E_F near a gap is only *flagged* "doping-sensitive" (Mg₂IrH₆ has this) | bulk electron–phonon pairing needs a Fermi surface | MG: `dft_band_gap` ≈ 0 for one budget share; PF: DFT |
| H2 | Non-magnetic ground state | moment ≈ 0 after spin-polarised relaxation | our Tc machinery (isotropic singlet Eliashberg, certified bounds) assumes a non-magnetic normal state, and magnetism suppresses singlet pairing. This is a modelling restriction, not a law of nature (coexistence exists at low Tc). No element is excluded: Mg₂FeH₆-type and NaNiH₃ parents are allowed | MG: `dft_mag_density` = 0, or `dft_mag_density_hhi_score` with HHI left free; PF |
| H3 | Spectral weight | single-band: A = ∫α²F ≥ 1300 K. Multiband (distinct coupled sheets): the band-resolved, κ-based certified bound ≥ 300 K replaces A | certified bound k_B·Tc ≤ A/4.33148 (T1, `rtt.core.certify`); κ for multiband (`rtt.core`, wMB) | PF (DFPT/EPW); BETE-NET only for triage |
| H4 | Coupling within the lattice-stability limit | λ ≤ 4.7; λ > 3.7 flagged "near-instability" | Semenok λ* = 3.69 (Einstein, the minimum) to 4.72 (Debye) | PF |
| H5 | Frequency scale | from H3 and H4: ⟨ω⟩ ≥ 2600 K/λ ≥ ~550 K. ω_max ≥ 940 K if metastable, ≥ 1500 K if stable | first-moment identity plus λ*; Semenok Tc < 0.32 or 0.20 ω_max | PF; MatterSim phonons as a pre-screen on ω_max |
| H6 | Dynamically and kinetically stable | no imaginary modes on the quantum energy surface (SSCHA for H-containing phases); survives MD or PIMD at ≥ 300 K without H₂ formation, collapse or decomposition by diffusion | the material has to exist; harmonic instabilities count as "undecided", not rejected | PF: MatterSim harmonic first, then SSCHA or PIMD for finalists |
| H7 | Thermodynamic outer window | E_hull ≤ 0.35 eV/atom (ZPE-corrected for H); synthesis route stated above 0.1 | a generous outer bound. Li₂AgH₆ (0.319) is the most metastable published candidate; beyond this window, persistence is implausible | PF; inside the window E_hull is soft (S5) |

H1–H7 do not mention hydrogen. H3 and H5 together force light atoms or very stiff bonds that
carry strong coupling.

### 5(b) SOFT priors (ranking and conditioning weights)

| # | prior | preferred direction / reference values | justification and caveat | use |
|---|---|---|---|---|
| S1 | H present and carrying the coupling | strong prior: ~80–90 % of the budget on H-containing systems | in every published dataset, only H reaches A ≥ 1300 K with λ ≤ 3. This is not a bound (Sec. 3.3) | MG: `chemical_system` or `chemical_system_energy_above_hull` over enumerated X–Y–H systems |
| S2 | H_f (H atom fraction) | score rises with H_f; soft floor 0.5 in the H channel only. Reference: A₂MH₆ 0.667, AMH₃ 0.60, H₃S 0.75, LaH₁₀ 0.91 | the Belli arithmetic (Sec. 2.1) is untested at 1 atm and contradicted by Mg₂IrH₆. Priority family 1 (H_f 0.60–0.67) stays as the best existing family and is not filtered out | RF; MG indirectly through the composition in `chemical_system` |
| S3 | H share of N(E_F), and η_H | H_DOS higher is better (≥ 0.4 good, ≥ 0.6 better). η_H ≳ 10 eV/Å² per coupled H, which is robust in the single-band model when H carries the coupling | Belli; Sec. 1.3. Becomes effectively hard only once the coupling is shown to be H-carried | RF (DFT projected DOS, Gaspari–Gyorffy) |
| S4 | Bonding motif | favour any of: weak-covalent H–H network (1.0–1.55 Å), M–H–M bridges, **σ-bonded [MHₙ] complexes with H-1s/M-e_g states at E_F**. Penalise intact H₂ (0.74 Å) and fully ionic closed-shell H⁻ | Belli's classes (high pressure) plus the Mg₂IrH₆ counterexample (isolated octahedra, 160–175 K) | RF: geometry is cheap; ELF φ needs DFT |
| S5 | E_hull (inside the H7 window) | Pareto axis against predicted Tc. Reference: Sun 2016 median 15 / 90th percentile 67 meV/atom; Mg₂IrH₆ 60, Li₂AuH₆ 172, Li₂AgH₆ 319 meV/atom | conditioning on low E_hull reproduces the 40–80 K near-hull region that has already been searched | MG: `chemical_system_energy_above_hull`, with the conditioning value **sampled across 0–0.3 eV/atom**; RF |
| S6 | Space group | spread the budget over ≥ 10 groups (191, 194, 139, 123, 166, 221, 223, 225, 227, 229) plus an unconditioned share. **No symmetry post-filter** | high symmetry puts degenerate states and vHs at E_F (MgB₂ 191, H₃S 229, LaH₁₀ 225, Mg₂IrH₆ 225, CeH₉/YH₉ 194, Li₂MgH₁₆ 227, YH₄ 139, LiMoN₂ 166). The evidence is biased by which structures were studied, and symmetric phases tend to distort | MG: `space_group`, one group per run |
| S7 | Electron count | a metallic composition next to a closed-shell insulating parent (doped Mg₂IrH₅ or Mg₂IrH₇, σ-hole MgB₂, doped LiBC). E_F near a gap edge is expected here and flagged doping-sensitive | every ambient high-Tc candidate is a doped insulator | RF |
| S8 | Volume per H | lower scores better, relative to an **ambient reference set**; no threshold. LiH ≈ 17 Å³/H and Mg₂FeH₆ ≈ 11 Å³/H (R02, from lattice constants); Mg₂IrH₆ similar, not computed. Compressed H₃S ≈ 4.9 for contrast | the old ≲ 5–6 Å³/H threshold was a megabar value and is withdrawn | RF only |
| S9 | Stiffness | bulk modulus, as a tie-breaker for covalent frameworks only | stiff sp³ frameworks give kinetic protection, but in hydrides the soft H sublattice is what produces λ | MG: `ml_bulk_modulus` for the non-hydride share only; RF |
| S10 | Element palette | favoured: H plus electropositive donors (Li, Na, Mg, Ca, Sr, Ba, Al, Sc, Y, La) plus late 4d/5d metals (Rh, Pd, Ir, Pt, Os, Au, Ag). No element is excluded; H2 handles magnetism. Do **not** condition on low HHI, because Ir, Pt and Os score poorly | the composition of every top ambient candidate | MG: enumeration of `chemical_system` |
| S11 | Predicted (λ, ω_log) position | distance to the 300 K iso-Tc curve (Sec. 3.4, ±15 % band) | this is the actual objective | RF, after DFPT |

**Priority families inside the H channel**
1. **Metallic A₂MH₆ / AMH₃ complexes** (the Mg₂IrH₆ / Li₂AuH₆ family). This is the best
   existing ambient family. By the Belli arithmetic it is probably a ≤ ~200 K calibration and
   extension family rather than a direct route to 300 K.
2. **Light-metal hydrides with H-network or clathrate fragments**, held together by chemical
   precompression. This family has the highest ceiling and the worst stability.

### 5(c) Minority non-hydride channel (10–20 % of the generation budget)

No proven bound excludes metastable Li- or Be-carried 300 K superconductivity (Semenok
metastable form: Li 359 K, Be 317 K). Nor does the ω_max form exclude stiff B, C or N frameworks
with ω_max ≥ 940 K.

The obstacle is coupling strength. Reaching A ≥ 1300 K needs λ ≈ 2 at ⟨ω⟩ ≈ 1300 K. The best
known values are λ 0.87 (MgB₂) and a predicted λ ≈ 1.3 (p-doped graphane). Isotropic screens
are biased low for these multiband materials (Sec. 1.4).

| family | examples / rationale | use |
|---|---|---|
| Multiband σ-hole layered borides and carbides, with stiffer bonds than MgB₂ | MgB₂ analogues, hole-doped LiBC, Li₂B₃C, hole-doped graphane-like sheets | MG: `chemical_system` (Li–B–C, Mg–B–C, Be–B, Li–Be–B); PF using **κ (multiband)**, not A |
| Metallized sp³ B–C–N frameworks and clathrates | SrB₃C₃-type, ClB₂C₈, heavily B-doped carbon | MG: `chemical_system`; `ml_bulk_modulus` (high); PF |
| Li- and Be-rich intermetallics and borocarbides | phonons carried by light cations, up to ~1000 K | MG: `chemical_system`; PF |

The hard conditions H1–H7 apply unchanged, with κ in place of A in H3 for multiband cases.

### 5(d) Generation plan implied by the available checkpoints

1. **~60 %:** `chemical_system_energy_above_hull` over enumerated X–Y–H systems (S1, S10). The
   E_hull conditioning value is sampled uniformly over 0–0.3 eV/atom (S5).
2. **~15 %:** `space_group` runs over the S6 list.
3. **~10 %:** separate `dft_band_gap` = 0 and `dft_mag_density` = 0 runs, to enrich for
   metallic and non-magnetic structures.
4. **10–20 %:** the non-hydride channel (5c).
5. **A small unconditioned share**, as a bias control.

Everything else is a post-filter or a ranking feature. Joint conditioning would need
fine-tuning; that prerequisite is flagged for the MatterGen audit agent.

---

## 6. Open problems and unverified items

**Not verified**
- Journal references: Semenok et al.; Moussa–Cohen.
- The Trachenko JPCM DOI.
- Page numbers for McMillan, Hopfield, Allen–Dynes, Nagamatsu and Ekimov.
- Details of the Duan review.
- The negative LiBC experiments.
- The Rb₀.₅Sr₀.₅B₃C₃ value of 78 K.
- The source of the Li₃B₄C₂ prediction.
- The "−10 to −40 % anharmonic" range, which is a rough rule.

**Taken from R02 without my own re-fetch**
- The Ding et al. DOI.
- Chen & Shi, PRB 109, L140505.
- The arXiv:2401.17024 threshold of 80 meV/atom.
- Savini et al. (graphane).
- Mg₂IrH₆ at 60 meV/atom.

**Not computed (need DFT and/or the structure file)**
- ⟨I_H²⟩ ≈ 25 eV²/Å² is an order-of-magnitude assumption. Absolute η_H for H₃S and Mg₂IrH₆
  should be computed to calibrate S3.
- φ and the volume per H for Mg₂IrH₆.

**Other open questions**
- The Belli fit is untested at 1 atm.
- The Gao dataset's E_hull cutoff is not stated.
- The multiband κ machinery is taken from the multiband notes of `rtt.core` and was not re-checked here.

---

## Appendix A — McMillan / Allen–Dynes grid and iso-Tc solver (computed, POSIX awk)

```sh
awk 'function ad(l,wl,w2,mu){L1=2.46*(1+3.8*mu);L2=1.82*(1+6.3*mu)*(w2/wl);
 f1=(1+(l/L1)^1.5)^(1/3);f2=1+(w2/wl-1)*l*l/(l*l+L2*L2);
 return f1*f2*wl/1.2*exp(-1.04*(1+l)/(l-mu*(1+0.62*l)))}
function mm(l,wl,mu){return wl/1.2*exp(-1.04*(1+l)/(l-mu*(1+0.62*l)))}
function solve(l,r,mu, lo,hi,m,i){lo=100;hi=5000;for(i=0;i<60;i++){m=(lo+hi)/2;
 if(ad(l,m,r*m,mu)<300) lo=m; else hi=m}; return m}
BEGIN{mu=0.1; for(l=1.5;l<=4.01;l+=0.5) for(wl=600;wl<=1800;wl+=300)
 printf "lam=%.1f wlog=%d McM=%.0f AD=%.0f\n",l,wl,mm(l,wl,mu),ad(l,wl,1.2*wl,mu);
 for(l=2;l<=4.01;l+=0.5) printf "lam=%.1f wlog300: %.0f (r1.5) %.0f (mu.15) %.0f\n",
  l,solve(l,1.2,0.1),solve(l,1.5,0.1),solve(l,1.2,0.15)}'
```

Selected outputs (K), McMillan / Allen–Dynes:
- λ = 2.0, ω_log = 1800: 259 / 307. This matches Gao's 260 K.
- λ = 3.0, ω_log = 1200: 216 / 286.
- λ = 4.0, ω_log = 900: 181 / 264.

The iso-Tc values in the Sec. 3.4 table come from the solver; R02 reproduced them
independently.

---

## Response to review R02

| finding | change made / response |
|---|---|
| **C1** (critical): "H required", "non-hydrides cannot reach 300 K", and the 325 K self-contradiction | **Accepted.** "H required" is no longer a condition; H is the soft prior S1. Sec. 3.3 now states the correct consequences: the metastable √A form allows Li and Be (359 and 317 K); the ω_max form allows any phase with ω_max ≥ 940 K; MgB₂ (325 K) and diamond (384/614 K) are not excluded; and Semenok et al. restrict their claim to stable phases. Added the non-hydride channel 5(c) with 10–20 % of the budget, as a real exploration share rather than controls. Sec. 3.5 separates expectation from impossibility. |
| **C2** (critical): row 8 rejects Mg₂IrH₆ and Li₂AuH₆ | **Accepted.** The bonding motif is now ranking feature S4, and σ-bonded [MHₙ] complexes are explicitly allowed. Sec. 1.2 and 2.1 state that Mg₂IrH₆ falsifies Belli's isolated-class cap at ambient pressure. |
| **C3** (critical): ≲ 5–6 Å³/H is a megabar number | **Accepted.** The threshold is deleted. S8 is a ranking feature against ambient reference values (LiH 17, Mg₂FeH₆ ≈ 11 Å³/H), with H₃S at 4.9 shown for contrast. |
| **C4** (critical): conditioning MatterGen on E_hull ≤ 50 meV reproduces the ceiling | **Accepted.** E_hull is now a Pareto axis (S5), with the conditioning value sampled over 0–0.3 eV/atom. The hard outer window is ≤ 0.35 eV/atom, ZPE-corrected (H7), with a synthesis route required above 0.1. The "Dolui 50 meV" attribution is withdrawn. ZPE and pressure-quench caveats are added to the Sun 2016 statistics. |
| **C5** (major): five space groups as a hard condition | **Accepted.** S6 spreads the budget over 10 groups plus an unconditioned share, with no symmetry post-filter. The section notes the bias in the evidence and that each checkpoint run takes one space group. |
| **C6** (major): rows 1 and 13 contradict each other | **Accepted.** The gap-edge clause is gone from H1 and is now a "doping-sensitive" flag. S7 keeps the doped-insulator prior. |
| **C7** (major): Belli used as a necessary condition; family 1 fails its own H_f target | **Accepted.** H_f and H_DOS are ranking features (S2, S3), with only a 0.5 soft floor inside the H channel. Added how the floor depends on φ (0.51 / 0.64 / 0.73) and the extrapolation and transferability caveats. Family 1 is now labelled a ≤ ~200 K calibration and extension family, which is consistent with the Belli arithmetic (A₂MH₆ needs φ ≳ 0.77). |
| **C8** (major): multiband caveat missing | **Accepted.** Added in Sec. 1.4 and in H3 (κ replaces A when there are distinct coupled sheets). Noted that the isotropic ceiling datasets are biased low. Added a multiband σ-hole family to 5(c). |
| **C9** (major): λ capped at 3.7 | **Accepted.** H4 is now λ ≤ 4.7, with λ > 3.7 flagged as near-instability. |
| **C10** (major): MatterGen conditioning column overstated | **Accepted.** Section 5 lists only the installed single-property checkpoints and the two pair checkpoints. It states that no joint band gap + magnetism + E_hull + space group model exists and that H-rich metals are out of distribution. Added a concrete generation plan, 5(d). Fine-tuning is flagged as a prerequisite for any joint conditioning. |
| **D1** (major): ≤ 150 K overstated as a ceiling | **Accepted.** Sec. 3.5 and summary item 6 keep the prior expectation but separate it from the statement "no fundamental bound excludes 300 K for metastable H-, Li- or Be-based compounds". They list four caveats: non-independent searches, sampling bias, isotropic and harmonic methods, and a moving frontier. Gao's 100–120 K is described as the Tc of "compounds one can design", not a limit. |
| **A8** (major): Li₂MgH₁₆ misattributed | **Accepted.** Sun et al., PRL 2019, is the 473 K harmonic prediction. Chen & Shi, PRB 109, L140505 (2024), arXiv:2305.04875, is the 277 K superionic path-integral result. Reclassified as a kinetic/superionic effect, not "anharmonicity" (Sec. 4.2, rule 3). |
| **A9** (major): arXiv:1307.1323 misdescribed | **Accepted.** It is now cited as Gao, Lu & Xiang, PRB 91, 045132 (2015), on Li₂B₃C only, claiming Tc "might be even higher than" MgB₂. No 50 K figure. Li₃B₄C₂ is marked as coming from a paper I did not locate. |
| A2: S formula and the "more stable → lower DOS" quote not verified by the reviewer | **Partly disagree, with evidence.** While writing revision 1, I extracted the text of the arXiv:1910.00385 PDF with pdftotext. Eq. 8 reads S = Tc/√(Tc,MgB₂² + P²), with S(MgB₂) = 1, S(390 K, 0 GPa) = 10, S(H₃S) = S(LaH₁₀) = 1.3 and S(HgBaCaCuO) = 3.5. The sentence "more stable structures or compositions will typically have lower densities of electronic states at the Fermi level" is in the same text. I keep both and state the extraction method. I added the recomputed LaH₁₀ value of 1.43. |
| A5: Dolui 50 meV | Withdrawn (Sec. 4.1); replaced by Dolui's own wording, "close to the convex hull". |
| A6, A7, A10, A11, C14 (small fixes) | Mg₂IrH₆ ≈ 60 meV/atom: UNVERIFIED tag dropped, now cited to Hansen and Sinha. Ding et al. DOI added. Sanna moved to "claimed stable, ≤ 80 K". Cerqueira (Adv. Funct. Mater.) reworded to "some even reaching above 70 K". S target corrected to 7.69. H_DOS threshold corrected to 0.935. "Avoid magnetic 3d metals" removed (S10); magnetism is handled by H2. Semenok journal ref marked UNVERIFIED. |
| C11 (minor): "≥ 50 % of λ in H bands above 1000 K" | Replaced by the A/κ condition (H3) plus ω_max flags (H5). |
| C12 (minor): wrong "corner" | Replaced by the iso-Tc table with μ* and ω₂ sensitivity (Sec. 3.4). The "λ ≳ 3 with low ω_log" observation is now labelled as two examples. |
| C13 (minor): mislabelled corrections | Sec. 4.2 rule 3 now separates anharmonic renormalisation from kinetic instability and notes that the kinetic results come from one group and one family of methods. |
| D2 (minor): non-hydride failures are synthesis failures, not physics | **Accepted.** Stated in Sec. 1.2 and used to justify 5(c). |
| A1 notes: Belli inconsistencies (1.55 Å as H–host vs H–H; Si₂H₆; BeH₂) | Noted. I use only the H–H reading, and only as a ranking feature, so these inconsistencies no longer affect any filter. |
