# 34. Campaign 2: MatterGen lanes G2, G3, G4-bulk (partial; closed 2026-11-01)

The worker that ran these lanes stopped (usage limit) after S3, before S4 and before its
analysis. This note records what finished. Driver: `scripts/campaign2_g_lanes.py`; report helper
`scripts/campaign2_g_lanes_report.py`; checkpoint check `scripts/g_lanes_ckpt_check.py`.

## Checkpoints

`space_group`, `dft_band_gap`, `dft_mag_density`, `ml_bulk_modulus` were downloaded from Hugging
Face `microsoft/mattergen` (user-approved). `checkpoint_sha256` was fixed: on the run machine
the HF blob is a symlink into a second store whose file name is not the content sha256, so the
first 64-hex name on the symlink chain is used (`src/rtt/source/mattergen.py`, test in
`tests/front/test_lanes_mattergen.py`).

## Funnel (own ledger, 577 candidates)

| Lane | Generated | S1 PASS |
|---|---|---|
| G2 (space group, S6 list) | 320 | 213 |
| G3 (band gap 0 / magnetic density 0) | 192 | 107 |
| G4-bulk (ml_bulk_modulus 300 GPa) | 64 | 64 |
| control (Mg₂IrH₆) | 1 | 1 |

S2: 380 PASS, 3 INVALID, 2 DUPLICATE. S3: 41 PASS (rank cut 40 + control), 1 SCREENED_OUT,
338 in the rank queue.

## Reading

These lanes condition on properties (space group, zero gap, zero moment, stiffness) that do not
select hydrogen-rich chemistry: in the last 64-structure batch only one cell contained H. They are
useful as metallicity and structure-type diversity lanes, not as a route to A ≥ 1300 K. No S4/S5
was run on them; their S3 survivors are in the ledger for any later run.

NOT CHECKED: novelty statistics, metallicity at DFT level, any phonon property.
