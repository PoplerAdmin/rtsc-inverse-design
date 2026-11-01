# Electron–phonon label sources (agent D1-labels, 2026-10-29)

Raw data are kept outside the repo in `~/rtsc_targeted_work/D1-labels/raw/`. Only manifests and
summaries are in `data/labels/`. Scripts are in `src/ml/labels/`. The design note is
`docs/design/14_property_conditioning.md`.

| # | Dataset | URL / DOI | Licence | Version / date | What it provides | Used? |
|---|---|---|---|---|---|---|
| S1 | Cerqueira, Sanna, Marques, "Sampling the materials space for conventional superconducting compounds" (arXiv:2307.10728; *Adv. Mater.* 36, e2307085 (2024)) | Materials Cloud Archive record `3kbt5-r3n56`, DOI 10.24435/materialscloud:qv-bq, MCID 2023.163 | CC BY 4.0 (record metadata) | published 2023-10-30, record revision 4 | 8,249 QE DFPT calculations, all at 0 GPa (PBEsol, vc-relax with no `press`): relaxed CIF, α²F(ω) for 10 smearings, λ, ω_log, McMillan/Allen–Dynes/Eliashberg Tc, phonon dispersion | **yes, primary** |
| S2 | JARVIS-DFT EPC, Choudhary & Garrity, *npj Comput. Mater.* 8, 244 (2022) | figshare file 38307921 (`jarvis_epc_data_figshare_1058.json.zip`), article DOI 10.6084/m9.figshare.21370572 ("JARVIS-SuperconDB"); jarvis-tools key `supercon_3d` | CC BY 4.0 (figshare metadata) | figshare article version 4; jarvis-tools 2026.6.12 | 1,058 entries, all `press = 0 GPa`: structure, 50-bin α²F in meV, λ, ω_log, Tc, stable/unstable flag | **yes** (stable entries only) |
| S3 | JARVIS `supercon_2d` (2D EPC set) | figshare file 38950433, same article as S2 | CC BY 4.0 | as S2 | 161 2D monolayers | downloaded, **not used** (2D slabs with vacuum are outside MatterGen's bulk domain) |
| S4 | JARVIS `alex_supercon` ("Alexandria_SuperConDB 8253") | figshare file 55551155; jarvis-tools cites doi:10.1002/adfm.202404043 | CC BY 4.0 per the JARVIS figshare listing (article 10.6084/m9.figshare.27174897; that this file belongs to that article is **UNVERIFIED**) | jarvis-tools 2026.6.12 | 8,253 entries: structure, λ, ω_log, Tc, DOS(E_F), Debye T. No α²F | **cross-check only.** 8,248 ids match S1 one-to-one; λ equals the S1 `McMillan.dat` λ in 100 % of matches. It is the same calculation set, so it adds no labels |
| S5 | Bercx, Poncé, …, Marzari, "Charting the landscape of BCS superconductors in experimentally known compounds" (doi:10.1103/sb28-fjc9) | Materials Cloud `0kqen-w8h90`, DOI 10.24435/materialscloud:x2-ap, MCID 2025.175 | CC BY 4.0 | 2025-11-04 | 4,533 QE screenings plus 250 EPW α²F. Labels exist **only inside AiiDA archives** (15 GB + 8.4 GB) | **not used.** Only the four Table_*.csv files were fetched; they hold no labels. Future work: worth adding, since these are experimentally known structures |

Not used:
- **High-pressure hydride sets** (e.g. the source of the 953-hydride ALIGNN set): out of scope by
  the ≤ 1 GPa training rule.
- **BEE-NET training data**: open release not located (**UNVERIFIED**).
- **SISSO α²F set** (Materials Cloud `skdzy-kyx54`): no crystal structures.

## Checksums

S1 was streamed (`curl | fetch_cerqueira.py`) and never stored whole. The md5 of each stream equals
the archive's published md5 for all 17 files. sha256 of each stream:

| file | bytes | md5 (= published) | sha256 | materials |
|---|---|---|---|---|
| batch-a.tar.bz2 | 3227483011 | b7f75b330e0b6e0573d15f172c322707 | d4f8df81446da10e8cd621bb115eeb48ab37665aab3571c07f36ddb4efda92c4 | 500 |
| batch-b.tar.bz2 | 3438986099 | 63ee110453a8d2588c5271d71134e5a3 | 55cdfae976f51a238520edede895e727d5cecdbebf809250b1531139cd62ab11 | 500 |
| batch-c.tar.bz2 | 3424761873 | 77a3927603fa4654f003742b9f37e740 | b3ab7a4e7f1cfca84ac5b1bff39753be9d0e78a01e5057986f5def7b39a9bdf9 | 500 |
| batch-d.tar.bz2 | 3494234906 | 00d933a8bf6041997c0242ff7dd53f5a | c326744a17fa6b3c3029aec621a19af45d63980257d80fa8f6dc26ae749587a1 | 500 |
| batch-e.tar.bz2 | 3378142502 | 5e4724b49246796a9b8125f8ac0b3467 | dc48314c74be71108d371abdba87bab294a8305c7c398977d54932d1a2ef7197 | 500 |
| batch-f.tar.bz2 | 3242369863 | 96801dd4489b6c3351edee326b58f7bb | dc054c540f46542532f5feabb5d846f6b0ee5ea06cad82119098ecc03367dbef | 500 |
| batch-g.tar.bz2 | 3121983436 | ff32a36bd5e7917b0167a803db073756 | fde1d56fb139b67c576ffe54eba9ba998600b1ae35ba87ad94caee75bfdcf496 | 500 |
| batch-h.tar.bz2 | 3221363829 | c9739508448664a4628116acc2c86912 | 77ef4d1551fc6446af66448f6c66d9f73d6084edc1272588788e1a2c724a459e | 500 |
| batch-i.tar.bz2 | 3028598664 | 0d10ed458c403bd9c31705cc823135d3 | 5279af33170559088c1edaf2125764d0fe3a68b1f153aa1ea37aa6d6508a2bf5 | 500 |
| batch-j.tar.bz2 | 3097509131 | b1ab88ba69dd299b99550700e628d46d | 6e3045b6773bae54d9b40dc0e8743345951d7f15653f61c8698d49da8e46b7c0 | 500 |
| batch-k.tar.bz2 | 3026215208 | d4f879b6d33864777b08380c61cdb8ab | 32c62289434d3346b523fbe3fd45674360b669108669f5d8678ddd8b7bd2c1ed | 500 |
| batch-l.tar.bz2 | 3023483754 | dda7c023c250a343f7de7753de6ed67c | 9d023f5f53bf33e3ba8cf23c0f49b1fb691f41362b0292b2b5335ceed7d731d4 | 500 |
| batch-m.tar.bz2 | 3025352999 | 0bbecf3bf187552fbd4a76e92f926ce3 | b77eb75bc8de897718f15f091f4b4177229244545225068346e0f08430cf5b17 | 500 |
| batch-n.tar.bz2 | 3120868222 | 94adedd49761a0a02d603f02f6b64869 | d7353d5cdfea3fb20a7f10a6907bbd55b1ff47dbd16ec5ffb9d0b4dc1b1c8891 | 500 |
| batch-o.tar.bz2 | 2992368798 | 7a5fe557d8f0aae014a5fad7e1050359 | a63a758536c80c37e6c4d77fb92acaf32a20949d6b6d2970ae83c22c2745ab32 | 500 |
| batch-p.tar.bz2 | 2623380238 | a2b8ca9b8116879ac41a984f738a1980 | cb28fe375f7500349e2112e939b1db463b8d760a09d982cf62a171b1e5b6ee74 | 500 |
| batch-q.tar.bz2 | 1130098806 | f8a9e53bcc73de5a893ffd8134ebc2fd | d51f067e2505df5dcd5e4d8acc1827639b150c3b59122391b7f6edc9bdaf8d65 | 249 |

JARVIS zips (downloaded whole):

| file | md5 | sha256 |
|---|---|---|
| f_38307921.zip (supercon_3d) | bb4a91e899458f34cbb3b9bf03e9bd56 (= figshare computed_md5) | 4213651cee9b4c13a376dd74ac1f8d3bd9efcd451f1b8da615255c8df0e1487f |
| f_38950433.zip (supercon_2d) | 721df1cc069296d02d6653bbed527229 (= figshare computed_md5) | 2f7edceaaf54bd026fc32c3d2985682cb924a2168b4b52cbddbb99416fcfea69 |
| f_55551155.zip (alex_supercon) | 980a8d38210c24925b526f68a059181f | b6807b7002b607eba163a042846563efea7bbb9c64fe8bc516b7a81dcd2e638a |

Derived files (scratch; sha256):

| file | sha256 |
|---|---|
| labels/cerq.jsonl | f5173155f93ce34c29a8b895abd470e22293a91c15e2a2cac958b2368049d95f |
| labels/jarvis.jsonl | 68d17308ce3c6baf05dc97a89222754e1fd644730445fe6c420f1abf23aa6dce |
| labels/alexsc.jsonl | 6bcc24b656f817b66f898d54d4cfc4de2ec95928b6dd905a9a0f975ee6852c76 |
| ds/csv/train.csv | 4f7fc41c067d27f80b8dc9febd0e07fbb4b251c30f1ba4f337d68f02e47ddf26 |
| ds/csv/val.csv | 5a7138c03293d6c0c407a2e3aa526d641be14e382ec0041b223e7b82e05890b1 |
| ds/csv/test.csv | d9c9eff8f25667273508c5672b9b0bdd1cc90259516b2dfaa7ec17a13646e6d1 |
| ds/manifest.tsv (= data/labels/manifest_log10A.tsv.gz, gunzipped) | 71ab4eb30c31fef63ff78e0c8e1fa26108a742463e511a2aab1e0dc2f73d220d |

## Files in the repo

- `data/labels/labels_all.tsv.gz`: all 9,307 S1+S2 entries (8,249 + 1,058), with labels and flags and no
  structures. Unstable entries are included.
- `data/labels/manifest_log10A.tsv.gz`: the 8,675 training rows with their split.
- `data/labels/summary_log10A.json`: counts, distributions, top-40 A.

Attribution (CC BY 4.0): S1 © T. F. T. Cerqueira, A. Sanna, M. A. L. Marques; S2/S3/S4 © JARVIS /
NIST (K. Choudhary et al.). The labels here are derived works (moments recomputed from their α²F).
