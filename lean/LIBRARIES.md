# Frozen Lean libraries

Libraries: `TcR1`, `PipelineLogic` (includes `PipelineLogic.Multiband`, which imports `MultiBand`),
`MultiBand`, plus the axiom-check files `Check.lean`, `CheckPipeline.lean`, `CheckMultiBand.lean`.
These files are frozen: a change to any of them means re-checking the proofs and the axiom report
and updating this table. `lakefile.toml` declares these libraries and `TargetLogic`;
`lean-toolchain` and `lake-manifest.json` pin the toolchain and dependencies.

| file (path relative to `lean/`) | sha256 |
|---|---|
| `Check.lean` | `3d02b3373c96b0963d662061b0398f56c7b390ba7ae5c9f7d0fb70be1a0234ef` |
| `CheckMultiBand.lean` | `7a22ab6f3ef9368465b29ff18dff0a1cecbb1cb1ec847dd85d44199fcfcd6d43` |
| `CheckPipeline.lean` | `94654ae6bb8560b7822ba71cc521b36eed1fafdd892ecbee2806e9db5491049b` |
| `MultiBand.lean` | `f84f4bbc42ee901f70bbfba708aca103c5e7f6b39369b9f31ab65871e74735da` |
| `MultiBand/Algebra.lean` | `7e42a9d1e3b20d3c56b0520e1f3315cbc48ed1c24b080836d50771993255fa67` |
| `MultiBand/Defs.lean` | `f720b66d2bf02ec8a07a69d093c2896e42a916c5199a410e2a1f291276970a63` |
| `MultiBand/Faithful.lean` | `5cee68435b8479d5a608c4da7f22db4fd17b186c8b80fc186479a2cbfd81f9c2` |
| `MultiBand/Identities.lean` | `88b94d7ce013c26530b062816584886cbfb771f661e0503d4c60653cf7cd1831` |
| `MultiBand/Main.lean` | `65a5ed50f637d1f6c7f2007101af1d0a515620ff0703be50331487af3a948096` |
| `MultiBand/Odd.lean` | `f359b9f849c9c113b020165f33f66ca14478e24bcb0e87c40f1ff75a5d67da04` |
| `MultiBand/Pipeline.lean` | `6dfe208383b17387ba124241a5245a9256ce6f5a0cb263ec54712ef324bfa29a` |
| `PipelineLogic.lean` | `0d2a49e789f459da6ed4225a167cf9f410be9a2ea0dd1636a14018905a6898e9` |
| `PipelineLogic/Assembly.lean` | `97e165bce083ce23fd29bd2d8bdeb0e6a8e64df57db57392a5c4a251cb0a30c7` |
| `PipelineLogic/Classify.lean` | `ce796fc8019f15ae63ab6a64a0b47304189dc038e1075ed44b1f9885ececac4f` |
| `PipelineLogic/Common.lean` | `990101b2419cf7d9590014fc593ccca32efe4f35b794d03f9b450a6b05abfc83` |
| `PipelineLogic/Cover.lean` | `6085e86ca7f85b43f8e6a616b21a65ea1584edb08dfb1462f163ae8e6b9f6a9d` |
| `PipelineLogic/Monotone.lean` | `9c58cbd1a475703692df4afef13c1aef7a34049c19dbab1e26793be6800f61c1` |
| `PipelineLogic/Multiband.lean` | `96ad41439272758d096b385cefb57beb61f4656786acffb852c20642bb30a612` |
| `PipelineLogic/Rayleigh.lean` | `dcde7393d54465d135c68bbdf6396c2af38181f627cdc2b3ab431c02ee59b0ef` |
| `PipelineLogic/Schur.lean` | `eb7a33ee7b99ba8a1d6d8aa3f0c9b6279717ab09fb53241fb5ec053a2b125c23` |
| `TcR1.lean` | `a14a871ba6f649cf1dc2225ae64c49029738229cce1e0206b3207bee25269c84` |
| `TcR1/Basic.lean` | `af693f241843be61f391c9f898333913a4fbd9a98bb147194b4b653b0b4d4e93` |
| `lake-manifest.json` | `188dd64d3fab820639a6604230313aa6f969c27337387c77a2957d4b4f625020` |
| `lean-toolchain` | `ca5f1f601d34c49079e38d7e7041eb0f2a9c227fbc1df0f404820069a052330c` |

Re-check: `cd lean && sha256sum -c` against the table.
