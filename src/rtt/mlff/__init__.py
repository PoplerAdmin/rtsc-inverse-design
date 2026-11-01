"""S3 (mlff): MatterSim relax, ML hulls, features, phonons (ranking only), ALIGNN (recorded only).

Heavy imports (torch, mattersim, phonopy) happen inside functions, so importing this package is
cheap and CPU-safe."""
from .screen import ORIGIN_S3, apply_rank_cut, compute_s3, decide_s3, run_s3, score_s3  # noqa: F401
