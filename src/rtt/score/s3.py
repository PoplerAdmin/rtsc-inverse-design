"""S3 ranking score (doc 11 T1a; architecture §3.1 S3, §3.1 notes). Ranking only (contract R2):
the score orders candidates for the S3 rank cut (``DEFERRED(rank_cut)``); it never removes.

Components (each clipped to [0, 1]):
  hull   = 1 - e_hull_ml_self / window, window = tiers [S3] ehull_self_max_eV (H7 + ML margin)
  H      = H_f / 0.75                                         (H fraction, doc 02 S2)
  mass   = sqrt(m_H / M_eff) / sqrt(1/2)                      (doc 11 T1a; light-atom weight)
  omega  = omega2_light / omega2_needed,  omega2_needed = 2 * floor / lambda_aim_hi   (λ_needed
           rule of §3.1 notes read backwards: the ω₂ that a λ at the top of the aim box needs)
  light  = light-element (H, Li, Be, B, C, N) atomic fraction
base  = 0.35 hull + 0.25 H + 0.20 omega + 0.10 mass + 0.10 light
score = (1 - w_m) base + w_m metal - 0.2 [H2 molecular motif]
  metal = composition-based metallicity prior (``rtt.score.metallicity``; doc 31 §2), w_m =
          tiers [S3] metallicity_weight (default 0.20). Campaign 1: without it the rank cut
          favoured charge-balanced insulating binaries (doc 30 §3, §7). Ranking only.
Missing components count 0.5 (neutral). ALIGNN numbers have weight 0 for H/light-element
candidates (doc 11 §2.3) and are not used at all here.
"""
from __future__ import annotations

import math
from typing import Any, Mapping, Optional

WEIGHTS = {"hull": 0.35, "H": 0.25, "omega": 0.20, "mass": 0.10, "light": 0.10}
H2_PENALTY = 0.2
METAL_WEIGHT = 0.20       # default; tiers.toml [S3] metallicity_weight is authoritative (doc 31 §2)
SCORE_VERSION = "s3-score/2"   # /2 adds the metallicity prior (doc 31 §2)
EHULL_WINDOW = 0.60       # default; matches tiers.toml [S3] ehull_self_max_eV (doc 21 §6, F15)
LIGHT = frozenset({"H", "Li", "Be", "B", "C", "N"})
FLOOR_300 = 433148145 / 1e8 * 300.0


def _clip(x: Optional[float]) -> float:
    if x is None or (isinstance(x, float) and math.isnan(x)):
        return 0.5
    return min(1.0, max(0.0, float(x)))


def components(features: Mapping[str, Any], phonons: Mapping[str, Any] | None = None, *,
               e_hull_ml_self: float | None = None, T_star_K: float = 300.0,
               lambda_aim_hi: float = 3.0) -> dict[str, float]:
    f = dict(features or {})
    ph = dict(phonons or {})
    hull = None if e_hull_ml_self is None else 1.0 - float(e_hull_ml_self) / EHULL_WINDOW
    hf = f.get("H_f")
    H = None if hf is None else float(hf) / 0.75
    mass = None
    if f.get("M_eff_amu"):
        mass = math.sqrt(1.008 / float(f["M_eff_amu"])) / math.sqrt(0.5)
    o2 = ph.get("omega2_light_K")
    need = 2.0 * FLOOR_300 * (T_star_K / 300.0) / lambda_aim_hi
    omega = None if o2 is None else float(o2) / need
    light = f.get("light_fraction")
    from .metallicity import prior_from_metrics
    return {"hull": _clip(hull), "H": _clip(H), "omega": _clip(omega), "mass": _clip(mass), "light": _clip(light),
            "metal": _clip(prior_from_metrics(f)["p_metal"]),
            "h2_molecular": 1.0 if f.get("motif") == "H2_molecular" else 0.0}


def metal_weight(tiers: Any = None) -> float:
    """tiers [S3] metallicity_weight (sourced threshold), else :data:`METAL_WEIGHT`."""
    if tiers is not None:
        try:
            return float(tiers.threshold("S3.metallicity_weight").value)
        except Exception:  # noqa: BLE001 - older tiers.toml without the key
            pass
    return METAL_WEIGHT


def score_s3(features_or_metrics: Mapping[str, Any], phonons: Any = None, targets: Any = None, *,
             tiers: Any = None) -> float:
    """Ranking score in [0, 1]. Accepts either (features, phonons, targets) as in architecture
    §5.1 or (s3_metrics, targets) as called by ``rtt.mlff.screen`` (metrics with keys
    ``features``, ``phonons``, ``e_hull_ml_self_eV``)."""
    m = dict(features_or_metrics or {})
    if "features" in m or "e_hull_ml_self_eV" in m:          # metrics form: (metrics, targets)
        if targets is None and phonons is not None and not isinstance(phonons, Mapping):
            targets, phonons = phonons, None
        feats, ph, eh = m.get("features") or {}, m.get("phonons") or {}, m.get("e_hull_ml_self_eV")
    else:
        feats, ph, eh = m, phonons if isinstance(phonons, Mapping) else (vars(phonons) if phonons else {}), None
    T = 300.0
    lam_hi = 3.0
    if targets is not None:
        try:
            T = float(targets.T_star)
            lam_hi = float(targets.get("aim.lambda")[1])
        except Exception:
            pass
    c = components(feats, ph, e_hull_ml_self=eh, T_star_K=T, lambda_aim_hi=lam_hi)
    wm = metal_weight(tiers)
    base = sum(WEIGHTS[k] * c[k] for k in WEIGHTS)
    s = (1.0 - wm) * base + wm * c["metal"] - H2_PENALTY * c["h2_molecular"]
    return round(max(0.0, min(1.0, s)), 6)
