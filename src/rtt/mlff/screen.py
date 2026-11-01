"""S3 (mlff, GPU) runner and decision (architecture §3.1 S3; §4 R2; doc 10 §3.3-3.4).

Order, as §3.1 fixes it: (a) batched MatterSim-5M relax (at the candidate's pressure: 0 GPa, or
the control override); (b) MatterSim-only hull (+ DFT-reference hull as a recorded second
opinion, 0 GPa only); (c) cheap features (+ ALIGNN numbers, weight 0 for H/light); (d) MatterSim
phonons for (b)-survivors -- ranking only.

Allowed outcomes: ``PASS | SCREENED_OUT | DEFERRED``. The only removal is

    SCREENED_OUT(H7-ml)  if  e_hull_ml_self > tiers[S3].ehull_self_max_eV (0.42 = 0.35 + 0.07)

MatterSim phonons, ALIGNN, the DFT-reference hull and every feature are **never** a removal
criterion (R2; doc 10 §3.4): :func:`decide_s3` does not read them. A missing hull (incomplete
references, or P above ``ml_hull_max_pressure_GPa``) cannot remove. An E_hull above the window
from an *unconverged* relax gives ``DEFERRED(relax_not_converged)``, not a removal. A relax failure gives ``DEFERRED(ml_failed)``. Survivors beyond the
rank cut get ``DEFERRED(rank_cut(rank, cut))``; the cut is applied once per pool (day or
campaign), with per-lane minimum shares, and recorded in ``rank_pools`` (doc 31 §1); controls
(lane C) are exempt from the
rank cut. A screen disabled through ``screen_status`` (R8) is recorded as PASS with tag
``SCREEN_DISABLED:S3:H7-ml``.
"""
from __future__ import annotations

import json
import math
import time
from pathlib import Path
from typing import Any, Iterable, Sequence

from ..records import DEFERRED, PASS, SCREENED_OUT, Outcome
from ..source.ingest import write_outcome

ORIGIN_S3 = "rtt.mlff.screen"
SCREEN_H7 = "S3:H7-ml"
DEFAULT_EHULL_MAX = (0.56, "RC1#F15 interim floor (0.35 + 2*sqrt(2)*0.074); tiers.toml [S3] is authoritative")
DEFAULT_RANK_CUT = (25, "doc12#8.3 (S4 capacity)")


def _th(tiers: Any, key: str, default: tuple[Any, str]) -> tuple[Any, str, bool]:
    if tiers is not None:
        try:
            t = tiers.threshold(key)
            return t.value, t.source_id, t.provisional
        except Exception:  # noqa: BLE001
            pass
    return default[0], default[1], False


# ----------------------------------------------------------------------------- score (ranking)
def rank_metrics(metrics: dict) -> dict:
    """Metrics as seen by the rank-cut score. For H/Li/Be compounds (or unknown composition) the
    MatterSim phonon block is removed, so its ω₂ term counts neutral (0.5): MatterSim phonons are
    unreliable exactly there (doc 10 §3.4) and, through the rank cut, would otherwise decide
    which candidates wait (RC1 F12)."""
    feats = metrics.get("features") or {}
    if feats.get("contains_HLiBe") is False:
        return metrics
    return {k: v for k, v in metrics.items() if k != "phonons"}


def score_s3(metrics: dict, targets: Any = None, tiers: Any = None) -> float:
    """Ranking score in [0, 1] (higher = earlier to S4). Uses ``rtt.score.s3.score_s3`` when that
    module exists (with the doc 31 §2 metallicity prior, weight from ``tiers``); otherwise this
    provisional default (doc 11 T1a-like, documented in doc 21):

        0.35·hull + 0.30·H-term + 0.25·ω₂-term + 0.10·light − 0.2·[H2 molecular]

    hull = clip(1 − e_hull_self/0.42); H-term = clip(H_f/0.75); ω₂-term = clip(ω₂_light/1500 K)
    (0.5 when phonons were not run); light = light-element fraction. ALIGNN has weight 0 here."""
    metrics = rank_metrics(metrics)
    try:
        from ..score.s3 import score_s3 as ext   # type: ignore
        return float(ext(metrics, targets, tiers=tiers))
    except ImportError:
        pass
    eh = metrics.get("e_hull_ml_self_eV")
    hull = 0.5 if eh is None else min(1.0, max(0.0, 1.0 - eh / DEFAULT_EHULL_MAX[0]))
    feats = metrics.get("features") or {}
    hterm = min(1.0, (feats.get("H_f") or 0.0) / 0.75)
    ph = metrics.get("phonons") or {}
    o2 = ph.get("omega2_light_K")
    o2t = 0.5 if o2 is None else min(1.0, max(0.0, o2 / 1500.0))
    light = feats.get("light_fraction") or 0.0
    pen = 0.2 if feats.get("motif") == "H2_molecular" else 0.0
    return round(max(0.0, 0.35 * hull + 0.30 * hterm + 0.25 * o2t + 0.10 * light - pen), 5)


# ----------------------------------------------------------------------------- decision (pure)
def window_calibrated_for(metrics: dict) -> bool:
    """True if the H7-ml window's calibration domain (H-containing compounds, doc 21) covers
    this candidate. Unknown composition -> False (cannot remove on an unverified domain)."""
    f = metrics.get("features") or {}
    if f.get("n_H") is not None:
        return int(f["n_H"]) > 0
    rs = metrics.get("relaxed_structure") or {}
    sp = rs.get("species") if isinstance(rs, dict) else None
    return bool(sp) and "H" in sp



def decide_s3(metrics: dict, tiers: Any = None, *, screen_enabled: bool = True) -> Outcome:
    """Per-candidate S3 outcome before the rank cut. Reads only ``e_hull_ml_self_eV`` and
    ``error``; never phonons (R2)."""
    if metrics.get("error"):
        return Outcome(stage="S3", outcome=DEFERRED, criterion="ml_failed",
                       reason=f"ml_failed: {str(metrics['error'])[:300]}")
    eh = metrics.get("e_hull_ml_self_eV")
    thr, src, prov = _th(tiers, "S3.ehull_self_max_eV", DEFAULT_EHULL_MAX)
    converged = (metrics.get("relax") or {}).get("converged", True)
    if eh is not None and eh > float(thr) and not converged:
        # an unconverged relax overestimates E_hull (and was not reproducible in the smoke run,
        # doc 21 §4): never remove on it
        return Outcome(stage="S3", outcome=DEFERRED, criterion="relax_not_converged", value=eh,
                       threshold=float(thr), source_id=src,
                       reason="relax_not_converged: E_hull above window but relax unconverged")
    if eh is not None and eh > float(thr) and not window_calibrated_for(metrics):
        # R30 F1(a), doc 31 §7: the 0.60 window is calibrated on H-containing entries only
        # (doc 21: 160 H-entries + H ZPE allowance); outside that domain an exceedance parks
        return Outcome(stage="S3", outcome=DEFERRED, criterion="uncalibrated_window", value=eh,
                       threshold=float(thr), source_id=f"{src}; doc31#7 (calibration domain: hydrides)",
                       provisional=True, tags=["H7_ML_WINDOW_UNCALIBRATED"],
                       reason=("uncalibrated_window: E_hull,ML above the H7-ml window, but the window is "
                               "calibrated on hydrides only; not a removal for this composition"))
    if eh is not None and eh > float(thr):
        if not screen_enabled:
            return Outcome(stage="S3", outcome=PASS, tags=[f"SCREEN_DISABLED:{SCREEN_H7}"],
                           value={"e_hull_ml_self_eV": eh})
        return Outcome(stage="S3", outcome=SCREENED_OUT, criterion="H7-ml:e_hull_ml_self_eV", value=eh,
                       threshold=float(thr), source_id=src, provisional=prov)
    tags = [] if eh is not None else ["NO_ML_HULL"]
    if not converged:
        tags.append("RELAX_NOT_CONVERGED")
    return Outcome(stage="S3", outcome=PASS, tags=tags)


def _lane_in(lane: str, key: str) -> bool:
    """Lane-share key match: an exact lane ('G1', 'D') or a one-letter family ('G' = G1..G5, GA)."""
    return lane == key or (len(key) == 1 and str(lane).startswith(key))


def lane_min_share(tiers: Any) -> tuple[dict[str, float], str | None]:
    """tiers [S3.rank_lane_min_share] -> ({lane key: share}, source id). Doc 31 §2: replaces the
    campaign-1 hand quota CAMPAIGN_GENERATED_QUOTA."""
    if tiers is None:
        return {}, None
    try:
        tab = tiers.get("S3.rank_lane_min_share", None)
    except Exception:  # noqa: BLE001
        tab = None
    if not isinstance(tab, dict):
        return {}, None
    shares = {str(k): float(v) for k, v in tab.items()
              if k not in ("source", "source_id", "provisional") and isinstance(v, (int, float))
              and not isinstance(v, bool) and float(v) > 0}
    return shares, tab.get("source") or tab.get("source_id")


def rank_pool(items: Sequence[tuple[str, str, Outcome, float]], cut: int, cut_src: str, *,
              min_share: dict[str, float] | None = None, share_src: str | None = None,
              pool_id: str | None = None) -> tuple[dict[str, Outcome], list[dict]]:
    """One deterministic rank cut over a whole pool (B9, doc 31 §1).

    items: (cid, lane, pre-cut outcome, score). Non-control PASS items are ranked by
    (score desc, cid asc). Per-lane minimum shares are reserved first: for each key (sorted),
    the best ``floor(share * cut)`` members of that lane family take a slot; the remaining slots
    go by global rank. Everything else becomes ``DEFERRED(rank_cut(rank, cut))``. A slot
    obtained only through a share (global rank > cut) is tagged ``LANE_MIN_SHARE:<key>``.
    Controls (lane C) are exempt. Returns (outcomes, ranking rows for the ledger)."""
    out = {cid: o for cid, _, o, _ in items}
    lanes = {cid: lane for cid, lane, _, _ in items}
    scores = {cid: float(s) for cid, _, _, s in items}
    cut = max(0, int(cut))
    ranked = [cid for _s, cid in sorted([(s, cid) for cid, lane, o, s in items if o.outcome == PASS and lane != "C"],
                                        key=lambda t: (-t[0], t[1]))]
    grank = {cid: i for i, cid in enumerate(ranked, start=1)}
    reserved: dict[str, str] = {}
    for key in sorted(min_share or {}):
        k = min(int(math.floor(float(min_share[key]) * cut + 1e-9)), cut - len(reserved))
        for cid in [c for c in ranked if _lane_in(lanes[c], key) and c not in reserved][:max(0, k)]:
            reserved[cid] = key
    rest = [c for c in ranked if c not in reserved][:max(0, cut - len(reserved))]
    selected = set(reserved) | set(rest)
    ptag = [f"RANK_POOL:{pool_id}"] if pool_id else []
    rows = []
    for cid in ranked:
        r = grank[cid]
        via = None
        if cid in selected:
            o = out[cid]
            tags = [*o.tags, *ptag]
            via = "rank"
            if cid in reserved and r > cut:
                via = f"lane_min_share:{reserved[cid]}"
                tags.append(f"LANE_MIN_SHARE:{reserved[cid]}")
            out[cid] = Outcome(stage="S3", outcome=PASS, tags=tags, value=o.value, tier_run_id=o.tier_run_id,
                               source_id=(share_src if via != "rank" else o.source_id))
        else:
            out[cid] = Outcome(stage="S3", outcome=DEFERRED, criterion="rank_cut", value=r, threshold=cut,
                               source_id=cut_src, reason=f"rank_cut({r}, {cut})", tags=list(ptag))
        rows.append({"cid": cid, "lane": lanes[cid], "score": scores[cid], "rank": r,
                     "selected": cid in selected, "via": via})
    for cid, lane, o, _ in items:
        if lane == "C" and o.outcome == PASS:
            out[cid] = Outcome(stage="S3", outcome=PASS, tags=[*o.tags, "CONTROL_EXEMPT_RANK_CUT"], value=o.value)
    return out, rows


def apply_rank_cut(items: Sequence[tuple[str, str, Outcome, float]], cut: int, cut_src: str, **kw: Any
                   ) -> dict[str, Outcome]:
    """items: (cid, lane, pre-cut outcome, score). PASS items of non-control lanes beyond rank
    ``cut`` (by score, descending) become DEFERRED(rank_cut). See :func:`rank_pool`."""
    return rank_pool(items, cut, cut_src, **kw)[0]


# ----------------------------------------------------------------------------- runner
def _versions() -> dict:
    v = {}
    for mod in ("mattersim", "phonopy", "pymatgen", "ase", "torch"):
        try:
            v[mod] = __import__(mod).__version__
        except Exception:  # noqa: BLE001
            pass
    return v


def compute_s3(structs: Sequence[dict], pressures: Sequence[float], *, model: str = "MatterSim-v1.0.0-5M",
               device: str | None = None, phonons: bool = True, dftref: bool = True, alignn: bool = True,
               tiers: Any = None, machine: Any = None, relax_fn: Any = None, hull_fn: Any = None,
               phonon_fn: Any = None, dftref_fn: Any = None, max_hull_pressure_GPa: float | None = None,
               max_steps: int | None = None) -> list[dict]:
    """Pure compute for a batch (no ledger): one metrics dict per structure.

    The MatterSim-only hull is built only for P <= ``max_hull_pressure_GPa`` (default
    ``tiers [S3] ml_hull_max_pressure_GPa`` or 1.0: the project scope is 0-1 GPa and H7 is an
    ambient window). Above it (the H3S calibration entry) E_hull is None (cannot remove): in the
    smoke run the 200 GPa H-S reference hull cost 493 s and hit CUDA OOM three times.

    The ``*_fn`` hooks exist for tests (CPU, no MatterSim); defaults are the real implementations."""
    from .features import features
    from .hull import build_self_hull, ehull_dftref, ehull_self
    from .phonons import ml_phonons
    from .relax import relax_batch
    relax_fn = relax_fn or relax_batch
    phonon_fn = phonon_fn or ml_phonons
    dftref_fn = dftref_fn or ehull_dftref
    thr, _, _ = _th(tiers, "S3.ehull_self_max_eV", DEFAULT_EHULL_MAX)
    if max_hull_pressure_GPa is None:
        try:
            max_hull_pressure_GPa = float(tiers.get("S3.ml_hull_max_pressure_GPa", 1.0)) if tiers is not None else 1.0
        except Exception:  # noqa: BLE001
            max_hull_pressure_GPa = 1.0

    t0 = time.time()
    kw = {} if max_steps is None else {"max_steps": int(max_steps)}
    rel = relax_fn(list(structs), model=model, pressure_GPa=list(pressures), device=device, **kw)
    t_relax = time.time() - t0
    mets: list[dict] = []
    hulls: dict = {}
    for r in rel:
        m: dict[str, Any] = {"model": model, "pressure_GPa": r.pressure_GPa,
                             "relax": {k: getattr(r, k) for k in ("energy_eV", "energy_per_atom_eV",
                                       "enthalpy_per_atom_eV", "volume_per_atom_A3", "max_force_eV_A",
                                       "internal_pressure_GPa", "n_steps", "converged")},
                             "relaxed_structure": r.structure}
        if r.error or not math.isfinite(r.enthalpy_per_atom_eV):
            m["error"] = r.error or "non-finite energy"
        mets.append(m)
    # (b) hulls
    for r, m in zip(rel, mets):
        if m.get("error"):
            continue
        cs = "-".join(sorted(set(r.structure["species"])))
        key = (cs, r.pressure_GPa)
        if r.pressure_GPa > max_hull_pressure_GPa:
            m["e_hull_ml_self_eV"] = None
            m["hull_self"] = {"note": f"no ML hull above {max_hull_pressure_GPa} GPa (out of scope)"}
            continue
        try:
            if key not in hulls:
                hulls[key] = (hull_fn or build_self_hull)(cs, pressure_GPa=r.pressure_GPa, model=model,
                                                          device=device, machine=machine, relax_fn=relax_fn)
            eh, info = ehull_self(r.structure, r.enthalpy_per_atom_eV, cs, hull=hulls[key])
        except Exception as ex:  # noqa: BLE001 - a missing hull cannot remove
            eh, info = None, {"error": f"{type(ex).__name__}: {ex}"[:300]}
        m["e_hull_ml_self_eV"] = eh
        m["hull_self"] = info
    if dftref:
        items = [(str(i), r.structure, r.energy_eV) for i, (r, m) in enumerate(zip(rel, mets))
                 if not m.get("error") and r.pressure_GPa == 0.0]
        try:
            res = dftref_fn(items, machine=machine) if items else {}
        except Exception as ex:  # noqa: BLE001
            res = {k: {"e_above_hull_eV": None, "note": f"{type(ex).__name__}: {ex}"[:300]} for k, _, _ in items}
        for i, m in enumerate(mets):
            row = res.get(str(i))
            m["e_hull_ml_dftref_eV"] = None if row is None else row.get("e_above_hull_eV")
            if row is not None and row.get("note"):
                m["dftref_note"] = row["note"]
    # (c) features (+ ALIGNN, recorded only)
    for r, m in zip(rel, mets):
        if m.get("error"):
            continue
        m["features"] = features(r.structure)
        if alignn:
            from .alignn import alignn_numbers
            m["alignn"] = alignn_numbers(r.structure)
    # (d) phonons for (b)-survivors, ranking only
    for r, m in zip(rel, mets):
        if not phonons or m.get("error"):
            continue
        eh = m.get("e_hull_ml_self_eV")
        if eh is not None and eh > float(thr):
            continue
        try:
            ps = phonon_fn(r.structure, model=model, device=device)
            m["phonons"] = ps.to_json() if hasattr(ps, "to_json") else dict(ps)
        except Exception as ex:  # noqa: BLE001
            m["phonons"] = {"error": f"{type(ex).__name__}: {ex}"[:300]}
    for m in mets:
        m["timing"] = {"batch_relax_s": round(t_relax, 2), "batch_total_s": round(time.time() - t0, 2),
                       "n_batch": len(mets)}
    return mets


def _score_version() -> str | None:
    try:
        from ..score.s3 import SCORE_VERSION
        return SCORE_VERSION
    except ImportError:
        return None


def _stored_score(ledger: Any, cid: str, targets: Any = None, tiers: Any = None) -> float:
    """Stored S3 score; recomputed from the stored metrics when they were scored by an older
    score version (so one pool is ranked by one score function, doc 31 §1)."""
    o = ledger.current_outcome(cid, "S3")
    tr = ledger.tier_run(o["tier_run_id"]) if o and o.get("tier_run_id") else None
    m = json.loads((tr or {}).get("metrics_json") or "{}") or {}
    ver = _score_version()
    if ver and m.get("score_s3_version") != ver and m.get("features") and not m.get("error"):
        return float(score_s3(m, targets, tiers))
    return float(m.get("score_s3") or 0.0)


def _is_rank_cut(o: dict | None) -> bool:
    return bool(o) and o["outcome"] == DEFERRED and o.get("criterion") == "rank_cut"


def _waiting_passes(ledger: Any, exclude: set) -> int:
    """Non-control S3 PASSes not yet picked up by S4: they already hold rank-cut slots."""
    n = 0
    for cid, lane in ledger.qa("SELECT id, lane FROM candidates"):
        if cid in exclude or lane == "C":
            continue
        o = ledger.current_outcome(cid, "S3")
        if o and o["outcome"] == PASS and ledger.current_outcome(cid, "S4") is None:
            n += 1
    return n


def _default_pool_id(ledger: Any) -> str:
    import datetime as _dt
    day = _dt.datetime.now(_dt.timezone.utc).strftime("%Y-%m-%d")
    try:
        n = int(ledger.q1("SELECT COUNT(*) FROM rank_pools WHERE stage='S3'")[0]) + 1
    except Exception:  # noqa: BLE001
        n = 1
    return f"S3:{day}#{n}"


def rerank_s3(ledger: Any, *, targets: Any, tiers: Any, rank_cut: int | None = None,
              new_items: Sequence[tuple[str, str, Outcome, float]] = (), pool_id: str | None = None
              ) -> dict[str, Outcome]:
    """Close the S3 rank pool: ONE rank cut over every queued candidate (B9, doc 31 §1; RC1 F12).

    The pool is the new batch's PASSes plus every candidate currently ``DEFERRED(rank_cut)``
    (including ``rank_pending`` entries queued by ``run_s3(close_pool=False)``), ranked by
    :func:`rank_pool` with the tiers per-lane minimum shares. Pool scores are the stored S3
    scores (recomputed from stored metrics only when an older score version made them). Free
    slots = cut − non-control PASSes still waiting for S4. The ranking (cid, lane, score, rank,
    selected, via) is recorded in the ledger table ``rank_pools``; every written outcome carries
    the tag ``RANK_POOL:<pool_id>``. Writes outcomes only for queued candidates that are promoted
    (PASS, tag ``PROMOTED_FROM_RANK_CUT``); returns the final outcome of every new item plus the
    promotions (new items are written by the caller)."""
    cut, cut_src, _ = _th(tiers, "S3.rank_cut_per_day", DEFAULT_RANK_CUT)
    cut = int(rank_cut if rank_cut is not None else cut)
    new_ids = {cid for cid, *_ in new_items}
    pool = []
    for cid, lane in sorted(ledger.qa("SELECT id, lane FROM candidates")):
        if cid not in new_ids and _is_rank_cut(ledger.current_outcome(cid, "S3")):
            pool.append((cid, lane, Outcome(stage="S3", outcome=PASS), _stored_score(ledger, cid, targets, tiers)))
    waiting = _waiting_passes(ledger, new_ids)
    free = max(0, cut - waiting)
    shares, share_src = lane_min_share(tiers)
    pool_id = pool_id or _default_pool_id(ledger)
    final, rows = rank_pool(list(new_items) + pool, free, cut_src, min_share=shares, share_src=share_src,
                            pool_id=pool_id)
    if rows or new_items:
        _record_pool(ledger, pool_id, cut=cut, free=free, waiting=waiting, shares=shares, share_src=share_src,
                     cut_src=cut_src, rows=rows, tiers=tiers)
    out = {cid: final[cid] for cid in new_ids}
    for cid, _lane, _o, _s in pool:
        if final[cid].outcome == PASS:
            rec = ledger.candidate_record(cid)
            o = Outcome(stage="S3", outcome=PASS, tags=["PROMOTED_FROM_RANK_CUT", *final[cid].tags],
                        source_id=final[cid].source_id,
                        tier_run_id=ledger.current_outcome(cid, "S3").get("tier_run_id"))
            write_outcome(ledger, cid, o, origin_obj=rec.origin, targets=targets, tiers=tiers, origin=ORIGIN_S3)
            out[cid] = o
    return out


#: public name: close the day's (or campaign's) pool once, after all S3 chunks
close_rank_pool = rerank_s3


def _record_pool(ledger: Any, pool_id: str, *, cut: int, free: int, waiting: int, shares: dict, share_src: Any,
                 cut_src: str, rows: list[dict], tiers: Any) -> None:
    import hashlib
    ranking = json.dumps(rows, sort_keys=True)
    rec = {"pool_id": pool_id, "stage": "S3", "cut": cut, "free_slots": free, "waiting_s4": waiting,
           "cut_source": cut_src, "lane_min_share": shares, "lane_min_share_source": share_src,
           "score_version": _score_version(), "ranking": rows,
           "ranking_sha256": hashlib.sha256(ranking.encode()).hexdigest(),
           "tiers_sha256": getattr(tiers, "sha256", None)}
    ins = getattr(ledger, "insert_rank_pool", None)
    if ins is not None:
        ins(rec)


def run_s3(ledger: Any, cids: Iterable[str] | None = None, *, targets: Any, tiers: Any, machine: Any = None,
           device: str | None = None, phonons: bool = True, dftref: bool = True, alignn: bool = True,
           rank_cut: int | None = None, workdir: str | Path | None = None, close_pool: bool = True,
           pool_id: str | None = None, chunk_size: int | None = None, **hooks: Any) -> dict[str, Outcome]:
    """Run S3 for every S2-PASS candidate that has no S3 outcome yet **or** whose current S3
    outcome is retryable (``DEFERRED(ml_failed | relax_not_converged)``) within the
    ``[retry] max_attempts`` budget of :mod:`rtt.engine.retry` (RC1 F12). Candidates in the rank
    queue (``DEFERRED(rank_cut)``) are not recomputed: :func:`rerank_s3` re-ranks them from their
    stored score together with the new batch. A retry doubles the relax ``max_steps`` per failed
    attempt; ``Outcome.retry`` carries the attempt count.

    Writes one ``tier_runs`` row per computed candidate (kind 'ml', state 'started' before
    compute), the metrics as an artifact, and one S3 outcome. Returns cid -> Outcome (including
    rank-queue promotions).

    Rank cut (B9, doc 31 §1): the cut is applied ONCE over the whole pool, never per chunk.
    ``chunk_size`` computes in chunks (memory); each chunk's survivors are written as
    ``DEFERRED(rank_cut)`` "rank_pending" queue entries and the pool is closed once at the end
    (:func:`close_rank_pool`). ``close_pool=False`` only queues: call :func:`close_rank_pool`
    once after the last call of the day/campaign. Each close is recorded in ``rank_pools``."""
    from ..engine.retry import failed_attempts, needs_run
    model = (tiers.get("S3.ml_model", "MatterSim-v1.0.0-5M") if tiers is not None else "MatterSim-v1.0.0-5M")
    todo = list(cids) if cids is not None else ledger.candidate_ids()
    recs, attempts = [], {}
    for cid in todo:
        s2 = ledger.current_outcome(cid, "S2")
        if not (s2 and s2["outcome"] == PASS):
            continue
        cur = ledger.current_outcome(cid, "S3")
        if _is_rank_cut(cur):
            continue                                   # rank queue: re-ranked, not recomputed
        if needs_run(ledger, cid, "S3", tiers):
            recs.append(ledger.candidate_record(cid))
            attempts[cid] = failed_attempts(ledger, cid, "S3") if cur is not None else 0
    if not recs:
        return rerank_s3(ledger, targets=targets, tiers=tiers, rank_cut=rank_cut, pool_id=pool_id) if close_pool else {}
    n_fail = max(attempts.values())
    max_steps = 1000 * 2 ** n_fail if n_fail else None
    params = {"model": model, "fmax": 0.02, "optimizer": "FIRE+FrechetCellFilter", "phonons": phonons,
              "dftref": dftref, "device": device, "max_steps": max_steps or 1000}
    vers = _versions()
    wd = Path(workdir or ".") / "s3_artifacts"
    wd.mkdir(parents=True, exist_ok=True)
    cs = int(chunk_size) if chunk_size else len(recs)
    chunks = [recs[i:i + cs] for i in range(0, len(recs), cs)]
    queued_mode = len(chunks) > 1 or not close_pool
    final: dict[str, Outcome] = {}
    trs: dict[str, int] = {}
    items: list = []
    for chunk in chunks:
        citems = _s3_chunk(ledger, chunk, trs, params, vers, wd, model=model, targets=targets, tiers=tiers,
                           machine=machine, device=device, phonons=phonons, dftref=dftref, alignn=alignn,
                           max_steps=max_steps, hooks=hooks)
        if queued_mode:          # B9: survivors enter the pool queue; the cut is applied once, below
            for cid, lane, o, _s in citems:
                q = _pending(o, lane, pool_id)
                q.tier_run_id, q.retry = trs[cid], attempts.get(cid, 0)
                rec = next(r for r in chunk if r.candidate_id == cid)
                write_outcome(ledger, cid, q, origin_obj=rec.origin, targets=targets, tiers=tiers, origin=ORIGIN_S3)
                final[cid] = q
        else:
            items += citems
    if queued_mode:
        if close_pool:
            final.update(rerank_s3(ledger, targets=targets, tiers=tiers, rank_cut=rank_cut, pool_id=pool_id))
        return final
    final = rerank_s3(ledger, targets=targets, tiers=tiers, rank_cut=rank_cut, new_items=items, pool_id=pool_id)
    for rec in recs:
        o = final[rec.candidate_id]
        o.tier_run_id = trs[rec.candidate_id]
        o.retry = attempts.get(rec.candidate_id, 0)
        write_outcome(ledger, rec.candidate_id, o, origin_obj=rec.origin, targets=targets, tiers=tiers,
                      origin=ORIGIN_S3)
    return final


def _pending(o: Outcome, lane: str, pool_id: str | None) -> Outcome:
    """Queue entry for a non-control S3 PASS awaiting the pool-wide cut (B9)."""
    if o.outcome != PASS:
        return o
    if lane == "C":
        return Outcome(stage="S3", outcome=PASS, tags=[*o.tags, "CONTROL_EXEMPT_RANK_CUT"], value=o.value)
    return Outcome(stage="S3", outcome=DEFERRED, criterion="rank_cut", tags=[*o.tags, "RANK_PENDING"],
                   reason=f"rank_pending(pool={pool_id or 'next close'}): the rank cut is applied once per pool")


def _s3_chunk(ledger: Any, recs: Sequence[Any], trs: dict, params: dict, vers: dict, wd: Path, *, model: str,
              targets: Any, tiers: Any, machine: Any, device: Any, phonons: bool, dftref: bool, alignn: bool,
              max_steps: Any, hooks: dict) -> list:
    """Compute one chunk: tier_runs, artifacts, per-candidate pre-cut outcome and score."""
    from ..ledger.artifacts import register
    for r in recs:
        trs[r.candidate_id] = ledger.start_tier_run(r.candidate_id, "S3", "ml", profile="mlff", params=params,
                                                    versions=vers, lane="G")
    t0 = time.time()
    try:
        mets = compute_s3([r.structure.to_json() for r in recs], [r.pressure_GPa for r in recs], model=model,
                          device=device, phonons=phonons, dftref=dftref, alignn=alignn, tiers=tiers,
                          machine=machine, max_steps=max_steps, **hooks)
    except Exception as ex:  # noqa: BLE001 - whole-batch failure: defer everyone, keep the reason
        mets = [{"error": f"batch: {type(ex).__name__}: {ex}"[:500]} for _ in recs]
    wall = time.time() - t0
    gpu_h_each = (wall / 3600 / len(recs)) if (device or "cuda") != "cpu" else 0.0
    _relaxed_links(ledger, recs, mets)
    items = []
    for rec, m in zip(recs, mets):
        m["score_s3"] = None if m.get("error") else score_s3(m, targets, tiers)
        m["score_s3_version"] = _score_version()
        m["score_s3_phonon_term"] = "used" if "phonons" in rank_metrics(m) else "excluded_or_absent"
        tr = trs[rec.candidate_id]
        art = wd / f"{rec.candidate_id}_s3_tr{tr}.json"
        art.write_text(json.dumps(m, sort_keys=True, default=str))
        register(ledger, tr, art)
        ledger.finish_tier_run(tr, "failed" if m.get("error") else "done", wall_s=wall / len(recs),
                               gpu_h=gpu_h_each, metrics=m)
        o = decide_s3(m, tiers, screen_enabled=ledger.screen_enabled(SCREEN_H7))
        o.tier_run_id = tr
        items.append((rec.candidate_id, rec.origin.lane, o, m["score_s3"] or 0.0))
    return items


def _relaxed_links(ledger: Any, recs: Sequence[Any], mets: Sequence[dict]) -> None:
    """Post-relaxation duplicate *link* (RC1 F16): record ``relaxed_match_of`` when a relaxed cell
    matches an earlier relaxed S3 cell at the same P. Informational only, never an outcome."""
    from ..structure.dedupe import DedupeIndex
    idx = DedupeIndex()
    try:
        for cid, mj in ledger.qa("SELECT candidate_id, metrics_json FROM tier_runs WHERE stage='S3' AND state='done'"):
            m = json.loads(mj or "{}") or {}
            if m.get("relaxed_structure"):
                idx.add(cid, m["relaxed_structure"], m.get("pressure_GPa", 0.0))
        for rec, m in zip(recs, mets):
            if m.get("error") or not m.get("relaxed_structure"):
                continue
            hit = idx.find(m["relaxed_structure"], m.get("pressure_GPa", 0.0))
            if hit is not None and hit != rec.candidate_id:
                m["relaxed_match_of"] = hit
            idx.add(rec.candidate_id, m["relaxed_structure"], m.get("pressure_GPa", 0.0))
    except Exception as ex:  # noqa: BLE001 - informational
        for m in mets:
            m.setdefault("relaxed_match_error", f"{type(ex).__name__}: {ex}"[:200])


__all__ = ["run_s3", "rerank_s3", "close_rank_pool", "rank_pool", "lane_min_share", "compute_s3", "decide_s3",
           "apply_rank_cut", "score_s3", "rank_metrics",
           "ORIGIN_S3", "SCREEN_H7"]
