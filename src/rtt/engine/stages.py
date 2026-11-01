"""S1..S9 wired to the sibling builders' functions (architecture §3.1, §5.1).

Batch stages (S1-S3) delegate to the sibling runners, which write their own tier_runs and
outcomes. S4-S8 are per-candidate here; their compute comes from ``rtt.dft`` (B3) and their
decisions from ``rtt.certify`` + ``rtt.decide`` (this builder). The soundness-critical glue
(core verdict -> recheck -> ledger; decide_t -> decision -> stage outcome; the §3.2
combination rule) is :func:`record_dfpt_decision`, which is pure ledger work and is tested
without QE in tests/test_engine.py.

Heavy imports (QE wrappers, MatterSim, MatterGen) are local to each stage so that importing this
module needs neither GPU nor QE.
"""
from __future__ import annotations

import json
import statistics
import traceback
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Optional, Sequence

from ..records import (D_PASS, D_REJECTED_PROVEN, D_SCREENED_OUT, D_UNDECIDED, ERROR, ORIGIN_CERTIFY,
                       ORIGIN_DECIDE, PROFILE_DECIDE_TIER, CandidateRecord, CostEstimate, DecideInputs, Outcome,
                       StageResult)
from .stage import Context, StageBase

ENGINE_ORIGIN = "rtt.engine.stages"


# ------------------------------------------------------------------------------------ helpers
def dft_structure(rec: CandidateRecord) -> dict:
    """DFT input cell: the candidate's cell, idealized to the space group it is labelled with
    (spglib at 0.01 Å) when it is only approximately symmetric (campaign 1: generated cells ran
    as P1 in QE). Exact cells are returned unchanged."""
    from ..structure.canonical import symmetrized_for_dft
    sd, _info = symmetrized_for_dft(rec.structure.to_json())
    return sd


def write_outcome(ctx: Context, cid: str, o: Outcome, decision_id: int | None = None, rec: CandidateRecord | None = None):
    tgt = ctx.targets_for(rec) if rec is not None else ctx.targets
    ctx.ledger.register_configs(tgt)
    return ctx.ledger.insert_outcome(cid, o, targets_sha256=tgt.sha256, tiers_sha256=ctx.tiers.sha256,
                                     origin=ENGINE_ORIGIN, decision_id=decision_id)


def combine(stage: str, decide_outcome: str, rule: Optional[str], core_status: Optional[str]) -> tuple[str, Optional[str]]:
    """Architecture §3.2 combination rule, mapped to the outcomes the stage may emit (§3.1).

    * rejectedProven by N1-ISO needs the core verdict REJECTED_PROVEN for the same tier_run,
      else ``undecided`` with reason CROSSCHECK_DISAGREE.
    * At S6-S8 a ``screenedOut`` (a floor holds on lower-tier inputs, Lean v2) is not a removal:
      stage outcome ``undecided`` with reason WOULD_FIRE_LOWER_TIER.
    * At S5 the decide outcome is evidence only; the stage outcome comes from the coarse margin
      rule, so this function is not used there."""
    if decide_outcome == D_REJECTED_PROVEN and rule == "N1-ISO" and core_status != "REJECTED_PROVEN":
        return D_UNDECIDED, "CROSSCHECK_DISAGREE"
    if decide_outcome == D_SCREENED_OUT and stage in ("S6", "S7", "S8"):
        return D_UNDECIDED, "WOULD_FIRE_LOWER_TIER"
    return decide_outcome, None


def verdict_spectra_rows(spectra: Sequence[Any]) -> list[dict]:
    """[{alt_index, spectrum_digest, source_artifact_id, sha256}] for the verdicts_bind trigger.
    Frequency-scaled copies point at their parent's source artifact."""
    from ..certify.recheck import digests
    rows = []
    for i, (sp, dg) in enumerate(zip(spectra, digests(spectra))):
        p = sp.provenance
        rows.append({"alt_index": i, "spectrum_digest": dg, "source_artifact_id": p.get("source_artifact_id"),
                     "sha256": p.get("source_sha256")})
    return rows


def store_verdict(ctx: Context, cid: str, tier_run_id: int, spectrum_tier: str, verdict: Any, spectra: Sequence[Any],
                  ccfg: dict, workdir: Path) -> int:
    """Recheck (R4) then insert the core verdict with its certificate artifact."""
    from ..certify.recheck import RecheckError, recheck_or_raise
    from ..ledger.artifacts import register
    ok = True
    try:
        recheck_or_raise(verdict, spectra, ccfg)
    except RecheckError:
        ok = False
    workdir.mkdir(parents=True, exist_ok=True)
    cert = workdir / f"verdict_{spectrum_tier}_{tier_run_id}_mu{ccfg['mu_hi']}.json"
    cert.write_text(verdict.to_json(sort_keys=True), encoding="utf-8")
    aref = register(ctx.ledger, tier_run_id, cert)
    return ctx.ledger.insert_verdict(
        cid, tier_run_id, spectrum_tier, status=verdict.status.value, T_low_K=verdict.T_low_K, T_up_K=verdict.T_up_K,
        T_star_K=float(ccfg["T_target_K"]), mu_hi=float(ccfg["mu_hi"]),
        rel_envelope=max(s.rel_envelope for s in spectra), rule_ids=verdict.rule_ids,
        assumptions=verdict.assumptions, flags=verdict.flags, certificate_artifact_id=aref.id, recheck_ok=ok,
        spectra=verdict_spectra_rows(spectra), origin=ORIGIN_CERTIFY)


def record_dfpt_decision(ctx: Context, rec: CandidateRecord, stage: str, tier_run_id: int, di: DecideInputs,
                         crosscheck_verdict_id: Optional[int], core_status: Optional[str],
                         emit_stage_outcome: bool = True) -> tuple[Any, Optional[Outcome]]:
    """decide_t -> decide_inputs + decisions rows -> stage outcome (combination rule)."""
    from ..decide.core import DecideConfig, decide_t
    from ..targets.spec import TargetSpec
    spec = TargetSpec.from_targets(ctx.targets_for(rec))
    cfg = DecideConfig.from_spec(spec)
    di.cfg = cfg.record()
    o = decide_t(cfg, di)
    L = ctx.ledger
    di_id = L.insert_decide_input(rec.candidate_id, di)
    did = L.insert_decision(rec.candidate_id, stage, di.tier, di_id, o.outcome, firing_rule=o.rule,
                            assumptions=o.assumptions,
                            crosscheck_verdict_id=crosscheck_verdict_id if o.rule == "N1-ISO" else None,
                            origin=ORIGIN_DECIDE)
    if not emit_stage_outcome:
        return o, None
    so, reason = combine(stage, o.outcome, o.rule, core_status)
    out = Outcome(stage=stage, outcome=so, tier_run_id=tier_run_id, decide=o.to_json(),
                  reason=reason or (";".join(o.reasons) or None), criterion=o.rule)
    write_outcome(ctx, rec.candidate_id, out, decision_id=did, rec=rec)
    return o, out


def ladder_from_runs(collected_by_mesh: Sequence[Any]) -> dict:
    """The S6 ``ladder`` metric from collected runs (coarser mesh first, finest last)."""
    def moments(sps):
        return {"A_K": [s.A() for s in sps], "lambda": [s.lambda_() for s in sps],
                "m_K2": [2.0 * sum(o * w for o, w in zip(s.omega_K, s.weight_K)) for s in sps]}
    fine = collected_by_mesh[-1].spectra
    coarse = collected_by_mesh[-2].spectra if len(collected_by_mesh) > 1 else fine
    idx_f = [s.provenance.get("smearing_index") for s in fine]
    by_idx_c = {s.provenance.get("smearing_index"): s for s in coarse}
    common = [i for i in idx_f if i in by_idx_c]
    fine_c = [s for s in fine if s.provenance.get("smearing_index") in common]
    coarse_c = [by_idx_c[i] for i in common]
    return {"q_meshes": [c.spectra[0].provenance.get("q") for c in collected_by_mesh if c.spectra],
            "fine": moments(fine_c), "coarse": moments(coarse_c),
            "omega_log_K_by_q": [statistics.median([s.omega_log() for s in c.spectra]) for c in collected_by_mesh]}


# ------------------------------------------------------------------------------------ S1-S3 (batch, siblings)
class S1Source(StageBase):
    id, kind, lane = "S1", "heuristic", "G"

    def run_batch(self, ctx: Context, *, lane_jobs: Sequence[Any] = (), database_systems: Sequence[str] = (),
                  controls: bool = False) -> list[str]:
        from ..source import controls as C
        from ..source.runner import run_database_lane, run_lane_job
        cids: list[str] = []
        for job in lane_jobs:
            cids += run_lane_job(ctx.ledger, job, targets=ctx.targets, tiers=ctx.tiers, machine=ctx.machine,
                                 workdir=ctx.workdir / "S1")
        if database_systems:
            cids += run_database_lane(ctx.ledger, list(database_systems), targets=ctx.targets, tiers=ctx.tiers,
                                      machine=ctx.machine)
        if controls:
            cids += C.inject_controls(ctx.ledger, targets=ctx.targets, tiers=ctx.tiers)
        return cids


class S2Structure(StageBase):
    id, kind, lane, requires = "S2", "heuristic", "G", ("S1",)

    def run_batch(self, ctx: Context, cids: Sequence[str] | None = None) -> dict:
        from ..structure.stage import run_s2
        return run_s2(ctx.ledger, cids, targets=ctx.targets, tiers=ctx.tiers, machine=ctx.machine)


class S3MLFF(StageBase):
    id, kind, lane, requires = "S3", "heuristic", "G", ("S2",)

    def run_batch(self, ctx: Context, cids: Sequence[str] | None = None, **kw) -> dict:
        from ..mlff.screen import run_s3
        from .retry import needs_run
        if cids is None:     # RC1 F12: retryable DEFERRED/ERROR re-enter, rank_cut is re-ranked
            cids = [c for c in ctx.ledger.candidate_ids() if needs_run(ctx.ledger, c, "S3", ctx.tiers)]
        cut = kw.pop("rank_cut", None)
        if cut is None:
            cut = int(ctx.tiers.threshold("S3.rank_cut_per_day").value)
        return run_s3(ctx.ledger, cids, targets=ctx.targets, tiers=ctx.tiers, machine=ctx.machine, rank_cut=cut,
                      workdir=ctx.workdir / "S3", **kw)


# ------------------------------------------------------------------------------------ S4 (B3 compute)
class S4Electronic(StageBase):
    id, kind, lane, requires = "S4", "physical", "S", ("S3",)

    def estimate_cost(self, rec, machine):
        return CostEstimate(core_h={"C1": 0.2, "C2": 2.0, "C3": 3.0}.get(rec.cost_class or "C2", 2.0), lane="S",
                            cost_class=rec.cost_class)

    def run(self, rec: CandidateRecord, ctx: Context) -> StageResult:
        from ..dft.electronic import run_electronic
        from ..dft.profiles import profiles_from_tiers
        L = ctx.ledger
        prof = profiles_from_tiers(ctx.tiers)["s4"]
        tr = L.start_tier_run(rec.candidate_id, "S4", "dft", profile="s4", params={"profile": "s4"}, lane="S")
        wd = ctx.workdir / rec.candidate_id / "S4"
        try:
            res = run_electronic(dft_structure(rec), wd, profile=prof, launcher=ctx.launcher,
                                 material_id=rec.candidate_id, pressure_GPa=rec.pressure_GPa)
            m = getattr(res, "metrics", None)
            metrics = (m() if callable(m) else m) or {}
            # H7-heur (campaign 1 fix): E_form vs elemental references at the same s4 profile
            metrics.update(self.formation_energy(res, prof, ctx))
            from ..dft.electronic import electronic_artifacts
            electronic_artifacts(res, ledger=L, tier_run_id=tr)       # RC1 N1: S4 removals re-verifiable
            run = getattr(res, "run", None)      # campaign 1: S4 tier_runs had no wall/core-h
            L.finish_tier_run(tr, "done", metrics=metrics, wall_s=getattr(run, "wall_s", None),
                              core_h=getattr(run, "core_h", None))
            return StageResult("S4", rec.candidate_id, True, metrics, tr, payload=res)
        except Exception as e:  # noqa: BLE001 - recorded as ERROR
            L.finish_tier_run(tr, "failed", metrics={"error": repr(e)})
            return StageResult("S4", rec.candidate_id, False, {}, tr, error=traceback.format_exc(limit=3))

    @staticmethod
    def formation_energy(res: Any, prof: Any, ctx: Context) -> dict:
        """{e_form_eV, element_refs} (e_form_eV None when any reference failed or the element
        is pure: then H7-heur stays NOT CHECKED). Hook ``ctx.extra['element_refs_fn']`` for tests."""
        from ..dft.element_refs import element_refs, formation_energy
        comp = dict(getattr(res, "composition", {}) or {})
        if len(comp) < 2:
            return {"e_form_eV": None, "element_refs": {"note": "elemental: no formation energy"}}
        fn = ctx.extra.get("element_refs_fn", element_refs)
        try:
            mu, info = fn(list(comp), prof, ctx.workdir / "element_refs", launcher=ctx.launcher,
                          cache_dir=ctx.workdir / "element_refs" / "cache", machine=ctx.machine)
        except Exception as ex:  # noqa: BLE001 - no E_form -> H7_HEUR_NOT_CHECKED
            return {"e_form_eV": None, "element_refs": {"error": f"{type(ex).__name__}: {ex}"[:300]}}
        return {"e_form_eV": formation_energy(getattr(res, "e_total_Ry", None), comp, mu), "element_refs": info}

    def decide(self, res: StageResult, th, tgt) -> Outcome:
        from ..dft.electronic import decide_s4
        if not res.ok:
            return Outcome("S4", ERROR, reason=(res.error or "error")[:500], tier_run_id=res.tier_run_id)
        # RC1 F1: sourced Threshold objects, never the raw TOML tables
        return decide_s4(res.payload, gap_max=tgt.threshold("hard.H1_metallic"),
                         moment_tol=tgt.threshold("hard.H2_nonmagnetic"),
                         ehull_heur_max=th.threshold("S4.ehull_heur_max_eV"),
                         e_form_eV=(res.metrics or {}).get("e_form_eV"), tier_run_id=res.tier_run_id)


# ------------------------------------------------------------------------------------ S5 coarse DFPT
class S5Coarse(StageBase):
    id, kind, lane, requires = "S5", "dfpt", "S", ("S4",)

    def run(self, rec: CandidateRecord, ctx: Context) -> StageResult:
        from ..dft.collect import collect_run
        from ..dft.coarse import coarse_metrics
        from ..dft.profiles import profiles_from_tiers
        from ..dft.workflow import run_coarse
        L = ctx.ledger
        from ..dft.profiles import apply_ecut_floor
        prof = profiles_from_tiers(ctx.tiers)["coarse"]
        # doc 16: per-element cutoff floor (H >= 80 Ry; SG15 H unconverged at 60 Ry, doc 15 §1)
        prof = apply_ecut_floor(prof, list(rec.structure.species),
                                ctx.tiers.get("S5.ecut_floor_Ry", None))
        # doc 13 T2-SCREEN: coarse runs on S5.q_meshes (q3^3 and q4^3), coarser first, on ONE relaxed
        # cell (run_coarse relaxes once), so delta_q compares one structure; one tier_run per mesh;
        # the finest mesh carries the S5 metrics (incl. per-mesh detail)
        meshes = [list(q) for q in (prof.q_meshes or ())] or [ctx.tiers.get("S5.q")]
        s4 = rec.stages.get("S4")
        s4m = s4.metrics if s4 else {}
        trs = [L.start_tier_run(rec.candidate_id, "S5", "dfpt", profile="coarse",
                                params={"q": q, "ecutwfc_Ry": prof.ecutwfc_Ry}, lane="S")
               for q in meshes]
        try:
            runs = run_coarse(dft_structure(rec), prof, ctx.workdir / rec.candidate_id / "S5",
                              ctx.launcher, material_id=rec.candidate_id, pressure_GPa=rec.pressure_GPa)
            if len(runs) != len(meshes):
                raise RuntimeError(f"run_coarse returned {len(runs)} runs for {len(meshes)} meshes")
            cols, mets = [], []
            for run, tr in zip(runs, trs):
                col = collect_run(run, ledger=L, tier_run_id=tr) if run.ok else None
                met = coarse_metrics(col, run, ctx.tiers.get("S5.plateau_sigma_min_Ry", None)) if col else None
                cols.append(col); mets.append(met)
            ok = all(r.ok for r in runs) and all(c is not None and c.ok for c in cols)
            metrics = mets[-1].to_metrics() if mets[-1] else {"error": runs[-1].error}
            metrics["by_mesh"] = [{"q": q, **(m.to_metrics() if m else {"error": r.error})}
                                  for q, m, r in zip(meshes, mets, runs)]
            metrics["n_sheets"], metrics["sheet_w_min"] = s4m.get("n_sheets"), s4m.get("sheet_w_min")
            for i, (tr, run) in enumerate(zip(trs, runs)):
                L.finish_tier_run(tr, "done" if run.ok else "failed", wall_s=run.wall_s, core_h=run.core_h,
                                  metrics=metrics if i == len(trs) - 1 else {"q": meshes[i]})
            err = None
            if not ok:      # campaign 1: the ERROR outcome said only "error"; name the failing mesh and why
                err = "; ".join(f"q{q}: " + (r.error or "; ".join(getattr(c, "errors", None) or []) or "collect failed")
                                for q, r, c in zip(meshes, runs, cols) if not (r.ok and c is not None and c.ok))[:1500]
            return StageResult("S5", rec.candidate_id, ok, metrics, trs[-1], error=err,
                               payload={"metrics": mets[-1], "metrics_by_mesh": mets, "collected": cols[-1],
                                        "collected_by_mesh": cols, "runs": runs})
        except Exception:  # noqa: BLE001
            for tr in trs:
                if (L.tier_run(tr) or {}).get("state") == "started":
                    L.finish_tier_run(tr, "failed")
            return StageResult("S5", rec.candidate_id, False, {}, trs[-1], error=traceback.format_exc(limit=3))

    def decide(self, res: StageResult, th, tgt) -> Outcome:
        """Heuristic removal only: SCREENED_OUT(coarse-A) with M_c (provisional) — B3's rule."""
        from ..dft.coarse import decide_coarse
        from ..dft.admission import admission_deferral
        if (o := admission_deferral(res)) is not None:     # doc 31 §4: one mesh failed admission only
            return o
        pay = res.payload or {}
        met = pay.get("metrics") if res.ok else None
        by_mesh = pay.get("metrics_by_mesh") or []
        coarser = by_mesh[-2] if res.ok and len(by_mesh) >= 2 else None      # doc 13 T2-SCREEN (q3^3 + q4^3)
        return decide_coarse(met, floor_K=tgt.floor, margin_M_c=th.threshold("S5.margin_M_c"),
                             e_c=th.threshold("S5.e_c"), metrics_coarser=coarser,
                             eps_sys=tgt.threshold("envelope.eps_sys"), tier_run_id=res.tier_run_id,
                             n_sheets=res.metrics.get("n_sheets"), w_min=res.metrics.get("sheet_w_min"),
                             error=None if res.ok else (res.error or "error"),
                             ecut_floors=(th.get("S5.ecut_floor_Ry", None)        # doc 16 §3
                                          if hasattr(th, "get") else None))

    def redecide(self, rec: CandidateRecord, ctx: Context) -> Optional[Outcome]:
        """Doc 16: re-run :func:`decide_coarse` (current rule and tiers) on the S5 metrics stored
        in the ledger (no compute). Returns the new Outcome, or None when the current S5 outcome
        has no stored per-mesh metrics. The caller appends it (outcomes are append-only; the
        latest row is current, R5). Records written before doc 16 lack the cutoff: it is read
        from the run's ``dfpt_run.json`` under the work dir when present, else left unknown
        (unknown cutoff blocks a removal when a floor applies)."""
        from ..dft.coarse import CoarseMetrics, decide_coarse
        L = ctx.ledger
        cur = L.current_outcome(rec.candidate_id, "S5")
        if not cur or not cur.get("tier_run_id"):
            return None
        tr = L.tier_run(cur["tier_run_id"]) or {}
        met = json.loads(tr.get("metrics_json") or "null") or {}
        by_mesh = [m for m in (met.get("by_mesh") or []) if "rows" in m]
        if not by_mesh:
            return None
        els = sorted(set(rec.structure.species))
        params = json.loads(tr.get("params_json") or "null") or {}
        mets = []
        for m in by_mesh:
            ec = m.get("ecutwfc_Ry") or params.get("ecutwfc_Ry")
            if ec is None and m.get("q"):
                f = ctx.workdir / rec.candidate_id / "S5" / ("q" + "x".join(str(x) for x in m["q"])) / "dfpt_run.json"
                try:
                    ec = (json.loads(f.read_text()).get("settings") or {}).get("ecutwfc_Ry")
                except Exception:  # noqa: BLE001 - missing/unreadable: cutoff stays unknown
                    ec = None
            mets.append(CoarseMetrics.from_metrics(m, elements=els, ecutwfc_Ry=ec))
        tgt = ctx.targets_for(rec)
        th = ctx.tiers
        o = decide_coarse(mets[-1], floor_K=tgt.floor, margin_M_c=th.threshold("S5.margin_M_c"),
                          e_c=th.threshold("S5.e_c"), metrics_coarser=mets[-2] if len(mets) >= 2 else None,
                          eps_sys=tgt.threshold("envelope.eps_sys"), tier_run_id=cur["tier_run_id"],
                          n_sheets=met.get("n_sheets"), w_min=met.get("sheet_w_min"),
                          ecut_floors=th.get("S5.ecut_floor_Ry", None))
        o.reason = f"redecide (doc 16 rule; supersedes outcome #{cur['id']} {cur['outcome']}): " + (o.reason or "")
        o.tags = list(o.tags) + ["REDECIDED_DOC16"]
        return o

    def evidence(self, rec: CandidateRecord, ctx: Context, res: StageResult) -> None:
        """Core certify in decide mode on the coarse spectra (evidence only, spectrum_tier
        'coarse') and decide_t at coarse provenance (recorded, never a stage outcome)."""
        from ..certify import certify, certify_config
        from ..decide.inputs import tags_from_spectra
        col = (res.payload or {}).get("collected")
        if not res.ok or col is None:
            return
        tgt = ctx.targets_for(rec)
        sp = [s if s.rel_envelope >= 0.10 else s.with_envelope(0.10) for s in col.spectra]
        cc = certify_config(float(tgt.T_star), 0.15, mode="decide")
        v = certify(rec.candidate_id, sp, cc)
        store_verdict(ctx, rec.candidate_id, res.tier_run_id, "coarse", v, sp, cc, ctx.workdir / rec.candidate_id / "S5")


# ------------------------------------------------------------------------------------ S6 full DFPT
class S6Full(StageBase):
    id, kind, lane, requires = "S6", "dfpt", "F", ("S5",)

    def run(self, rec: CandidateRecord, ctx: Context) -> StageResult:
        from ..dft.collect import collect_run
        from ..dft.profiles import profiles_from_tiers
        from ..dft.workflow import run_full
        L = ctx.ledger
        prof = profiles_from_tiers(ctx.tiers)["full"]
        meshes = ctx.tiers.get("S6.q_meshes")
        trs = [L.start_tier_run(rec.candidate_id, "S6", "dfpt", profile="full", params={"q": q}, lane="F") for q in meshes]
        try:
            runs = run_full(dft_structure(rec), prof, ctx.workdir / rec.candidate_id / "S6", ctx.launcher,
                            material_id=rec.candidate_id, pressure_GPa=rec.pressure_GPa)
            cols = [collect_run(run, ledger=L, tier_run_id=tr) if run.ok else None for tr, run in zip(trs, runs)]
            ok = len(cols) == len(trs) and all(c is not None and c.ok for c in cols)
            metrics = {}
            if ok:
                # full-grid stability of the final structure: min over both meshes (dyn files + matdyn)
                from ..dft.coarse import coarse_metrics
                mf = [coarse_metrics(c, r).min_freq_K for c, r in zip(cols, runs)]
                mf = [v for v in mf if v is not None]
                mfk = min(mf) if len(mf) == len(cols) else None      # missing on any mesh -> MISSING
                metrics = {"ladder": ladder_from_runs(cols), "run_id": trs[-1],
                           "omega_max_K": max(s.omega_max() for s in cols[-1].spectra),
                           "min_freq_K": mfk, "max_imag_mode_K": None if mfk is None else max(0.0, -mfk),
                           "mu_star_negative_suspected": False}
            # the finest mesh run carries the S6 metrics (tier_runs are insert-then-finish-once)
            for i, (tr, run) in enumerate(zip(trs, runs)):
                m = metrics if i == len(trs) - 1 else {"q": run.q, "error": run.error or None}
                L.finish_tier_run(tr, "done" if run.ok else "failed", wall_s=run.wall_s, core_h=run.core_h, metrics=m)
            for tr in trs[len(runs):]:
                L.finish_tier_run(tr, "failed")
            return StageResult("S6", rec.candidate_id, ok, metrics, trs[-1], payload={"collected": cols, "runs": runs})
        except Exception:  # noqa: BLE001
            for tr in trs:
                if (L.tier_run(tr) or {}).get("state") == "started":
                    L.finish_tier_run(tr, "failed")
            return StageResult("S6", rec.candidate_id, False, {}, trs[-1], error=traceback.format_exc(limit=3))

    def certify_and_decide(self, rec: CandidateRecord, ctx: Context, res: StageResult) -> Outcome:
        from ..certify import certify, certify_config
        from ..certify.alternatives import envelope_e, freq_scaled
        from ..decide.inputs import build_inputs, tags_from_spectra
        if not res.ok:
            o = Outcome("S6", ERROR, reason=(res.error or "error")[:500], tier_run_id=res.tier_run_id)
            write_outcome(ctx, rec.candidate_id, o, rec=rec)
            return o
        cols = res.payload["collected"]
        tgt = ctx.targets_for(rec)
        lad = res.metrics["ladder"]
        A_f, A_c = statistics.median(lad["fine"]["A_K"]), statistics.median(lad["coarse"]["A_K"])
        dq = abs(A_f - A_c) / A_f
        eps_sys = float(tgt.threshold("envelope.eps_sys").value)
        e = envelope_e(dq, eps_sys)
        eps_omega = max(float(tgt.threshold("envelope.eps_omega_min").value),
                        abs(lad["omega_log_K_by_q"][-1] - lad["omega_log_K_by_q"][-2]) / lad["omega_log_K_by_q"][-1]
                        if len(lad["omega_log_K_by_q"]) > 1 else 0.0)
        base = cols[-1].spectra
        vid, vstatus = None, None
        if e is not None:
            alts = freq_scaled([s.with_envelope(e) for s in base], eps_omega)
            cc_cfg = ctx.tiers.get("S6.certify")
            cc = certify_config(float(tgt.T_star), float(cc_cfg["mu_hi"]), mode=cc_cfg.get("mode", "bracket"),
                                omega_c_factor=float(cc_cfg.get("omega_c_factor", 10.0)))
            v = certify(rec.candidate_id, alts, cc)
            vid = store_verdict(ctx, rec.candidate_id, res.tier_run_id, "full", v, alts, cc,
                                ctx.workdir / rec.candidate_id / "S6")
            vstatus = v.status.value
        rec2 = ctx.ledger.candidate_record(rec.candidate_id)
        from ..records import StageEntry
        rec2.stages["S6"] = StageEntry(outcome="pending", tier_run_id=res.tier_run_id, metrics=res.metrics)
        di = build_inputs(rec2, "fullDFPT", stage="S6", tags=tags_from_spectra(base), eps_sys=eps_sys)
        di.provenance["tier_run_ids"] = sorted(set(di.provenance.get("tier_run_ids", []) + [res.tier_run_id]))
        _, out = record_dfpt_decision(ctx, rec, "S6", res.tier_run_id, di, vid, vstatus)
        return out


# ------------------------------------------------------------------------------------ S7-S9
class S7Validity(StageBase):
    """Validity closure: own-DFT hull + ZPE (H7), anharmonic check, k-resolved data, then decideT
    is re-run. Compute pieces are B3/B2 work in progress; this stage only re-decides from records."""
    id, kind, lane, requires = "S7", "physical", "F", ("S6",)

    def redecide(self, rec: CandidateRecord, ctx: Context) -> Optional[Outcome]:
        from ..decide.inputs import build_inputs
        s6 = rec.stages.get("S6")
        if s6 is None or s6.tier_run_id is None:
            return None
        di = build_inputs(rec, "fullDFPT", stage="S6")
        v = ctx.ledger.q1("SELECT id, status FROM verdicts WHERE candidate_id=? AND tier_run_id=? ORDER BY id DESC",
                          (rec.candidate_id, s6.tier_run_id))
        _, out = record_dfpt_decision(ctx, rec, "S7", s6.tier_run_id, di, v[0] if v else None, v[1] if v else None)
        return out


class S8Converged(StageBase):
    id, kind, lane, requires = "S8", "dfpt", "F", ("S7",)       # EPW protocol: not implemented (adapt later)


class S9Evidence(StageBase):
    id, kind, lane, requires = "S9", "report", "none", ("S6",)

    def statement(self, cid: str, ctx: Context) -> str:
        from ..report.evidence import evidence_statement
        return evidence_statement(cid, ctx.ledger)


STAGES: dict[str, StageBase] = {s.id: s for s in (S1Source(), S2Structure(), S3MLFF(), S4Electronic(), S5Coarse(),
                                                  S6Full(), S7Validity(), S8Converged(), S9Evidence())}
