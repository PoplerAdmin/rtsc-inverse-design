"""SQLite ledger (architecture §6.2). Append-only; the triggers in ``schema.sql`` enforce the
soundness contract R1b/R1c/R2-R5 at insert time, independent of this Python layer."""
from __future__ import annotations

import contextlib
import datetime as _dt
import hashlib
import json
import sqlite3
from pathlib import Path
from typing import Any, Iterable, Iterator, Optional, Sequence

from ..records import (CandidateRecord, Context, DecideInputs, Origin, Outcome, StageEntry,
                       StructureData, ORIGIN_CERTIFY, ORIGIN_DECIDE)

SCHEMA_PATH = Path(__file__).with_name("schema.sql")
SCHEMA_ID = "rtt.ledger/1"


def now() -> str:
    return _dt.datetime.now(_dt.timezone.utc).isoformat(timespec="seconds")


def _j(x: Any) -> Optional[str]:
    if x is None:
        return None
    return json.dumps(x, sort_keys=True, default=str)


class LedgerError(RuntimeError):
    pass


class Ledger:
    def __init__(self, path: str | Path, *, run_id: str | None = None):
        self.path = str(path)
        if self.path != ":memory:":
            Path(self.path).parent.mkdir(parents=True, exist_ok=True)
        self.conn = sqlite3.connect(self.path, isolation_level=None, timeout=60)
        self.conn.execute("PRAGMA foreign_keys = ON")
        if self.path != ":memory:":
            self.conn.execute("PRAGMA journal_mode = WAL")
        self.conn.executescript(SCHEMA_PATH.read_text(encoding="utf-8"))
        self._migrate_stage_outcomes()
        self._meta_default("schema", SCHEMA_ID)
        self._meta_default("schema_sql_sha256", hashlib.sha256(SCHEMA_PATH.read_bytes()).hexdigest())
        if run_id:
            self._meta_default("run_id", run_id)

    #: (stage, outcome) pairs added to the fixed §3.1 table after ledgers already existed
    MIGRATED_PAIRS = (("S1", "DEFERRED"),)            # R30 F1(b), doc 31 §7

    def _migrate_stage_outcomes(self) -> None:
        """One-time, recorded migration of an older ledger's fixed ``stage_outcomes`` table: add the
        missing pairs of :data:`MIGRATED_PAIRS` and re-create the insert guard for the new size.
        The guard is the schema's own trigger text; a ``meta`` row records the migration."""
        missing = [p for p in self.MIGRATED_PAIRS
                   if self.q1("SELECT 1 FROM stage_outcomes WHERE stage=? AND outcome=?", p) is None]
        if not missing:
            return
        import re
        sql = SCHEMA_PATH.read_text(encoding="utf-8")
        m = re.search(r"CREATE TRIGGER IF NOT EXISTS stage_outcomes_no_insert .*?END;", sql, re.S)
        with self.tx() as c:
            c.execute("DROP TRIGGER IF EXISTS stage_outcomes_no_insert")
            for p in missing:
                c.execute("INSERT INTO stage_outcomes(stage, outcome) VALUES (?, ?)", p)
            c.execute(m.group(0))
            c.execute("INSERT OR IGNORE INTO meta(key, value) VALUES (?, ?)",
                      ("migration:stage_outcomes:" + ",".join("/".join(p) for p in missing), now()))

    # ---------------------------------------------------------------- basics
    def close(self) -> None:
        self.conn.close()

    def __enter__(self):
        return self

    def __exit__(self, *a):
        self.close()

    @contextlib.contextmanager
    def tx(self) -> Iterator[sqlite3.Connection]:
        self.conn.execute("BEGIN IMMEDIATE")
        try:
            yield self.conn
        except BaseException:
            self.conn.execute("ROLLBACK")
            raise
        else:
            self.conn.execute("COMMIT")

    def x(self, sql: str, args: Sequence = ()) -> sqlite3.Cursor:
        return self.conn.execute(sql, args)

    def q1(self, sql: str, args: Sequence = ()):
        return self.conn.execute(sql, args).fetchone()

    def qa(self, sql: str, args: Sequence = ()) -> list:
        return self.conn.execute(sql, args).fetchall()

    def _meta_default(self, key: str, value: str) -> None:
        self.x("INSERT OR IGNORE INTO meta(key, value) VALUES (?, ?)", (key, value))

    def meta(self, key: str) -> Optional[str]:
        r = self.q1("SELECT value FROM meta WHERE key = ?", (key,))
        return r[0] if r else None

    def set_meta_once(self, key: str, value: str) -> None:
        cur = self.meta(key)
        if cur is not None and cur != value:
            raise LedgerError(f"meta {key!r} already set to {cur!r}")
        self._meta_default(key, value)

    # ---------------------------------------------------------------- inputs / structures / candidates
    def register_input(self, kind: str, sha256: str, content: str) -> None:
        self.x("INSERT OR IGNORE INTO inputs(kind, sha256, content, created) VALUES (?,?,?,?)",
               (kind, sha256, content, now()))

    def register_configs(self, *cfgs) -> None:
        for c in cfgs:
            self.register_input(c.kind, c.sha256, c.text)

    def insert_structure(self, sha256: str, s: StructureData | dict) -> None:
        js = s.dumps() if isinstance(s, StructureData) else json.dumps(s, sort_keys=True)
        self.x("INSERT OR IGNORE INTO structures(sha256, json) VALUES (?, ?)", (sha256, js))

    def insert_candidate(self, cid: str, structure_sha256: str, structure: StructureData | dict,
                         origin: Origin, pressure_GPa: float = 0.0, cost_class: str | None = None) -> str:
        with self.tx():
            self.insert_structure(structure_sha256, structure)
            self.x("INSERT INTO candidates(id, structure_sha256, lane, origin_json, cost_class, pressure_GPa, created)"
                   " VALUES (?,?,?,?,?,?,?)",
                   (cid, structure_sha256, origin.lane, origin.dumps(), cost_class, float(pressure_GPa), now()))
        return cid

    def has_candidate(self, cid: str) -> bool:
        return self.q1("SELECT 1 FROM candidates WHERE id = ?", (cid,)) is not None

    def candidate_ids(self, lane: str | None = None) -> list[str]:
        if lane:
            return [r[0] for r in self.qa("SELECT id FROM candidates WHERE lane = ? ORDER BY created, id", (lane,))]
        return [r[0] for r in self.qa("SELECT id FROM candidates ORDER BY created, id")]

    def insert_duplicate(self, cid: str, dup_of: str, method: str) -> None:
        self.x("INSERT INTO duplicates VALUES (?,?,?)", (cid, dup_of, method))

    # ---------------------------------------------------------------- tier runs and artifacts
    def start_tier_run(self, candidate_id: str, stage: str, kind: str, *, profile: str | None = None,
                       params: dict | None = None, versions: dict | None = None,
                       lane: str | None = None) -> int:
        """Written BEFORE launch (W12)."""
        attempt = 1 + (self.q1("SELECT COUNT(*) FROM tier_runs WHERE candidate_id=? AND stage=? AND kind=?",
                               (candidate_id, stage, kind))[0])
        cur = self.x("INSERT INTO tier_runs(candidate_id, stage, kind, profile, attempt, state, params_json,"
                     " versions_json, lane, started) VALUES (?,?,?,?,?,'started',?,?,?,?)",
                     (candidate_id, stage, kind, profile, attempt, _j(params or {}), _j(versions or {}),
                      lane, now()))
        return int(cur.lastrowid)

    def finish_tier_run(self, tier_run_id: int, state: str = "done", *, wall_s: float | None = None,
                        core_h: float | None = None, gpu_h: float | None = None,
                        metrics: dict | None = None) -> None:
        self.x("UPDATE tier_runs SET state=?, finished=?, wall_s=?, core_h=?, gpu_h=?, metrics_json=? WHERE id=?",
               (state, now(), wall_s, core_h, gpu_h, _j(metrics), tier_run_id))

    def tier_run(self, tier_run_id: int) -> Optional[dict]:
        cur = self.x("SELECT * FROM tier_runs WHERE id = ?", (tier_run_id,))
        row = cur.fetchone()
        if row is None:
            return None
        return dict(zip([d[0] for d in cur.description], row))

    def insert_artifact(self, tier_run_id: int, path: str, sha256: str, nbytes: int, retained: bool = True) -> int:
        cur = self.x("INSERT INTO artifacts(tier_run_id, path, sha256, bytes, retained) VALUES (?,?,?,?,?)",
                     (tier_run_id, path, sha256, int(nbytes), 1 if retained else 0))
        return int(cur.lastrowid)

    def insert_artifact_event(self, artifact_id: int, event: str, note: str | None) -> None:
        self.x("INSERT INTO artifact_events(artifact_id, event, note, created) VALUES (?,?,?,?)",
               (artifact_id, event, note, now()))

    def artifact_pruned(self, artifact_id: int) -> bool:
        return self.q1("SELECT 1 FROM artifact_events WHERE artifact_id=? AND event='pruned'", (artifact_id,)) is not None

    # ---------------------------------------------------------------- outcomes
    def insert_outcome(self, candidate_id: str, o: Outcome, *, targets_sha256: str, tiers_sha256: str,
                       origin: str, decision_id: int | None = None) -> int:
        cur = self.x(
            "INSERT INTO outcomes(candidate_id, stage, outcome, criterion, value, threshold, source_id, provisional,"
            " tier_run_id, targets_sha256, tiers_sha256, origin, reason, decision_id, tags, created)"
            " VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?)",
            (candidate_id, o.stage, o.outcome, o.criterion, _j(o.value), _j(o.threshold), o.source_id,
             1 if o.provisional else 0, o.tier_run_id, targets_sha256, tiers_sha256, origin,
             o.reason, decision_id, _j(o.tags) if o.tags else None, now()))
        return int(cur.lastrowid)

    def current_outcome(self, candidate_id: str, stage: str) -> Optional[dict]:
        """Latest outcome row of a stage (re-decisions append; R5: exactly one *current*)."""
        cur = self.x("SELECT * FROM outcomes WHERE candidate_id=? AND stage=? ORDER BY id DESC LIMIT 1",
                     (candidate_id, stage))
        row = cur.fetchone()
        if row is None:
            return None
        d = dict(zip([c[0] for c in cur.description], row))
        for k in ("value", "threshold", "tags"):
            if d.get(k) is not None:
                d[k] = json.loads(d[k])
        return d

    def current_outcomes(self, candidate_id: str) -> dict[str, dict]:
        stages = [r[0] for r in self.qa("SELECT DISTINCT stage FROM outcomes WHERE candidate_id=?", (candidate_id,))]
        return {s: self.current_outcome(candidate_id, s) for s in stages}

    # ---------------------------------------------------------------- decide inputs / decisions / verdicts
    def insert_decide_input(self, candidate_id: str, di: DecideInputs) -> int:
        cur = self.x("INSERT INTO decide_inputs(candidate_id, json, created) VALUES (?,?,?)",
                     (candidate_id, di.dumps(), now()))
        return int(cur.lastrowid)

    def decide_input(self, decide_input_id: int) -> DecideInputs:
        r = self.q1("SELECT json FROM decide_inputs WHERE id=?", (decide_input_id,))
        return DecideInputs.from_json(json.loads(r[0]))

    def insert_decision(self, candidate_id: str, stage: str, tier: str, decide_input_id: int, outcome: str, *,
                        firing_rule: str | None = None, assumptions: Iterable[str] = (),
                        crosscheck_verdict_id: int | None = None, origin: str = ORIGIN_DECIDE) -> int:
        cur = self.x("INSERT INTO decisions(candidate_id, stage, tier, decide_input_id, outcome, firing_rule,"
                     " assumptions, crosscheck_verdict_id, origin, created) VALUES (?,?,?,?,?,?,?,?,?,?)",
                     (candidate_id, stage, tier, decide_input_id, outcome, firing_rule,
                      _j(list(assumptions)), crosscheck_verdict_id, origin, now()))
        return int(cur.lastrowid)

    def insert_verdict(self, candidate_id: str, tier_run_id: int, spectrum_tier: str, *, status: str,
                       T_low_K: float | None, T_up_K: float | None, T_star_K: float, mu_hi: float,
                       rel_envelope: float, rule_ids: Iterable[str] = (), assumptions: Iterable[str] = (),
                       flags: Iterable[str] = (), certificate_artifact_id: int, recheck_ok: bool,
                       spectra: Sequence[dict], origin: str = ORIGIN_CERTIFY) -> int:
        """``spectra``: [{alt_index, spectrum_digest, source_artifact_id, sha256}] (R4)."""
        spectra = [dict(s) for s in spectra]
        with self.tx():
            cur = self.x(
                "INSERT INTO verdicts(candidate_id, tier_run_id, spectrum_tier, status, T_low_K, T_up_K, T_star_K,"
                " mu_hi, rel_envelope, rule_ids, assumptions, flags, certificate_artifact_id, recheck_ok,"
                " spectra_json, origin, created) VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?)",
                (candidate_id, tier_run_id, spectrum_tier, status, T_low_K, T_up_K, float(T_star_K), float(mu_hi),
                 float(rel_envelope), _j(list(rule_ids)), _j(list(assumptions)), _j(list(flags)),
                 certificate_artifact_id, 1 if recheck_ok else 0, json.dumps(spectra, sort_keys=True),
                 origin, now()))
            vid = int(cur.lastrowid)
            for s in spectra:
                self.x("INSERT INTO verdict_spectra(verdict_id, alt_index, spectrum_digest, source_artifact_id)"
                       " VALUES (?,?,?,?)", (vid, s["alt_index"], s["spectrum_digest"], s["source_artifact_id"]))
        return vid

    def verdict(self, verdict_id: int) -> Optional[dict]:
        cur = self.x("SELECT * FROM verdicts WHERE id=?", (verdict_id,))
        row = cur.fetchone()
        return None if row is None else dict(zip([c[0] for c in cur.description], row))

    def verdicts_for(self, candidate_id: str) -> list[dict]:
        cur = self.x("SELECT * FROM verdicts WHERE candidate_id=? ORDER BY id", (candidate_id,))
        cols = [c[0] for c in cur.description]
        return [dict(zip(cols, r)) for r in cur.fetchall()]

    # ---------------------------------------------------------------- misc tables
    def insert_label(self, candidate_id: str, name: str, value: float, stage: str, tier_run_id: int | None = None) -> None:
        self.x("INSERT INTO labels VALUES (?,?,?,?,?,?)", (candidate_id, name, float(value), stage, tier_run_id, now()))

    def insert_calibration(self, material: str, coarse_run: int | None, ref_run: int | None, sigma: int,
                           A_coarse: float, A_ref: float, note: str = "") -> int:
        cur = self.x("INSERT INTO calibration(material, coarse_run, ref_run, sigma, A_coarse, A_ref, note)"
                     " VALUES (?,?,?,?,?,?,?)", (material, coarse_run, ref_run, sigma, A_coarse, A_ref, note))
        return int(cur.lastrowid)

    def insert_audit(self, candidate_id: str, audit: str, inclusion_prob: float) -> None:
        self.x("INSERT INTO audits VALUES (?,?,?,?)", (candidate_id, audit, float(inclusion_prob), now()))

    def set_screen_enabled(self, screen: str, enabled: bool, reason: str, actor: str) -> None:
        self.x("INSERT INTO screen_status(screen, enabled, reason, actor, created) VALUES (?,?,?,?,?)",
               (screen, 1 if enabled else 0, reason, actor, now()))

    def screen_enabled(self, screen: str) -> bool:
        r = self.q1("SELECT enabled FROM screen_status WHERE screen=? ORDER BY id DESC LIMIT 1", (screen,))
        return True if r is None else bool(r[0])

    # rank pools (B9, doc 31 §1): insert-only record of each pool-wide rank cut
    def insert_rank_pool(self, rec: dict) -> int:
        cur = self.x("INSERT INTO rank_pools(pool_id, stage, cut, free_slots, json, ranking_sha256, tiers_sha256,"
                     " created) VALUES (?,?,?,?,?,?,?,?)",
                     (rec["pool_id"], rec.get("stage", "S3"), int(rec["cut"]), int(rec["free_slots"]),
                      _j(rec), rec["ranking_sha256"], rec.get("tiers_sha256"), now()))
        return int(cur.lastrowid)

    def rank_pools(self, stage: str = "S3") -> list[dict]:
        return [json.loads(r[0]) for r in self.qa("SELECT json FROM rank_pools WHERE stage=? ORDER BY id", (stage,))]

    # queue: derived and mutable
    def queue_put(self, candidate_id: str, stage: str, priority: float, state: str = "QUEUED") -> None:
        self.x("INSERT INTO queue VALUES (?,?,?,?) ON CONFLICT(candidate_id, stage) DO UPDATE SET"
               " priority=excluded.priority, state=excluded.state", (candidate_id, stage, priority, state))

    def queue_items(self, stage: str | None = None, state: str = "QUEUED") -> list[tuple[str, str, float]]:
        if stage:
            return self.qa("SELECT candidate_id, stage, priority FROM queue WHERE stage=? AND state=?"
                           " ORDER BY priority DESC, candidate_id", (stage, state))
        return self.qa("SELECT candidate_id, stage, priority FROM queue WHERE state=?"
                       " ORDER BY priority DESC, candidate_id", (state,))

    def queue_clear(self) -> None:
        self.x("DELETE FROM queue")

    # ---------------------------------------------------------------- record assembly / export
    def candidate_record(self, candidate_id: str) -> CandidateRecord:
        r = self.q1("SELECT c.structure_sha256, c.origin_json, c.cost_class, c.pressure_GPa, s.json"
                    " FROM candidates c JOIN structures s ON s.sha256 = c.structure_sha256 WHERE c.id=?",
                    (candidate_id,))
        if r is None:
            raise KeyError(candidate_id)
        sha, origin_json, cost_class, p, sjson = r
        st = StructureData.from_json(json.loads(sjson))
        rec = CandidateRecord(candidate_id=candidate_id, structure_sha256=sha, formula=st.formula, structure=st,
                              origin=Origin.from_json(json.loads(origin_json)), pressure_GPa=p,
                              cost_class=cost_class)
        ctx = None
        for stage, o in sorted(self.current_outcomes(candidate_id).items()):
            if ctx is None:
                ctx = Context(targets_sha256=o["targets_sha256"], tiers_sha256=o["tiers_sha256"],
                              code_git_sha=self.meta("code_git_sha"))
            metrics = {}
            profile = None
            if o.get("tier_run_id"):
                tr = self.tier_run(o["tier_run_id"]) or {}
                metrics = json.loads(tr.get("metrics_json") or "null") or {}
                profile = tr.get("profile")
            decide = None
            if o.get("decision_id"):
                d = self.q1("SELECT outcome, firing_rule, decide_input_id, assumptions FROM decisions WHERE id=?",
                            (o["decision_id"],))
                decide = {"outcome": d[0], "rule": d[1], "decide_input_id": d[2],
                          "assumptions": json.loads(d[3]) if d[3] else []}
            vids = [v[0] for v in self.qa("SELECT id FROM verdicts WHERE candidate_id=? AND tier_run_id=?",
                                          (candidate_id, o.get("tier_run_id")))] if o.get("tier_run_id") else []
            rec.stages[stage] = StageEntry(outcome=o["outcome"], tier_run_id=o.get("tier_run_id"), metrics=metrics,
                                           profile=profile, decide=decide, verdict_ids=vids,
                                           flags=o.get("tags") or [], criterion=o.get("criterion"),
                                           reason=o.get("reason"))
        rec.context = ctx
        return rec

    TABLES = ("meta", "inputs", "structures", "candidates", "duplicates", "tier_runs", "artifacts",
              "artifact_events", "outcomes", "decide_inputs", "decisions", "verdicts", "verdict_spectra",
              "labels", "calibration", "audits", "screen_status", "rank_pools")

    def export_jsonl(self, outdir: str | Path, verify_artifacts: bool = True) -> Path:
        """One JSONL file per table + ``manifest.sha256``. Re-verifies every retained artifact
        it can still read (appends verified/missing/mismatch events) before exporting."""
        from .artifacts import verify
        outdir = Path(outdir)
        outdir.mkdir(parents=True, exist_ok=True)
        if verify_artifacts:
            for (aid,) in self.qa("SELECT id FROM artifacts"):
                verify(self, aid)
        lines = []
        for t in self.TABLES:
            cur = self.x(f"SELECT * FROM {t} ORDER BY rowid")
            cols = [c[0] for c in cur.description]
            p = outdir / f"{t}.jsonl"
            with open(p, "w", encoding="utf-8", newline="\n") as fh:
                for row in cur.fetchall():
                    fh.write(json.dumps(dict(zip(cols, row)), sort_keys=True) + "\n")
            lines.append(f"{hashlib.sha256(p.read_bytes()).hexdigest()}  {p.name}")
        man = outdir / "manifest.sha256"
        man.write_text("\n".join(lines) + "\n", encoding="utf-8")
        return man

    def summary_counts(self) -> dict:
        """R5: current outcome per (candidate, stage); REJECTED_PROVEN split by spectrum tier."""
        out: dict[str, dict[str, int]] = {}
        for (cid,) in self.qa("SELECT id FROM candidates"):
            for stage, o in self.current_outcomes(cid).items():
                out.setdefault(stage, {}).setdefault(o["outcome"], 0)
                out[stage][o["outcome"]] += 1
        vt = {}
        for tier, n in self.qa("SELECT spectrum_tier, COUNT(*) FROM verdicts WHERE status='REJECTED_PROVEN'"
                               " GROUP BY spectrum_tier"):
            vt[tier] = n
        return {"outcomes": out, "core_rejected_proven_by_tier": vt,
                "candidates": self.q1("SELECT COUNT(*) FROM candidates")[0]}
