-- rtt ledger schema (architecture §6.2). SQLite >= 3.38 (JSON1 built in), WAL.
-- All tables except `queue` are insert-only. `tier_runs` allows exactly one state transition
-- (started -> done|failed) that fills the completion columns; nothing else changes.
-- Deviations from §6.2 are listed in docs/design/20_build_core.md.
PRAGMA foreign_keys = ON;

CREATE TABLE IF NOT EXISTS meta (key TEXT PRIMARY KEY, value TEXT NOT NULL);

CREATE TABLE IF NOT EXISTS inputs (
  id INTEGER PRIMARY KEY,
  kind TEXT NOT NULL CHECK (kind IN ('targets','tiers','machine')),
  sha256 TEXT UNIQUE NOT NULL, content TEXT NOT NULL, created TEXT NOT NULL);

CREATE TABLE IF NOT EXISTS structures (sha256 TEXT PRIMARY KEY, json TEXT NOT NULL);

CREATE TABLE IF NOT EXISTS candidates (
  id TEXT PRIMARY KEY,
  structure_sha256 TEXT NOT NULL REFERENCES structures,
  lane TEXT NOT NULL CHECK (lane IN ('G1','G2','G3','G4','G5','GA','D','C','P')),   -- P: doc 32
  origin_json TEXT NOT NULL, cost_class TEXT CHECK (cost_class IS NULL OR cost_class IN ('C1','C2','C3')),
  pressure_GPa REAL NOT NULL, created TEXT NOT NULL);

CREATE TABLE IF NOT EXISTS duplicates (
  candidate_id TEXT NOT NULL REFERENCES candidates, duplicate_of TEXT NOT NULL REFERENCES candidates,
  method TEXT NOT NULL);

CREATE TABLE IF NOT EXISTS tier_runs (
  id INTEGER PRIMARY KEY,
  candidate_id TEXT NOT NULL REFERENCES candidates,
  stage TEXT NOT NULL,
  kind TEXT NOT NULL CHECK (kind IN ('ml','dft','dfpt','epw','md','report')),
  profile TEXT, attempt INTEGER NOT NULL,
  state TEXT NOT NULL CHECK (state IN ('started','done','failed')),
  params_json TEXT NOT NULL, versions_json TEXT NOT NULL, lane TEXT,
  started TEXT, finished TEXT, wall_s REAL, core_h REAL, gpu_h REAL, metrics_json TEXT);

CREATE TABLE IF NOT EXISTS artifacts (
  id INTEGER PRIMARY KEY,
  tier_run_id INTEGER NOT NULL REFERENCES tier_runs,
  path TEXT NOT NULL, sha256 TEXT NOT NULL CHECK (length(sha256) = 64),
  bytes INTEGER NOT NULL CHECK (bytes >= 0), retained INTEGER NOT NULL);

-- prune events (artifacts are insert-only; pruning is recorded, never an UPDATE)
CREATE TABLE IF NOT EXISTS artifact_events (
  id INTEGER PRIMARY KEY, artifact_id INTEGER NOT NULL REFERENCES artifacts,
  event TEXT NOT NULL CHECK (event IN ('pruned','verified','missing','mismatch')),
  note TEXT, created TEXT NOT NULL);

-- allowed (stage, outcome) pairs: the §3.1 table (immutable; see triggers)
CREATE TABLE IF NOT EXISTS stage_outcomes (stage TEXT NOT NULL, outcome TEXT NOT NULL,
  PRIMARY KEY (stage, outcome));
INSERT INTO stage_outcomes (stage, outcome) SELECT column1, column2 FROM (VALUES
  ('S1','PASS'),('S1','SCREENED_OUT'),('S1','DEFERRED'),
  ('S2','PASS'),('S2','INVALID'),('S2','DUPLICATE'),
  ('S3','PASS'),('S3','SCREENED_OUT'),('S3','DEFERRED'),
  ('S4','PASS'),('S4','EXCLUDED_ASSUMPTION'),('S4','SCREENED_OUT'),('S4','DEFERRED'),('S4','ERROR'),
  ('S5','PASS'),('S5','SCREENED_OUT'),('S5','DEFERRED'),('S5','ERROR'),
  ('S5','screenedOut'),('S5','undecided'),('S5','pass'),
  ('S6','rejectedProven'),('S6','rejectedAssuming'),('S6','undecided'),('S6','pass'),('S6','ERROR'),
  ('S7','EXCLUDED_ASSUMPTION'),('S7','rejectedProven'),('S7','rejectedAssuming'),('S7','undecided'),
  ('S7','pass'),('S7','ERROR'),
  ('S8','rejectedProven'),('S8','rejectedAssuming'),('S8','undecided'),('S8','pass'),('S8','ERROR'))
  WHERE NOT EXISTS (SELECT 1 FROM stage_outcomes);

CREATE TABLE IF NOT EXISTS outcomes (
  id INTEGER PRIMARY KEY,
  candidate_id TEXT NOT NULL REFERENCES candidates, stage TEXT NOT NULL,
  outcome TEXT NOT NULL, criterion TEXT, value TEXT, threshold TEXT, source_id TEXT,
  provisional INTEGER NOT NULL DEFAULT 0, tier_run_id INTEGER REFERENCES tier_runs,
  targets_sha256 TEXT NOT NULL, tiers_sha256 TEXT NOT NULL, origin TEXT NOT NULL,
  reason TEXT, decision_id INTEGER REFERENCES decisions, tags TEXT,
  created TEXT NOT NULL);

CREATE TABLE IF NOT EXISTS decide_inputs (
  id INTEGER PRIMARY KEY, candidate_id TEXT NOT NULL REFERENCES candidates,
  json TEXT NOT NULL CHECK (json_valid(json)), created TEXT NOT NULL);

CREATE TABLE IF NOT EXISTS decisions (
  id INTEGER PRIMARY KEY, candidate_id TEXT NOT NULL REFERENCES candidates, stage TEXT NOT NULL,
  tier TEXT NOT NULL CHECK (tier IN ('ml','coarseDFPT','fullDFPT','converged')),
  decide_input_id INTEGER NOT NULL REFERENCES decide_inputs,
  outcome TEXT NOT NULL CHECK (outcome IN ('rejectedProven','rejectedAssuming','screenedOut','undecided','pass')),
  firing_rule TEXT, assumptions TEXT, crosscheck_verdict_id INTEGER REFERENCES verdicts,
  origin TEXT NOT NULL, created TEXT NOT NULL);

CREATE TABLE IF NOT EXISTS verdicts (
  id INTEGER PRIMARY KEY, candidate_id TEXT NOT NULL REFERENCES candidates,
  tier_run_id INTEGER NOT NULL REFERENCES tier_runs,
  spectrum_tier TEXT NOT NULL CHECK (spectrum_tier IN ('coarse','full','converged')),
  status TEXT NOT NULL CHECK (status IN ('REJECTED_PROVEN','CERTIFIED_CANDIDATE','UNDECIDED')),
  T_low_K REAL, T_up_K REAL, T_star_K REAL NOT NULL, mu_hi REAL NOT NULL,
  rel_envelope REAL NOT NULL, rule_ids TEXT, assumptions TEXT, flags TEXT,
  certificate_artifact_id INTEGER NOT NULL REFERENCES artifacts,
  recheck_ok INTEGER NOT NULL CHECK (recheck_ok IN (0,1)),
  -- [{alt_index, spectrum_digest, source_artifact_id, sha256}], checked by verdicts_bind
  spectra_json TEXT NOT NULL CHECK (json_valid(spectra_json)),
  origin TEXT NOT NULL, created TEXT NOT NULL);

CREATE TABLE IF NOT EXISTS verdict_spectra (
  verdict_id INTEGER NOT NULL REFERENCES verdicts, alt_index INTEGER NOT NULL,
  spectrum_digest TEXT NOT NULL, source_artifact_id INTEGER NOT NULL REFERENCES artifacts);

CREATE TABLE IF NOT EXISTS labels (
  candidate_id TEXT NOT NULL REFERENCES candidates, name TEXT NOT NULL, value REAL NOT NULL,
  stage TEXT NOT NULL, tier_run_id INTEGER REFERENCES tier_runs, created TEXT);

CREATE TABLE IF NOT EXISTS calibration (
  id INTEGER PRIMARY KEY, material TEXT, coarse_run INTEGER, ref_run INTEGER,
  sigma INTEGER, A_coarse REAL, A_ref REAL, note TEXT);

CREATE TABLE IF NOT EXISTS audits (
  candidate_id TEXT NOT NULL REFERENCES candidates, audit TEXT NOT NULL,
  inclusion_prob REAL NOT NULL CHECK (inclusion_prob > 0 AND inclusion_prob <= 1), created TEXT NOT NULL);

CREATE TABLE IF NOT EXISTS screen_status (   -- R8: control-recall auto-disable log (insert-only)
  id INTEGER PRIMARY KEY, screen TEXT NOT NULL, enabled INTEGER NOT NULL, reason TEXT NOT NULL,
  actor TEXT NOT NULL, created TEXT NOT NULL);

-- B9 (doc 31 §1): one row per closed rank pool (insert-only): the whole ranking, the cut, the
-- per-lane minimum shares and the score version, so every rank-cut decision is reproducible
CREATE TABLE IF NOT EXISTS rank_pools (
  id INTEGER PRIMARY KEY, pool_id TEXT NOT NULL, stage TEXT NOT NULL,
  cut INTEGER NOT NULL, free_slots INTEGER NOT NULL,
  json TEXT NOT NULL CHECK (json_valid(json)), ranking_sha256 TEXT NOT NULL,
  tiers_sha256 TEXT, created TEXT NOT NULL);
CREATE TRIGGER IF NOT EXISTS rank_pools_no_update BEFORE UPDATE ON rank_pools BEGIN SELECT RAISE(ABORT, 'append-only: rank_pools'); END;
CREATE TRIGGER IF NOT EXISTS rank_pools_no_delete BEFORE DELETE ON rank_pools BEGIN SELECT RAISE(ABORT, 'append-only: rank_pools'); END;

CREATE TABLE IF NOT EXISTS queue (candidate_id TEXT, stage TEXT, priority REAL, state TEXT,
  PRIMARY KEY (candidate_id, stage));   -- derived, mutable, rebuildable

CREATE INDEX IF NOT EXISTS ix_outcomes_cs ON outcomes(candidate_id, stage, id);
CREATE INDEX IF NOT EXISTS ix_tier_runs_c ON tier_runs(candidate_id);
CREATE INDEX IF NOT EXISTS ix_artifacts_tr ON artifacts(tier_run_id);
CREATE INDEX IF NOT EXISTS ix_decisions_c ON decisions(candidate_id, stage);
CREATE INDEX IF NOT EXISTS ix_verdicts_c ON verdicts(candidate_id);

-- ============================================================== append-only triggers
CREATE TRIGGER IF NOT EXISTS meta_no_update BEFORE UPDATE ON meta BEGIN SELECT RAISE(ABORT, 'append-only: meta'); END;
CREATE TRIGGER IF NOT EXISTS meta_no_delete BEFORE DELETE ON meta BEGIN SELECT RAISE(ABORT, 'append-only: meta'); END;
CREATE TRIGGER IF NOT EXISTS inputs_no_update BEFORE UPDATE ON inputs BEGIN SELECT RAISE(ABORT, 'append-only: inputs'); END;
CREATE TRIGGER IF NOT EXISTS inputs_no_delete BEFORE DELETE ON inputs BEGIN SELECT RAISE(ABORT, 'append-only: inputs'); END;
CREATE TRIGGER IF NOT EXISTS structures_no_update BEFORE UPDATE ON structures BEGIN SELECT RAISE(ABORT, 'append-only: structures'); END;
CREATE TRIGGER IF NOT EXISTS structures_no_delete BEFORE DELETE ON structures BEGIN SELECT RAISE(ABORT, 'append-only: structures'); END;
CREATE TRIGGER IF NOT EXISTS candidates_no_update BEFORE UPDATE ON candidates BEGIN SELECT RAISE(ABORT, 'append-only: candidates'); END;
CREATE TRIGGER IF NOT EXISTS candidates_no_delete BEFORE DELETE ON candidates BEGIN SELECT RAISE(ABORT, 'append-only: candidates'); END;
CREATE TRIGGER IF NOT EXISTS duplicates_no_update BEFORE UPDATE ON duplicates BEGIN SELECT RAISE(ABORT, 'append-only: duplicates'); END;
CREATE TRIGGER IF NOT EXISTS duplicates_no_delete BEFORE DELETE ON duplicates BEGIN SELECT RAISE(ABORT, 'append-only: duplicates'); END;
CREATE TRIGGER IF NOT EXISTS tier_runs_no_delete BEFORE DELETE ON tier_runs BEGIN SELECT RAISE(ABORT, 'append-only: tier_runs'); END;
CREATE TRIGGER IF NOT EXISTS artifacts_no_update BEFORE UPDATE ON artifacts BEGIN SELECT RAISE(ABORT, 'append-only: artifacts'); END;
CREATE TRIGGER IF NOT EXISTS artifacts_no_delete BEFORE DELETE ON artifacts BEGIN SELECT RAISE(ABORT, 'append-only: artifacts'); END;
CREATE TRIGGER IF NOT EXISTS artifact_events_no_update BEFORE UPDATE ON artifact_events BEGIN SELECT RAISE(ABORT, 'append-only: artifact_events'); END;
CREATE TRIGGER IF NOT EXISTS artifact_events_no_delete BEFORE DELETE ON artifact_events BEGIN SELECT RAISE(ABORT, 'append-only: artifact_events'); END;
CREATE TRIGGER IF NOT EXISTS stage_outcomes_no_insert BEFORE INSERT ON stage_outcomes
  WHEN (SELECT COUNT(*) FROM stage_outcomes) >= 37 BEGIN SELECT RAISE(ABORT, 'stage_outcomes is fixed (§3.1)'); END;
CREATE TRIGGER IF NOT EXISTS stage_outcomes_no_update BEFORE UPDATE ON stage_outcomes BEGIN SELECT RAISE(ABORT, 'stage_outcomes is fixed (§3.1)'); END;
CREATE TRIGGER IF NOT EXISTS stage_outcomes_no_delete BEFORE DELETE ON stage_outcomes BEGIN SELECT RAISE(ABORT, 'stage_outcomes is fixed (§3.1)'); END;
CREATE TRIGGER IF NOT EXISTS outcomes_no_update BEFORE UPDATE ON outcomes BEGIN SELECT RAISE(ABORT, 'append-only: outcomes'); END;
CREATE TRIGGER IF NOT EXISTS outcomes_no_delete BEFORE DELETE ON outcomes BEGIN SELECT RAISE(ABORT, 'append-only: outcomes'); END;
CREATE TRIGGER IF NOT EXISTS decide_inputs_no_update BEFORE UPDATE ON decide_inputs BEGIN SELECT RAISE(ABORT, 'append-only: decide_inputs'); END;
CREATE TRIGGER IF NOT EXISTS decide_inputs_no_delete BEFORE DELETE ON decide_inputs BEGIN SELECT RAISE(ABORT, 'append-only: decide_inputs'); END;
CREATE TRIGGER IF NOT EXISTS decisions_no_update BEFORE UPDATE ON decisions BEGIN SELECT RAISE(ABORT, 'append-only: decisions'); END;
CREATE TRIGGER IF NOT EXISTS decisions_no_delete BEFORE DELETE ON decisions BEGIN SELECT RAISE(ABORT, 'append-only: decisions'); END;
CREATE TRIGGER IF NOT EXISTS verdicts_no_update BEFORE UPDATE ON verdicts BEGIN SELECT RAISE(ABORT, 'append-only: verdicts'); END;
CREATE TRIGGER IF NOT EXISTS verdicts_no_delete BEFORE DELETE ON verdicts BEGIN SELECT RAISE(ABORT, 'append-only: verdicts'); END;
CREATE TRIGGER IF NOT EXISTS verdict_spectra_no_update BEFORE UPDATE ON verdict_spectra BEGIN SELECT RAISE(ABORT, 'append-only: verdict_spectra'); END;
CREATE TRIGGER IF NOT EXISTS verdict_spectra_no_delete BEFORE DELETE ON verdict_spectra BEGIN SELECT RAISE(ABORT, 'append-only: verdict_spectra'); END;
CREATE TRIGGER IF NOT EXISTS labels_no_update BEFORE UPDATE ON labels BEGIN SELECT RAISE(ABORT, 'append-only: labels'); END;
CREATE TRIGGER IF NOT EXISTS labels_no_delete BEFORE DELETE ON labels BEGIN SELECT RAISE(ABORT, 'append-only: labels'); END;
CREATE TRIGGER IF NOT EXISTS calibration_no_update BEFORE UPDATE ON calibration BEGIN SELECT RAISE(ABORT, 'append-only: calibration'); END;
CREATE TRIGGER IF NOT EXISTS calibration_no_delete BEFORE DELETE ON calibration BEGIN SELECT RAISE(ABORT, 'append-only: calibration'); END;
CREATE TRIGGER IF NOT EXISTS audits_no_update BEFORE UPDATE ON audits BEGIN SELECT RAISE(ABORT, 'append-only: audits'); END;
CREATE TRIGGER IF NOT EXISTS audits_no_delete BEFORE DELETE ON audits BEGIN SELECT RAISE(ABORT, 'append-only: audits'); END;
CREATE TRIGGER IF NOT EXISTS screen_status_no_update BEFORE UPDATE ON screen_status BEGIN SELECT RAISE(ABORT, 'append-only: screen_status'); END;
CREATE TRIGGER IF NOT EXISTS screen_status_no_delete BEFORE DELETE ON screen_status BEGIN SELECT RAISE(ABORT, 'append-only: screen_status'); END;

-- ============================================================== tier_runs life cycle
CREATE TRIGGER IF NOT EXISTS tier_runs_insert_started BEFORE INSERT ON tier_runs
  WHEN NEW.state <> 'started'
BEGIN SELECT RAISE(ABORT, 'tier_runs: a run is inserted as started (written before launch)'); END;

CREATE TRIGGER IF NOT EXISTS tier_runs_transition BEFORE UPDATE ON tier_runs
  WHEN OLD.state <> 'started' OR NEW.state NOT IN ('done','failed')
    OR NEW.id IS NOT OLD.id OR NEW.candidate_id IS NOT OLD.candidate_id OR NEW.stage IS NOT OLD.stage
    OR NEW.kind IS NOT OLD.kind OR NEW.profile IS NOT OLD.profile OR NEW.attempt IS NOT OLD.attempt
    OR NEW.params_json IS NOT OLD.params_json OR NEW.versions_json IS NOT OLD.versions_json
    OR NEW.lane IS NOT OLD.lane OR NEW.started IS NOT OLD.started
BEGIN SELECT RAISE(ABORT, 'tier_runs: only started -> done|failed with completion fields'); END;

-- artifacts may only be attached to a run that is not failed
CREATE TRIGGER IF NOT EXISTS artifacts_bind BEFORE INSERT ON artifacts
  WHEN (SELECT state FROM tier_runs WHERE id = NEW.tier_run_id) IS NOT 'started'
   AND (SELECT state FROM tier_runs WHERE id = NEW.tier_run_id) IS NOT 'done'
BEGIN SELECT RAISE(ABORT, 'artifacts: tier_run missing or failed'); END;

-- ============================================================== outcomes_allowed
CREATE TRIGGER IF NOT EXISTS outcomes_allowed BEFORE INSERT ON outcomes
  WHEN NOT EXISTS (SELECT 1 FROM stage_outcomes s WHERE s.stage = NEW.stage AND s.outcome = NEW.outcome)
BEGIN SELECT RAISE(ABORT, 'outcomes_allowed: (stage, outcome) not in the §3.1 table'); END;

CREATE TRIGGER IF NOT EXISTS outcomes_inputs_known BEFORE INSERT ON outcomes
  WHEN NOT EXISTS (SELECT 1 FROM inputs WHERE sha256 = NEW.targets_sha256 AND kind = 'targets')
    OR NOT EXISTS (SELECT 1 FROM inputs WHERE sha256 = NEW.tiers_sha256 AND kind = 'tiers')
BEGIN SELECT RAISE(ABORT, 'outcomes: targets/tiers sha256 not registered in inputs (R7)'); END;

-- R2/R3/R7: removals carry criterion and source id; DEFERRED and ERROR carry a reason
CREATE TRIGGER IF NOT EXISTS outcomes_criterion BEFORE INSERT ON outcomes
  WHEN (NEW.outcome IN ('SCREENED_OUT','EXCLUDED_ASSUMPTION')
        AND (NEW.criterion IS NULL OR NEW.source_id IS NULL OR NEW.value IS NULL OR NEW.threshold IS NULL))
    OR (NEW.outcome IN ('DEFERRED','ERROR') AND NEW.reason IS NULL AND NEW.criterion IS NULL)
BEGIN SELECT RAISE(ABORT, 'outcomes: removal without criterion/value/threshold/source_id, or DEFERRED/ERROR without reason'); END;

-- R3: EXCLUDED_ASSUMPTION only from S4 (H1, H2) and S7 (H6, H7)
CREATE TRIGGER IF NOT EXISTS outcomes_excluded_criteria BEFORE INSERT ON outcomes
  WHEN NEW.outcome = 'EXCLUDED_ASSUMPTION' AND NOT (
       (NEW.stage = 'S4' AND (NEW.criterion LIKE 'H1%' OR NEW.criterion LIKE 'H2%'))
    OR (NEW.stage = 'S7' AND (NEW.criterion LIKE 'H6%' OR NEW.criterion LIKE 'H7%')))
BEGIN SELECT RAISE(ABORT, 'outcomes: EXCLUDED_ASSUMPTION criterion not allowed at this stage (R3)'); END;

-- decide outcomes in `outcomes` must point to a matching decisions row (same candidate, stage);
-- the outcome equals the decision, or is 'undecided' (combination rule downgrade, §3.2)
CREATE TRIGGER IF NOT EXISTS outcomes_decide_bound BEFORE INSERT ON outcomes
  WHEN NEW.outcome IN ('rejectedProven','rejectedAssuming','screenedOut','undecided','pass')
   AND NOT EXISTS (SELECT 1 FROM decisions d WHERE d.id = NEW.decision_id
                   AND d.candidate_id = NEW.candidate_id AND d.stage = NEW.stage
                   AND (d.outcome = NEW.outcome OR NEW.outcome = 'undecided'))
BEGIN SELECT RAISE(ABORT, 'outcomes: decide outcome without a matching decisions row'); END;

-- removal by decide outcome also needs the tier_run (R4) to be the bound DFPT/EPW run
CREATE TRIGGER IF NOT EXISTS outcomes_rejection_bind BEFORE INSERT ON outcomes
  WHEN NEW.outcome IN ('rejectedProven','rejectedAssuming')
   AND NOT EXISTS (SELECT 1 FROM tier_runs t WHERE t.id = NEW.tier_run_id
                   AND t.candidate_id = NEW.candidate_id AND t.state = 'done'
                   AND ((t.kind = 'dfpt' AND t.profile = 'full') OR (t.kind = 'epw' AND t.profile = 'epw')))
BEGIN SELECT RAISE(ABORT, 'outcomes: rejection without a done full-DFPT/EPW tier_run of this candidate (R4)'); END;

-- S5 coarse removal is SCREENED_OUT(coarse-A) only
CREATE TRIGGER IF NOT EXISTS outcomes_s5_screen BEFORE INSERT ON outcomes
  WHEN NEW.stage = 'S5' AND NEW.outcome = 'SCREENED_OUT' AND NEW.criterion NOT LIKE 'coarse-A%'
BEGIN SELECT RAISE(ABORT, 'outcomes: S5 removes only by SCREENED_OUT(coarse-A)'); END;

-- ============================================================== verdicts_bind (R4)
CREATE TRIGGER IF NOT EXISTS verdicts_bind_origin BEFORE INSERT ON verdicts
  WHEN NEW.origin <> 'rtt.core.certify_policy.certify_spectra'
BEGIN SELECT RAISE(ABORT, 'verdicts_bind: origin must be the core certify_spectra (R1b)'); END;

CREATE TRIGGER IF NOT EXISTS verdicts_bind_recheck BEFORE INSERT ON verdicts
  WHEN NEW.status IN ('REJECTED_PROVEN','CERTIFIED_CANDIDATE') AND NEW.recheck_ok <> 1
BEGIN SELECT RAISE(ABORT, 'verdicts_bind: proven status without recheck_ok = 1'); END;

CREATE TRIGGER IF NOT EXISTS verdicts_bind_run BEFORE INSERT ON verdicts
  WHEN NOT EXISTS (SELECT 1 FROM tier_runs t WHERE t.id = NEW.tier_run_id
                   AND t.candidate_id = NEW.candidate_id AND t.state = 'done'
                   AND ((t.kind = 'dfpt' AND t.profile = 'coarse' AND NEW.spectrum_tier = 'coarse')
                     OR (t.kind = 'dfpt' AND t.profile = 'full' AND NEW.spectrum_tier = 'full')
                     OR (t.kind = 'epw' AND t.profile = 'epw' AND NEW.spectrum_tier = 'converged')))
BEGIN SELECT RAISE(ABORT, 'verdicts_bind: tier_run not a done dfpt/epw run of this candidate with matching profile'); END;

CREATE TRIGGER IF NOT EXISTS verdicts_bind_certificate BEFORE INSERT ON verdicts
  WHEN NOT EXISTS (SELECT 1 FROM artifacts a WHERE a.id = NEW.certificate_artifact_id
                   AND a.tier_run_id = NEW.tier_run_id)
BEGIN SELECT RAISE(ABORT, 'verdicts_bind: certificate artifact not of the same tier_run'); END;

-- every alternative names a source artifact of a done dfpt/epw run of this candidate whose
-- registered sha256 equals the sha256 the certify step hashed; at least one alternative
CREATE TRIGGER IF NOT EXISTS verdicts_bind_spectra BEFORE INSERT ON verdicts
  WHEN json_type(NEW.spectra_json) IS NOT 'array' OR json_array_length(NEW.spectra_json) < 1
    OR EXISTS (SELECT 1 FROM json_each(NEW.spectra_json) j
               WHERE json_extract(j.value, '$.spectrum_digest') IS NULL
                  OR json_extract(j.value, '$.alt_index') IS NULL
                  OR NOT EXISTS (SELECT 1 FROM artifacts a JOIN tier_runs t ON t.id = a.tier_run_id
                                 WHERE a.id = json_extract(j.value, '$.source_artifact_id')
                                   AND a.sha256 = json_extract(j.value, '$.sha256')
                                   AND t.candidate_id = NEW.candidate_id AND t.state = 'done'
                                   AND t.kind IN ('dfpt','epw')))
BEGIN SELECT RAISE(ABORT, 'verdicts_bind: spectra without hashed DFPT source artifacts (R4)'); END;

CREATE TRIGGER IF NOT EXISTS verdict_spectra_bind BEFORE INSERT ON verdict_spectra
  WHEN NOT EXISTS (SELECT 1 FROM verdicts v, json_each(v.spectra_json) j
                   WHERE v.id = NEW.verdict_id
                     AND json_extract(j.value, '$.alt_index') = NEW.alt_index
                     AND json_extract(j.value, '$.spectrum_digest') = NEW.spectrum_digest
                     AND json_extract(j.value, '$.source_artifact_id') = NEW.source_artifact_id)
BEGIN SELECT RAISE(ABORT, 'verdict_spectra_bind: row does not match the verdict''s bound spectra'); END;

-- ============================================================== decisions_guard (R1c, R4)
CREATE TRIGGER IF NOT EXISTS decisions_guard_origin BEFORE INSERT ON decisions
  WHEN NEW.origin <> 'rtt.decide.core.decide_t'
BEGIN SELECT RAISE(ABORT, 'decisions_guard: origin must be rtt.decide.core.decide_t (R1c)'); END;

CREATE TRIGGER IF NOT EXISTS decisions_guard_input BEFORE INSERT ON decisions
  WHEN NOT EXISTS (SELECT 1 FROM decide_inputs di WHERE di.id = NEW.decide_input_id
                   AND di.candidate_id = NEW.candidate_id
                   AND json_extract(di.json, '$.tier') = NEW.tier)
BEGIN SELECT RAISE(ABORT, 'decisions_guard: decide input missing, of another candidate, or of another tier'); END;

-- Lean decideT v2: screenedOut = 'a floor holds on lower-tier inputs' and can occur at any stage
-- tier (mixed provenance at S6); it never removes at S6-S8 (stage outcome becomes 'undecided').
-- rejectedAssuming is not produced by decideT v2 (NM-G gives rejectedProven with NUM-* assumptions).
CREATE TRIGGER IF NOT EXISTS decisions_guard_no_assuming BEFORE INSERT ON decisions
  WHEN NEW.outcome = 'rejectedAssuming'
BEGIN SELECT RAISE(ABORT, 'decisions_guard: rejectedAssuming is not an outcome of decideT v2'); END;

CREATE TRIGGER IF NOT EXISTS decisions_guard_reject_tier BEFORE INSERT ON decisions
  WHEN NEW.outcome IN ('rejectedProven','rejectedAssuming') AND NEW.tier NOT IN ('fullDFPT','converged')
BEGIN SELECT RAISE(ABORT, 'decisions_guard: rejection requires tier fullDFPT or converged'); END;

-- the decide input's tier must equal the tier of the DFPT runs it references: for a
-- proof-grade decision at least one referenced tier_run is a done run of this candidate with the
-- matching profile, and no referenced tier_run is of another candidate
CREATE TRIGGER IF NOT EXISTS decisions_guard_provenance BEFORE INSERT ON decisions
  WHEN EXISTS (SELECT 1 FROM decide_inputs di, json_each(di.json, '$.provenance.tier_run_ids') j
               WHERE di.id = NEW.decide_input_id
                 AND NOT EXISTS (SELECT 1 FROM tier_runs t WHERE t.id = j.value AND t.candidate_id = NEW.candidate_id))
    OR (NEW.tier IN ('fullDFPT','converged') AND NOT EXISTS (
          SELECT 1 FROM decide_inputs di, json_each(di.json, '$.provenance.tier_run_ids') j, tier_runs t
          WHERE di.id = NEW.decide_input_id AND t.id = j.value AND t.state = 'done'
            AND ((NEW.tier = 'fullDFPT' AND t.kind = 'dfpt' AND t.profile = 'full')
              OR (NEW.tier = 'converged' AND t.kind = 'epw' AND t.profile = 'epw'))))
    OR (NEW.tier = 'coarseDFPT' AND NOT EXISTS (
          SELECT 1 FROM decide_inputs di, json_each(di.json, '$.provenance.tier_run_ids') j, tier_runs t
          WHERE di.id = NEW.decide_input_id AND t.id = j.value AND t.kind = 'dfpt' AND t.profile = 'coarse'))
BEGIN SELECT RAISE(ABORT, 'decisions_guard: decide-input tier does not match its referenced tier_runs'); END;

-- removal decisions (R4): a rechecked core verdict bound to hashed spectra of a referenced
-- proof-grade run must exist; N1-ISO additionally needs a REJECTED_PROVEN crosscheck verdict
CREATE TRIGGER IF NOT EXISTS decisions_guard_evidence BEFORE INSERT ON decisions
  WHEN NEW.outcome IN ('rejectedProven','rejectedAssuming') AND NOT EXISTS (
          SELECT 1 FROM verdicts v, decide_inputs di, json_each(di.json, '$.provenance.tier_run_ids') j
          WHERE di.id = NEW.decide_input_id AND v.tier_run_id = j.value
            AND v.candidate_id = NEW.candidate_id AND v.recheck_ok = 1
            AND v.spectrum_tier = CASE NEW.tier WHEN 'fullDFPT' THEN 'full' ELSE 'converged' END)
BEGIN SELECT RAISE(ABORT, 'decisions_guard: rejection without a rechecked verdict on hashed DFPT spectra (R4)'); END;

CREATE TRIGGER IF NOT EXISTS decisions_guard_crosscheck BEFORE INSERT ON decisions
  WHEN NEW.outcome = 'rejectedProven' AND NEW.firing_rule = 'N1-ISO' AND NOT EXISTS (
          SELECT 1 FROM verdicts v WHERE v.id = NEW.crosscheck_verdict_id
            AND v.candidate_id = NEW.candidate_id AND v.status = 'REJECTED_PROVEN'
            AND v.recheck_ok = 1 AND v.spectrum_tier IN ('full','converged'))
BEGIN SELECT RAISE(ABORT, 'decisions_guard: N1-ISO rejection needs a REJECTED_PROVEN crosscheck verdict'); END;

CREATE TRIGGER IF NOT EXISTS decisions_guard_rule BEFORE INSERT ON decisions
  WHEN (NEW.outcome = 'rejectedProven' AND (NEW.firing_rule IS NULL
          OR NEW.firing_rule NOT IN ('N1-KAPPA','N1-ROWMAX','N1-DOS','N1-ISO','NM-CERT','NM-G')
          OR NEW.assumptions IS NULL OR NEW.assumptions NOT LIKE '%ME-MODEL%'))
BEGIN SELECT RAISE(ABORT, 'decisions_guard: rejection without a matching firing rule'); END;
