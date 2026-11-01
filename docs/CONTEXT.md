# Agent and contributor brief (read first)

## Goal
Build a pipeline that
1. derives, from theory, the properties a conventional (electron–phonon) superconductor with
   Tc ≥ 300 K at ~1 atm (0–1 GPa) must have,
2. turns those properties into conditioning targets / filters for MatterGen (generation) and
   MatterSim (relaxation, stability, phonons), producing meaningful outputs,
3. ends in one or more concrete chemical compounds that are the best-supported candidates,
   with an honest statement of how strong the evidence is.

## Theory first
Every claim in `docs/theory/` must be either derived (show it), computed (give the
script), or cited (give the paper: authors, year, journal/arXiv id). Mark anything recalled from
memory and not checked as **UNVERIFIED**. Adversarial reviewers will attack every document.

## Where things are
- `docs/theory/`: theoretical framework; `docs/reviews/`: adversarial reviews; `docs/design/`:
  pipeline architecture and build notes.
- `src/rtt/core/`: the certified engine. Isotropic Migdal–Eliashberg Tc bounds with interval
  arithmetic (python-flint/Arb) and Lean-proved lemmas, e.g. the universal bound
  k_B·Tc ≤ A/4.33148 where A = ∫ α²F(ω) dω (exact definition in `rtt/core/certify`), and the
  QE 7.6 DFPT workflow (`rtt/core/dfpt/`). Statuses: REJECTED_PROVEN, CERTIFIED_CANDIDATE,
  UNDECIDED, … ML never issues verdicts. The core is frozen: sha256 of every file in
  `src/rtt/core/PINS.json`, checked by `tests/test_core_pins.py`; its tests are in `tests/core/`.
- `lean/`: Lean libraries `TcR1`, `PipelineLogic`, `MultiBand` (frozen) and `TargetLogic`.
- `src/`, `tests/`: the rest of the code; `runs/`: run records.
- Agents write files; commits are made by whoever coordinates the work.

## Compute etiquette
- CPU work: `nice -n 19`, at most 2 threads unless your brief says otherwise: export
  `OMP_NUM_THREADS=2 OPENBLAS_NUM_THREADS=2 MKL_NUM_THREADS=2`.
- Keep RAM use modest (~3 GB per agent). Do not run large test suites in parallel with other jobs.
- GPU: only the agent whose brief says it owns the GPU may use it.
- Do not stop or renice jobs you did not start.
- Anything expected to run > 15 min: launch in the background with
  `setsid nohup ... > log 2>&1 < /dev/null & sleep 5` and report the log path instead of
  waiting. The trailing `sleep` matters when the launching shell exits immediately.
- Scratch dir: `~/rtsc_targeted_work/<your-agent-name>/` (create it).

## Privacy
No personal names, emails, hostnames, account names, or home paths in anything written to the
repo. Use `~/...` or `/home/user/...` in docs.

## Output
Write your deliverable(s) to the path named in your brief. End with a short summary: key
conclusions, numbers, open problems, and what you could not verify.
