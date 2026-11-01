# rtsc-targeted

Theory-targeted inverse design of conventional (phonon-mediated) superconductors aimed at
Tc ≥ 300 K at 1 atm.

The pipeline first derives from theory which properties a room-temperature, ambient-pressure
superconductor must have, then conditions generation (MatterGen) and screening (MatterSim) on
those properties, and only then spends DFPT time on the survivors. It ends in a certified
Eliashberg Tc bound (the frozen engine in `src/rtt/core/`).

## Quick start

On Linux or WSL2 with an NVIDIA GPU (requirements and details: [docs/INSTALL.md](docs/INSTALL.md)):

```bash
./install.sh      # once: QE, Python env, MatterGen + checkpoints, pseudopotentials -> ~/.rtt
./start.sh        # dashboard: s = start, d = drain, x = stop now, q = quit
```

How it works, with flow diagrams and the reasons behind each design choice:
[docs/RUNNER.md](docs/RUNNER.md).

## Layout

Work in progress.

- `docs/theory/`: the theoretical framework (necessary conditions, target property set)
- `docs/reviews/`: adversarial reviews of the theory and design, with responses
- `docs/design/`: pipeline architecture
- `src/`, `tests/`: code
- `runs/`: run records (large outputs are not committed; SHA-256 manifests are)
