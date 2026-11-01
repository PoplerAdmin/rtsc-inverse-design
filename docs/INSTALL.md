# Installing and running rtsc-targeted on your own machine

This page is the how-to. For what each part does and why it was designed that way, with block
diagrams of the flow, see [RUNNER.md](RUNNER.md).

## Requirements

| | Minimum (reference machine) | Notes |
|---|---|---|
| OS | Linux x86-64, or Windows 11 with WSL2 Ubuntu | macOS is not supported (QE + CUDA) |
| CPU | 8 physical cores | parallelism scales up automatically on larger machines |
| RAM | 16 GB (15 GB visible to Linux) | ~1.5 GB per QE rank |
| GPU | NVIDIA, ≥ 12 GB (reference: RTX 4090, 24 GB) | without one, MatterGen lanes are off and MatterSim runs on CPU |
| Disk | 200 GB free under `~/.rtt` | DFPT scratch is large; the installer needs ~25 GB |
| Tools | `bash curl git tar gzip` | no root, nothing installed system-wide |

On Windows: install WSL2 (`wsl --install -d Ubuntu`), the NVIDIA Windows driver (CUDA works
inside WSL2 without a Linux driver), then clone this repository **inside** WSL (for example
`~/rtsc-targeted`) and run everything from a WSL terminal. A clone on a Windows drive
(`/mnt/c/...`) also works, but is slower.

## Install (once, ~20–40 min, mostly downloads)

```bash
git clone <this repository> rtsc-targeted
cd rtsc-targeted
./install.sh                    # add --with-reference for lane D and the DFT-reference hull (+5 GB)
```

Everything goes into `~/.rtt` (override with `RTT_HOME=/path ./install.sh`):

- Quantum ESPRESSO 7.6 (conda-forge, OpenMPI) via a private micromamba;
- a Python 3.12 environment with `rtt`, MatterSim 1.2.5, the dashboard (`textual`);
- MatterGen at commit 9242366 in its own Python 3.10 environment, plus the seven shipped
  checkpoints the lanes use (Hugging Face `microsoft/mattergen`, ~0.5 GB each);
- the SG15 ONCV pseudopotentials;
- `~/.rtt/env.sh`, the environment every `rtt` process uses;
- `~/.rtt/machine.toml`, written by `rtt init` from the detected cores, RAM and GPU.

Versions are pinned in `tools/versions.env`. Each download's sha256 goes to
`~/.rtt/install_manifest.txt`, so two installs can be compared. The installer ends with
`rtt doctor`, which lists every component as `ok`, `FAIL` (required) or `off` (optional; it only
disables a lane).

## Check the install (about an hour)

```bash
source ~/.rtt/env.sh
rtt campaign start --smoke --run-id smoke
./start.sh --run-id smoke          # watch it; it ends by itself
```

The smoke run takes the Al control through S1–S5: MatterSim, a DFT electronic screen and a
coarse DFPT run on 4 ranks. It should end with S4 and S5 PASS for Al. Delete `~/.rtt/runs/smoke`
afterwards.

## Run

```bash
./start.sh
```

This opens the dashboard. Keys:

| Key | Action |
|---|---|
| `s` | start the campaign (a background controller; closing the dashboard does not stop it) |
| `d` | drain: running jobs finish, nothing new starts |
| `x` | stop now: every job is terminated (asks first); interrupted jobs are retried on the next start |
| `l` | switch the log pane between the event log and each running job's log |
| `q` | quit the dashboard (if a campaign is running it asks: keep running, drain, or stop) |

The dashboard shows the machine load, the stage funnel (current outcome per candidate), the
running jobs with elapsed vs estimated time, the candidates that got furthest, and jobs that were
skipped as **compute-limited** (estimated longer than `max_job_wall_h`, 48 h by default).

Without the dashboard:

```bash
source ~/.rtt/env.sh
rtt campaign start [--rounds N | --rounds 0] [--run-id NAME] [--no-front]
rtt campaign status
./stop.sh            # stop now      (./stop.sh --drain to drain)
```

`--rounds 0` keeps generating and verifying until stopped. `--no-front` skips generation and only
continues verifying what is already in the ledger. Each run lives in `~/.rtt/runs/<run-id>/`:
`ledger.sqlite` (every decision, append-only), `work/` (QE/MatterGen scratch), `events.log`,
`jobs/*.log`, `status.json`.

## What one round does

1. **Front end** (GPU): MatterGen lanes G1–G5 for `gpu_h_per_round` GPU hours, the prototype lane P
   and controls (round 1), the database lane D (with `--with-reference`), then S2 (validity,
   dedupe) and S3 (MatterSim relax, hull, ranking; ML only ranks or defers, never rejects).
2. **Back end** (CPU): S4 (DFT electronic screen), S5 (coarse DFPT screen), S6 (full DFPT +
   certified Eliashberg bounds + decision), as many in parallel as `max_ranks` allows.

## Sizing

`rtt init` sizes the parallelism from the hardware: one physical core is kept free; QE gets
`max_ranks = cores − 1` (fewer if RAM is short); S4/S5 jobs use 4 ranks each (the setting all
measured costs refer to), S6 up to 8. Edit `~/.rtt/machine.toml` to change any of this
(`rtt init --force` rewrites it). On the reference machine one 9-atom hydride S5 takes about
1 day on 4 ranks and an S6 1.5–9 days, so useful campaigns run for days to weeks; machines with
more cores run more candidates at once, not each one faster.

## Known limits

- The core mode-positive admission tolerance (1e-6) rejected every smearing of the first
  hydride S6 (NaIrH₃; doc 33): S6 on hydrides currently ends in ERROR until it is recalibrated.
- MatterGen needs an NVIDIA GPU; on CPU-only machines only lanes P, D and C feed the pipeline.
- Measured costs come from one machine (docs/design/30, 31, 33).
