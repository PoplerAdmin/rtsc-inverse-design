"""Hardware detection and a ``machine.toml`` sized to this machine (``rtt init``).

The reference machine (where every cost in docs/design was measured) has 8 physical cores /
16 threads, 15 GB RAM visible to Linux and one 24 GB NVIDIA GPU. ``minimum_report`` compares
the detected machine against it; the sizing rules below scale the parallelism up on larger
machines and never past the reference per-job settings (ranks per job), so per-job costs stay
comparable with the measured ones.
"""
from __future__ import annotations

import os
import shutil
import subprocess
from dataclasses import asdict, dataclass, field
from pathlib import Path

REFERENCE = {"physical_cores": 8, "ram_GB": 15.0, "gpu_mem_GB": 24.0}
MIN_GPU_GB_MATTERGEN = 12.0        # MatterGen batch 16 + MatterSim measured peak ~ 9 GB (doc 10)
RAM_GB_PER_RANK = 1.5              # Mg2IrH6 S5 at 80 Ry on 4 ranks peaked ~ 5 GB (doc 15/16)
S_RANKS = 4                        # per S4/S5 job: the setting the S5 cost fit was measured at
F_RANKS_MAX = 8                    # per S6 job


@dataclass
class Hardware:
    logical_cpus: int
    physical_cores: int
    ram_GB: float
    gpus: list = field(default_factory=list)      # [{"name":..., "mem_GB":...}]
    disk_free_GB: float | None = None
    os: str = ""

    @property
    def gpu_mem_GB(self) -> float:
        return max((g["mem_GB"] for g in self.gpus), default=0.0)

    def as_dict(self) -> dict:
        d = asdict(self)
        d["gpu_mem_GB"] = self.gpu_mem_GB
        return d


def _physical_cores() -> int:
    try:
        out = subprocess.run(["lscpu", "-p=core,socket"], capture_output=True, text=True, timeout=10).stdout
        pairs = {ln for ln in out.splitlines() if ln and not ln.startswith("#")}
        if pairs:
            return len(pairs)
    except (OSError, subprocess.SubprocessError):
        pass
    n = os.cpu_count() or 1
    return max(1, n // 2)


def _ram_GB() -> float:
    try:
        for ln in Path("/proc/meminfo").read_text().splitlines():
            if ln.startswith("MemTotal:"):
                return round(int(ln.split()[1]) / 1024 ** 2, 1)
    except OSError:
        pass
    return 0.0


def _gpus() -> list:
    if not shutil.which("nvidia-smi"):
        return []
    try:
        out = subprocess.run(["nvidia-smi", "--query-gpu=name,memory.total", "--format=csv,noheader,nounits"],
                             capture_output=True, text=True, timeout=30).stdout
    except (OSError, subprocess.SubprocessError):
        return []
    gpus = []
    for ln in out.splitlines():
        parts = [p.strip() for p in ln.split(",")]
        if len(parts) == 2:
            try:
                gpus.append({"name": parts[0], "mem_GB": round(float(parts[1]) / 1024, 1)})
            except ValueError:
                pass
    return gpus


def detect(path_for_disk: str | os.PathLike | None = None) -> Hardware:
    disk = None
    try:
        p = Path(path_for_disk or os.path.expanduser("~"))
        while not p.exists():
            p = p.parent
        disk = round(shutil.disk_usage(p).free / 1024 ** 3, 1)
    except OSError:
        pass
    osname = ""
    try:
        osname = Path("/proc/version").read_text().split("(")[0].strip()
    except OSError:
        osname = os.name
    return Hardware(logical_cpus=os.cpu_count() or 1, physical_cores=_physical_cores(), ram_GB=_ram_GB(),
                    gpus=_gpus(), disk_free_GB=disk, os=osname)


def sizing(hw: Hardware) -> dict:
    """Parallelism for this machine. One physical core is left for the OS, the controller and the
    UI on machines with at least 4 cores. QE runs one rank per physical core."""
    usable = hw.physical_cores - (1 if hw.physical_cores >= 4 else 0)
    by_ram = int(hw.ram_GB // RAM_GB_PER_RANK) if hw.ram_GB else usable
    max_ranks = max(1, min(usable, by_ram))
    s_ranks = min(S_RANKS, max_ranks)
    f_ranks = min(F_RANKS_MAX, max_ranks)
    gpu_ok = hw.gpu_mem_GB >= MIN_GPU_GB_MATTERGEN
    return {"max_ranks": max_ranks, "s_ranks": s_ranks, "s_jobs": max(1, max_ranks // s_ranks),
            "f_ranks": f_ranks, "gpu": gpu_ok, "mattergen_batch": 32 if hw.gpu_mem_GB >= 40 else 16}


def minimum_report(hw: Hardware) -> list[str]:
    """Warnings for a machine smaller than the reference one."""
    w = []
    if hw.physical_cores < REFERENCE["physical_cores"]:
        w.append(f"{hw.physical_cores} physical cores < reference {REFERENCE['physical_cores']}: "
                 "DFPT stages will take proportionally longer")
    if hw.ram_GB and hw.ram_GB < REFERENCE["ram_GB"] - 1:
        w.append(f"{hw.ram_GB} GB RAM < reference {REFERENCE['ram_GB']} GB: fewer QE ranks will run at once")
    if not hw.gpus:
        w.append("no NVIDIA GPU found: MatterGen generation and MatterSim screening are disabled; "
                 "only database (D), prototype (P) and control (C) lanes can feed the pipeline, and S3 "
                 "will be very slow on CPU")
    elif hw.gpu_mem_GB < MIN_GPU_GB_MATTERGEN:
        w.append(f"GPU memory {hw.gpu_mem_GB} GB < {MIN_GPU_GB_MATTERGEN} GB: MatterGen lanes disabled")
    if hw.disk_free_GB is not None and hw.disk_free_GB < 100:
        w.append(f"{hw.disk_free_GB} GB free disk < 100 GB: DFPT scratch can need tens of GB per candidate")
    return w


def machine_toml(hw: Hardware, *, scratch: str = "~/.rtt/scratch", max_job_wall_h: float = 48.0,
                 gpu_h_per_round: float = 2.0) -> str:
    """``machine.toml`` text (same schema as ``src/rtt/configs/machine.toml``, plus [budget])."""
    s = sizing(hw)
    gpu_names = ", ".join(g["name"] for g in hw.gpus) or "none"
    return f"""schema = "rtt.machine/1"
# Written by `rtt init` from the detected hardware. Edit freely; `rtt init --force` rewrites it.
# detected: {hw.physical_cores} physical cores / {hw.logical_cpus} logical, {hw.ram_GB} GB RAM, GPU: {gpu_names}

[budget]
max_ranks = {s['max_ranks']}            # total QE MPI ranks running at once (all jobs)
max_job_wall_h = {max_job_wall_h}       # jobs estimated longer than this are skipped and reported
gpu_h_per_round = {gpu_h_per_round}     # MatterGen GPU hours per generation round

[lanes.G]
gpu = {str(s['gpu']).lower()}
cpu_threads = 2
mattergen_batch = {s['mattergen_batch']}
mattergen_gpu_h = {gpu_h_per_round}

[lanes.S]
mpi_ranks = {s['s_ranks']}
jobs = {s['s_jobs']}

[lanes.F]
mpi_ranks = {s['f_ranks']}
jobs = 1

[env]
OMP_NUM_THREADS = 1
TORCH_DISABLE_NATIVE_JIT = 1
nice = 19

[disk]
scratch = "{scratch}"
min_free_GB = 80
prune_per_q = true

[derate]
contention = 0.85

[capacity]
physical_cores = {hw.physical_cores}
core_h_per_day = {s['max_ranks'] * 24}
gpu_h_per_day = {24 if s['gpu'] else 0}
"""
