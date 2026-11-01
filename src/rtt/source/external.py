"""Locations of external tools and data used by the front end (S1-S3), and the subprocess
bridge into the MatterGen venv.

Defaults follow docs/design/10_ml_capabilities.md §1 and §3.3. Each can be overridden by an
environment variable or by a ``[paths]`` table in ``machine.toml``:

=====================  ==========================  ==============================================
key                    env var                     default
=====================  ==========================  ==============================================
mattergen_venv         RTT_MATTERGEN_VENV          ~/mattergen/.venv
alexmp_lmdb            RTT_ALEXMP_LMDB             ~/rtsc_targeted_work/T3-mlaudit/ref/reference_MP2020correction.lmdb
hf_hub                 RTT_HF_HUB                  ~/.cache/huggingface/hub
pseudo_dir             RTT_PSEUDO_DIR              tiers.toml [S1] pseudo_dir, else ~/pseudo/sg15
front_cache            RTT_FRONT_CACHE             ~/rtt_scratch/front_cache
=====================  ==========================  ==============================================
"""
from __future__ import annotations

import json
import os
import subprocess
import tempfile
import time
from pathlib import Path
from typing import Any

DEFAULTS = {
    "mattergen_venv": ("RTT_MATTERGEN_VENV", "~/mattergen/.venv"),
    "alexmp_lmdb": ("RTT_ALEXMP_LMDB", "~/rtsc_targeted_work/T3-mlaudit/ref/reference_MP2020correction.lmdb"),
    "hf_hub": ("RTT_HF_HUB", "~/.cache/huggingface/hub"),
    "pseudo_dir": ("RTT_PSEUDO_DIR", "~/pseudo/sg15"),
    "front_cache": ("RTT_FRONT_CACHE", "~/rtt_scratch/front_cache"),
}

TOOL = Path(__file__).resolve().parent / "mgvenv_tool.py"


class ExternalToolError(RuntimeError):
    pass


def path_of(key: str, machine: Any = None) -> Path:
    env, default = DEFAULTS[key]
    v = os.environ.get(env)
    if not v and machine is not None:
        try:
            v = machine.get(f"paths.{key}", None)
        except Exception:  # noqa: BLE001 - machine may be a plain dict
            v = (machine.get("paths", {}) or {}).get(key) if isinstance(machine, dict) else None
    return Path(os.path.expanduser(v or default))


def child_env(extra: dict[str, str] | None = None) -> dict[str, str]:
    """CONTEXT.md compute rules: at most 2 threads for every child."""
    env = dict(os.environ)
    env.update({"OMP_NUM_THREADS": "2", "OPENBLAS_NUM_THREADS": "2", "MKL_NUM_THREADS": "2",
                "TORCH_DISABLE_NATIVE_JIT": "1"})
    env.update(extra or {})
    return env


def gpu_used_mib() -> int | None:
    """GPU memory in use (MiB) from nvidia-smi, or None if unavailable."""
    try:
        cp = subprocess.run(["nvidia-smi", "--query-gpu=memory.used", "--format=csv,noheader,nounits"],
                            capture_output=True, text=True, timeout=30)
        return int(cp.stdout.strip().splitlines()[0])
    except Exception:  # noqa: BLE001
        return None


def gpu_free_for_us(max_used_mib: int = 4000) -> bool:
    """Shared-GPU etiquette: only start GPU work when other jobs use < ``max_used_mib``."""
    u = gpu_used_mib()
    return u is not None and u < max_used_mib


def run_mgvenv_tool(args: list[str], *, machine: Any = None, timeout_s: float = 3600) -> None:
    """Run ``mgvenv_tool.py`` with the MatterGen venv's python (isolated mode, nice 19)."""
    py = path_of("mattergen_venv", machine) / "bin" / "python"
    if not py.exists():
        raise ExternalToolError(f"MatterGen venv python not found: {py}")
    cmd = ["nice", "-n", "19", str(py), "-I", str(TOOL), *args]
    t0 = time.time()
    cp = subprocess.run(cmd, env=child_env(), capture_output=True, text=True, timeout=timeout_s,
                        cwd=tempfile.gettempdir())
    if cp.returncode != 0:
        raise ExternalToolError(f"mgvenv_tool {args[0]} failed (rc={cp.returncode}, "
                                f"{time.time() - t0:.0f} s): {(cp.stderr or cp.stdout)[-2000:]}")


def refquery(chemsys: str, *, machine: Any = None, cache_dir: Path | None = None) -> dict:
    """All Alex-MP reference entries of ``chemsys`` and its subsystems (cached as JSON)."""
    cs = "-".join(sorted(set(chemsys.split("-"))))
    cache_dir = Path(cache_dir or path_of("front_cache", machine)) / "alexmp_ref"
    cache_dir.mkdir(parents=True, exist_ok=True)
    f = cache_dir / f"ref_{cs}.json"
    if not f.exists():
        tmp = f.with_suffix(".tmp")
        run_mgvenv_tool(["refquery", "--lmdb", str(path_of("alexmp_lmdb", machine)),
                         "--chemsys", cs, "--out", str(tmp)], machine=machine)
        tmp.replace(f)
    return json.loads(f.read_text())
