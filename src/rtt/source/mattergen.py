"""MatterGen adapter for lanes G1-G5 (architecture §3.3, §5.1 ``source/mattergen.py``).

Built on the measured script ``src/ml/mg_generate.sh``
(doc 10 §3.1). MatterGen lives in its own venv (Python 3.10, torch 2.2.1) and is never imported:
``mattergen-generate`` is run as a subprocess, under ``nice -n 19`` with 2 threads.

Only shipped checkpoints are used (fetched from Hugging Face ``microsoft/mattergen`` into the HF
cache on first use). :class:`~rtt.source.lanes.LaneSpec` has already refused any conditioning a
checkpoint was not trained for. The checkpoint's sha256 is read from the HF blob name (HF stores
LFS files under their sha256) and recorded in every structure's ``Origin``.

Every structure is flagged ``HYPOTHETICAL_STRUCTURE``. MatterGen has no pressure conditioning
(doc 10 §1): structures are 0 GPa proposals.
"""
from __future__ import annotations

import hashlib
import json
import math
import re
import shlex
import subprocess
import time
from pathlib import Path
from typing import Any

from ..records import Origin, RawStructure, StructureData
from .external import child_env, path_of
from .lanes import LaneJob, LaneSpec, norm_chemsys

EXTXYZ = "generated_crystals.extxyz"
GENERATOR = "mattergen@9242366"


class MatterGenRunError(RuntimeError):
    pass


def _literal(d: dict[str, Any]) -> str:
    """python-fire dict literal, no whitespace (MatterGen README)."""
    parts = []
    for k, v in d.items():
        parts.append(f"'{k}':" + (f"'{v}'" if isinstance(v, str) else repr(v)))
    return "{" + ",".join(parts) + "}"


def build_argv(spec: LaneSpec, n: int, seed: int, outdir: str | Path, *, machine: Any = None,
               executable: str | None = None) -> list[str]:
    exe = executable or str(path_of("mattergen_venv", machine) / "bin" / "mattergen-generate")
    bs = max(1, min(int(spec.batch_size), int(n)))
    nb = int(math.ceil(int(n) / bs))
    argv = ["nice", "-n", "19", exe, str(outdir)]
    if spec.model_path:
        argv += [f"--model_path={spec.model_path}", "--checkpoint_epoch=last"]
    else:
        argv.append(f"--pretrained_name={spec.checkpoint}")
    argv += [f"--batch_size={bs}", f"--num_batches={nb}", f"--seed={int(seed)}", "--record_trajectories=False"]
    if spec.conditioned:
        argv += [f"--properties_to_condition_on={_literal(dict(spec.conditioning))}",
                 f"--diffusion_guidance_factor={spec.effective_guidance}"]
    return argv


def checkpoint_sha256(name: str, machine: Any = None) -> str | None:
    """sha256 of ``checkpoints/<name>/checkpoints/last.ckpt`` in the HF cache (blob name).

    HF names an LFS blob by its content sha256. The symlink chain is followed one hop at a
    time and the *first* 64-hex name is returned: on this machine the HF blob is itself a
    symlink into a second store (``blobs/<xx>/<other-hash>``) whose name is NOT the content
    sha256, so ``Path.resolve()`` returned the wrong hash (doc 34 §1; content hashes checked
    with sha256sum for all six local checkpoints)."""
    import os
    hub = path_of("hf_hub", machine) / "models--microsoft--mattergen" / "snapshots"
    for p in sorted(hub.glob(f"*/checkpoints/{name}/checkpoints/last.ckpt")):
        cur = p
        for _ in range(16):
            if re.fullmatch(r"[0-9a-f]{64}", cur.name):
                return cur.name
            if not cur.is_symlink():
                break
            cur = cur.parent / os.readlink(cur)
        tgt = p.resolve()
        h = hashlib.sha256()
        with open(tgt, "rb") as fh:
            for chunk in iter(lambda: fh.read(1 << 22), b""):
                h.update(chunk)
        return h.hexdigest()
    return None


def read_extxyz(path: str | Path) -> list[StructureData]:
    from ase.io import read
    out = []
    for at in read(str(path), ":"):
        out.append(StructureData(lattice=[[float(x) for x in r] for r in at.cell.array],
                                 species=list(at.get_chemical_symbols()),
                                 frac_coords=[[float(x) % 1.0 for x in r] for r in at.get_scaled_positions(wrap=True)]))
    return out


def _tilde(p: Path) -> str:
    """Path with the home directory written as ~ (no user name in ledger rows / exports)."""
    s, home = str(p), str(Path.home())
    return "~" + s[len(home):] if s == home or s.startswith(home + "/") else s


def _run_key(spec: LaneSpec, n: int, seed: int) -> str:
    return hashlib.sha256(json.dumps({"s": spec.to_json(), "n": n, "seed": seed}, sort_keys=True)
                          .encode()).hexdigest()[:12]


def generate(lane: LaneSpec | LaneJob, n: int | None = None, seed: int | None = None, *,
             workdir: str | Path | None = None, batch_id: str | None = None, machine: Any = None,
             timeout_s: float = 12 * 3600, executable: str | None = None,
             reuse: bool = True) -> list[RawStructure]:
    """Run MatterGen for one lane job; returns the structures as :class:`RawStructure`.

    ``lane`` may be a :class:`LaneJob` (then ``n``/``seed``/``batch_id`` come from it).
    Output goes to ``workdir/<batch_id>_<key>/``; a finished run with the same argv is reused.
    """
    if isinstance(lane, LaneJob):
        spec, n, seed, batch_id = lane.spec, lane.n, lane.seed, batch_id or lane.batch_id
    else:
        spec = lane
    if n is None or seed is None:
        raise ValueError("n and seed are required")
    key = _run_key(spec, n, seed)
    batch_id = batch_id or f"{spec.lane}-{key}"
    base = Path(workdir or path_of("front_cache", machine) / "mattergen_runs")
    outdir = (base / f"{re.sub(r'[^A-Za-z0-9_.=-]+', '_', batch_id)}_{key}").resolve()
    outdir.mkdir(parents=True, exist_ok=True)
    argv = build_argv(spec, n, seed, outdir, machine=machine, executable=executable)
    stamp = outdir / "rtt_invocation.json"
    argv_s = shlex.join(argv)
    info: dict[str, Any] = {}
    prev = json.loads(stamp.read_text()) if stamp.exists() else {}
    if reuse and prev.get("argv") == argv_s and prev.get("finished") and (outdir / EXTXYZ).exists():
        info = {k: prev[k] for k in ("wall_s", "started") if k in prev}
        info["reused"] = True
    else:
        (outdir / EXTXYZ).unlink(missing_ok=True)
        stamp.write_text(json.dumps({"argv": argv_s, "finished": False}))
        t0 = time.time()
        with open(outdir / "mattergen.log", "w") as log:
            log.write(f"$ {argv_s}\n")
            log.flush()
            try:
                cp = subprocess.run(argv, env=child_env(), stdout=log, stderr=subprocess.STDOUT,
                                    timeout=timeout_s, cwd=str(outdir))
            except FileNotFoundError as ex:
                raise MatterGenRunError(f"mattergen-generate not found ({ex})") from ex
            except subprocess.TimeoutExpired as ex:
                raise MatterGenRunError(f"mattergen-generate timed out after {timeout_s} s") from ex
        info = {"wall_s": round(time.time() - t0, 1), "started": time.strftime("%Y-%m-%dT%H:%M:%S", time.localtime(t0)),
                "returncode": cp.returncode}
        if cp.returncode != 0 or not (outdir / EXTXYZ).exists():
            tail = (outdir / "mattergen.log").read_text()[-2000:]
            raise MatterGenRunError(f"mattergen-generate rc={cp.returncode}: {tail}")
        stamp.write_text(json.dumps({"argv": argv_s, "finished": True, **info}))
    structs = read_extxyz(outdir / EXTXYZ)
    ck_sha = checkpoint_sha256(spec.checkpoint, machine) if spec.checkpoint else None
    want = set(spec.conditioning.get("chemical_system", "").split("-")) - {""}
    out: list[RawStructure] = []
    for i, sd in enumerate(structs):
        flags = ["HYPOTHETICAL_STRUCTURE"]
        if want:
            got = set(sd.species)
            if not got <= want:
                flags.append("OFF_REQUESTED_SYSTEM")
            elif got != want:
                flags.append("SUBSYSTEM_OF_REQUESTED")
        origin = Origin(lane=spec.lane, generator=GENERATOR, checkpoint=spec.checkpoint or spec.model_path,
                        checkpoint_sha256=ck_sha, conditioning=dict(spec.conditioning) or None,
                        guidance=spec.effective_guidance, batch_id=batch_id, seed=int(seed), flags=flags)
        out.append(RawStructure(structure=sd, origin=origin, pressure_GPa=0.0,
                                source_ref=f"{_tilde(outdir / EXTXYZ)}#{i}",
                                meta={"index_in_batch": i, "lane_source_id": spec.source_id,
                                      "requested_chemsys": norm_chemsys(want) if want else None,
                                      "run": info, "n_requested": int(n)}))
    return out[: int(n)]


__all__ = ["generate", "build_argv", "checkpoint_sha256", "read_extxyz", "MatterGenRunError"]
