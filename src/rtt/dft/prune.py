"""Hash-then-delete scratch pruning and a disk guard for QE runs (architecture W5, §8.3).

Evidence policy: **every file is SHA-256 hashed and its manifest line
written and fsync'ed to every manifest before it is deleted.** Manifests are the
outside-verification method: :func:`verify_manifest` re-hashes the listed files (of a re-run, or
of a restored copy) and reports MATCH / MISMATCH / MISSING.

Manifests: ``<root>/PRUNED_MANIFEST.tsv`` plus a central copy in ``central_dir`` (default from
the environment variable ``RTT_PRUNE_MANIFEST_DIR``; none if unset). Columns::

    utc_time  machine  path  bytes  sha256  mtime_utc  reason

Paths are written relative to the home directory (``~/...``) and the machine by a neutral label
(``RTT_MACHINE_LABEL``, default ``run-machine-1``); never a hostname or user name.

What is regenerable and may be pruned (QE 7.6, ``ldisp``):

* wavefunction files ``<prefix>.wfc<N>`` and ``wfc<N>.dat``/``wfcup|dw<N>.dat`` outside
  ``<prefix>.save/``;
* per-q (``_ph0/<prefix>.q_<N>/``): the wavefunctions of a *finished* q (its dyn file is complete
  and ph.x has written the next one), and optionally its ``dvscf`` files
  (``prune_dvscf_per_q``; this prevents a later EPW run from reusing them);
* after a successful run, ``tmp/`` as a whole (``prune_final = "tmp"``).

Never touched by the per-q or wfc modes: dyn files, ``elph_dir``, ``*.xml``, the charge density,
anything inside ``<prefix>.save/``, ``<prefix>.a2Fsave``.

Disk guard: :func:`disk_guard` raises :class:`DiskSpaceLow` before a
step starts; :class:`ScratchWatcher` runs beside a long ph.x step, prunes finished q dirs and,
below ``min_free_gb``, pauses (SIGSTOP) **only the descendants of this Python process**, resuming
(SIGCONT) above ``resume_free_gb``. Other users' QE processes are never signalled.
"""
from __future__ import annotations

import hashlib
import os
import re
import shutil
import signal
import threading
import time
from dataclasses import dataclass, field
from pathlib import Path
from typing import Callable

HEADER = "utc_time\tmachine\tpath\tbytes\tsha256\tmtime_utc\treason\n"
WFC_RE = re.compile(r"(^.+\.wfc\d+$)|(^wfc(up|dw)?\d+\.(dat|hdf5)$)")
DVSCF_RE = re.compile(r"^.+\.dvscf\d*(_paw)?\d*$")
QDIR_RE = re.compile(r"\.q_(\d+)$")


def machine_label() -> str:
    return os.environ.get("RTT_MACHINE_LABEL", "run-machine-1")


def central_dir_default() -> Path | None:
    v = os.environ.get("RTT_PRUNE_MANIFEST_DIR")
    return Path(v).expanduser() if v else None


def disp(p: Path | str) -> str:
    """Path for publication: '~/...' under the home directory, else as is."""
    s, home = str(Path(p).resolve()), str(Path.home())
    return "~" + s[len(home):] if s == home or s.startswith(home + "/") else s


def tag_of(root: Path) -> str:
    return re.sub(r"[^A-Za-z0-9_.-]+", "_", disp(root).lstrip("~/"))


def sha256_file(path: Path | str) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 22), b""):
            h.update(chunk)
    return h.hexdigest()


def _utc(t: float | None = None) -> str:
    return time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime(t))


def manifest_paths(root: Path, central_dir: Path | None) -> list[Path]:
    out = [Path(root) / "PRUNED_MANIFEST.tsv"]
    if central_dir is not None:
        Path(central_dir).mkdir(parents=True, exist_ok=True)
        out.append(Path(central_dir) / f"{tag_of(Path(root))}.tsv")
    return out


def _append(paths: list[Path], line: str) -> None:
    for m in paths:
        new = not m.exists()
        with open(m, "a") as f:
            if new:
                f.write(HEADER)
            f.write(line)
            f.flush()
            os.fsync(f.fileno())


@dataclass
class PruneReport:
    files: int = 0
    bytes: int = 0
    manifests: list[str] = field(default_factory=list)


def prune_files(files: list[Path], root: Path, reason: str, *, central_dir: Path | None = None,
                dry_run: bool = False) -> PruneReport:
    """Hash, record (all manifests, fsync), then delete each file. Symlinks, non-files and the
    manifests themselves are skipped."""
    if not reason:
        raise ValueError("a reason is required")
    root = Path(root)
    mans = manifest_paths(root, central_dir)
    rep = PruneReport(manifests=[disp(m) for m in mans])
    for p in sorted(set(Path(x) for x in files)):
        if p.is_symlink() or not p.is_file() or p.name in ("PRUNED_MANIFEST.tsv", "RUN_MANIFEST.tsv"):
            continue
        st = p.stat()
        digest = sha256_file(p)
        line = "\t".join([_utc(), machine_label(), disp(p), str(st.st_size), digest,
                          _utc(st.st_mtime), reason]) + "\n"
        if dry_run:
            continue
        _append(mans, line)
        p.unlink()
        rep.files += 1
        rep.bytes += st.st_size
    return rep


def wfc_files(root: Path) -> list[Path]:
    """Wavefunction files under root, excluding anything inside a ``*.save`` directory."""
    root = Path(root)
    if not root.is_dir():
        return []
    return [p for p in root.rglob("*") if WFC_RE.match(p.name)
            and not any(part.endswith(".save") for part in p.relative_to(root).parts[:-1])]


def dvscf_files(root: Path) -> list[Path]:
    root = Path(root)
    if not root.is_dir():
        return []
    return [p for p in root.rglob("*") if DVSCF_RE.match(p.name)
            and not any(part.endswith(".save") for part in p.relative_to(root).parts[:-1])]


def _dyn_complete(dyn: Path) -> bool:
    return dyn.is_file() and dyn.stat().st_size > 0 and bool(dyn.read_text(errors="replace").strip())


def finished_q_numbers(dyndir: Path, fildyn: str = "dyn") -> set[int]:
    """q indices N whose dyn file is complete AND whose successor dyn exists (ph.x moved on).
    ph.x creates every ``_ph0/<prefix>.q_N`` directory up front, so directory presence says
    nothing about progress; the dyn files do. The last q is never "finished" here."""
    out = set()
    for f in Path(dyndir).glob(f"{fildyn}[0-9]*"):
        m = re.fullmatch(re.escape(fildyn) + r"(\d+)", f.name)
        if not m or m.group(1) == "0":
            continue
        n = int(m.group(1))
        if _dyn_complete(f) and (Path(dyndir) / f"{fildyn}{n + 1}").exists():
            out.add(n)
    return out


def q_dirs(ph0: Path) -> dict[int, Path]:
    out = {}
    if Path(ph0).is_dir():
        for d in Path(ph0).iterdir():
            m = QDIR_RE.search(d.name)
            if d.is_dir() and m:
                out[int(m.group(1))] = d
    return out


def prune_q_scratch(workdir: Path | str, q_done: set[int] | None = None, *,
                    with_dvscf: bool = False, central_dir: Path | None = None,
                    reason: str = "per-q prune: finished q of a running/finished ph.x") -> int:
    """Hash-then-delete the wavefunctions (and optionally dvscf) of finished q points of a
    ph.x ``ldisp`` run in ``workdir`` (scratch ``workdir/tmp/_ph0``). ``q_done`` defaults to
    :func:`finished_q_numbers`. q = 1 (Gamma) lives in ``_ph0`` itself, not in a q dir, and is
    left for the final prune. Returns bytes freed."""
    w = Path(workdir)
    done = finished_q_numbers(w) if q_done is None else set(q_done)
    files: list[Path] = []
    for n, d in q_dirs(w / "tmp" / "_ph0").items():
        if n in done:
            files += wfc_files(d)
            if with_dvscf:
                files += dvscf_files(d)
    if not files:
        return 0
    return prune_files(files, w, reason, central_dir=central_dir).bytes


def prune_final(workdir: Path | str, mode: str, *, central_dir: Path | None = None,
                reason: str = "final prune after successful run") -> int:
    """mode 'none' | 'wfc' (all wfc files under tmp/) | 'tmp' (every file under tmp/)."""
    w = Path(workdir)
    tmp = w / "tmp"
    if mode == "none" or not tmp.is_dir():
        return 0
    if mode == "wfc":
        files = wfc_files(tmp)
    elif mode == "tmp":
        files = [p for p in tmp.rglob("*") if p.is_file()]
    else:
        raise ValueError(f"unknown prune mode {mode!r}")
    rep = prune_files(files, w, reason, central_dir=central_dir)
    if mode == "tmp":
        for d in sorted((p for p in tmp.rglob("*") if p.is_dir()), key=lambda p: -len(p.parts)):
            try:
                d.rmdir()
            except OSError:
                pass
    return rep.bytes


def inventory(root: Path | str, *, central_dir: Path | None = None) -> Path:
    """Hash every kept file under root into ``<root>/RUN_MANIFEST.tsv`` (+ central copy)."""
    root = Path(root)
    rows = [HEADER.replace("\treason", "\tstatus")]
    stamp = _utc()
    for p in sorted(root.rglob("*")):
        if p.is_file() and not p.is_symlink() and p.name not in ("RUN_MANIFEST.tsv", "PRUNED_MANIFEST.tsv"):
            st = p.stat()
            rows.append("\t".join([stamp, machine_label(), disp(p), str(st.st_size), sha256_file(p),
                                   _utc(st.st_mtime), "kept"]) + "\n")
    out = root / "RUN_MANIFEST.tsv"
    targets = [out] + ([Path(central_dir) / f"{tag_of(root)}.RUN.tsv"] if central_dir else [])
    for m in targets:
        m.parent.mkdir(parents=True, exist_ok=True)
        m.write_text("".join(rows))
    return out


def verify_manifest(manifest: Path | str, base: Path | str | None = None) -> dict[str, int]:
    """Re-hash the files a manifest lists ('~' resolved to ``base``, default the home dir) and
    count MATCH / MISMATCH / MISSING."""
    base = Path(base) if base is not None else Path.home()
    counts = {"MATCH": 0, "MISMATCH": 0, "MISSING": 0}
    for row in Path(manifest).read_text().splitlines()[1:]:
        f = row.split("\t")
        if len(f) < 5:
            continue
        p = base / f[2][2:] if f[2].startswith("~/") else Path(f[2])
        st = "MISSING" if not p.is_file() else ("MATCH" if sha256_file(p) == f[4] else "MISMATCH")
        counts[st] += 1
    return counts


# ---------------------------------------------------------------------- disk guard
class DiskSpaceLow(RuntimeError):
    pass


def disk_free_gb(path: Path | str) -> float:
    p = Path(path)
    while not p.exists() and p != p.parent:
        p = p.parent
    return shutil.disk_usage(p).free / 1e9


def disk_guard(path: Path | str, min_free_gb: float) -> None:
    """Raise :class:`DiskSpaceLow` if the file system holding ``path`` has less than
    ``min_free_gb`` free. Under WSL the scratch lives in the VM disk image on the host drive, so
    pass a path on the host drive (e.g. /mnt/c) as ``path`` when that is the real limit."""
    free = disk_free_gb(path)
    if free < float(min_free_gb):
        raise DiskSpaceLow(f"only {free:.1f} GB free at {disp(path)} (< {min_free_gb} GB)")


def descendants(pid: int) -> list[int]:
    """PIDs of all descendants of ``pid`` (Linux /proc)."""
    children: dict[int, list[int]] = {}
    for d in Path("/proc").iterdir() if Path("/proc").is_dir() else []:
        if not d.name.isdigit():
            continue
        try:
            stat = (d / "stat").read_text()
        except OSError:
            continue
        rest = stat[stat.rfind(")") + 2:].split()
        try:
            ppid = int(rest[1])
        except (IndexError, ValueError):
            continue
        children.setdefault(ppid, []).append(int(d.name))
    out, stack = [], [pid]
    while stack:
        for c in children.get(stack.pop(), []):
            out.append(c)
            stack.append(c)
    return out


class ScratchWatcher:
    """Background thread for a running QE step: every ``interval_s`` it prunes finished q dirs
    (if ``prune``), and pauses/resumes this process's descendants on low disk space."""

    def __init__(self, workdir: Path, *, interval_s: float = 60.0, prune: bool = True,
                 with_dvscf: bool = False, central_dir: Path | None = None,
                 disk_path: Path | None = None, min_free_gb: float | None = None,
                 resume_free_gb: float | None = None, log: Callable[[str], None] | None = None):
        self.workdir = Path(workdir)
        self.interval_s = interval_s
        self.prune = prune
        self.with_dvscf = with_dvscf
        self.central_dir = central_dir
        self.disk_path = Path(disk_path) if disk_path else self.workdir
        self.min_free_gb = min_free_gb
        self.resume_free_gb = resume_free_gb if resume_free_gb is not None else (
            None if min_free_gb is None else min_free_gb + 10.0)
        self.log = log or (lambda s: None)
        self.pruned_bytes = 0
        self.paused_s = 0.0
        self._paused_at: float | None = None
        self._stop = threading.Event()
        self._thread = threading.Thread(target=self._loop, daemon=True)

    def __enter__(self) -> "ScratchWatcher":
        self._thread.start()
        return self

    def __exit__(self, *exc) -> None:
        self._stop.set()
        self._thread.join(timeout=30)
        if self._paused_at is not None:
            self._signal(signal.SIGCONT)
            self.paused_s += time.time() - self._paused_at
            self._paused_at = None
        # a last pass after the step: q dirs finished in the final interval
        if self.prune:
            self._prune_once()

    def _signal(self, sig: int) -> None:
        for pid in descendants(os.getpid()):
            try:
                os.kill(pid, sig)
            except OSError:
                pass

    def _prune_once(self) -> None:
        try:
            self.pruned_bytes += prune_q_scratch(self.workdir, with_dvscf=self.with_dvscf,
                                                 central_dir=self.central_dir)
        except Exception as exc:  # pruning must never kill the run
            self.log(f"prune error: {exc}")

    def _loop(self) -> None:
        while not self._stop.wait(self.interval_s):
            if self.prune:
                self._prune_once()
            if self.min_free_gb is None:
                continue
            free = disk_free_gb(self.disk_path)
            if self._paused_at is None and free < self.min_free_gb:
                self._signal(signal.SIGSTOP)
                self._paused_at = time.time()
                self.log(f"PAUSED own QE processes: {free:.1f} GB free < {self.min_free_gb}")
            elif self._paused_at is not None and free >= (self.resume_free_gb or 0):
                self._signal(signal.SIGCONT)
                self.paused_s += time.time() - self._paused_at
                self._paused_at = None
                self.log(f"RESUMED own QE processes: {free:.1f} GB free")
