"""Configuration loading (architecture §7, contract R7).

* ``load_targets()``, ``load_tiers()``, ``load_machine()`` read the TOML files in
  ``src/rtt/configs/`` (or an explicit path) and return a :class:`LoadedConfig` carrying the raw
  text, its sha256 and the parsed dict.
* Every numeric threshold in ``targets.toml`` and ``tiers.toml`` must have a source id (R7):
  either as an inline table ``{ value = ..., source = "...", provisional = ... }`` or through a
  ``source`` key of the enclosing table (section-level source). The loader raises
  :class:`ConfigError` otherwise.
* :class:`Threshold` is the value object handed to stages: ``(value, source_id, provisional)``.
"""
from __future__ import annotations

import hashlib
import os
import tomllib
from dataclasses import dataclass, field
from fractions import Fraction
from pathlib import Path
from typing import Any, Iterator

CONFIG_DIR = Path(__file__).resolve().parent / "configs"
REPO_ROOT = Path(__file__).resolve().parents[2]

SOURCE_KEYS = ("source", "source_id")


class ConfigError(ValueError):
    """Invalid or unsourced configuration."""


@dataclass(frozen=True)
class Threshold:
    value: Any
    source_id: str
    provisional: bool = False

    def as_dict(self) -> dict:
        return {"value": self.value, "source_id": self.source_id, "provisional": self.provisional}


def sha256_of(obj: str | bytes | os.PathLike) -> str:
    """sha256 hex of bytes, of a str (utf-8), or of a file given as a Path."""
    if isinstance(obj, (bytes, bytearray)):
        return hashlib.sha256(obj).hexdigest()
    if isinstance(obj, os.PathLike):
        h = hashlib.sha256()
        with open(obj, "rb") as fh:
            for chunk in iter(lambda: fh.read(1 << 20), b""):
                h.update(chunk)
        return h.hexdigest()
    return hashlib.sha256(obj.encode("utf-8")).hexdigest()


def _is_number(v: Any) -> bool:
    return isinstance(v, (int, float)) and not isinstance(v, bool)


def _is_numeric_leaf(v: Any) -> bool:
    if _is_number(v):
        return True
    if isinstance(v, list) and v and all(_is_number(x) or (isinstance(x, list) and all(_is_number(y) for y in x)) for x in v):
        return True
    return False


def _is_threshold_table(v: Any) -> bool:
    return isinstance(v, dict) and "value" in v


def _source_of(d: dict) -> str | None:
    for k in SOURCE_KEYS:
        if isinstance(d.get(k), str) and d[k].strip():
            return d[k]
    return None


def _walk_unsourced(d: dict, inherited: str | None, prefix: str) -> Iterator[str]:
    here = _source_of(d) or inherited
    for k, v in d.items():
        if k in SOURCE_KEYS or k in ("schema", "provisional"):
            continue
        path = f"{prefix}{k}"
        if isinstance(v, dict):
            if _is_threshold_table(v):
                if not (_source_of(v) or here):
                    yield path
            else:
                # inline tables of thresholds (e.g. H7_ehull = {window_eV=..., source=...})
                yield from _walk_unsourced(v, here, path + ".")
        elif _is_numeric_leaf(v) and here is None:
            yield path


def validate_sources(raw: dict) -> None:
    missing = list(_walk_unsourced(raw, None, ""))
    if missing:
        raise ConfigError("thresholds without a source id (contract R7): " + ", ".join(missing))


@dataclass(frozen=True)
class LoadedConfig:
    kind: str
    path: str
    text: str
    sha256: str
    raw: dict = field(repr=False)

    def get(self, dotted: str, default: Any = ...) -> Any:
        cur: Any = self.raw
        for part in dotted.split("."):
            if not isinstance(cur, dict) or part not in cur:
                if default is ...:
                    raise ConfigError(f"{self.kind}: missing key {dotted!r}")
                return default
            cur = cur[part]
        return cur

    def threshold(self, dotted: str) -> Threshold:
        """Return a sourced threshold. ``dotted`` may name an inline threshold table
        (``{value, source, provisional}``) or a plain leaf covered by a section source."""
        parts = dotted.split(".")
        node: Any = self.raw
        source = _source_of(self.raw)
        provisional = bool(self.raw.get("provisional", False)) if isinstance(self.raw, dict) else False
        for part in parts:
            if not isinstance(node, dict) or part not in node:
                raise ConfigError(f"{self.kind}: missing threshold {dotted!r}")
            if isinstance(node, dict):
                source = _source_of(node) or source
                provisional = bool(node.get("provisional", provisional))
            node = node[part]
        if isinstance(node, dict):
            if "value" not in node:
                raise ConfigError(f"{self.kind}: {dotted!r} is a table without 'value'")
            src = _source_of(node) or source
            prov = bool(node.get("provisional", provisional))
            value = node["value"]
        else:
            src, prov, value = source, provisional, node
        if not src:
            raise ConfigError(f"{self.kind}: threshold {dotted!r} has no source id (R7)")
        return Threshold(value=value, source_id=src, provisional=prov)


class Targets(LoadedConfig):
    """targets.toml with derived quantities (doc 01 §0.4)."""

    @property
    def T_star(self) -> Fraction:
        return Fraction(str(self.get("tc.T_star_K")))

    @property
    def a0(self) -> Fraction:
        return Fraction(self.get("tc.a0"))

    @property
    def s1(self) -> Fraction:
        return self.T_star / 300

    @property
    def s2(self) -> Fraction:
        return self.s1 ** 2

    @property
    def floor(self) -> Fraction:
        """a0 * T* (= 1299.444435 K at T* = 300 K)."""
        return self.a0 * self.T_star

    def g_table_path(self) -> Path:
        p = Path(self.get("nm.g_table"))
        return p if p.is_absolute() else REPO_ROOT / p

    def with_T_star(self, T_star_K: float) -> "Targets":
        """Surrogate target (control lane C): same file, different T*. The sha256 is that of the
        file plus the override, so outcomes record the surrogate."""
        raw = {**self.raw, "tc": {**self.raw["tc"], "T_star_K": float(T_star_K)}}
        text = self.text + f"\n# surrogate override: tc.T_star_K = {float(T_star_K)!r}\n"
        return Targets(kind=self.kind, path=self.path, text=text, sha256=sha256_of(text), raw=raw)


def _resolve(path: str | os.PathLike | None, default_name: str) -> Path:
    if path is None:
        env = os.environ.get("RTT_CONFIG_DIR")
        base = Path(env) if env else CONFIG_DIR
        return base / default_name
    return Path(path)


def _load(kind: str, path: str | os.PathLike | None, default_name: str, cls=LoadedConfig, check_sources=True):
    p = _resolve(path, default_name)
    text = p.read_text(encoding="utf-8")
    try:
        raw = tomllib.loads(text)
    except tomllib.TOMLDecodeError as e:
        raise ConfigError(f"{kind}: {p}: {e}") from e
    expected = {"targets": "rtt.targets/1", "tiers": "rtt.tiers/1", "machine": "rtt.machine/1"}[kind]
    if raw.get("schema") != expected:
        raise ConfigError(f"{kind}: schema must be {expected!r}, got {raw.get('schema')!r}")
    if check_sources:
        validate_sources(raw)
    return cls(kind=kind, path=str(p), text=text, sha256=sha256_of(text), raw=raw)


def load_targets(path: str | os.PathLike | None = None) -> Targets:
    t = _load("targets", path, "targets.toml", Targets)
    a0 = t.a0
    if a0 < Fraction(433148145, 100000000):
        raise ConfigError("targets: a0 below the certified 4.33148145 (would lower every floor)")
    if t.T_star <= 0:
        raise ConfigError("targets: T_star_K must be > 0")
    return t


def load_tiers(path: str | os.PathLike | None = None) -> LoadedConfig:
    return _load("tiers", path, "tiers.toml")


def load_machine(path: str | os.PathLike | None = None) -> LoadedConfig:
    # machine.toml holds resources, not physics thresholds: no source ids required.
    return _load("machine", path, "machine.toml", check_sources=False)


def thread_env(machine: LoadedConfig | None = None) -> dict[str, str]:
    """Environment for child processes (CONTEXT.md compute rules: <= 2 threads)."""
    env = {"OMP_NUM_THREADS": "2", "OPENBLAS_NUM_THREADS": "2", "MKL_NUM_THREADS": "2"}
    if machine is not None:
        for k, v in machine.get("env", {}).items():
            if k != "nice":
                env[k] = str(v)
    return env
