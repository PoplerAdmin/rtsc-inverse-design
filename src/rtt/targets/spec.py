"""TargetSpec (architecture §5.1): everything S0 derives from targets.toml."""
from __future__ import annotations

from dataclasses import dataclass, field
from fractions import Fraction
from pathlib import Path
from typing import Any, Optional

from ..config import REPO_ROOT, Targets, Threshold, load_targets
from .gtable import GTable, load_gtable


@dataclass
class TargetSpec:
    T_star: Fraction
    a0: Fraction
    floor: Fraction            # a0 * T*
    s1: Fraction
    s2: Fraction
    g_table: GTable
    eps_sys_min: Fraction
    eps_omega_min: Fraction
    aim_box: dict
    hard: dict
    soft: dict
    chemsys_list: list[str] = field(default_factory=list)
    targets_sha256: str = ""
    source: Optional[Targets] = None
    decide_flags: dict = field(default_factory=dict)

    @classmethod
    def from_targets(cls, t: Targets) -> "TargetSpec":
        cap = t.get("nm.lambda_cap", None)
        g = load_gtable(t.g_table_path(), None if cap is None else float(cap))
        cs: list[str] = []
        p = Path(t.get("chemsys.list", ""))
        if not p.is_absolute():
            p = REPO_ROOT / p
        if p.is_file():
            cs = [ln.split("#")[0].strip() for ln in p.read_text(encoding="utf-8").splitlines()]
            cs = [c for c in cs if c]
        return cls(T_star=t.T_star, a0=t.a0, floor=t.floor, s1=t.s1, s2=t.s2, g_table=g,
                   eps_sys_min=Fraction(str(t.threshold("envelope.eps_sys").value)),
                   eps_omega_min=Fraction(str(t.threshold("envelope.eps_omega_min").value)),
                   aim_box=dict(t.get("aim")), hard=dict(t.get("hard")), soft=dict(t.get("soft")),
                   chemsys_list=cs, targets_sha256=t.sha256, source=t,
                   decide_flags={k: bool(v) for k, v in t.get("decide", {}).items()
                                 if k in ("w_lo_sum_check", "moment_consistency", "nmg_certified")})

    def surrogate(self, T_star_K: float) -> "TargetSpec":
        assert self.source is not None
        return TargetSpec.from_targets(self.source.with_T_star(T_star_K))

    def hard_threshold(self, name: str, key: str) -> Threshold:
        d = self.hard[name]
        return Threshold(d[key], d["source"], bool(d.get("provisional", False)))


def load_spec(path=None) -> TargetSpec:
    return TargetSpec.from_targets(load_targets(path))
