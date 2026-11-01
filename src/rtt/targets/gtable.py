"""Certified g table for NM-G (doc 03 §0, §2.5; doc 01 §0.3; NM-G is [THM-C], assumption CERT-NM).

Read exactly as lean/scripts/gen_nm_table.py does: every ``docs/theory/certificates/nm_*.json``
must have verdict CERTIFIED; its exact decimal ``m0_over_T2`` (never the float ``g_lo_K2``) gives
g_lo(λ_k) = 90000·m0 K² at 300 K, and g scales as T*², so NM-G uses ``m0(λ_up)·T*²`` exactly.
Lookup ('up' rule): λ_up = smallest tabulated λ_k >= λ_hi; beyond the table (or λ_hi unknown) the
λ = ∞ certificate. g is decreasing, so this is conservative (doc 03 §2.5).
A ``SUMMARY.json`` path is still accepted (same numbers, older format).
"""
from __future__ import annotations

import glob
import hashlib
import json
from dataclasses import dataclass, field
from fractions import Fraction
from pathlib import Path
from typing import Optional

from ..config import REPO_ROOT

DEFAULT_PATH = REPO_ROOT / "docs" / "theory" / "certificates"


@dataclass(frozen=True)
class GTable:
    lams: tuple[Fraction, ...]          # ascending, finite
    m0: tuple[Fraction, ...]            # exact m0_over_T2 per finite lambda
    m0_inf: Fraction
    lambda_cap: Optional[Fraction] = None
    source: str = str(DEFAULT_PATH)
    files_sha256: dict = field(default_factory=dict, compare=False, hash=False)

    def lambda_up(self, lam_hi: Optional[Fraction]) -> Optional[Fraction]:
        """None means lambda = infinity."""
        if lam_hi is None or (self.lambda_cap is not None and lam_hi > self.lambda_cap):
            return None
        for lk in self.lams:
            if lk >= lam_hi:
                return lk
        return None

    def m0_up(self, lam_hi: Optional[Fraction]) -> Fraction:
        lu = self.lambda_up(lam_hi)
        return self.m0_inf if lu is None else self.m0[self.lams.index(lu)]

    def g_up(self, lam_hi: Optional[Fraction], T_star: Fraction) -> Fraction:
        """g_lo(lambda_up)·s2 in K² (exact rational)."""
        return self.m0_up(lam_hi) * Fraction(T_star) ** 2


def _rows_from_certificates(d: Path) -> tuple[list, Fraction, dict]:
    rows, inf, shas = [], None, {}
    files = sorted(glob.glob(str(d / "nm_*.json")))
    if not files:
        raise FileNotFoundError(f"no nm_*.json in {d}")
    for f in files:
        raw = Path(f).read_bytes()
        c = json.loads(raw)
        if c.get("verdict", "CERTIFIED") != "CERTIFIED":
            raise ValueError(f"{f}: not CERTIFIED")
        shas[Path(f).name] = hashlib.sha256(raw).hexdigest()
        if str(c["lambda"]).lower() == "inf":
            inf = Fraction(c["m0_over_T2"])
        else:
            rows.append((Fraction(c["lambda"]), Fraction(c["m0_over_T2"])))
    if inf is None:
        raise ValueError(f"{d}: no lambda = inf certificate")
    return rows, inf, shas


def _rows_from_summary(p: Path) -> tuple[list, Fraction, dict]:
    d = json.loads(p.read_text(encoding="utf-8"))
    rows, inf = [], None
    for r in d["NM_G"]:
        if r.get("verdict") != "CERTIFIED":
            continue
        if str(r["lam"]).lower() == "inf":
            inf = Fraction(r["m0_over_T2"])
        else:
            rows.append((Fraction(r["lam"]), Fraction(r["m0_over_T2"])))
    if inf is None:
        raise ValueError(f"{p}: no certified lambda = inf row")
    return rows, inf, {p.name: hashlib.sha256(p.read_bytes()).hexdigest()}


def load_gtable(path: str | Path | None = None, lambda_cap: float | None = None) -> GTable:
    p = Path(path) if path else DEFAULT_PATH
    rows, inf, shas = _rows_from_certificates(p) if p.is_dir() else _rows_from_summary(p)
    rows.sort()
    for (l1, a), (l2, b) in zip(rows, rows[1:]):
        if not b < a:
            raise ValueError(f"g table not decreasing at lambda {l1} -> {l2}")
    if rows and not inf < rows[-1][1]:
        raise ValueError("g(inf) not below the last row")
    cap = None if lambda_cap is None else Fraction(str(lambda_cap))
    return GTable(tuple(r[0] for r in rows), tuple(r[1] for r in rows), inf, cap, str(p), shas)
