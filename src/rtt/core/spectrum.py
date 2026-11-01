"""Spectrum: alpha^2F as a finite atomic measure on (0, inf), energies in kelvin (ARCH.md section 2).

The model (SPEC.md) uses k_B = hbar = 1, so frequencies and temperatures are both in kelvin.
``omega_K[j] > 0`` are the atom positions, ``weight_K[j] >= 0`` the masses (in kelvin, because
alpha^2F is dimensionless and the masses are alpha^2F * d omega).  ``A = sum(weight_K)``.

``rel_envelope = e`` is an ASSUMPTION recorded in every verdict: the true masses lie in
[(1-e) w_j, (1+e) w_j] (same positions).  The rigorous engine treats every spectrum in that box.

All stored numbers are Python floats (IEEE doubles), i.e. exact dyadic rationals; the rigorous
engine converts them to Arb balls exactly, so JSON round trips (repr) preserve the exact input.
"""
from __future__ import annotations

import json
import math
from dataclasses import dataclass, field
from typing import Any, Sequence

# --- unit conversions to kelvin (energy E -> E / k_B).  CODATA 2018 exact SI values. ---
_KB_EV = 8.617333262e-5            # eV / K (exact in SI 2019: k_B = 1.380649e-23 J/K, e exact)
_H_J = 6.62607015e-34              # J s (exact)
_KB_J = 1.380649e-23               # J / K (exact)
_C_CM = 2.99792458e10              # cm / s (exact)
_RY_EV = 13.605693122994           # eV (CODATA 2018)

UNIT_TO_K: dict[str, float] = {
    "K": 1.0,
    "meV": 1e-3 / _KB_EV,                      # 11.6045181 K
    "eV": 1.0 / _KB_EV,
    "THz": _H_J * 1e12 / _KB_J,                # h * 1 THz / k_B = 47.9924 K  (frequency nu, not angular)
    "cm-1": _H_J * _C_CM / _KB_J,              # 1.43877688 K
    "Ry": _RY_EV / _KB_EV,                     # 157887.51 K
    "mRy": 1e-3 * _RY_EV / _KB_EV,
}
UNIT_ALIASES = {"kelvin": "K", "k": "K", "mev": "meV", "ev": "eV", "thz": "THz",
                "cm^-1": "cm-1", "cm-1": "cm-1", "1/cm": "cm-1", "ry": "Ry", "mry": "mRy"}


def unit_factor(unit: str) -> float:
    """Multiplicative factor converting an energy in ``unit`` to kelvin."""
    u = UNIT_ALIASES.get(unit.strip().lower(), unit.strip())
    if u not in UNIT_TO_K:
        raise ValueError(f"unknown energy unit {unit!r}; known: {sorted(UNIT_TO_K)}")
    return UNIT_TO_K[u]


@dataclass(frozen=True)
class Spectrum:
    omega_K: tuple[float, ...]
    weight_K: tuple[float, ...]
    rel_envelope: float = 0.0
    label: str = ""
    provenance: dict = field(default_factory=dict, compare=False, hash=False)

    def __post_init__(self) -> None:
        om = tuple(float(x) for x in self.omega_K)
        w = tuple(float(x) for x in self.weight_K)
        object.__setattr__(self, "omega_K", om)
        object.__setattr__(self, "weight_K", w)
        object.__setattr__(self, "rel_envelope", float(self.rel_envelope))
        if len(om) != len(w):
            raise ValueError("omega_K and weight_K must have the same length")
        if len(om) == 0:
            raise ValueError("empty spectrum")
        for x in om:
            if not (math.isfinite(x) and x > 0):
                raise ValueError(f"atom position must be finite and > 0, got {x}")
        for x in w:
            if not (math.isfinite(x) and x >= 0):
                raise ValueError(f"atom weight must be finite and >= 0, got {x}")
        e = self.rel_envelope
        if not (math.isfinite(e) and 0.0 <= e < 1.0):
            raise ValueError(f"rel_envelope must lie in [0, 1), got {e}")

    # ---------------- float helpers (informational; NOT used for decisions) ----------------
    def A(self) -> float:
        """Total weight A = int alpha^2F d omega (kelvin), nominal masses."""
        return math.fsum(self.weight_K)

    def A_hi(self) -> float:
        """Envelope maximum (1+e) A (float; the rigorous engine recomputes it in Arb)."""
        return (1.0 + self.rel_envelope) * self.A()

    def lambda_(self) -> float:
        """lambda = 2 int alpha^2F(w)/w dw."""
        return 2.0 * math.fsum(w / o for o, w in zip(self.omega_K, self.weight_K))

    def omega_log(self) -> float:
        """omega_log = exp[(2/lambda) int alpha^2F(w) ln(w)/w dw] (kelvin)."""
        lam = self.lambda_()
        if lam <= 0:
            return float("nan")
        s = math.fsum(w / o * math.log(o) for o, w in zip(self.omega_K, self.weight_K))
        return math.exp(2.0 * s / lam)

    def omega2(self) -> float:
        """sqrt(<w^2>) with <f> = (2/lambda) int f(w) alpha^2F(w)/w dw (kelvin)."""
        lam = self.lambda_()
        s = math.fsum(w * o for o, w in zip(self.omega_K, self.weight_K))
        return math.sqrt(2.0 * s / lam) if lam > 0 else float("nan")

    def omega_max(self) -> float:
        return max(o for o, w in zip(self.omega_K, self.weight_K)) if self.omega_K else 0.0

    def n_atoms(self) -> int:
        return len(self.omega_K)

    def with_envelope(self, e: float) -> "Spectrum":
        return Spectrum(self.omega_K, self.weight_K, e, self.label, dict(self.provenance))

    # ---------------- constructors ----------------
    @classmethod
    def from_atoms(cls, omegas: Sequence[float], weights: Sequence[float], unit: str = "K",
                   rel_envelope: float = 0.0, label: str = "", provenance: dict | None = None,
                   drop_zero: bool = True) -> "Spectrum":
        f = unit_factor(unit)
        om, w = [], []
        for o, p in zip(omegas, weights):
            if drop_zero and p == 0:
                continue
            om.append(float(o) * f)
            w.append(float(p) * f)
        return cls(tuple(om), tuple(w), rel_envelope, label, dict(provenance or {}))

    @classmethod
    def from_grid(cls, omega: Sequence[float], a2f: Sequence[float], unit: str = "K",
                  rel_envelope: float = 0.0, label: str = "", provenance: dict | None = None,
                  ) -> "Spectrum":
        """Atomic measure from a sampled density alpha^2F(omega) on a grid.

        Midpoint weights: grid point omega_i (converted to kelvin) carries mass
        a2F(omega_i) * (b_i - a_i), where [a_i, b_i] is bounded by the midpoints to its
        neighbours (the first/last cells extend by half a spacing, the first is clipped at 0).
        Points with omega <= 0 are dropped; negative a2F values (smearing artefacts) are clipped
        to 0.  Both are recorded in ``provenance`` ("dropped_nonpositive_omega_mass",
        "clipped_negative_mass", in kelvin) -- this discretisation is an input convention, not part
        of any proof.
        """
        f = unit_factor(unit)
        pts = sorted((float(o) * f, float(a)) for o, a in zip(omega, a2f))
        if len(pts) < 2:
            raise ValueError("need at least two grid points")
        xs = [p[0] for p in pts]
        n = len(xs)
        edges = [0.0] * (n + 1)
        for i in range(1, n):
            edges[i] = 0.5 * (xs[i - 1] + xs[i])
        edges[0] = xs[0] - 0.5 * (xs[1] - xs[0])
        edges[n] = xs[-1] + 0.5 * (xs[-1] - xs[-2])
        om, w = [], []
        dropped = 0.0
        clipped = 0.0
        for i, (x, a) in enumerate(pts):
            lo = max(edges[i], 0.0)
            width = max(edges[i + 1] - lo, 0.0)
            mass = a * width
            if x <= 0:
                dropped += abs(mass)
                continue
            if mass < 0:
                clipped += -mass
                continue
            if mass == 0:
                continue
            om.append(x)
            w.append(mass)
        prov = dict(provenance or {})
        prov.setdefault("discretisation", "midpoint weights (rtt.core.spectrum.Spectrum.from_grid)")
        prov["input_unit"] = unit
        prov["dropped_nonpositive_omega_mass_K"] = dropped
        prov["clipped_negative_mass_K"] = clipped
        return cls(tuple(om), tuple(w), rel_envelope, label, prov)

    # ---------------- serialisation (exact: floats round-trip through repr) ----------------
    def to_dict(self) -> dict[str, Any]:
        return {"omega_K": list(self.omega_K), "weight_K": list(self.weight_K),
                "rel_envelope": self.rel_envelope, "label": self.label,
                "provenance": dict(self.provenance)}

    def to_json(self, **kw) -> str:
        return json.dumps(self.to_dict(), **kw)

    @classmethod
    def from_dict(cls, d: dict[str, Any]) -> "Spectrum":
        return cls(tuple(d["omega_K"]), tuple(d["weight_K"]), d.get("rel_envelope", 0.0),
                   d.get("label", ""), dict(d.get("provenance", {})))

    @classmethod
    def from_json(cls, s: str) -> "Spectrum":
        return cls.from_dict(json.loads(s))

    def summary(self) -> str:
        return (f"Spectrum({self.label or '?'}: {self.n_atoms()} atoms, A={self.A():.4g} K, "
                f"lambda={self.lambda_():.4g}, omega_log={self.omega_log():.4g} K, e={self.rel_envelope})")
