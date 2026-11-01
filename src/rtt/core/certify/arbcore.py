"""Rigorous (Arb ball arithmetic) primitives for the certify engine.

Everything returned here is an Arb ball (python-flint ``arb``) that CONTAINS the exact real
quantity, given exact inputs (IEEE doubles convert to Arb exactly).  Lemma numbers refer to
rtt/core/certify/PROOFS.md.

Coupling (dimensionless, k >= 1, all energies in kelvin):
    lambda(k; T) = sum_j p_j * 2 Omega_j / (Omega_j^2 + (2 pi k T)^2)
"""
from __future__ import annotations

import math
from fractions import Fraction

from flint import arb, arb_poly, ctx, fmpq

from ..constants import ARB_PREC_BITS

ctx.prec = ARB_PREC_BITS

PI = arb.pi()
TWO_PI = 2 * PI


def set_prec(bits: int = ARB_PREC_BITS) -> None:
    ctx.prec = bits


# --------------------------------------------------------------------------- exact conversions
def q_of_float(x: float) -> fmpq:
    """Exact rational value of a double."""
    n, d = float(x).as_integer_ratio()
    return fmpq(n, d)


def q_of_fraction(f: Fraction) -> fmpq:
    return fmpq(f.numerator, f.denominator)


def arb_exact(x) -> arb:
    """Exact Arb ball of a double / int / fmpq (radius 0 for doubles and ints)."""
    if isinstance(x, float):
        return arb(x)
    if isinstance(x, int):
        return arb(x)
    if isinstance(x, fmpq):
        return arb(x)  # rounded at working precision, ball contains it
    if isinstance(x, Fraction):
        return arb(q_of_fraction(x))
    return arb(x)


def certainly_lt(a: arb, b) -> bool:
    """True only if a < b holds for every point of the balls (rigorous)."""
    return bool(a < b)


def certainly_le(a: arb, b) -> bool:
    return bool(a <= b)


def upper_float(x: arb) -> float:
    """A double >= every point of x (for display / proposals only)."""
    m = float(x.mid())
    r = float(x.rad())
    u = m + r
    return math.nextafter(math.nextafter(u, math.inf), math.inf) if math.isfinite(u) else u


def lower_float(x: arb) -> float:
    m = float(x.mid())
    r = float(x.rad())
    lo = m - r
    return math.nextafter(math.nextafter(lo, -math.inf), -math.inf) if math.isfinite(lo) else lo


# --------------------------------------------------------------------------- spectrum in Arb
SERIES_RMAX = 60            # max number of moment-series terms
SERIES_RATIO = 2.0          # series used when 2 pi k T >= SERIES_RATIO * Omega_max


class ArbSpectrum:
    """Exact Arb data of a Spectrum's nominal atoms (positions and masses are doubles)."""

    def __init__(self, omega_K, weight_K, rel_envelope: float):
        pairs = [(float(o), float(w)) for o, w in zip(omega_K, weight_K) if w > 0]
        if not pairs:
            raise ValueError("spectrum has no positive weight")
        self.n_atoms = len(pairs)
        self.om_f = [o for o, _ in pairs]
        self.p_f = [w for _, w in pairs]
        self.e_f = float(rel_envelope)
        self.e = arb(self.e_f)                       # exact
        self.e_q = q_of_float(self.e_f)
        self.om = [arb(o) for o in self.om_f]        # exact
        self.p = [arb(w) for w in self.p_f]          # exact
        # a_j = 2 p_j Omega_j and o2_j = Omega_j^2 (exact at 128 bits: 53x53 -> 106 bits)
        self.a = [2 * p * o for p, o in zip(self.p, self.om)]
        self.o2 = [o * o for o in self.om]
        self.om_max_f = max(self.om_f)
        self.om_max = arb(self.om_max_f)
        # exact total weight A (rational) and A_hi = (1+e) A
        A_q = fmpq(0)
        for w in self.p_f:
            A_q += q_of_float(w)
        self.A_q = A_q
        self.A_hi_q = (1 + self.e_q) * A_q
        # sum_j p_j Omega_j (for lambda(k) <= b / k^2, Lemma T0)
        s = arb(0)
        for p, o in zip(self.p, self.om):
            s += p * o
        self.sum_pom = s
        # scaled odd moments mu_s = sum_j p_j (Omega_j / Omega_max)^(2s+1), s = 0..RMAX
        self._mu = None

    @property
    def mu(self):
        if self._mu is None:
            mu = [arb(0)] * (SERIES_RMAX + 1)
            for p, o in zip(self.p, self.om):
                t = o / self.om_max
                t2 = t * t
                c = p * t
                for s in range(SERIES_RMAX + 1):
                    mu[s] = mu[s] + c
                    c = c * t2
            self._mu = mu
        return self._mu

    # ------------------------------------------------------------------ lambda(k; T)
    def k_series_start(self, T: float) -> int:
        """Smallest k with 2 pi k T >= SERIES_RATIO * Omega_max (float proposal; the series
        routine re-checks the ratio condition rigorously)."""
        return max(1, int(math.ceil(SERIES_RATIO * self.om_max_f / (2 * math.pi * T) * (1 + 1e-12))) + 1)

    def lam(self, T: float, kmax: int, factor: arb | None = None) -> list:
        """[0, lambda(1;T), ..., lambda(kmax;T)] as balls, nominal masses times ``factor``.

        Direct atom sums for small k; for 2 pi k T >= 2 Omega_max the alternating moment series
        of Lemma S1 with a rigorous remainder.
        """
        T_a = arb(float(T))
        twopiT = TWO_PI * T_a
        ks = min(self.k_series_start(T), kmax + 1)
        out = [arb(0)] * (kmax + 1)
        a, o2 = self.a, self.o2
        for k in range(1, ks):
            x = twopiT * k
            x2 = x * x
            s = arb(0)
            for aj, oj in zip(a, o2):
                s += aj / (oj + x2)
            out[k] = s
        if ks <= kmax:
            # Lemma S1: lambda(k) = (2 Om/x^2) [P_R(r^2) +- mu_R r^(2R)], r = Om/x (r^2 checked < 1)
            mu = self.mu
            om_max = self.om_max
            polys = {}
            for k in range(ks, kmax + 1):
                x = twopiT * k
                r2 = (om_max / x) ** 2
                rf = float(r2.mid())
                if not rf < 0.3:
                    raise ArithmeticError("series used outside r^2 < 0.3")
                # R terms with r2^R <= 2^-110:  r2 <= 1/4 -> 55,  <= 1/64 -> 19, <= 1/1024 -> 11
                R = 55 if rf > 1 / 64 else (19 if rf > 1 / 1024 else 11)
                P = polys.get(R)
                if P is None:
                    P = arb_poly([mu[s] if s % 2 == 0 else -mu[s] for s in range(R)])
                    polys[R] = (P, R)
                    P = polys[R]
                if not bool(r2 < 1):
                    raise ArithmeticError("series ratio not < 1")
                s = P[0](r2)
                rem = mu[P[1]] * r2 ** P[1]
                s = (s - rem).union(s + rem)
                out[k] = (2 * om_max / (x * x)) * s
        if factor is not None:
            out = [factor * v for v in out]
        out[0] = arb(0)
        return out

    def b_const(self, T: float, factor: arb) -> arb:
        """b with lambda(k; T') * factor/(...)  <= b / k^2 for all T' >= T  (Lemma T0):
        b = factor * sum_j p_j Omega_j / (2 pi^2 T^2)."""
        T_a = arb(float(T))
        return factor * self.sum_pom / (2 * PI * PI * T_a * T_a)


# --------------------------------------------------------------------------- poly helpers
def poly_coeffs(p: arb_poly, n: int) -> list:
    """First n coefficients of p (zero-padded)."""
    c = p.coeffs()
    if len(c) < n:
        c = list(c) + [arb(0)] * (n - len(c))
    return c


def tail_f(m: int, p: int = 3) -> arb:
    """f(m) = (2m+1)^-p as a ball."""
    return arb(1) / ((2 * m + 1) ** p)
