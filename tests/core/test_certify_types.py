"""Tests for the shared types owned by team CERT: rtt.core.spectrum, rtt.core.status, rtt.core.constants."""
from __future__ import annotations

import math
import sys
from fractions import Fraction
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from rtt.core.constants import A0, A0_FRACTION, PROVENANCE  # noqa: E402
from rtt.core.spectrum import Spectrum, unit_factor  # noqa: E402
from rtt.core.status import Status, Verdict  # noqa: E402


def test_spectrum_validation():
    with pytest.raises(ValueError):
        Spectrum((0.0,), (1.0,))
    with pytest.raises(ValueError):
        Spectrum((1.0,), (-1.0,))
    with pytest.raises(ValueError):
        Spectrum((1.0,), (1.0,), rel_envelope=1.0)     # ERRATA D5: keep e < 1
    with pytest.raises(ValueError):
        Spectrum((1.0, 2.0), (1.0,))


def test_spectrum_helpers_einstein():
    s = Spectrum((100.0,), (50.0,), label="synthetic")
    assert s.A() == 50.0
    assert s.lambda_() == pytest.approx(1.0)
    assert s.omega_log() == pytest.approx(100.0)
    assert s.A_hi() == 50.0


def test_spectrum_json_exact_roundtrip():
    s = Spectrum((1 / 3, math.pi), (0.1, 1e-300), 0.07, "synthetic", {"k": 1})
    t = Spectrum.from_json(s.to_json())
    assert t == s and t.omega_K == s.omega_K and t.weight_K == s.weight_K
    hash(s)


def test_from_grid_units_and_mass():
    # constant a2F = 0.5 on [1, 10] meV, midpoint weights
    om = [1 + 0.01 * i for i in range(901)]
    s = Spectrum.from_grid(om, [0.5] * len(om), unit="meV")
    assert unit_factor("meV") == pytest.approx(11.604518, rel=1e-6)
    # total mass ~ 0.5 * (10.005 - 0.995) meV in K
    assert s.A() == pytest.approx(0.5 * 9.01 * unit_factor("meV"), rel=1e-9)
    assert unit_factor("cm-1") == pytest.approx(1.4387769, rel=1e-6)
    assert unit_factor("THz") == pytest.approx(47.992430, rel=1e-6)


def test_verdict_roundtrip():
    v = Verdict("m1", Status.UNDECIDED, 1.0, 2.0, ["r"], ["R-UNIV"], ["a"], ["f"], {"x": 1})
    w = Verdict.from_json(v.to_json())
    assert w == v


def test_constants():
    assert A0_FRACTION == Fraction(433148145, 10**8)
    assert A0 == 4.33148145
    assert "A0" in PROVENANCE
