"""Unit conversions: re-export of the core ``rtt.core.dfpt.units`` (energies -> kelvin)."""
from .core.dfpt.units import (CMM1_TO_K, EV_TO_K, MEV_TO_K, RY_TO_K, THZ_TO_K,  # noqa: F401
                                      canonical_unit, to_kelvin_factor)

K_PER_MEV = MEV_TO_K
K_PER_RY = RY_TO_K
K_PER_CM1 = CMM1_TO_K


def to_K(value: float, unit: str) -> float:
    return float(value) * to_kelvin_factor(unit)
