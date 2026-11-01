"""Energy-unit conversions to kelvin, built from the same SI constants Quantum ESPRESSO uses.

Source checked (fetched from GitHub, tags ``qe-6.7.0`` and ``develop``,
``upflib/constants.f90`` / ``Modules/constants.f90``)::

    H_PLANCK_SI     = 6.62607015E-34      J s     (exact, SI 2019)
    K_BOLTZMANN_SI  = 1.380649E-23        J/K     (exact)
    ELECTRONVOLT_SI = 1.602176634E-19     J       (exact)
    HARTREE_SI      = 4.3597447222071E-18 J       (CODATA 2018)
    C_SI            = 2.99792458E+8       m/s     (exact)
    RYDBERG_SI      = HARTREE_SI / 2
    K_BOLTZMANN_RY  = K_BOLTZMANN_SI / RYDBERG_SI      (used by matdyn for omega_ln [K])
    RY_TO_THZ       = 1 / AU_TERAHERTZ / FPI  ( = Ry / h in THz, ordinary frequency nu)
    RY_TO_CMM1      = 1E10 * RY_TO_THZ / C_SI

Using QE's own constants means omega_ln [K] printed by matdyn.x (``E / K_BOLTZMANN_RY``) and our
kelvin conversion of the a2F.dos frequency axis (Rydberg) use the identical factor.

Conventions: "THz" is ordinary frequency nu (E = h nu), as in QE's RY_TO_THZ; "cm-1" is
wavenumber (E = h c nu~). Every factor converts an ENERGY in that unit to kelvin (E / k_B).
"""
from __future__ import annotations

import math

H_PLANCK_SI = 6.62607015e-34
K_BOLTZMANN_SI = 1.380649e-23
ELECTRONVOLT_SI = 1.602176634e-19
HARTREE_SI = 4.3597447222071e-18
C_SI = 2.99792458e8
RYDBERG_SI = HARTREE_SI / 2.0

RY_TO_K = RYDBERG_SI / K_BOLTZMANN_SI                 # 157887.512402 K
EV_TO_K = ELECTRONVOLT_SI / K_BOLTZMANN_SI            # 11604.5181216 K
MEV_TO_K = 1e-3 * EV_TO_K                             # 11.6045181216 K
THZ_TO_K = H_PLANCK_SI * 1e12 / K_BOLTZMANN_SI        # 47.9924307337 K
CMM1_TO_K = H_PLANCK_SI * C_SI * 100.0 / K_BOLTZMANN_SI   # 1.43877687750 K
_AU_SEC = H_PLANCK_SI / (2.0 * math.pi) / HARTREE_SI
RY_TO_THZ = 1.0 / (_AU_SEC * 1e12) / (4.0 * math.pi)  # 3289.84196025 THz (QE RY_TO_THZ)
RY_TO_CMM1 = 1e10 * RY_TO_THZ / C_SI                  # 109737.315682 cm^-1 (QE RY_TO_CMM1)
RY_TO_GHZ = RY_TO_THZ * 1000.0

#: canonical unit name -> factor to kelvin
TO_KELVIN: dict[str, float] = {
    "K": 1.0,
    "meV": MEV_TO_K,
    "eV": EV_TO_K,
    "THz": THZ_TO_K,
    "cm-1": CMM1_TO_K,
    "Ry": RY_TO_K,
    "mRy": 1e-3 * RY_TO_K,
    "Ha": 2.0 * RY_TO_K,
}

_ALIASES = {
    "k": "K", "kelvin": "K",
    "mev": "meV",
    "ev": "eV",
    "thz": "THz",
    "cm-1": "cm-1", "cm^-1": "cm-1", "1/cm": "cm-1", "cm**-1": "cm-1", "invcm": "cm-1",
    "ry": "Ry", "rydberg": "Ry",
    "mry": "mRy",
    "ha": "Ha", "hartree": "Ha",
}


def canonical_unit(unit: str) -> str:
    u = unit.strip()
    if u in TO_KELVIN:
        return u
    c = _ALIASES.get(u.lower())
    if c is None:
        raise ValueError(f"unknown energy unit {unit!r}; known: {sorted(TO_KELVIN)} "
                         f"(aliases {sorted(_ALIASES)})")
    return c


def to_kelvin_factor(unit: str) -> float:
    """Multiplicative factor converting an energy expressed in ``unit`` to kelvin."""
    return TO_KELVIN[canonical_unit(unit)]


GPA_TO_KBAR = 10.0
