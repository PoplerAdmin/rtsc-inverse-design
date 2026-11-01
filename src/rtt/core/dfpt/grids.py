"""k- and q-grid selection for the interpolated electron-phonon workflow.

Density convention: ``spacing`` in 1/angstrom INCLUDING the 2*pi factor (VASP KSPACING
convention): n_i = max(1, ceil(|b_i| / spacing)).

Constraints enforced (from QE ``PHonon/PH/elphon.f90``, SUBROUTINE elphsum, checked in
``develop`` and ``qe-6.7.0``):
  * every q of the ph.x grid must be a vector of the dense ("fit") k grid, otherwise
    ``errore('elphsum','q is not a vector in the dense grid')`` -> nk_fit_i % nq_i == 0;
  * the dense grid is mapped onto the coarse scf grid with integer ratios
    ``nti = nk1fit/nk1`` -> nk_fit_i % nk_i == 0;
  * both grids are generated unshifted (``lint(..., 0,0,0, ...)``) -> no k-point offsets.
We additionally require nk_i % nq_i == 0 (k+q stays on the coarse grid), the usual practice.
"""
from __future__ import annotations

import math
from dataclasses import dataclass

import numpy as np


def grid_from_spacing(recip_lengths_2pi: np.ndarray, spacing: float) -> tuple[int, int, int]:
    if spacing <= 0:
        raise ValueError("spacing must be > 0")
    return tuple(max(1, int(math.ceil(float(b) / spacing - 1e-9))) for b in recip_lengths_2pi)  # type: ignore[return-value]


def _ceil_multiple(n: int, m: int) -> int:
    return int(math.ceil(n / m)) * m


@dataclass(frozen=True)
class EPhGrids:
    k_scf: tuple[int, int, int]
    k_fit: tuple[int, int, int]
    q: tuple[int, int, int]
    k_matdyn: tuple[int, int, int]   # q-interpolation grid used by matdyn.x for a2F / DOS

    def validate(self) -> None:
        for name, g in (("k_scf", self.k_scf), ("k_fit", self.k_fit), ("q", self.q),
                        ("k_matdyn", self.k_matdyn)):
            if len(g) != 3 or any(int(x) < 1 for x in g):
                raise ValueError(f"{name} grid must be three positive integers, got {g}")
        for i in range(3):
            if self.k_fit[i] % self.q[i]:
                raise ValueError(f"dense k grid {self.k_fit} must be a multiple of q grid {self.q} "
                                 "(QE elphsum: 'q is not a vector in the dense grid')")
            if self.k_fit[i] % self.k_scf[i]:
                raise ValueError(f"dense k grid {self.k_fit} must be a multiple of scf k grid "
                                 f"{self.k_scf} (QE elphsum: nti = nk1fit/nk1)")
            if self.k_scf[i] % self.q[i]:
                raise ValueError(f"scf k grid {self.k_scf} must be a multiple of q grid {self.q}")

    def to_dict(self) -> dict:
        return {"k_scf": list(self.k_scf), "k_fit": list(self.k_fit), "q": list(self.q),
                "k_matdyn": list(self.k_matdyn)}


def select_grids(recip_lengths_2pi: np.ndarray, *, k_spacing: float = 0.25,
                 q_spacing: float = 0.8, fit_factor: int = 2,
                 k_fit_spacing: float | None = None,
                 matdyn_spacing: float | None = None,
                 k_scf: tuple[int, int, int] | None = None,
                 k_fit: tuple[int, int, int] | None = None,
                 q: tuple[int, int, int] | None = None,
                 k_matdyn: tuple[int, int, int] | None = None) -> EPhGrids:
    """Choose commensurate (k_scf, k_fit, q, k_matdyn) grids.

    Explicit grids override the spacing rules; the result is validated either way (an
    explicit inconsistent choice raises instead of being silently changed).
    """
    b = np.asarray(recip_lengths_2pi, dtype=float)
    qq = tuple(q) if q is not None else grid_from_spacing(b, q_spacing)
    if k_scf is not None:
        kk = tuple(k_scf)
    else:
        k0 = grid_from_spacing(b, k_spacing)
        kk = tuple(_ceil_multiple(k0[i], qq[i]) for i in range(3))
    if k_fit is not None:
        kf = tuple(k_fit)
    else:
        if k_fit_spacing is not None:
            f0 = grid_from_spacing(b, k_fit_spacing)
        else:
            f0 = tuple(kk[i] * max(1, int(fit_factor)) for i in range(3))
        kf = tuple(_ceil_multiple(max(f0[i], kk[i]), kk[i]) for i in range(3))
    if k_matdyn is not None:
        km = tuple(k_matdyn)
    elif matdyn_spacing is not None:
        km = grid_from_spacing(b, matdyn_spacing)
    else:
        km = tuple(4 * x for x in qq)
    g = EPhGrids(tuple(int(x) for x in kk), tuple(int(x) for x in kf),  # type: ignore[arg-type]
                 tuple(int(x) for x in qq), tuple(int(x) for x in km))  # type: ignore[arg-type]
    g.validate()
    return g
