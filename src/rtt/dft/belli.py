"""Belli networking value phi and H_DOS (architecture S4; doc 11 T1c, §3.3).

Ported (not imported) from ``src/scoring/phi_elf.py``.

phi = the highest ELF value t such that {ELF >= t} has a connected component percolating
through the periodic crystal in all three directions (rank-3 wrap lattice). Belli, Novoa,
Contreras-Garcia, Errea, Nat. Commun. 12, 5381 (2021), doi:10.1038/s41467-021-25687-0 fit
(hydrogen compounds only, scatter ~ +-60 K):

    Tc_phi = 750 * phi * H_f * H_DOS**(1/3) - 85   [K]

**Ranking proxy only; never a stage decision.** Deviations from the paper (doc 11 §3.3):
6-connectivity on the FFT grid and bisection (the paper used Critic2); H_DOS is the fraction of
the integrated LDOS in E_F +- window inside H Voronoi cells, because SG15 pseudopotentials have no
atomic wavefunctions for projwfc.x (UNVERIFIED equivalence to the paper's projected DOS).
For structures without H, H_f = 0 and Tc_phi is None (the fit does not apply).
"""
from __future__ import annotations

from collections import deque
from dataclasses import dataclass
from pathlib import Path

import numpy as np

TC_PHI_SOURCE = "Belli et al., Nat. Commun. 12, 5381 (2021); doc11#3.3"


def _wrap_rank(mask: np.ndarray) -> int:
    """Max over connected components (periodic) of the rank of the component's wrap lattice."""
    from scipy import ndimage
    lab, n = ndimage.label(mask)
    if n == 0:
        return 0
    edges: dict[int, list[tuple[int, np.ndarray]]] = {i: [] for i in range(1, n + 1)}
    for ax in range(3):
        a = np.take(lab, -1, axis=ax)
        b = np.take(lab, 0, axis=ax)
        sel = (a > 0) & (b > 0)
        e = np.zeros(3, int)
        e[ax] = 1
        for u, v in set(zip(a[sel].tolist(), b[sel].tolist())):
            edges[u].append((v, e))
            edges[v].append((u, -e))
    pos: dict[int, np.ndarray] = {}
    best = 0
    for start in range(1, n + 1):
        if start in pos:
            continue
        pos[start] = np.zeros(3, int)
        wraps = []
        dq = deque([start])
        while dq:
            u = dq.popleft()
            for v, e in edges[u]:
                p = pos[u] + e
                if v not in pos:
                    pos[v] = p
                    dq.append(v)
                elif not np.array_equal(pos[v], p):
                    wraps.append(p - pos[v])
        if wraps:
            best = max(best, int(np.linalg.matrix_rank(np.array(wraps))))
            if best == 3:
                return 3
    return best


def networking_value(elf: np.ndarray, tol: float = 1e-4) -> float:
    elf = np.asarray(elf, dtype=float)
    lo, hi = 0.0, float(elf.max())
    if _wrap_rank(elf >= lo) < 3:
        return 0.0
    while hi - lo > tol:
        mid = 0.5 * (lo + hi)
        if _wrap_rank(elf >= mid) == 3:
            lo = mid
        else:
            hi = mid
    return lo


def h_fraction_voronoi(density: np.ndarray, cell: np.ndarray, positions: np.ndarray,
                       symbols: list[str]) -> float:
    """Fraction of a real-space density (grid covering one cell) inside H Voronoi cells
    (nearest atom under the 27 nearest periodic images)."""
    n = np.array(density.shape)
    cell = np.asarray(cell, float)
    g = np.stack(np.meshgrid(*[np.arange(k) / k for k in n], indexing="ij"), -1).reshape(-1, 3)
    cart = g @ cell
    best = np.full(len(g), np.inf)
    owner = np.zeros(len(g), int)
    shifts = np.array([[i, j, k] for i in (-1, 0, 1) for j in (-1, 0, 1) for k in (-1, 0, 1)])
    for ia, p in enumerate(np.asarray(positions, float)):
        for s in shifts:
            d = np.sum((cart - (p + s @ cell)) ** 2, axis=1)
            m = d < best
            best[m] = d[m]
            owner[m] = ia
    isH = np.array([sym == "H" for sym in symbols])
    w = np.clip(np.asarray(density, float).reshape(-1), 0, None)
    tot = w.sum()
    return float(w[isH[owner]].sum() / tot) if tot > 0 else 0.0


def tc_phi(phi: float, h_f: float, h_dos: float) -> float | None:
    if h_f <= 0:
        return None
    return 750.0 * phi * h_f * max(h_dos, 0.0) ** (1.0 / 3.0) - 85.0


@dataclass
class BelliResult:
    phi: float
    H_f: float
    H_DOS: float | None
    Tc_phi_K: float | None
    grid: tuple[int, int, int]


def phi_and_hdos(workdir: str | Path, elf_cube: str = "elf.cube",
                 ildos_cube: str = "ildos.cube") -> BelliResult:
    """phi from ``<workdir>/elf.cube`` (pp.x plot_num=8); H_DOS from ``ildos.cube``
    (pp.x plot_num=10 around E_F) when H is present."""
    from ase.io.cube import read_cube_data
    w = Path(workdir)
    data, atoms = read_cube_data(str(w / elf_cube))
    phi = networking_value(np.asarray(data))
    syms = atoms.get_chemical_symbols()
    hf = syms.count("H") / len(syms)
    hdos = None
    if hf > 0 and (w / ildos_cube).is_file():
        ild, at2 = read_cube_data(str(w / ildos_cube))
        hdos = h_fraction_voronoi(np.asarray(ild), at2.cell.array, at2.get_positions(),
                                  at2.get_chemical_symbols())
    elif hf == 0:
        hdos = 0.0
    tc = tc_phi(phi, hf, hdos) if hdos is not None else None
    return BelliResult(phi=phi, H_f=hf, H_DOS=hdos, Tc_phi_K=tc,
                       grid=tuple(int(x) for x in np.asarray(data).shape))  # type: ignore[arg-type]
