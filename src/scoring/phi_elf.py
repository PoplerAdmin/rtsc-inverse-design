"""Networking value phi (Belli, Novoa, Contreras-Garcia, Errea, Nat. Commun. 12, 5381 (2021)).

phi = the highest ELF value t such that the region {ELF >= t} contains a connected component
that percolates through the periodic crystal in all three directions (i.e. the component's set of
wrap vectors has rank 3).  The paper's fit (hydrogen-based compounds only, ~135 structures,
literature Tc, scatter ~ +/-60 K):

    Tc_Belli = (750 * phi * H_f * H_DOS**(1/3) - 85) K,
    H_f = N_H / N_atoms,  H_DOS = fraction of N(E_F) projected on H.

TIER-1 PHYSICS PROXY: ordering only, never a verdict.  This implementation is ours (6-connectivity
on the FFT grid, bisection on t); the paper used Critic2.  Numbers from it are NOT validated
against the paper's table (UNVERIFIED agreement).

Input: an ELF cube file covering one full unit cell, e.g. from QE pp.x
    &inputpp plot_num = 8 /  &plot iflag = 3, output_format = 6, fileout = 'elf.cube' /
Usage: python phi_elf.py elf.cube [--hdos 0.5 | --ildos-cube ildos.cube]
"""
from __future__ import annotations

import argparse
from collections import deque

import numpy as np
from scipy import ndimage


def _wrap_rank(mask: np.ndarray) -> int:
    """Max over connected components (periodic) of the rank of the component's wrap lattice."""
    lab, n = ndimage.label(mask)  # 6-connectivity inside the box
    if n == 0:
        return 0
    edges: dict[int, list[tuple[int, np.ndarray]]] = {i: [] for i in range(1, n + 1)}
    for ax in range(3):
        a = np.take(lab, -1, axis=ax)  # last layer
        b = np.take(lab, 0, axis=ax)   # first layer (periodic neighbour, offset +e_ax)
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


def h_fraction_voronoi(ildos: np.ndarray, atoms) -> float:
    """Fraction of a real-space density (e.g. ILDOS in a window around E_F) inside H Voronoi cells.

    Substitute for the projwfc H_DOS when the pseudopotentials carry no atomic wavefunctions
    (SG15 ONCV files have none).  Voronoi partition, nearest atom under periodic images.
    This is NOT the same partition as Belli et al.'s projected DOS (UNVERIFIED equivalence).
    """
    n = np.array(ildos.shape)
    g = np.stack(np.meshgrid(*[np.arange(k) / k for k in n], indexing="ij"), -1).reshape(-1, 3)
    cart = g @ atoms.cell.array
    best = np.full(len(g), np.inf)
    owner = np.zeros(len(g), int)
    shifts = np.array([[i, j, k] for i in (-1, 0, 1) for j in (-1, 0, 1) for k in (-1, 0, 1)])
    for ia, p in enumerate(atoms.get_positions()):
        for s in shifts:
            d = np.sum((cart - (p + s @ atoms.cell.array)) ** 2, axis=1)
            m = d < best
            best[m] = d[m]
            owner[m] = ia
    isH = np.array([sym == "H" for sym in atoms.get_chemical_symbols()])
    w = np.clip(ildos.reshape(-1), 0, None)
    return float(w[isH[owner]].sum() / w.sum())


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("cube")
    ap.add_argument("--hdos", type=float, default=None, help="H fraction of N(E_F) (from projwfc)")
    ap.add_argument("--ildos-cube", default=None,
                    help="ILDOS cube around E_F (pp.x plot_num=10); H_DOS from H Voronoi cells")
    a = ap.parse_args()
    from ase.io.cube import read_cube_data
    data, atoms = read_cube_data(a.cube)
    phi = networking_value(np.asarray(data))
    syms = atoms.get_chemical_symbols()
    hf = syms.count("H") / len(syms)
    out = f"phi={phi:.4f} H_f={hf:.4f} grid={data.shape}"
    if a.hdos is None and a.ildos_cube:
        ild, at2 = read_cube_data(a.ildos_cube)
        a.hdos = h_fraction_voronoi(np.asarray(ild), at2)
    if a.hdos is not None:
        out += f" H_DOS={a.hdos:.4f} Tc_Belli={750 * phi * hf * a.hdos ** (1 / 3) - 85:.1f} K (+/-60 K, ordering only)"
    print(out)


if __name__ == "__main__":
    main()
