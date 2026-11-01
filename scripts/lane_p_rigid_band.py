"""Rigid-band N(E) scan from a pw.x/dos.x DOS file (doc 32 §1.5; R30 §3.3 S4 note). Ranking only.

    ~/.rtt/envs/py/bin/python scripts/lane_p_rigid_band.py <pwscf.dos> <n_atoms> [N_ref=19] [n_complex=1]

dos.x columns: E (eV), DOS (states/eV/cell), integrated DOS. The electron count of the cell at
energy E is the integrated DOS; the per-complex count N(E) = N_ref + (IDOS(E) - IDOS(E_F)) /
n_complex. Prints DOS at N = 18.5 ... 20.5 in steps of 0.25 (states/eV/atom), smoothed over
+-0.05 eV.
"""
import re
import sys

import numpy as np


def scan(path, nat, nref=19.0, ncx=1):
    txt = open(path).read()
    ef = float(re.search(r"EFermi\s*=\s*(-?[\d.]+)", txt).group(1))
    d = np.loadtxt(path)
    E, D, I = d[:, 0], d[:, 1], d[:, 2]
    i_ef = np.interp(ef, E, I)
    N = nref + (I - i_ef) / ncx
    dE = E[1] - E[0]
    w = max(1, int(round(0.05 / dE)))
    Ds = np.convolve(D, np.ones(2 * w + 1) / (2 * w + 1), mode="same")
    rows = []
    for n in np.arange(18.25, 20.76, 0.25):
        j = int(np.argmin(np.abs(N - n)))
        rows.append((round(float(n), 2), round(float(E[j] - ef), 3), round(float(Ds[j]) / nat, 4)))
    return ef, rows


if __name__ == "__main__":
    a = sys.argv
    ef, rows = scan(a[1], int(a[2]), float(a[3]) if len(a) > 3 else 19.0, int(a[4]) if len(a) > 4 else 1)
    print(f"E_F = {ef} eV")
    print("N_per_complex  E-E_F(eV)  DOS(states/eV/atom, 0.1 eV box)")
    for r in rows:
        print(*r, sep="  ")
