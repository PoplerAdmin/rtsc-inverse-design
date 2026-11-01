"""Band-level check of the lane-P electron-count rule (doc 32 §1.3).

    ~/.rtt/envs/py/bin/python scripts/lane_p_band_count.py <pw.x nscf/scf output with eigenvalues> [n_closed_bands]

Reads all k-point eigenvalues of a non-spin-polarised pw.x output (verbosity high), the number of
electrons and the Fermi energy. Reports, for the band index n_c that the 18-electron closed shell
fills (default: n_electrons // 2 rounded down to the closed-shell count passed in), the top of
band n_c, the bottom of band n_c + 1, the gap between them, and where E_F sits. For Mg2IrH6
(43 SG15 electrons: 24 semicore + 19 valence) the closed shell is 42 electrons = 21 bands.
"""
import re
import sys

import numpy as np


def main(path, n_closed=None):
    t = open(path, errors="replace").read()
    nel = float(re.search(r"number of electrons\s+=\s+([\d.]+)", t).group(1))
    ef = float(re.findall(r"the Fermi energy is\s+(-?[\d.]+) ev", t)[-1])
    blocks = re.split(r"\n\s+k =.*?bands \(ev\):\s*\n", t)
    eig = []
    for b in blocks[1:]:
        b = b.split("occupation numbers")[0]
        b = re.split(r"\n\s*\n\s*(?:k =|the Fermi|highest)", b)[0]
        vals = [float(x) for x in re.findall(r"-?\d+\.\d+", b.replace("-", " -"))]
        if vals:
            eig.append(vals)
    nb = min(len(e) for e in eig)
    E = np.array([e[:nb] for e in eig])
    nc = int(n_closed) if n_closed else int(nel // 2)
    top, bot = E[:, nc - 1].max(), E[:, nc].min()
    occ_above = None
    print(f"{path}: n_el={nel} n_k={len(E)} n_bands={nb} E_F={ef:.3f} eV")
    print(f"closed-shell bands n_c={nc}: top(band {nc})={top:.3f}  bottom(band {nc+1})={bot:.3f}  "
          f"gap={bot - top:.3f} eV;  E_F - top = {ef - top:.3f}, E_F - bottom = {ef - bot:.3f}")
    for j in range(max(0, nc - 3), min(nb, nc + 4)):
        print(f"  band {j+1}: [{E[:, j].min():.3f}, {E[:, j].max():.3f}]")
    return E, ef


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2] if len(sys.argv) > 2 else None)
