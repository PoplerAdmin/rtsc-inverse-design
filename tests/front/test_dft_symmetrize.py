"""Campaign 1: noisy generated cells must reach QE in their labelled space group."""
import numpy as np
import spglib

from rtt.structure.canonical import symmetrized_for_dft


def _perov(noise):
    rng = np.random.default_rng(0)
    fc = np.array([[0, 0, 0], [0.5, 0.5, 0.5], [0.5, 0.5, 0], [0.5, 0, 0.5], [0, 0.5, 0.5]], float)
    fc = (fc + noise * rng.standard_normal(fc.shape)) % 1.0
    lat = np.eye(3) * 4.0 + noise * rng.standard_normal((3, 3))
    return {"lattice": lat.tolist(), "species": ["Ba", "Pt", "H", "H", "H"], "frac_coords": fc.tolist()}


def _sg(d, sp):
    els = sorted(set(d["species"]))
    return spglib.get_symmetry_dataset((d["lattice"], d["frac_coords"], [els.index(x) + 1 for x in d["species"]]),
                                       symprec=sp).number


def test_noisy_cell_is_idealized():
    d = _perov(3e-4)
    assert _sg(d, 1e-5) == 1 and _sg(d, 0.01) == 221
    out, info = symmetrized_for_dft(d)
    assert info["symmetrized"] and _sg(out, 1e-6) == 221 and len(out["species"]) == 5


def test_exact_cell_unchanged():
    d = _perov(0.0)
    out, info = symmetrized_for_dft(d)
    assert not info["symmetrized"] and out["frac_coords"] == d["frac_coords"]
