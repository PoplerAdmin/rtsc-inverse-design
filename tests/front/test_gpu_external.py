"""Tests that need real MatterSim / the MatterGen venv. Skipped unless enabled:

    RTT_GPU_TESTS=1        (CUDA + MatterSim; obey the shared-GPU rule: < 4000 MiB in use)
    RTT_EXTERNAL_TESTS=1   (~/mattergen/.venv + Alex-MP LMDB)
"""
from __future__ import annotations

import pytest


@pytest.mark.gpu
def test_mattersim_relax_and_phonons_al():
    from rtt.mlff.phonons import ml_phonons
    from rtt.mlff.relax import relax_batch
    from rtt.source.controls import control
    sd = control("Al").raw().structure.to_json()
    r = relax_batch([sd, sd], device="cuda")
    assert len(r) == 2 and all(x.converged for x in r)
    assert 15.0 < r[0].volume_per_atom_A3 < 18.5            # fcc Al ~16.5 A^3/atom
    ph = ml_phonons(r[0].structure, device="cuda")
    assert ph.min_freq_cm1 > -15.0 and 250 < ph.omega_max_cm1 < 450


@pytest.mark.external
def test_refquery_and_database_lane():
    from rtt.source.database import query_alexmp
    rs = query_alexmp("H-Ir-Mg")
    forms = {r.meta["formula"] for r in rs}
    assert "Mg2H6Ir" in forms or "Mg2IrH6" in forms
    assert all(r.origin.lane == "D" and "GAP_UNKNOWN" in r.origin.flags for r in rs)
    assert all(r.meta["e_above_hull_dft_eV"] <= 0.1 for r in rs)


@pytest.mark.external
def test_dftref_hull_roundtrip():
    from rtt.mlff.hull import ehull_dftref
    from rtt.source.controls import control
    sd = control("Mg2IrH6").raw().structure.to_json()
    # DFT energy of the reference entry itself (-34.304 eV / 9 atoms): MatterSim-free check of the bridge
    n = len(sd["species"])                                # conventional cell: 36 atoms
    res = ehull_dftref([("x", sd, -34.30429345 * n / 9)])
    print(res)
    assert res["x"]["e_above_hull_eV"] is not None and abs(res["x"]["e_above_hull_eV"]) < 0.25
