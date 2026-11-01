"""Fixtures for the front-end tests (S1-S3). CPU only unless marked ``gpu``/``external``.

Markers
* ``gpu``: needs CUDA + MatterSim; runs only with ``RTT_GPU_TESTS=1`` (shared GPU etiquette).
* ``external``: needs the MatterGen venv and the Alex-MP LMDB; runs only with ``RTT_EXTERNAL_TESTS=1``.
"""
from __future__ import annotations

import os
from pathlib import Path

import pytest

PSEUDO_ELEMENTS = ("H", "Li", "Be", "B", "C", "N", "O", "Na", "Mg", "Al", "S", "Ca", "Nb", "Pd", "Ir",
                   "Pt", "Pb", "La", "Y", "Xe")


def pytest_configure(config):
    config.addinivalue_line("markers", "gpu: needs CUDA + MatterSim (RTT_GPU_TESTS=1)")
    config.addinivalue_line("markers", "external: needs MatterGen venv + Alex-MP LMDB (RTT_EXTERNAL_TESTS=1)")


def pytest_collection_modifyitems(config, items):
    for it in items:
        if "gpu" in it.keywords and os.environ.get("RTT_GPU_TESTS") != "1":
            it.add_marker(pytest.mark.skip(reason="GPU test: set RTT_GPU_TESTS=1"))
        if "external" in it.keywords and os.environ.get("RTT_EXTERNAL_TESTS") != "1":
            it.add_marker(pytest.mark.skip(reason="external test: set RTT_EXTERNAL_TESTS=1"))


@pytest.fixture(autouse=True)
def _front_cache(tmp_path, monkeypatch):
    if os.environ.get("RTT_EXTERNAL_TESTS") != "1":
        monkeypatch.setenv("RTT_FRONT_CACHE", str(tmp_path / "front_cache"))
    monkeypatch.setenv("OMP_NUM_THREADS", "2")


@pytest.fixture
def pseudo_dir(tmp_path) -> Path:
    d = tmp_path / "sg15"
    d.mkdir()
    for el in PSEUDO_ELEMENTS:
        (d / f"{el}_ONCV_PBE-1.2.upf").write_text("fake")
        (d / f"{el}_ONCV_PBE_FR-1.0.upf").write_text("fake")
    return d


@pytest.fixture
def targets():
    from rtt.config import load_targets
    return load_targets()


@pytest.fixture
def tiers(tmp_path, pseudo_dir):
    """Real tiers.toml with ``pseudo_dir`` pointed at the fake SG15 directory."""
    from rtt.config import CONFIG_DIR, load_tiers
    txt = (CONFIG_DIR / "tiers.toml").read_text()
    lines = [(f'pseudo_dir = "{pseudo_dir.as_posix()}"' if ln.strip().startswith("pseudo_dir") else ln)
             for ln in txt.splitlines()]
    p = tmp_path / "tiers.toml"
    p.write_text("\n".join(lines) + "\n")
    return load_tiers(p)


@pytest.fixture
def ledger(targets, tiers):
    from rtt.config import load_machine
    from rtt.ledger.store import Ledger
    L = Ledger(":memory:")
    L.register_configs(targets, tiers, load_machine())
    yield L
    L.close()
