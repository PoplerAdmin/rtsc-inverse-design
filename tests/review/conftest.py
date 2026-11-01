"""Fixtures for the RC1 code-review tests (docs/reviews/RC1_code_review.md). CPU only, no QE,
no MatterSim. Each test in this directory demonstrates one finding and is EXPECTED TO FAIL until
the finding is fixed."""
from __future__ import annotations

import os
from pathlib import Path

import pytest

PSEUDO_ELEMENTS = ("H", "Li", "Be", "B", "C", "N", "O", "Na", "Mg", "Al", "S", "Ca", "Nb", "Pd", "Ir",
                   "Pt", "Pb", "La", "Y")


@pytest.fixture(autouse=True)
def _env(tmp_path, monkeypatch):
    monkeypatch.setenv("RTT_FRONT_CACHE", str(tmp_path / "front_cache"))
    monkeypatch.setenv("OMP_NUM_THREADS", "2")


@pytest.fixture
def pseudo_dir(tmp_path) -> Path:
    d = tmp_path / "sg15"
    d.mkdir()
    for el in PSEUDO_ELEMENTS:
        (d / f"{el}_ONCV_PBE-1.2.upf").write_text("fake")
    return d


@pytest.fixture
def targets():
    from rtt.config import load_targets
    return load_targets()


@pytest.fixture
def tiers(tmp_path, pseudo_dir):
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
