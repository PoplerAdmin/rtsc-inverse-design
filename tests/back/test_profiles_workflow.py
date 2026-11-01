"""rtt.dft.profiles and the generic Chain runner of rtt.dft.workflow (fake QE executables)."""
from __future__ import annotations

import json
import stat
from pathlib import Path

import pytest

from rtt.dft import profiles as PR
from rtt.dft import workflow as W


def test_profiles_present_and_dfpt_settings():
    assert set(PR.PROFILES) == {"s4", "coarse", "full", "ref_hull"}
    c = PR.PROFILES["coarse"]
    assert c.spectrum_tier == "coarse" and c.smearing_admission == "mode_positive"
    assert c.ecutrho_Ry == 4 * c.ecutwfc_Ry
    s = c.eph_settings(q=(3, 3, 3))
    assert s.q == (3, 3, 3) and s.k_matdyn == (12, 12, 12)
    assert s.smearing_admission == "mode_positive" and s.el_ph_nsigma == 10
    with pytest.raises(ValueError):
        PR.PROFILES["s4"].eph_settings()


def test_coarse_grids_h3s_match_doc11():
    """H3S 200 GPa primitive bcc cell -> k_scf 12^3, k_fit 24^3 at q 4^3 (doc 11 §3.1 c3)."""
    from rtt.core.dfpt.grids import select_grids
    from rtt.core.dfpt.structure import CellData
    a = 1.494124044
    cell = CellData([[-a, a, a], [a, -a, a], [a, a, -a]], ["S", "H", "H", "H"],
                    [[0, 0, 0], [0, .5, .5], [.5, 0, .5], [.5, .5, 0]], {"S": 32.06, "H": 1.008})
    s = PR.PROFILES["coarse"].eph_settings(q=(4, 4, 4))
    g = select_grids(cell.reciprocal_lengths_2pi(), k_spacing=s.k_spacing, q_spacing=s.q_spacing,
                     fit_factor=s.fit_factor, k_fit_spacing=s.k_fit_spacing, q=s.q)
    assert g.k_scf == (12, 12, 12) and g.k_fit == (24, 24, 24)


def test_profiles_from_tiers_overrides():
    tiers = {"S5": {"q": [3, 3, 3], "plateau_sigma_min_Ry": 0.04},
             "S6": {"q_meshes": [[3, 3, 3], [5, 5, 5]]},
             "profiles": {"s4": {"ecutwfc_Ry": 50.0, "source": "test"}}}
    P = PR.profiles_from_tiers(tiers)
    assert P["coarse"].q == (3, 3, 3) and P["coarse"].plateau_sigma_min_Ry == 0.04
    assert P["full"].q_meshes == ((3, 3, 3), (5, 5, 5)) and P["full"].q == (5, 5, 5)
    assert P["s4"].ecutwfc_Ry == 50.0 and ("ecutwfc_Ry", 50.0) in P["s4"].overrides
    assert P["s4"].sha256() != PR.PROFILES["s4"].sha256()
    with pytest.raises(ValueError):
        PR.profiles_from_tiers({"profiles": {"nope": {}}})
    with pytest.raises(ValueError):
        PR.PROFILES["s4"].with_overrides(bogus=1)


def test_real_tiers_file_loads():
    from rtt.config import load_tiers
    P = PR.profiles_from_tiers(load_tiers())
    assert P["coarse"].q == (4, 4, 4) and P["coarse"].plateau_sigma_min_Ry == 0.03


def _upf(path: Path, ptype: str, rel: str = "scalar") -> None:
    path.write_text(f'<UPF version="2.0.1">\n<PP_HEADER\n element="X"\n pseudo_type="{ptype}"\n'
                    f' relativistic="{rel}"\n functional="PBE"\n/>\n')


def test_nc_only_pseudos(tmp_path):
    _upf(tmp_path / "Al_ONCV_PBE-1.2.upf", "NC")
    _upf(tmp_path / "Mg_ONCV_PBE-1.2.upf", "US")
    _upf(tmp_path / "Pb_ONCV_PBE-1.2.upf", "NC", "full")
    t = PR.sg15_pseudos(["Al"], tmp_path)
    assert t.files == {"Al": "Al_ONCV_PBE-1.2.upf"}
    with pytest.raises(PR.PseudoError, match="norm-conserving"):
        PR.sg15_pseudos(["Mg"], tmp_path)
    with pytest.raises(PR.PseudoError, match="relativistic"):
        PR.sg15_pseudos(["Pb"], tmp_path)
    with pytest.raises(PR.PseudoError, match="no SG15"):
        PR.sg15_pseudos(["Zr"], tmp_path)
    assert PR.available_elements(tmp_path) == {"Al", "Mg", "Pb"}


def test_wait_for_ranks():
    seq = iter([11, 11, 10, 10])
    waited = W.wait_for_ranks(2, 12, poll_s=0.01, counter=lambda: next(seq))
    assert waited > 0
    with pytest.raises(TimeoutError):
        W.wait_for_ranks(2, 12, poll_s=0.01, timeout_s=0.03, counter=lambda: 12)


@pytest.fixture
def fake_bin(tmp_path):
    b = tmp_path / "bin"
    b.mkdir()
    exe = b / "fake.x"
    exe.write_text("#!/bin/sh\necho \"ran $*\" >> ../calls.log\necho 'convergence has been achieved'\n"
                   "echo 'JOB DONE.'\n")
    exe.chmod(exe.stat().st_mode | stat.S_IEXEC)
    from rtt.core.dfpt.launcher import Launcher
    return Launcher(kind="plain", bin_dir=str(b), nice=19)


def _chain(w: Path, launcher, texts: dict) -> W.Chain:
    def chk(name):
        def f(wd):
            o = wd / f"{name}.out"
            return (o.exists() and "JOB DONE" in o.read_text()), "x"
        return f
    steps = [W.ChainStep("a", "fake.x", lambda: texts["a"], chk("a")),
             W.ChainStep("b", "fake.x", lambda: texts["b"], chk("b"), upstream=("a",))]
    return W.Chain(w, steps, launcher, material_id="m", profile="t",
                   policy=W.RunPolicy(max_machine_ranks=None))


def test_chain_tokens_adoption_and_timing(tmp_path, fake_bin):
    w = tmp_path / "run"
    texts = {"a": "input a\n", "b": "input b\n"}
    ok, err = _chain(w, fake_bin, texts).run()
    assert ok, err
    log = lambda: (tmp_path / "calls.log").read_text().count("ran")  # noqa: E731
    assert log() == 2
    c2 = _chain(w, fake_bin, texts)
    assert c2.run()[0] and log() == 2
    assert [t["state"] for t in c2.timings] == ["done (state)", "done (state)"]
    # change the first input: both steps rerun (token chain invalidates downstream)
    texts["a"] = "input a v2\n"
    assert _chain(w, fake_bin, texts).run()[0] and log() == 4
    # lost state file: identical inputs are adopted, not rerun
    (w / "chain_state.json").unlink()
    c4 = _chain(w, fake_bin, texts)
    assert c4.run()[0] and log() == 4
    assert [t["state"] for t in c4.timings] == ["adopted", "adopted"]
    rows = [json.loads(x) for x in (w / "timing.jsonl").read_text().splitlines()]
    ran = [r for r in rows if r["state"] == "ran"]
    assert len(ran) == 4 and all(r["wall_s"] >= 0 and "core_h" in r for r in ran)


def test_chain_reports_failure(tmp_path, fake_bin):
    bad = Path(fake_bin.bin_dir) / "bad.x"
    bad.write_text("#!/bin/sh\necho oops\nexit 3\n")
    bad.chmod(bad.stat().st_mode | stat.S_IEXEC)
    steps = [W.ChainStep("a", "bad.x", lambda: "x\n", lambda wd: (False, "no"))]
    ok, err = W.Chain(tmp_path / "r", steps, fake_bin, material_id="m", profile="t",
                      policy=W.RunPolicy(max_machine_ranks=None)).run()
    assert not ok and "step a failed" in err


def test_disk_guard_blocks_step(tmp_path, fake_bin):
    steps = [W.ChainStep("a", "fake.x", lambda: "x\n", lambda wd: (False, "no"))]
    ok, err = W.Chain(tmp_path / "r", steps, fake_bin, material_id="m", profile="t",
                      policy=W.RunPolicy(max_machine_ranks=None, min_free_gb=1e12)).run()
    assert not ok and "not started" in err


def test_dfptrun_roundtrip():
    r = W.DFPTRun(material_id="x", profile=PR.PROFILES["coarse"].to_dict(), workdir="/tmp/x",
                  status="done", steps=[{"wall_s": 3.0, "cpu_s": 7200.0}])
    d = r.to_dict()
    assert d["core_h"] == 2.0 and W.DFPTRun.from_dict(d).spectrum_tier == "coarse"
