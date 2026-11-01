"""S1 lanes and the MatterGen adapter (no GPU: a fake mattergen-generate is used)."""
from __future__ import annotations

import json
import os
import stat
import textwrap

import pytest

from rtt.source.lanes import (CHECKPOINT_PROPERTIES, DEFAULT_SHARES, LaneSpec, UnsupportedConditioningError,
                              default_chemsys_list, plan_day, plan_summary)
from rtt.source.mattergen import build_argv, generate


# ------------------------------------------------------------------ LaneSpec refuses
@pytest.mark.parametrize("ckpt,cond", [
    ("chemical_system", {"chemical_system": "H-Mg", "energy_above_hull": 0.1}),   # extra property
    ("chemical_system_energy_above_hull", {"chemical_system": "H-Mg"}),            # missing property
    ("space_group", {"dft_band_gap": 0.0}),
    ("mattergen_base", {"space_group": 225}),
    ("not_a_checkpoint", {}),
    ("chemical_system", {"chemical_system": "H-Mg-Xe"}),                          # Xe outside domain
    ("chemical_system", {"chemical_system": "Mg-H"}),                             # not sorted
    ("space_group", {"space_group": 231}),
])
def test_lanespec_refuses_unsupported(ckpt, cond):
    with pytest.raises(UnsupportedConditioningError):
        LaneSpec("G1", ckpt, cond)


def test_lanespec_adapter_refused_unless_explicit():
    with pytest.raises(UnsupportedConditioningError):
        LaneSpec("GA", None, {"rtt_A_K": 1000.0}, model_path="/x")
    s = LaneSpec("GA", None, {"rtt_A_K": 1000.0}, model_path="/x", allow_adapter=True)
    assert s.conditioned


def test_every_shipped_checkpoint_accepts_exactly_its_set():
    vals = {"chemical_system": "H-Mg", "energy_above_hull": 0.1, "space_group": 225, "dft_band_gap": 0.0,
            "dft_mag_density": 0.0, "ml_bulk_modulus": 300.0, "hhi_score": 1000.0}
    for ck, props in CHECKPOINT_PROPERTIES.items():
        s = LaneSpec("G5", ck, {p: vals[p] for p in props})
        assert s.effective_guidance == (2.0 if props else 0.0)


def test_plan_day_shares_and_specs(targets):
    jobs = plan_day(None, targets, gpu_h=15.0, day="d1", seed=3)
    summ = plan_summary(jobs)
    assert set(summ) == set(DEFAULT_SHARES)
    tot_gpu = sum(v["gpu_h_est"] for v in summ.values())
    assert 13.0 < tot_gpu <= 15.0 + 1e-9
    assert abs(summ["G1"]["gpu_h_est"] / tot_gpu - 0.55) < 0.05
    for j in jobs:
        assert j.n == 64 and j.spec.checkpoint in CHECKPOINT_PROPERTIES
        if j.spec.lane == "G1":
            eh = j.spec.conditioning["energy_above_hull"]
            assert 0.0 <= eh <= 0.3 and "H" in j.spec.conditioning["chemical_system"].split("-")
        if j.spec.lane == "G2":
            assert j.spec.conditioning["space_group"] in (191, 194, 139, 123, 166, 221, 223, 225, 227, 229)
        if j.spec.lane == "G4" and j.spec.checkpoint == "chemical_system":
            assert "H" not in j.spec.conditioning["chemical_system"].split("-")
    # deterministic
    assert [j.batch_id for j in plan_day(None, targets, gpu_h=15.0, day="d1", seed=3)] == [j.batch_id for j in jobs]


def test_plan_day_weights_floor(targets):
    systems = default_chemsys_list()
    w = {s: 0.0 for s in systems}
    w[systems[0]] = 100.0
    jobs = [j for j in plan_day(None, targets, weights=w, gpu_h=60.0, day="w") if j.spec.lane == "G1"]
    chosen = {j.spec.conditioning["chemical_system"] for j in jobs}
    assert systems[0] in chosen and len(chosen) > 1     # floor weight 0.2 keeps others alive


def test_build_argv():
    s = LaneSpec("G1", "chemical_system_energy_above_hull", {"chemical_system": "H-Ir-Mg", "energy_above_hull": 0.05},
                 batch_size=16)
    argv = build_argv(s, 20, 7, "/out", executable="mg")
    assert argv[:5] == ["nice", "-n", "19", "mg", "/out"]
    assert "--pretrained_name=chemical_system_energy_above_hull" in argv
    assert "--batch_size=16" in argv and "--num_batches=2" in argv and "--seed=7" in argv
    assert "--properties_to_condition_on={'chemical_system':'H-Ir-Mg','energy_above_hull':0.05}" in argv
    assert "--diffusion_guidance_factor=2.0" in argv
    u = build_argv(LaneSpec("G5", "mattergen_base", {}), 4, 1, "/o", executable="mg")
    assert not any(a.startswith("--properties_to_condition_on") for a in u)


FAKE = textwrap.dedent("""\
    #!/usr/bin/env python3
    import sys, os
    out = sys.argv[1]
    os.makedirs(out, exist_ok=True)
    with open(os.path.join(out, "generated_crystals.extxyz"), "w") as f:
        for k in range(3):
            f.write("3\\n")
            f.write('Lattice="3.0 0.0 0.0 0.0 3.0 0.0 0.0 0.0 3.0" Properties=species:S:1:pos:R:3 pbc="T T T"\\n')
            f.write("Mg 0.0 0.0 0.0\\n")
            f.write("H 1.5 1.5 0.0\\n")
            f.write(("Ir" if k < 2 else "Li") + " 1.5 0.0 1.5\\n")
    with open(os.path.join(out, "calls.txt"), "a") as f:
        f.write(" ".join(sys.argv[1:]) + "\\n")
""")


def test_generate_with_fake_executable(tmp_path):
    exe = tmp_path / "mattergen-generate"
    exe.write_text(FAKE)
    exe.chmod(exe.stat().st_mode | stat.S_IXUSR)
    spec = LaneSpec("G1", "chemical_system_energy_above_hull", {"chemical_system": "H-Ir-Mg", "energy_above_hull": 0.1},
                    batch_size=16)
    raws = generate(spec, 3, 11, workdir=tmp_path / "runs", batch_id="t-G1-000", executable=str(exe))
    assert len(raws) == 3
    r0 = raws[0]
    assert r0.origin.lane == "G1" and r0.origin.seed == 11 and r0.origin.guidance == 2.0
    assert r0.origin.conditioning == {"chemical_system": "H-Ir-Mg", "energy_above_hull": 0.1}
    assert "HYPOTHETICAL_STRUCTURE" in r0.origin.flags
    assert "OFF_REQUESTED_SYSTEM" in raws[2].origin.flags
    assert sorted(r0.structure.species) == ["H", "Ir", "Mg"] and r0.pressure_GPa == 0.0
    # same invocation is reused, not rerun
    raws2 = generate(spec, 3, 11, workdir=tmp_path / "runs", batch_id="t-G1-000", executable=str(exe))
    calls = list((tmp_path / "runs").glob("*/calls.txt"))[0].read_text().splitlines()
    assert len(calls) == 1 and raws2[0].meta["run"].get("reused")
    st = json.loads(list((tmp_path / "runs").glob("*/rtt_invocation.json"))[0].read_text())
    assert st["finished"]


def test_generate_failure_raises(tmp_path):
    from rtt.source.mattergen import MatterGenRunError
    exe = tmp_path / "bad"
    exe.write_text("#!/bin/sh\nexit 3\n")
    exe.chmod(exe.stat().st_mode | stat.S_IXUSR)
    with pytest.raises(MatterGenRunError):
        generate(LaneSpec("G5", "mattergen_base", {}), 2, 1, workdir=tmp_path, executable=str(exe))


def test_source_ref_has_no_home_dir():
    """Campaign 1: source_ref flags carried the absolute home path into ledger exports."""
    from pathlib import Path
    from rtt.source.mattergen import _tilde
    assert _tilde(Path.home() / "x" / "generated_crystals.extxyz") == "~/x/generated_crystals.extxyz"
    assert _tilde(Path("/opt/run/a.extxyz")) == "/opt/run/a.extxyz"


def test_checkpoint_sha256_first_hex_hop(tmp_path, monkeypatch):
    """Doc 34 §1: the HF blob (content-sha256 name) may itself be a symlink into a second store
    whose name is a different hash; the first 64-hex name in the chain is the content hash."""
    import hashlib
    from rtt.source.mattergen import checkpoint_sha256
    data = b"fake checkpoint"
    good = hashlib.sha256(data).hexdigest()
    other = "f" * 64
    blobs = tmp_path / "models--microsoft--mattergen" / "blobs"
    (blobs / "ff").mkdir(parents=True)
    (blobs / "ff" / other).write_bytes(data)
    (blobs / good).symlink_to(blobs / "ff" / other)
    ck = tmp_path / "models--microsoft--mattergen" / "snapshots" / "abc" / "checkpoints" / "x" / "checkpoints"
    ck.mkdir(parents=True)
    (ck / "last.ckpt").symlink_to(os.path.relpath(blobs / good, ck))
    monkeypatch.setenv("RTT_HF_HUB", str(tmp_path))
    assert checkpoint_sha256("x") == good
