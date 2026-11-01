"""Real QE 7.6 end-to-end smoke: fcc Al, profile s4 and profile coarse at q 3^3 + q 4^3, 2 MPI ranks,
through collect.py and the S5 margin rule, with artifacts registered in a ledger.

Opt-in (other agents run QE on the same machine): set RTT_QE_SMOKE=1. Uses QE from
``RTT_QE_BIN`` (default ~/micromamba/envs/qe/bin) and SG15 from ``RTT_PSEUDO_DIR``. Work dir:
``RTT_SMOKE_DIR`` (default: pytest tmp dir). Cheap settings (ecut 60 Ry -- SG15 Al-1.2 has 11 valence electrons, 25 Ry gives P = -6800 kbar --, k spacing 0.5 / k_fit 0.25 1/A)
are smoke-test settings, not physics. Wall time ~1-3 min.
"""
from __future__ import annotations

import os
from fractions import Fraction
from pathlib import Path

import pytest

from rtt.dft.profiles import PROFILES, qe_bin

pytestmark = pytest.mark.skipif(
    os.environ.get("RTT_QE_SMOKE") != "1" or not (qe_bin() / "ph.x").exists(),
    reason="set RTT_QE_SMOKE=1 (needs QE 7.6 + SG15) to run the real QE smoke test")

A_AL = 4.04


def al_fcc() -> dict:
    a = A_AL / 2
    return {"lattice": [[0, a, a], [a, 0, a], [a, a, 0]], "species": ["Al"],
            "frac_coords": [[0, 0, 0]], "masses_amu": {"Al": 26.9815385}}


@pytest.fixture(scope="module")
def workroot(tmp_path_factory):
    d = os.environ.get("RTT_SMOKE_DIR")
    return Path(d).expanduser() if d else tmp_path_factory.mktemp("qe_smoke")


def test_al_s4_coarse_end_to_end(workroot):
    from rtt.config import Threshold, load_machine, load_targets, load_tiers
    from rtt.dft.coarse import coarse_envelope, coarse_metrics, decide_coarse
    from rtt.dft.collect import collect_run, spectra_from_run
    from rtt.dft.electronic import decide_s4, electronic_artifacts, run_electronic
    from rtt.dft.workflow import RunPolicy, default_launcher, run_coarse
    from rtt.ledger.store import Ledger
    from rtt.records import Origin

    launcher = default_launcher(2)
    policy = RunPolicy(max_machine_ranks=12, rank_poll_s=20, rank_wait_timeout_s=3600,
                       min_free_gb=20, watcher_interval_s=2)

    # ---------------- S4
    s4 = PROFILES["s4"].with_overrides(ecutwfc_Ry=60.0, k_spacing=0.5, nscf_k_spacing=0.15)
    er = run_electronic(al_fcc(), workroot / "s4", profile=s4, launcher=launcher,
                        material_id="Al", policy=policy)
    assert er.ok, er.error
    assert er.n_sheets and er.n_sheets >= 1 and er.gap_eV == 0.0
    assert er.moment_muB is not None and er.moment_muB < 0.05
    assert er.N_EF_per_eV_atom and 0.2 < er.N_EF_per_eV_atom < 0.6   # free-electron Al: 0.39 /eV/atom
    # G-VALID E_F = shallowest crossing band (Al: small 3rd-band pockets, ~0.1 eV);
    # the deepest crossing band is the free-electron-like one (~11 eV)
    assert er.E_F_eV is not None and er.E_F_eV > 0
    assert max(b["depth_eV"] for b in er.crossing_bands) > 2.0
    assert er.phi is not None and er.H_f == 0.0 and er.Tc_phi_K is None
    assert er.scf_pressure_kbar is not None and abs(er.scf_pressure_kbar) < 30.0
    assert not any(f.startswith("RELAX_PRESSURE_MISMATCH") for f in er.flags)
    o4 = decide_s4(er, gap_max=Threshold(0.0, "doc02#5a:H1"),
                   moment_tol=Threshold(0.05, "doc02#5a:H2", True))
    assert o4.outcome == "PASS"
    steps = [t["step"] for t in er.run.steps]
    assert steps == ["relax", "scf", "elf", "nscf", "nscf_tet", "dos", "spin"]
    # tetrahedron N(E_F) (dos.x, windowed around E_F); free-electron Al: 0.39 /eV/atom
    assert 0.3 < er.N_EF_per_eV_atom < 0.5, er.N_EF_per_eV_atom
    assert not any(f.startswith("DOS_EF_MISMATCH") for f in er.flags), er.flags
    dos_t = [t for t in er.run.steps if t["step"] == "dos"][0]
    assert dos_t["wall_s"] < 120, dos_t                 # full-range dos.x took 470 s
    assert er.sheet_w_min is not None and 0 < er.sheet_w_min < 1
    arts = electronic_artifacts(er)
    assert {Path(a.path).name for a in arts} >= {"spin.out", "scf.out", "pwscf.dos", "data-file-schema.xml"}
    assert all(t["wall_s"] > 0 for t in er.run.steps if t["state"] == "ran")


    # ---------------- S5 coarse: q 3^3 and q 4^3 on the S4-relaxed cell (doc 13 T2-SCREEN)
    coarse = PROFILES["coarse"].with_overrides(ecutwfc_Ry=60.0, k_spacing=0.5,
                                               k_fit_spacing=0.25, q_meshes=((3, 3, 3), (4, 4, 4)),
                                               q=(4, 4, 4), matdyn_factor=4,
                                               plateau_sigma_min_Ry=0.03)
    runs = run_coarse(er.relaxed_cell, coarse, workroot / "coarse", launcher, material_id="Al",
                      policy=policy)
    assert len(runs) == 2 and all(r.ok for r in runs), [r.error for r in runs]
    r3, r4 = runs
    assert [s["step"] for s in r3.steps] == ["relax", "scf_fit", "scf", "ph", "q2r", "matdyn"]
    assert [s["step"] for s in r4.steps] == ["scf_fit", "scf", "ph", "q2r", "matdyn"]  # same cell
    assert r3.grids["q"] == [3, 3, 3] and r4.grids["q"] == [4, 4, 4]
    for r in runs:
        assert r.grids["k_fit"][0] % r.grids["k_scf"][0] == 0 and r.grids["k_scf"][0] % r.q[0] == 0
        assert r.core_h > 0 and r.qe_version and "7.6" in r.qe_version
        assert (Path(r.workdir) / "profile_done.json").is_file()

    # ---------------- collect into a ledger
    t, ti, m = load_targets(), load_tiers(), load_machine()
    L = Ledger(":memory:")
    L.register_configs(t, ti, m)
    L.insert_candidate("c-al", "a" * 64, {"lattice": al_fcc()["lattice"], "species": ["Al"],
                                          "frac_coords": [[0, 0, 0]]},
                       Origin(lane="C", control="Al", T_star_surrogate_K=0.6))
    tr = L.start_tier_run("c-al", "S5", "dfpt", profile="coarse", params=coarse.to_dict())
    mets, cols = [], []
    for r in runs:
        c = collect_run(r, ledger=L, tier_run_id=tr)
        assert c.spectra and all(sp.provenance["source_sha256"] for sp in c.spectra), c.errors
        assert all(a.id for a in c.artifacts)
        assert all(0.05 < sp.lambda_() < 2.0 for sp in c.spectra)
        met = coarse_metrics(c, r)
        assert met.min_freq_K is not None and "min_freq_K" in met.to_metrics()
        assert not any(f.startswith("RELAX_PRESSURE_MISMATCH") for f in met.flags), met.flags
        mets.append(met)
        cols.append(c)
    m3, m4 = mets
    L.finish_tier_run(tr, "done", wall_s=sum(r.wall_s for r in runs),
                      core_h=sum(r.core_h for r in runs),
                      metrics={"q3": m3.to_metrics(), "q4": m4.to_metrics()})
    env = coarse_envelope(m4, m3, 0.10)
    assert env["ok"], env

    # T2-SCREEN at the real target. Al has n_sheets = 2 at S4, so the screen DEFERs (multiband);
    # treated as single-sheet, Al's A (~tens of K) is far below A_target -> SCREENED_OUT unless a
    # precondition (admission, >= 3 admitted sigma on q 4^3, phonons) blocks it.
    kw = dict(floor_K=t.floor, eps_sys=t.threshold("envelope.eps_sys"),
              e_c=ti.threshold("S5.e_c"), metrics_coarser=m3, tier_run_id=tr)
    o_mb = decide_coarse(m4, n_sheets=er.n_sheets, w_min=er.sheet_w_min, **kw)
    assert o_mb.outcome == "DEFERRED" and "BLOCKED:multiband" in o_mb.tags
    o5 = decide_coarse(m4, n_sheets=1, **kw)
    assert o5.outcome in ("SCREENED_OUT", "DEFERRED"), o5
    if o5.outcome == "SCREENED_OUT":
        assert o5.criterion == "coarse-A envelope" and "PROVISIONAL_MARGIN" in o5.tags
        L.insert_outcome("c-al", o5, targets_sha256=t.sha256, tiers_sha256=ti.sha256,
                         origin="rtt.dft.coarse.decide_coarse")
    print("T2-SCREEN Al:", o5.outcome, o5.criterion, o5.tags, o5.reason)
    # a tiny surrogate floor (T* = 1 K): never removed
    o5b = decide_coarse(m4, n_sheets=1, **dict(kw, floor_K=t.with_T_star(1.0).floor))
    assert o5b.outcome == "PASS"

    # ---------------- adoption: a second call re-runs nothing
    runs2 = run_coarse(er.relaxed_cell, coarse, workroot / "coarse", launcher, material_id="Al",
                       policy=policy)
    assert all(r.ok and r.adopted for r in runs2)
    sp2, _ = spectra_from_run(runs2[1])
    assert [s.provenance["source_sha256"] for s in sp2] == \
        [s.provenance["source_sha256"] for s in cols[1].spectra]
    er2 = run_electronic(al_fcc(), workroot / "s4", profile=s4, launcher=launcher,
                         material_id="Al", policy=policy)
    assert er2.ok and all(t["state"] != "ran" for t in er2.run.steps)
