"""rtt.dft.collect (core parsers + mode-positive admission, sha256 refs) and
rtt.dft.coarse (metrics + the S5 margin rule), on real QE 7.6 Pb bench output."""
from __future__ import annotations

import shutil
from fractions import Fraction
from pathlib import Path

import pytest

from rtt.dft import coarse as C
from rtt.dft.collect import CollectError, collect_run, spectra_from_run, spectrum_digest
from rtt.dft.profiles import PROFILES
from rtt.dft.workflow import DFPTRun

DATA = Path(__file__).parent / "data" / "qe76_pb_bench"
A = 4.95 / 2  # fcc Pb, angstrom (only the cell's existence matters for collection)


def _run(w: Path, status="done") -> DFPTRun:
    prof = PROFILES["coarse"].with_overrides(q=(8, 8, 8), el_ph_sigma_Ry=0.01,
                                             plateau_sigma_min_Ry=0.03)
    s = prof.eph_settings(do_relax=False)
    return DFPTRun(material_id="Pb", profile=prof.to_dict(), workdir=str(w), status=status,
                   q=[8, 8, 8], settings=s.to_dict(),
                   cell0={"lattice_A": [[0, A, A], [A, 0, A], [A, A, 0]], "species": ["Pb"],
                          "frac_coords": [[0, 0, 0]], "masses_amu": {"Pb": 207.2}},
                   pseudo_dir="/nonexistent", pseudo_files={"Pb": "Pb_ONCV_PBE-1.2.upf"})


@pytest.fixture
def pbdir(tmp_path):
    w = tmp_path / "pb"
    shutil.copytree(DATA, w)
    return w


def test_collect_admission_and_hashes(pbdir):
    from rtt.dft.prune import sha256_file
    c = collect_run(_run(pbdir))
    assert c.ok, c.errors
    adm = c.admission
    assert adm["rule"] == "mode_positive" and adm["excluded"], adm   # Pb: smallest sigma(s) out
    idx = [sp.provenance["smearing_index"] for sp in c.spectra]
    assert idx == adm["admitted"] and min(idx) > 1
    names = {Path(a.path).name: a for a in c.artifacts}
    for sp in c.spectra:
        src = Path(sp.provenance["source_file"])
        assert sp.provenance["source_sha256"] == sha256_file(src) == names[src.name].sha256
        assert sp.provenance["spectrum_tier"] == "coarse" and sp.provenance["profile"] == "coarse"
        assert sp.provenance["spectrum_digest"] == spectrum_digest(sp)
        assert any(a.startswith("A-SIGMA-EXCLUDED") for a in sp.provenance.get("assumptions", []))
    assert "lambda" in names and "a2F.dos1" in names
    sp, arts = spectra_from_run(_run(pbdir))
    assert len(sp) == len(c.spectra) and len(arts) == len(c.artifacts)


def test_collect_errors(pbdir):
    (pbdir / "a2F.dos5").write_text("garbage\n")
    with pytest.raises(CollectError):
        spectra_from_run(_run(pbdir))
    c = collect_run(_run(pbdir, status="failed"))
    assert not c.ok and "not done" in c.errors[0]


def test_collect_registers_in_ledger(pbdir):
    from rtt.ledger.store import Ledger
    from rtt.records import Origin
    L = Ledger(":memory:")
    L.insert_candidate("c-1", "s" * 64, {"lattice": [[1, 0, 0], [0, 1, 0], [0, 0, 1]],
                                         "species": ["Pb"], "frac_coords": [[0, 0, 0]]},
                       Origin(lane="C", control="Pb"))
    tr = L.start_tier_run("c-1", "S5", "dfpt", profile="coarse")
    c = collect_run(_run(pbdir), ledger=L, tier_run_id=tr)
    assert all(a.id is not None for a in c.artifacts)
    assert all(isinstance(sp.provenance.get("source_artifact_id"), int) for sp in c.spectra)
    assert L.q1("SELECT COUNT(*) FROM artifacts WHERE tier_run_id=?", (tr,))[0] == len(c.artifacts)


def test_coarse_metrics_plateau(pbdir):
    run = _run(pbdir)
    m = C.coarse_metrics(collect_run(run), run)
    sig = [r.sigma_Ry for r in m.rows]
    assert sig[0] == pytest.approx(0.01) and len(sig) == 10
    assert [r.index for r in m.plateau] == [r.index for r in m.rows
                                             if r.admitted and r.sigma_Ry >= 0.03 - 1e-12]
    mt = m.to_metrics()
    assert mt["plateau_sigmas"] and all(a > 0 for a in mt["A_K_by_sigma"])
    assert mt["A_mid_K"] is not None
    # excluded smearings are still reported (informational)
    ex = [r for r in m.rows if not r.admitted]
    assert ex and all(r.A_K is not None for r in ex)


def test_zpe_from_dos(tmp_path):
    # single-frequency DOS at 1000 cm^-1 (normalised later) -> ZPE = 3/2 * hbar*omega per atom
    p = tmp_path / "phonon.dos"
    p.write_text("# freq dos\n" + "\n".join(f"{w} {1.0 if abs(w - 1000) < 1e-9 else 0.0}"
                                           for w in (998.0, 999.0, 1000.0, 1001.0, 1002.0)))
    z = C.zpe_meV_per_atom(p, nat=2)
    assert z == pytest.approx(1.5 * 1000 * C.CM1_TO_MEV)


def test_min_freq_K_for_s6(pbdir):
    """min_freq_K (signed, K) from dyn files of every irreducible q + the matdyn mesh."""
    run = _run(pbdir)
    (pbdir / "dyn1").write_text("     freq (    1) =      -1.000000 [THz] =     -50.000000 [cm-1]\n"
                                "     freq (    2) =       2.000000 [THz] =      66.700000 [cm-1]\n")
    (pbdir / "dyn2").write_text("     freq (    1) =       1.000000 [THz] =      33.000000 [cm-1]\n")
    m = C.coarse_metrics(collect_run(run), run)
    mt = m.to_metrics()
    assert mt["min_freq_K"] == pytest.approx(-50.0 * C.CM1_TO_K)
    assert mt["max_imag_mode_K"] == pytest.approx(50.0 * C.CM1_TO_K)
    assert mt["min_freq_sources"]["dyn_irreducible_q"] == 2
    assert mt["omega_max_K"] and mt["A_lo_K"] <= mt["A_mid_K"]
    assert m.needs_phonon_check


# ------------------------------------------------------------------ S5 rule (doc 13 §5 T2-SCREEN)
def _ladder(As, wlogs=None, q=(4, 4, 4), admitted=None, minf=None, flags=()):
    wlogs = wlogs or [300.0] * len(As)
    admitted = admitted or [True] * len(As)
    rows = [C.SigmaRow(index=i + 1, sigma_Ry=0.01 * (i + 1), admitted=ad, A_K=a,
                       omega_log_K=wl, on_plateau=ad and 0.01 * (i + 1) >= 0.03 - 1e-12)
            for i, (a, wl, ad) in enumerate(zip(As, wlogs, admitted))]
    return C.CoarseMetrics(rows=rows, plateau_sigma_min_Ry=0.03, min_freq_matdyn_cm1=minf,
                           q=list(q), flags=list(flags))


FLOOR = Fraction(433148145, 100000000) * 300             # A_target at T* = 300 K: 1299.44 K


def _th():
    from rtt.config import Threshold
    return dict(eps_sys=Threshold(0.10, "doc01#0.2b"),
                e_c=Threshold(0.5, "doc13#5", provisional=True),
                margin_M_c=Threshold(2.0, "doc13#5", provisional=True))


def test_coarse_envelope_formula():
    fine = _ladder([400.0, 420.0, 440.0, 430.0], wlogs=[300.0, 300, 300, 300])
    coarse = _ladder([380.0, 400.0, 420.0, 410.0], q=(3, 3, 3), wlogs=[330.0, 330, 330, 330])
    env = C.coarse_envelope(fine, coarse, 0.10)
    assert env["sigma_nom_index"] == 1 and env["A_nom_K"] == 400.0
    assert env["d_sig_A"] == pytest.approx(0.10) and env["d_q_A"] == pytest.approx(0.05)
    assert env["eps_omega"] == pytest.approx(0.10)
    assert env["A_hi_K"] == pytest.approx(400 * 1.15 * 1.10 / 0.90)
    # sigma admitted on the fine mesh only -> s_nom is the smallest common one
    coarse2 = _ladder([380.0, 400.0, 420.0, 410.0], q=(3, 3, 3), admitted=[False, True, True, True])
    assert C.coarse_envelope(fine, coarse2, 0.1)["sigma_nom_index"] == 2
    assert not C.coarse_envelope(fine, _ladder([1.0], admitted=[False]), 0.1)["ok"]


def test_t2_screen_two_meshes():
    fine = _ladder([400.0, 420.0, 440.0, 430.0])
    coarse = _ladder([380.0, 400.0, 420.0, 410.0], q=(3, 3, 3))
    # doc 16 statistic: A_nom = max(a_top) = 440 (fine max), d_sig = 40/440, d_q = 20/400 at
    # sigma 1, eps_w floor 0.02 -> 1.5 * 440 (1 + 1/11 + 0.05) 1.1/0.98 = 845 K < 1299 K
    o = C.decide_coarse(fine, metrics_coarser=coarse, floor_K=FLOOR, n_sheets=1, **_th())
    assert (o.outcome, o.criterion) == ("SCREENED_OUT", "coarse-A envelope (all-sigma)")
    assert o.provisional and "PROVISIONAL_MARGIN" in o.tags
    assert o.value == pytest.approx(1.5 * 440 * (1 + 40 / 440 + 0.05) * 1.10 / 0.98)
    # above target -> PASS
    big = _ladder([900.0, 920.0, 940.0])
    o = C.decide_coarse(big, metrics_coarser=_ladder([880.0, 900.0, 920.0], q=(3, 3, 3)),
                        floor_K=FLOOR, n_sheets=1, **_th())
    assert o.outcome == "PASS"


@pytest.mark.parametrize("mod,blocker", [
    (dict(n_sheets=2), "multiband"),
    (dict(n_sheets=None), "multiband"),
    (dict(fine=dict(minf=-100.0)), "imaginary_modes"),
    (dict(coarse=dict(flags=["A_DROP_NONPOS"])), "imaginary_modes"),
    (dict(fine=dict(As=[400.0, 420.0, 440.0], admitted=[False, True, True])), "ladder_short"),
    (dict(coarse=dict(q=(2, 2, 2))), "q2_not_allowed"),
])
def test_t2_screen_preconditions_defer(mod, blocker):
    fk = dict(As=[400.0, 420.0, 440.0, 430.0])
    fk.update(mod.get("fine", {}))
    ck = dict(As=[380.0, 400.0, 420.0, 410.0], q=(3, 3, 3))
    ck.update(mod.get("coarse", {}))
    o = C.decide_coarse(_ladder(**fk), metrics_coarser=_ladder(**ck), floor_K=FLOOR,
                        n_sheets=mod.get("n_sheets", 1), **_th())
    assert o.outcome == "DEFERRED", o
    assert f"BLOCKED:{blocker}" in o.tags
    if blocker == "imaginary_modes":
        assert "NEEDS_PHONON_CHECK" in o.tags


def test_t2_admission_failure_defers():
    fine = _ladder([400.0, 420.0, 440.0, 430.0])
    coarse = _ladder([380.0, 400.0, 420.0, 410.0], q=(3, 3, 3))
    coarse.admission_ok = False
    o = C.decide_coarse(fine, metrics_coarser=coarse, floor_K=FLOOR, n_sheets=1, **_th())
    assert o.outcome == "DEFERRED" and "BLOCKED:admission_failed" in o.tags


def test_single_mesh_fallback_only_at_q4():
    th = _th()
    m4 = _ladder([300.0, 310.0, 320.0])
    # doc 16 one-mesh A_hi = 320 (1 + 20/320) 1.1 / 0.98; times M_c = 2 -> 763 K < 1299
    o = C.decide_coarse(m4, floor_K=FLOOR, n_sheets=1, **th)
    assert (o.outcome, o.criterion) == ("SCREENED_OUT", "coarse-A max (all-sigma, single mesh)")
    assert o.value == pytest.approx(2 * 320 * (1 + 20 / 320) * 1.1 / 0.98)
    m3 = _ladder([300.0, 310.0, 320.0], q=(3, 3, 3))
    o = C.decide_coarse(m3, floor_K=FLOOR, n_sheets=1, **th)
    assert o.outcome == "DEFERRED" and "BLOCKED:single_mesh_q_below_4" in o.tags
    # one row: A_hi = A (1.1/0.98); 2*1.1*579/0.98 = 1299.8 >= 1299.44 -> PASS
    edge = _ladder([579.0])
    edge.rows[0].sigma_Ry, edge.rows[0].on_plateau = 0.03, True
    assert C.decide_coarse(edge, floor_K=FLOOR, n_sheets=1, **th).outcome == "PASS"


def test_coarse_rule_required_inputs_and_errors():
    from rtt.config import Threshold
    m = _ladder([300.0, 310.0, 320.0])
    with pytest.raises(ValueError):
        C.decide_coarse(m, floor_K=FLOOR, eps_sys=None, margin_M_c=Threshold(2.0, "x"), n_sheets=1)
    with pytest.raises(ValueError):
        C.decide_coarse(m, metrics_coarser=m, floor_K=FLOOR, eps_sys=Threshold(0.1, "x"),
                        n_sheets=1)                                 # e_c missing
    assert C.decide_coarse(None, floor_K=FLOOR, error="ph.x crashed", **_th()).outcome == "ERROR"
    # never a rejection outcome
    for As in ([1.0, 1.0, 1.0], [1e5, 1e5, 1e5]):
        o = C.decide_coarse(_ladder(As), metrics_coarser=_ladder(As, q=(3, 3, 3)),
                            floor_K=FLOOR, n_sheets=1, **_th())
        assert o.outcome in ("PASS", "SCREENED_OUT")
