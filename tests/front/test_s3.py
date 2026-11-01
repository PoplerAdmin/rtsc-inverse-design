"""S3 decision rules and the S3 runner with fake MatterSim hooks (CPU; no torch needed).

The central contract (R2, doc 10 §3.4): MatterSim phonons never remove a candidate. Mg2IrH6 has
imaginary MatterSim modes and must still PASS S3.
"""
from __future__ import annotations

import math

import numpy as np
import pytest

from rtt.mlff.hull import SelfHull, ehull_self
from rtt.mlff.phonons import PhononSummary, auto_supercell, summarize
from rtt.mlff.relax import RelaxResult
from rtt.mlff.screen import apply_rank_cut, compute_s3, decide_s3, run_s3, score_s3
from rtt.records import ALLOWED_OUTCOMES, Outcome
from rtt.source.controls import inject_controls
from rtt.structure.stage import run_s2


# ------------------------------------------------------------------ pure decision
def test_decide_s3_hull_rule(tiers):
    assert decide_s3({"e_hull_ml_self_eV": 0.10}, tiers).outcome == "PASS"
    thr = float(tiers.threshold("S3.ehull_self_max_eV").value)
    assert thr >= 0.55                                                               # RC1 F15 floor
    assert decide_s3({"e_hull_ml_self_eV": thr}, tiers).outcome == "PASS"             # boundary kept
    hyd = {"features": {"n_H": 2}}                                                   # doc 31 §7: hydride domain
    o = decide_s3({"e_hull_ml_self_eV": thr + 0.01, **hyd}, tiers)
    assert o.outcome == "SCREENED_OUT" and o.criterion.startswith("H7-ml") and o.threshold == thr
    assert o.source_id and o.value == thr + 0.01
    o = decide_s3({"e_hull_ml_self_eV": None}, tiers)
    assert o.outcome == "PASS" and "NO_ML_HULL" in o.tags
    o = decide_s3({"error": "boom"}, tiers)
    assert o.outcome == "DEFERRED" and "ml_failed" in o.reason
    o = decide_s3({"e_hull_ml_self_eV": 0.9, "relax": {"converged": False}}, tiers)
    assert o.outcome == "DEFERRED" and o.criterion == "relax_not_converged"
    o = decide_s3({"e_hull_ml_self_eV": 0.1, "relax": {"converged": False}}, tiers)
    assert o.outcome == "PASS" and "RELAX_NOT_CONVERGED" in o.tags
    o = decide_s3({"e_hull_ml_self_eV": 0.9, **hyd}, tiers, screen_enabled=False)
    assert o.outcome == "PASS" and any(t.startswith("SCREEN_DISABLED") for t in o.tags)


@pytest.mark.parametrize("minf", [-363.0, -1000.0, -1e6])
def test_phonons_never_remove(tiers, minf):
    m = {"e_hull_ml_self_eV": 0.004, "e_hull_ml_dftref_eV": 5.0,             # dftref is second opinion only
         "phonons": {"min_freq_cm1": minf, "imag_mode_fraction": 0.5, "stable_ml": False,
                     "needs_phonon_check": True},
         "alignn": {"Tc": 0.0, "weight": 0.0}, "features": {"H_f": 0.667, "motif": "isolated_H"}}
    assert decide_s3(m, tiers).outcome == "PASS"


def test_rank_cut_and_control_exemption():
    P = lambda: Outcome(stage="S3", outcome="PASS")  # noqa: E731
    items = [("a", "G1", P(), 0.9), ("b", "G1", P(), 0.8), ("c", "G2", P(), 0.7), ("k", "C", P(), 0.0),
             ("s", "G1", Outcome(stage="S3", outcome="SCREENED_OUT", criterion="x", value=1, threshold=0,
                                 source_id="t"), 0.95)]
    out = apply_rank_cut(items, 2, "doc12#8.3")
    assert out["a"].outcome == "PASS" and out["b"].outcome == "PASS"
    assert out["c"].outcome == "DEFERRED" and out["c"].reason == "rank_cut(3, 2)"
    assert out["k"].outcome == "PASS" and "CONTROL_EXEMPT_RANK_CUT" in out["k"].tags
    assert out["s"].outcome == "SCREENED_OUT"
    for o in out.values():
        assert o.outcome in ALLOWED_OUTCOMES["S3"]


def test_score_monotone_in_hull():
    base = {"features": {"H_f": 0.6, "light_fraction": 0.6, "motif": "isolated_H"},
            "phonons": {"omega2_light_K": 1200.0}}
    s1 = score_s3({**base, "e_hull_ml_self_eV": 0.0})
    s2 = score_s3({**base, "e_hull_ml_self_eV": 0.3})
    assert 0 <= s2 < s1 <= 1
    assert score_s3({**base, "e_hull_ml_self_eV": 0.0, "features": {**base["features"], "motif": "H2_molecular"}}) < s1


# ------------------------------------------------------------------ hull, phonon summary, features
def test_ehull_self_on_toy_hull():
    h = SelfHull(chemsys="H-Mg", pressure_GPa=0.0, model="m",
                 entries=[{"entry_id": "Mg", "formula": "Mg", "enthalpy_per_atom_eV": -1.5, "converged": True},
                          {"entry_id": "H", "formula": "H2", "enthalpy_per_atom_eV": -3.4, "converged": True},
                          {"entry_id": "MgH2", "formula": "MgH2", "enthalpy_per_atom_eV": -3.0, "converged": True}],
                 elements_covered=["H", "Mg"])
    mgh2 = {"lattice": [[3, 0, 0], [0, 3, 0], [0, 0, 3]], "species": ["Mg", "H", "H"],
            "frac_coords": [[0, 0, 0], [0.5, 0.5, 0], [0, 0.5, 0.5]]}
    eh, info = ehull_self(mgh2, -2.9, hull=h)
    assert abs(eh - 0.1) < 1e-6
    eh, _ = ehull_self(mgh2, -3.1, hull=h)
    assert eh < 0                                                         # below the ML hull
    h2 = SelfHull(chemsys="H-Mg", pressure_GPa=0.0, model="m", entries=h.entries[:1], elements_covered=["Mg"])
    eh, info = ehull_self(mgh2, -3.0, hull=h2)
    assert eh is None and "incomplete" in info["note"]


def test_phonon_summary_pure():
    f = np.array([[0.0, 0.0, 0.0, 500.0, 500.0, 1500.0], [-50.0, 10.0, 20.0, 600.0, 700.0, 1400.0]])
    w = np.array([1.0, 1.0])
    lw = np.array([[0, 0, 0, 0.1, 0.1, 1.0], [0, 0, 0, 0.1, 0.1, 1.0]])
    s = summarize(f, w, lw, n_atoms=2, imag_tol_cm1=15.0)
    assert s["min_freq_cm1"] == -50.0 and s["needs_phonon_check"] and not s["stable_ml"]
    assert 0 < s["imag_mode_fraction"] <= 1 / 6 + 1e-9
    assert s["omega2_light_K"] > 1400 * 1.4387768775 * 0.8
    assert s["zpe_ml_meV_atom"] > 0
    assert auto_supercell(np.eye(3) * 3.0, 9.0, 250, 2) == [3, 3, 3]
    assert np.prod(auto_supercell(np.eye(3) * 3.0, 9.0, 30, 2)) * 2 <= 30


def test_features_motifs():
    from rtt.mlff.features import features
    h2 = {"lattice": [[5, 0, 0], [0, 5, 0], [0, 0, 5]], "species": ["Mg", "H", "H"],
          "frac_coords": [[0, 0, 0], [0.5, 0.5, 0.5], [0.648, 0.5, 0.5]]}
    f = features(h2)
    assert f["motif"] == "H2_molecular" and abs(f["min_HH_A"] - 0.74) < 1e-6 and f["contains_HLiBe"]
    al = {"lattice": [[4, 0, 0], [0, 4, 0], [0, 0, 4]], "species": ["Al"], "frac_coords": [[0, 0, 0]]}
    f = features(al)
    assert f["motif"] == "no_H" and f["vol_per_H_A3"] is None and not f["contains_HLiBe"]


# ------------------------------------------------------------------ runner with fakes
def _fake_relax(structs, model, pressure_GPa, device=None, **kw):
    ps = pressure_GPa if isinstance(pressure_GPa, list) else [pressure_GPa] * len(structs)
    out = []
    for i, (s, p) in enumerate(zip(structs, ps)):
        n = len(s["species"])
        nH = s["species"].count("H")
        e = -2.0 * n - 0.5 * nH
        out.append(RelaxResult(index=i, structure=s, formula="x", n_atoms=n, energy_eV=e, energy_per_atom_eV=e / n,
                               enthalpy_per_atom_eV=e / n, volume_per_atom_A3=10.0, max_force_eV_A=0.01,
                               internal_pressure_GPa=0.0, n_steps=10, converged=True, pressure_GPa=p, model=model))
    return out


def _fake_hull(cs, pressure_GPa, **kw):
    els = cs.split("-")
    return SelfHull(chemsys=cs, pressure_GPa=pressure_GPa, model="fake",
                    entries=[{"entry_id": e, "formula": e, "enthalpy_per_atom_eV": -2.0, "converged": True} for e in els],
                    elements_covered=els)


def _fake_phonons(struct, **kw):
    # every structure gets a strong ML instability, like Mg2IrH6 in doc 10 §3.4
    return PhononSummary(min_freq_cm1=-300.0, omega_max_cm1=2000.0, omega2_light_K=1500.0, zpe_ml_meV_atom=100.0,
                         imag_mode_fraction=0.05, omega_log_phDOS_K=500.0, light_weight_fraction=0.5,
                         stable_ml=False, needs_phonon_check=True)


def _fake_dftref(items, machine=None):
    return {k: {"e_above_hull_eV": 9.9} for k, _, _ in items}


def test_compute_s3_with_fakes(tiers):
    s = {"lattice": [[3, 0, 0], [0, 3, 0], [0, 0, 3]], "species": ["Mg", "H"], "frac_coords": [[0, 0, 0], [.5, .5, .5]]}
    m = compute_s3([s], [0.0], tiers=tiers, relax_fn=_fake_relax, hull_fn=_fake_hull, phonon_fn=_fake_phonons,
                   dftref_fn=_fake_dftref, alignn=False)[0]
    assert m["e_hull_ml_self_eV"] < 0 and m["e_hull_ml_dftref_eV"] == 9.9
    assert m["phonons"]["needs_phonon_check"] and m["features"]["H_f"] == 0.5


def test_run_s3_end_to_end_with_fakes(ledger, targets, tiers, tmp_path):
    ids = inject_controls(ledger, targets=targets, tiers=tiers)
    run_s2(ledger, targets=targets, tiers=tiers, novelty=False)
    out = run_s3(ledger, targets=targets, tiers=tiers, workdir=tmp_path, device="cpu", alignn=False,
                 relax_fn=_fake_relax, hull_fn=_fake_hull, phonon_fn=_fake_phonons, dftref_fn=_fake_dftref,
                 rank_cut=0)
    assert set(out) == set(ids)
    for cid in ids:
        o = ledger.current_outcome(cid, "S3")
        assert o["outcome"] == "PASS"                                   # controls: imaginary ML modes, still PASS
        tr = ledger.tier_run(o["tier_run_id"])
        assert tr["kind"] == "ml" and tr["state"] == "done" and tr["stage"] == "S3"
        assert ledger.q1("SELECT COUNT(*) FROM artifacts WHERE tier_run_id=?", (tr["id"],))[0] == 1
    # H3S (200 GPa) has no DFT-reference hull (0 GPa set)
    import json
    h3s = ids[1]
    met = json.loads(ledger.tier_run(ledger.current_outcome(h3s, "S3")["tier_run_id"])["metrics_json"])
    assert met["pressure_GPa"] == 200.0 and met["e_hull_ml_dftref_eV"] is None
    assert met["e_hull_ml_self_eV"] is None and "out of scope" in met["hull_self"]["note"]   # no 200 GPa hull
    assert run_s3(ledger, targets=targets, tiers=tiers, workdir=tmp_path, relax_fn=_fake_relax) == {}


def test_run_s3_screens_far_above_hull_and_defers_on_error(ledger, targets, tiers, tmp_path):
    from rtt.records import Origin, RawStructure, StructureData
    from rtt.source.ingest import ingest_raw
    a = RawStructure(structure=StructureData(lattice=[[3, 0, 0], [0, 3, 0], [0, 0, 3]], species=["Mg", "H"],
                                             frac_coords=[[0, 0, 0], [.5, .5, .5]]), origin=Origin(lane="G1"))
    b = RawStructure(structure=StructureData(lattice=[[3.3, 0, 0], [0, 3.3, 0], [0, 0, 3.3]], species=["Mg", "Mg", "H"],
                                             frac_coords=[[0, 0, 0], [.5, .5, 0], [.5, .5, .5]]), origin=Origin(lane="G1"))
    ca = ingest_raw(ledger, a, targets=targets, tiers=tiers)
    cb = ingest_raw(ledger, b, targets=targets, tiers=tiers)
    run_s2(ledger, targets=targets, tiers=tiers, novelty=False)

    def relax_high(structs, model, pressure_GPa, device=None, **kw):
        out = _fake_relax(structs, model, pressure_GPa)
        for r in out:
            if r.n_atoms == 2:
                r.enthalpy_per_atom_eV = -2.0 + 1.0                     # 1 eV/atom above the toy hull
            else:
                r.error, r.enthalpy_per_atom_eV = "relax diverged", math.nan
        return out

    out = run_s3(ledger, targets=targets, tiers=tiers, workdir=tmp_path, alignn=False, relax_fn=relax_high,
                 hull_fn=_fake_hull, phonon_fn=_fake_phonons, dftref_fn=_fake_dftref)
    assert out[ca].outcome == "SCREENED_OUT" and out[ca].criterion == "H7-ml:e_hull_ml_self_eV"
    assert out[cb].outcome == "DEFERRED" and out[cb].criterion == "ml_failed"
    assert ledger.tier_run(ledger.current_outcome(cb, "S3")["tier_run_id"])["state"] == "failed"


# ------------------------------------------------------------------ RC1 F12: retry + rank queue + score
def test_score_ignores_mattersim_phonons_for_HLiBe():
    base = {"e_hull_ml_self_eV": 0.05, "features": {"H_f": 0.6, "light_fraction": 0.6, "motif": "isolated_H",
                                                    "contains_HLiBe": True}}
    lo = score_s3({**base, "phonons": {"omega2_light_K": 10.0}})
    hi = score_s3({**base, "phonons": {"omega2_light_K": 5000.0}})
    assert lo == hi
    nb = {"e_hull_ml_self_eV": 0.05, "features": {"H_f": 0.0, "light_fraction": 1.0, "motif": "no_H",
                                                  "contains_HLiBe": False}}
    assert score_s3({**nb, "phonons": {"omega2_light_K": 10.0}}) < score_s3({**nb, "phonons": {"omega2_light_K": 5000.0}})


def _gens(ledger, targets, tiers, n):
    from rtt.records import Origin, RawStructure, StructureData
    from rtt.source.ingest import ingest_raw
    out = []
    for k in range(n):
        a = 3.0 + 0.3 * k
        out.append(ingest_raw(ledger, RawStructure(structure=StructureData(
            lattice=[[a, 0, 0], [0, a, 0], [0, 0, a]], species=["Mg"] + ["H"] * (k + 1),
            frac_coords=[[0, 0, 0]] + [[0.5, 0.5, (0.5 + 0.3 * j) % 1] for j in range(k + 1)]),
            origin=Origin(lane="G1")), targets=targets, tiers=tiers))
    return out


def test_rank_queue_promotes_on_rerank(ledger, targets, tiers, tmp_path):
    from rtt.mlff.screen import rerank_s3
    cids = _gens(ledger, targets, tiers, 3)
    run_s2(ledger, targets=targets, tiers=tiers, novelty=False)
    out = run_s3(ledger, targets=targets, tiers=tiers, workdir=tmp_path, alignn=False, dftref=False, phonons=False,
                 relax_fn=_fake_relax, hull_fn=_fake_hull, rank_cut=1)
    assert sum(o.outcome == "PASS" for o in out.values()) == 1
    queued = [c for c in cids if ledger.current_outcome(c, "S3")["criterion"] == "rank_cut"]
    assert len(queued) == 2
    # not recomputed by a second run_s3; nothing promoted while the PASS still waits for S4
    assert run_s3(ledger, targets=targets, tiers=tiers, workdir=tmp_path, relax_fn=_fake_relax, rank_cut=1) == {}
    # next day, bigger capacity: the queue is re-ranked from stored scores
    prom = rerank_s3(ledger, targets=targets, tiers=tiers, rank_cut=3)
    assert set(prom) == set(queued) and all("PROMOTED_FROM_RANK_CUT" in o.tags for o in prom.values())


def test_ml_failure_retried_with_more_steps(ledger, targets, tiers, tmp_path):
    seen = {}

    def boom(structs, model, pressure_GPa, device=None, **kw):
        raise RuntimeError("CUDA out of memory")

    def ok(structs, model, pressure_GPa, device=None, **kw):
        seen["max_steps"] = kw.get("max_steps")
        return _fake_relax(structs, model, pressure_GPa)

    c = _gens(ledger, targets, tiers, 1)[0]
    run_s2(ledger, targets=targets, tiers=tiers, novelty=False)
    assert run_s3(ledger, targets=targets, tiers=tiers, workdir=tmp_path, alignn=False, dftref=False, phonons=False,
                  relax_fn=boom, hull_fn=_fake_hull)[c].outcome == "DEFERRED"
    o = run_s3(ledger, targets=targets, tiers=tiers, workdir=tmp_path, alignn=False, dftref=False, phonons=False,
               relax_fn=ok, hull_fn=_fake_hull)[c]
    assert o.outcome == "PASS" and o.retry == 1 and seen["max_steps"] == 2000


def test_auto_supercell_terminates_with_long_axis():
    """Campaign 1 hang: 20-atom cell, short a/b, long c (rep [2,2,1]) over the atom cap looped forever."""
    from rtt.mlff.phonons import auto_supercell
    rep = auto_supercell([[5.0, 0, 0], [0, 5.0, 0], [0, 0, 20.0]], 9.0, 60, 20)
    assert rep == [1, 1, 1] or 20 * rep[0] * rep[1] * rep[2] <= 60
    assert all(r >= 1 for r in rep)
    assert auto_supercell([[3.0, 0, 0], [0, 3.0, 0], [0, 0, 3.0]], 9.0, 250, 1) == [3, 3, 3]
