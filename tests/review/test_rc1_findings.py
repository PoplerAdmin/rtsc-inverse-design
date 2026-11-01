"""RC1 adversarial code review: one test per finding (docs/reviews/RC1_code_review.md).

Every test here encodes the behaviour the spec (architecture rev 2 §3.1/§4, doc 02 §5, doc 01
§0.4/§5.1) requires and FAILS on the code as reviewed (2026-10-29). Test names carry the finding
number (F<n>). CPU only; no QE, no MatterSim, no GPU.
"""
from __future__ import annotations

import hashlib
import math
from fractions import Fraction

import pytest

from rtt.config import Threshold
from rtt.dft import electronic as E
from rtt.dft.coarse import CoarseMetrics, SigmaRow, decide_coarse
from rtt.dft.workflow import DFPTRun

HA = E.HA_EV


# ----------------------------------------------------------------------------- helpers
def _res(tmp_path, **kw) -> E.ElectronicResult:
    run = DFPTRun(material_id="x", profile={"name": "s4", "kind": "electronic"}, workdir=str(tmp_path),
                  status="done")
    r = E.ElectronicResult(run=run, nat=kw.pop("nat", 4))
    for k, v in kw.items():
        setattr(r, k, v)
    return r


def _xml(path, eigs_ha, weights, ef_ha, nat=2):
    ks = "".join(f'<ks_energies><k_point weight="{w}">0 0 0</k_point><npw>10</npw>'
                 f'<eigenvalues size="{len(e)}">{" ".join(map(str, e))}</eigenvalues></ks_energies>'
                 for e, w in zip(eigs_ha, weights))
    path.write_text(f'<?xml version="1.0"?><qes:espresso xmlns:qes="http://www.quantum-espresso.org/ns/qes/qes-1.0">'
                    f'<output><atomic_structure nat="{nat}"/><band_structure><lsda>false</lsda>'
                    f'<fermi_energy>{ef_ha}</fermi_energy>{ks}</band_structure></output></qes:espresso>')


H1 = Threshold(0.0, "doc02#5a:H1")
H2 = Threshold(0.05, "doc02#5a:H2", provisional=True)


def _metrics(A_list, *, min_freq=None, flags=()):
    rows = [SigmaRow(index=i + 3, sigma_Ry=0.01 * (i + 3), admitted=True, A_K=a, on_plateau=True)
            for i, a in enumerate(A_list)]
    return CoarseMetrics(rows=rows, plateau_sigma_min_Ry=0.03, min_freq_matdyn_cm1=min_freq,
                         flags=list(flags))


# ============================================================================ S4 (H1/H2/H7-heur)
def test_F1_s4_engine_threshold_tables_must_not_silently_disable_h1_h2(tmp_path, targets):
    """engine/stages.py S4Electronic.decide passes targets['hard']['H1_metallic'] (a table
    {gap_eV_max, source}) and ['H2_nonmagnetic'] ({moment_muB_tol, source, provisional}).
    electronic._th reads only 'value' -> None -> H1/H2 skipped, and no *_NOT_CHECKED tag is
    added, so evidence.py then prints 'H1 metallic: checked'. An insulating ferromagnet PASSES."""
    hard = targets.get("hard")
    r = _res(tmp_path, gap_eV=2.0, moment_muB=3.0, abs_moment_muB=3.0)
    o = E.decide_s4(r, gap_max=hard["H1_metallic"], moment_tol=hard["H2_nonmagnetic"])
    assert o.outcome != "PASS" or {"H1_NOT_CHECKED", "H2_NOT_CHECKED"} <= set(o.tags), \
        f"insulating magnet passed S4 silently: {o.outcome} tags={o.tags}"


def test_F2_s4_h1_doped_insulator_with_states_at_EF_is_not_excluded(tmp_path):
    """doc 02 H1 is 'N(E_F) > 0 at the PBE level' and explicitly keeps E_F near a gap edge
    (doped Mg2IrH6-type, S7). analyse_bands instead calls a system gapped when no band crosses
    E_Fermi on the nscf mesh. With smearing (mv 0.02 Ry) E_Fermi of a lightly electron-doped
    insulator sits a few 10 meV below the CB minimum: N(E_F) is clearly > 0 but gap_eV = 0.55 eV
    -> EXCLUDED_ASSUMPTION(H1)."""
    ef = 0.30
    vb, cb0, cb1 = ef - 0.50 / HA, ef + 0.05 / HA, ef + 1.0 / HA
    _xml(tmp_path / "d.xml", [[vb, cb0], [vb - 0.01, cb1]], [1.0, 1.0], ef)
    a = E.analyse_bands(E.parse_pw_xml(tmp_path / "d.xml"), dos_sigma_eV=0.136)
    assert a.N_EF_per_eV_atom > 0.1                     # states at E_F by doc 02's own criterion
    r = _res(tmp_path, gap_eV=a.gap_eV, moment_muB=0.0, abs_moment_muB=0.0,
             N_EF_per_eV_atom=a.N_EF_per_eV_atom)
    o = E.decide_s4(r, gap_max=H1, moment_tol=H2)
    assert o.outcome != "EXCLUDED_ASSUMPTION", f"N(E_F)={a.N_EF_per_eV_atom:.2f}/eV/atom but H1 fired (gap {a.gap_eV:.2f} eV)"


def test_F3_s4_mesh_limited_marginal_gap_is_deferred_not_excluded(tmp_path):
    """gap_eV_max = 0.0 with a strict '>' and no tolerance: a 10 meV mesh-limited gap (the code
    itself flags GAP_SMALL_MESH_LIMITED) is EXCLUDED. doc 02 H1: 'check with a hybrid functional
    if the gap is marginal' -> must be DEFERRED (or PASS with a flag), never excluded."""
    r = _res(tmp_path, gap_eV=0.01, moment_muB=0.0, abs_moment_muB=0.0, flags=["GAP_SMALL_MESH_LIMITED"])
    o = E.decide_s4(r, gap_max=H1, moment_tol=H2)
    assert o.outcome != "EXCLUDED_ASSUMPTION"


def test_F4_s4_h2_metastable_fm_solution_above_nm_is_not_excluded(tmp_path):
    """H2 is 'non-magnetic GROUND state'. The spin scf starts ferromagnetic and can converge to a
    metastable moment-bearing solution. Its energy (e_total_spin_Ry) is recorded but never
    compared with the NM energy (e_total_Ry): a candidate whose NM state is LOWER is excluded."""
    r = _res(tmp_path, gap_eV=0.0, moment_muB=0.8, abs_moment_muB=0.8,
             e_total_Ry=-100.000, e_total_spin_Ry=-99.990)            # FM solution 0.136 eV higher
    o = E.decide_s4(r, gap_max=H1, moment_tol=H2)
    assert o.outcome != "EXCLUDED_ASSUMPTION"


def test_F5_s4_h2_tolerance_is_per_cell_not_per_atom(tmp_path):
    """moment_muB is |total magnetization| per CELL; the tolerance 0.05 muB (doc 02 gives no
    number, 'moment ~ 0'; provisional) is applied per cell, so it shrinks per atom with cell size.
    A 20-atom cell with 0.003 muB/atom of smearing-level residual moment is excluded."""
    r = _res(tmp_path, nat=20, gap_eV=0.0, moment_muB=0.06, abs_moment_muB=0.06)
    o = E.decide_s4(r, gap_max=H1, moment_tol=H2)
    assert o.outcome != "EXCLUDED_ASSUMPTION"


def test_F6_s4_ferrimagnet_with_large_absolute_moment_is_not_a_plain_pass(tmp_path):
    """Total moment ~0 but absolute magnetization 4 muB/cell: H2 (non-magnetic normal state, the
    assumption of the whole Tc machinery) is violated, yet decide_s4 returns PASS with only the
    flag ANTIFERRO_OR_FERRI_SUSPECTED, and evidence.py prints 'H2 non-magnetic: checked'."""
    r = _res(tmp_path, gap_eV=0.0, moment_muB=0.0, abs_moment_muB=4.0,
             flags=["ANTIFERRO_OR_FERRI_SUSPECTED"])
    o = E.decide_s4(r, gap_max=H1, moment_tol=H2)
    assert o.outcome != "PASS" or "H2_NOT_CHECKED" in o.tags


def test_F7_s4_provisional_h2_tolerance_is_recorded_on_pass(tmp_path):
    """moment_muB_tol is provisional (targets.toml). A PASS decided with it carries neither
    provisional=True nor a PROVISIONAL_* tag, so S9 ('provisional margins used upstream') reads
    'checked (none used)'."""
    r = _res(tmp_path, gap_eV=0.0, moment_muB=0.01, abs_moment_muB=0.01)
    o = E.decide_s4(r, gap_max=H1, moment_tol=H2)
    assert o.outcome == "PASS"
    assert o.provisional or any(t.startswith("PROVISIONAL") for t in o.tags)


def test_F8_s4_h7heur_does_not_remove_on_unconverged_cutoff(tmp_path):
    """The code flags RELAX_PRESSURE_MISMATCH when the fresh-basis scf is far from the target
    pressure (cutoff unconverged, doc 22 §2: SG15 semicore). The formation-energy lower-bound
    argument assumes converged total energies; H7-heur still removes on such a run."""
    r = _res(tmp_path, gap_eV=0.0, moment_muB=0.0, abs_moment_muB=0.0,
             flags=["RELAX_PRESSURE_MISMATCH(P=-374.0kbar,target=0.0)"])
    o = E.decide_s4(r, gap_max=H1, moment_tol=H2, ehull_heur_max=Threshold(0.42, "doc02#5a:H7", True),
                    e_form_eV=0.45)
    assert o.outcome != "SCREENED_OUT"


# ============================================================================ S5 coarse
@pytest.mark.parametrize("min_freq,flags", [(-500.0, []),
                                            (None, ["IMAGINARY_MODES(min=-500.00cm-1)"]),
                                            (None, ["NEGATIVE_FREQ_A2F_WEIGHT"])])
def test_F9_s5_imaginary_modes_cannot_drive_coarse_removal(targets, min_freq, flags):
    """R2: 'coarse imaginary modes never remove'. But the core spectrum DROPS all alpha^2F
    weight at omega <= 0 (assumption A-DROP-NONPOS, 'lowers A'), and matdyn's a2F is
    unreliable/zero for unstable structures; decide_coarse still applies the A margin and returns
    SCREENED_OUT (only adding the tag NEEDS_PHONON_CHECK). A soft-mode-driven high-lambda
    candidate is removed through its imaginary modes."""
    m = _metrics([400.0, 420.0, 410.0], min_freq=min_freq, flags=flags)
    o = decide_coarse(m, floor_K=targets.floor, margin_M_c=Threshold(2.0, "doc11#4", True),
                      eps_sys=Threshold(0.10, "doc01#0.2b"))
    assert o.outcome != "SCREENED_OUT", f"{o.outcome} {o.tags}"


def test_F10_s5_multiband_isotropic_A_must_not_remove(targets):
    """doc 01 §5.1: for n_sheets > 1 the isotropic A falls by the DOS fraction w WITHOUT LIMIT
    (decoupled two-band example: A_iso 806 K for a true-300 K material); only kappa,
    max_i A_i or A_iso/w_min are valid N1 inputs. decide_coarse has no multiband input at all
    (S4 n_sheets is never consulted), so a multiband candidate is screened on A_iso."""
    m = _metrics([500.0, 520.0, 510.0])                     # A_iso/w_min = 500/0.3 = 1667 K > floor
    try:
        o = decide_coarse(m, floor_K=targets.floor, margin_M_c=Threshold(2.0, "doc11#4", True),
                          eps_sys=Threshold(0.10, "doc01#0.2b"), n_sheets=2, w_min=0.3)
    except TypeError:
        pytest.fail("decide_coarse cannot be told the candidate is multiband (no n_sheets / w_min input)")
    assert o.outcome != "SCREENED_OUT"


def test_F11_s5_missing_eps_sys_is_not_silently_zero(targets):
    """eps_sys=None is treated as 0 (one_eps = 1). Missing envelope data must not make the
    removal less conservative: A = 600 K is SCREENED_OUT with eps None, PASS with eps 0.10."""
    m = _metrics([600.0, 610.0, 605.0])
    with_eps = decide_coarse(m, floor_K=targets.floor, margin_M_c=Threshold(2.0, "doc11#4", True),
                             eps_sys=Threshold(0.10, "doc01#0.2b"))
    assert with_eps.outcome == "PASS"
    try:
        o = decide_coarse(m, floor_K=targets.floor, margin_M_c=Threshold(2.0, "doc11#4", True), eps_sys=None)
    except ValueError:
        return
    assert o.outcome != "SCREENED_OUT"


# ============================================================================ S3 / engine: non-terminal outcomes
def _relax_ok(structs, model, pressure_GPa, device=None, **kw):
    from rtt.mlff.relax import RelaxResult
    ps = pressure_GPa if isinstance(pressure_GPa, list) else [pressure_GPa] * len(structs)
    out = []
    for i, (s, p) in enumerate(zip(structs, ps)):
        n = len(s["species"])
        e = -2.0 * n
        out.append(RelaxResult(index=i, structure=s, formula="x", n_atoms=n, energy_eV=e, energy_per_atom_eV=e / n,
                               enthalpy_per_atom_eV=e / n, volume_per_atom_A3=10.0, max_force_eV_A=0.01,
                               internal_pressure_GPa=0.0, n_steps=10, converged=True, pressure_GPa=p, model=model))
    return out


def _relax_oom(structs, model, pressure_GPa, device=None, **kw):
    raise RuntimeError("CUDA out of memory")                      # transient, whole batch


def _hull(cs, pressure_GPa, **kw):
    from rtt.mlff.hull import SelfHull
    els = cs.split("-")
    return SelfHull(chemsys=cs, pressure_GPa=pressure_GPa, model="fake",
                    entries=[{"entry_id": e, "formula": e, "enthalpy_per_atom_eV": -2.0, "converged": True} for e in els],
                    elements_covered=els)


def _gen(ledger, targets, tiers):
    from rtt.records import Origin, RawStructure, StructureData
    from rtt.source.ingest import ingest_raw
    a = RawStructure(structure=StructureData(lattice=[[3, 0, 0], [0, 3, 0], [0, 0, 3]], species=["Mg", "H"],
                                             frac_coords=[[0, 0, 0], [.5, .5, .5]]), origin=Origin(lane="G1"))
    return ingest_raw(ledger, a, targets=targets, tiers=tiers)


def test_F12_s3_transient_ml_failure_is_retried_not_terminal(ledger, targets, tiers, tmp_path):
    """A whole-batch failure (CUDA OOM) makes every candidate DEFERRED(ml_failed). run_s3 only
    processes candidates WITHOUT an S3 outcome and the scheduler only advances PASS-like
    outcomes, so nothing ever re-runs them: DEFERRED is a silent permanent removal (same for
    rank_cut, relax_not_converged; and ERROR at S4/S5, e.g. matdyn NaN from imaginary modes)."""
    from rtt.mlff.screen import run_s3
    from rtt.structure.stage import run_s2
    cid = _gen(ledger, targets, tiers)
    run_s2(ledger, targets=targets, tiers=tiers, novelty=False)
    out = run_s3(ledger, targets=targets, tiers=tiers, workdir=tmp_path, alignn=False, phonons=False, dftref=False,
                 relax_fn=_relax_oom, hull_fn=_hull)
    assert out[cid].outcome == "DEFERRED"
    again = run_s3(ledger, targets=targets, tiers=tiers, workdir=tmp_path, alignn=False, phonons=False,
                   dftref=False, relax_fn=_relax_ok, hull_fn=_hull)
    assert cid in again, "DEFERRED(ml_failed) candidate is never re-processed"


def test_F13_scheduler_error_is_terminal(ledger, targets, tiers, tmp_path):
    """S4 ERROR (QE non-convergence, e.g. a frustrated spin scf) -> next_stage() returns None and
    nothing retries; the Outcome.retry field is never incremented anywhere."""
    from rtt.engine.scheduler import next_stage
    from rtt.mlff.screen import run_s3
    from rtt.records import Outcome
    from rtt.source.ingest import write_outcome
    from rtt.structure.stage import run_s2
    cid = _gen(ledger, targets, tiers)
    run_s2(ledger, targets=targets, tiers=tiers, novelty=False)
    run_s3(ledger, targets=targets, tiers=tiers, workdir=tmp_path, alignn=False, phonons=False, dftref=False,
           relax_fn=_relax_ok, hull_fn=_hull)
    assert ledger.current_outcome(cid, "S3")["outcome"] == "PASS"
    rec = ledger.candidate_record(cid)
    write_outcome(ledger, cid, Outcome("S4", "ERROR", reason="spin scf not converged"), origin_obj=rec.origin,
                  targets=targets, tiers=tiers, origin="test")
    assert next_stage(ledger, cid) == "S4", "ERROR at S4 is never retried (terminal)"


def test_F14_control_recall_counts_deferred_controls_as_passed(ledger, targets, tiers, tmp_path):
    """summary.control_recall treats every outcome outside REMOVING as 'passed'. Controls that
    are DEFERRED (ml_failed) or ERROR did not pass the screen, yet recall reads 1, so the R8
    auto-disable can never trigger from them and S9 overstates recall."""
    from rtt.mlff.screen import run_s3
    from rtt.report.summary import control_recall
    from rtt.source.controls import inject_controls
    from rtt.structure.stage import run_s2
    ids = inject_controls(ledger, targets=targets, tiers=tiers)
    run_s2(ledger, targets=targets, tiers=tiers, novelty=False)
    run_s3(ledger, targets=targets, tiers=tiers, workdir=tmp_path, alignn=False, phonons=False, dftref=False,
           relax_fn=_relax_oom, hull_fn=_hull)
    assert all(ledger.current_outcome(c, "S3")["outcome"] == "DEFERRED" for c in ids)
    assert control_recall(ledger)["S3"]["passed"] == 0


# ============================================================================ S3 threshold arithmetic
def test_F15_h7ml_margin_covers_measured_mattersim_error(tiers):
    """tiers S3.ehull_self_max_eV = 0.42 = 0.35 + 0.07. The 0.07 is doc 10 §3.3's resolution
    near zero ('indistinguishable from 0'), not an error bound: the measured per-entry MAE for
    H-containing H-Ir-Mg entries is 0.074 eV/atom (a hull energy is a difference, so its error is
    ~sqrt(2) larger), and H7 is ZPE-corrected while the ML hull is not (doc 02 §4.1: ZPE
    differences 'tens of meV/atom'). The margin is below even ONE per-entry MAE."""
    thr = float(tiers.threshold("S3.ehull_self_max_eV").value)
    window_H7, mae_H = 0.35, 0.074
    assert thr - window_H7 >= 2 * math.sqrt(2) * mae_H, f"margin {thr - window_H7:.3f} eV/atom"


# ============================================================================ S2 dedupe
def test_F16_dedupe_does_not_merge_molecular_and_dissociated_hydrogen():
    """StructureMatcher(ltol 0.2, stol 0.3, scale=True) on UNRELAXED MatterGen cells: stol is
    relative to (V/N)^(1/3) (~2.4 A here), so H-H = 0.74 A (intact H2, penalised by S4) and
    H-H = 1.60 A (dissociated, a different bonding class) are 'the same'. The later one becomes
    DUPLICATE (terminal) and its fate follows the first representative."""
    from rtt.structure.dedupe import DedupeIndex

    def mgh2(d, a=3.5):
        h = d / (2 * a)
        return {"lattice": [[a, 0, 0], [0, a, 0], [0, 0, a]], "species": ["Mg", "H", "H"],
                "frac_coords": [[0, 0, 0], [0.5 - h, 0.5, 0.5], [0.5 + h, 0.5, 0.5]]}
    idx = DedupeIndex()
    idx.add("molecular", mgh2(0.74))
    assert idx.find(mgh2(1.60)) is None


# ============================================================================ collect provenance
def test_F17_collect_hash_is_of_the_bytes_that_were_parsed(tmp_path, monkeypatch):
    """collect_run parses the spectra (wf.collect()) and only afterwards hashes the files on disk
    (run_artifact_paths -> sha256_file). A rewrite in between (adopted re-run / second matdyn in
    the same workdir) attaches a source_sha256 that is not the parsed file. The hash must be taken
    of the same bytes the parser read."""
    import types
    from rtt.dft import collect as C
    f = tmp_path / "a2F.dos1"
    f.write_text("ORIGINAL\n")
    original = hashlib.sha256(b"ORIGINAL\n").hexdigest()

    class Sp:
        omega_K, weight_K, rel_envelope = [100.0], [1.0], 0.0

        def __init__(self):
            self.provenance = {"source_file": str(f)}

    class WF:
        def collect(self):
            sp = Sp()                                   # parsed ORIGINAL
            f.write_text("REWRITTEN\n")                 # concurrent rewrite after parsing
            return types.SimpleNamespace(spectra=[sp], provenance={}, errors=[], flags=[], warnings=[])

    monkeypatch.setattr(C, "_workflow_for", lambda run: WF())
    run = DFPTRun(material_id="x", profile={"name": "coarse", "kind": "dfpt", "spectrum_tier": "coarse"},
                  workdir=str(tmp_path), status="done")
    col = C.collect_run(run)
    assert col.spectra[0].provenance["source_sha256"] == original
