"""Smearing admission 'mode_positive' (rtt.core.dfpt.admit_smearings).

alpha^2F_nu(omega) >= 0 for every phonon mode nu (|g|^2 times non-negative double-delta weights), so
per-mode negative alpha^2F marks a smearing the k grid does not resolve. These tests pin:
  * the per-mode measure sees ringing that cancels in the total,
  * the admitted set is the contiguous block from the smallest passing smearing,
  * a failure above a passing smearing, or no passing smearing, is an error,
  * the blocking tag is direction-correct (CERTIFIED always, REJECTED only if an excluded
    smearing had the larger lambda),
  * on the REAL QE 7.6 Pb and Al benchmark outputs the rule excludes exactly the smearings with
    negative per-mode alpha^2F (Pb: sigma1, sigma2; Al: sigma1).
"""
from __future__ import annotations

import re
import sys
from pathlib import Path

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
sys.path.insert(0, str(Path(__file__).resolve().parent))

from rtt.core.dfpt import (admit_smearings, collect_qe_eph, mode_negative_fractions,  # noqa: E402
                       parse_qe_a2f_dos)
from rtt.core.certify_policy import blocking_tags  # noqa: E402
from test_dfpt_parsers import toy_a2f, write_qe_a2f_dos  # noqa: E402

DATA = Path(__file__).parent / "data_dfpt"
DELTA = 2e-5
E = np.arange(200) * DELTA + 0.5 * DELTA


def modes(ripple: float = 0.0, scale: float = 1.0) -> np.ndarray:
    """Two SYNTHETIC modes; `ripple` adds a low-omega negative lobe to mode 2 and the same positive
    lobe to mode 1, so the TOTAL stays exactly the ripple-free total."""
    base = toy_a2f(E) * scale
    lobe = ripple * np.exp(-0.5 * ((E - 3e-4) / 1e-4) ** 2)
    return np.column_stack([0.5 * base + lobe, 0.5 * base - lobe])


def write_ladder(d: Path, ripples, scales=None) -> None:
    d.mkdir(parents=True, exist_ok=True)
    scales = scales or [1.0] * len(ripples)
    lines = []
    for i, (r, s) in enumerate(zip(ripples, scales), start=1):
        lam = write_qe_a2f_dos(d / f"a2F.dos{i}", 0.0, DELTA, modes(r, s))
        lines.append(f" Broadening   {0.01 * i:.4f} lambda       {lam:.4f} dos(Ef)  5.0000 "
                     f"omega_ln [K]   100.0000")
    (d / "lambda").write_text("\n".join(lines) + "\n")


def test_mode_measure_sees_cancelling_ringing(tmp_path):
    write_qe_a2f_dos(tmp_path / "a2F.dos1", 0.0, DELTA, modes(0.05))
    tab = parse_qe_a2f_dos(tmp_path / "a2F.dos1")
    assert np.all(tab.a2f_total >= 0)                 # the total alone hides it
    f = mode_negative_fractions(tab)
    assert f[0] == 0.0 and f[1] > 1e-3


def test_ladder_excludes_low_block(tmp_path):
    write_ladder(tmp_path, [0.05, 0.02, 0.0, 0.0])
    tabs = {i: parse_qe_a2f_dos(tmp_path / f"a2F.dos{i}") for i in range(1, 5)}
    adm, rep, err = admit_smearings(tabs)
    assert err is None and adm == [3, 4] and rep["excluded"] == [1, 2]


@pytest.mark.parametrize("ripples,msg", [([0.0, 0.05, 0.0], "non-contiguous"),
                                         ([0.05, 0.05], "no smearing")])
def test_admission_errors(tmp_path, ripples, msg):
    write_ladder(tmp_path, ripples)
    tabs = {i: parse_qe_a2f_dos(tmp_path / f"a2F.dos{i}") for i in range(1, len(ripples) + 1)}
    adm, _, err = admit_smearings(tabs)
    assert adm == [] and msg in err


@pytest.mark.parametrize("scales,blocks", [([0.8, 1.0, 1.1], (False, True)),    # below the boxes
                                           ([1.3, 1.0, 1.1], (True, False)),    # above the boxes
                                           ([0.95, 1.0, 1.1], (False, False))])  # inside
def test_collect_tags_and_flags(tmp_path, scales, blocks):
    write_ladder(tmp_path, [0.05, 0.0, 0.0], scales)
    res = collect_qe_eph(tmp_path, "SYNTHETIC-adm", el_ph_nsigma=3, smearing_admission="mode_positive")
    assert res.ok, res.errors
    assert [s.provenance["smearing_index"] for s in res.spectra] == [2, 3]
    assert any("only smearing indices [2, 3]" in w for w in res.warnings)
    tags = [a for s in res.spectra for a in s.provenance["assumptions"] if a.startswith("A-SIGMA-EXCLUDED(sigma1:lam=")]
    assert len(tags) == len(res.spectra) and len(set(tags)) == 1
    assert res.provenance.get("smearing_admission", res.spectra[0].provenance
                              ["smearing_admission"])["excluded"] == [1]
    rej, cert = blocking_tags([s.with_envelope(0.10) for s in res.spectra])
    assert (bool(rej), bool(cert)) == blocks


def test_collect_admission_failure_is_error(tmp_path):
    write_ladder(tmp_path, [0.0, 0.05])
    res = collect_qe_eph(tmp_path, "SYNTHETIC-adm", el_ph_nsigma=2, smearing_admission="mode_positive")
    assert not res.ok and res.spectra == [] and "SMEARING_ADMISSION_FAIL" in res.flags


def test_default_rule_unchanged(tmp_path):
    write_ladder(tmp_path, [0.0, 0.0])
    res = collect_qe_eph(tmp_path, "SYNTHETIC-adm", el_ph_nsigma=2)
    assert [s.provenance["smearing_index"] for s in res.spectra] == [1, 2]
    assert "smearing_admission" not in res.spectra[0].provenance
    with pytest.raises(ValueError, match="mutually exclusive"):
        collect_qe_eph(tmp_path, "x", el_ph_nsigma=2, use_smearings=[2],
                       smearing_admission="mode_positive")


# Without an envelope the excluded smearings point outside the admitted set (Pb: lower lambda,
# blocks CERTIFIED; Al: sigma1's lambda 0.428 > 0.415, blocks REJECTED). With the project's
# envelope e = 0.10 both are inside the admitted boxes, so nothing is blocked.
@pytest.mark.parametrize("d,excluded,blocks_e0", [("qe76_pb_bench", [1, 2], (False, True)),
                                                  ("qe76_al_bench", [1], (True, False))])
def test_real_qe76_benchmarks(d, excluded, blocks_e0):
    tabs = {i: parse_qe_a2f_dos(DATA / d / f"a2F.dos{i}") for i in range(1, 11)}
    neg = {i: max(mode_negative_fractions(t)) for i, t in tabs.items()}
    assert [i for i in neg if neg[i] > 0] == excluded       # negatives only in the narrowest
    res = collect_qe_eph(DATA / d, f"REAL-{d}", el_ph_nsigma=10, smearing_admission="mode_positive")
    assert res.ok, res.errors
    assert [s.provenance["smearing_index"] for s in res.spectra] == \
        [i for i in range(1, 11) if i not in excluded]
    tag = res.spectra[0].provenance["assumptions"][0]
    assert tag.startswith("A-SIGMA-EXCLUDED(") and \
        [int(x) for x in re.findall(r"sigma(\d+):lam", tag)] == excluded
    lam1 = float(re.search(r"sigma\d+:lam=([-+0-9.eE]+)", tag).group(1))
    assert lam1 == pytest.approx(tabs[excluded[0]].lambda_reported, rel=1e-4)   # QE's own footer
    assert res.checks_ok          # QE footer and lambda-file consistency hold on every admitted one
    rej, cert = blocking_tags(res.spectra)
    assert (bool(rej), bool(cert)) == blocks_e0
    assert blocking_tags([s.with_envelope(0.10) for s in res.spectra]) == ([], [])
