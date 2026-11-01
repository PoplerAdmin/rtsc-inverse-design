"""SYNTHETIC demo candidates (no DFT).  These are invented spectra, NOT data for any compound.

Every candidate here has ``synthetic=True``, an id starting with ``SYNTHETIC-``, the formula
``(synthetic)``, an empty composition, and spectra labelled ``SYNTHETIC``.  Reports and summaries
print a banner for them.  They exist to exercise the pipeline end to end: one should be rejected by
the universal ceiling, one by the per-material upper bound, one should be strong enough to be
certified, one borderline, two excluded by (synthetic) existence criteria, one deferred by budget.
"""
from __future__ import annotations

import math

from .base import Candidate

SYNTH_NOTE = ("SYNTHETIC demo spectrum invented for testing the pipeline; it is not data for any "
              "real material")


def gaussian_peak_spectrum(center_K: float, width_K: float, A_K: float, n: int = 24,
                           label: str = "SYNTHETIC", rel_envelope: float = 0.0) -> dict:
    """Spectrum dict (Spectrum.to_dict layout) of a discretised Gaussian bump with total mass A."""
    lo = max(center_K - 3 * width_K, 1.0)
    hi = center_K + 3 * width_K
    xs = [lo + (hi - lo) * (i + 0.5) / n for i in range(n)]
    ws = [math.exp(-0.5 * ((x - center_K) / width_K) ** 2) for x in xs]
    s = sum(ws)
    ws = [A_K * w / s for w in ws]
    return {"omega_K": xs, "weight_K": ws, "rel_envelope": rel_envelope, "label": label,
            "provenance": {"synthetic": True, "note": SYNTH_NOTE,
                           "shape": f"gaussian(center={center_K} K, width={width_K} K, A={A_K} K, n={n})"}}


def demo_candidates() -> list[Candidate]:
    def cand(name, spectra, props, note):
        return Candidate(material_id=f"SYNTHETIC-{name}", composition={}, formula="(synthetic)",
                         source="synthetic_demo", spectra=spectra, properties=props,
                         synthetic=True,
                         provenance={"source": "synthetic_demo", "note": SYNTH_NOTE,
                                     "purpose": note})

    return [
        cand("weak-einstein",
             [gaussian_peak_spectrum(500.0, 30.0, 200.0, label="SYNTHETIC weak")],
             {"is_metal": True, "energy_above_hull": 0.0, "dos_at_fermi": 0.5},
             "A far below a0*T*: universal ceiling should reject"),
        cand("high-freq-moderate-coupling",
             [gaussian_peak_spectrum(6000.0, 300.0, 1400.0, label="SYNTHETIC high-freq")],
             {"is_metal": True, "energy_above_hull": 0.0, "dos_at_fermi": 0.6},
             "A above a0*T* but lambda ~ 0.47: per-material upper bound (R-CELLS) should reject"),
        cand("strong-coupling",
             [gaussian_peak_spectrum(2000.0, 120.0, 3000.0, label="SYNTHETIC strong, smearing a"),
              gaussian_peak_spectrum(2050.0, 140.0, 2900.0, label="SYNTHETIC strong, smearing b")],
             {"is_metal": True, "energy_above_hull": 0.0, "dos_at_fermi": 1.5},
             "lambda ~ 3 at 2000 K: expected certifiable; two alternative spectra (S5)"),
        cand("borderline",
             [gaussian_peak_spectrum(4000.0, 200.0, 1700.0, label="SYNTHETIC borderline")],
             {"is_metal": True, "energy_above_hull": 0.0, "dos_at_fermi": 1.0},
             "Tc_model near T* depending on mu*: expected UNDECIDED"),
        cand("above-hull",
             [gaussian_peak_spectrum(2000.0, 120.0, 3000.0, label="SYNTHETIC above-hull")],
             {"is_metal": True, "energy_above_hull": 0.35, "dos_at_fermi": 1.2},
             "synthetic energy above hull 0.35 eV/atom: hull filter should exclude"),
        cand("insulator",
             [gaussian_peak_spectrum(1500.0, 100.0, 1000.0, label="SYNTHETIC insulator")],
             {"is_metal": False, "band_gap": 1.2, "energy_above_hull": 0.0},
             "synthetic band gap 1.2 eV: metallicity filter should exclude"),
        cand("low-priority",
             [gaussian_peak_spectrum(800.0, 50.0, 300.0, label="SYNTHETIC low-priority")],
             {"is_metal": True, "energy_above_hull": 0.0, "dos_at_fermi": 0.01},
             "lowest heuristic priority: budget should defer it (never computed)"),
    ]
