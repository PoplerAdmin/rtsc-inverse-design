"""Consistency checks between parsed Spectra and the code's own reported numbers, and the
spread report across alternative spectra (smearings / q-grids), ARCH section 5.

Checks (each a ``ConsistencyCheck``; failures become flags, they never alter a Spectrum):

* ``lambda_vs_a2F_footer``  lambda recomputed from the Spectrum vs the ``lambda =`` footer of
  a2F.dos<i>. Same formula (2 sum a2F/E Delta), so they agree up to E16.6 print rounding and
  the bins at E <= 0 that we drop (we compare against QE's sum over the kept bins as well).
  A failure here means a parse/units error.
* ``lambda_vs_lambda_file`` lambda from the Spectrum vs the ``lambda`` file value (direct
  q-point sum, no histogram, printed F12.4). Differences measure histogram discretisation
  (ndos too small) -- tolerance is looser.
* ``omega_log_vs_lambda_file`` omega_log [K] likewise (lambda file prints omega_ln [K]).
* ``lambda_vs_epw``          EPW's 'Integrated el-ph coupling' vs the Spectrum.
"""
from __future__ import annotations

import math
from dataclasses import asdict, dataclass
from typing import Any, Sequence

import numpy as np

from .parsers import QEA2FDos, QELambdaEntry, spectrum_lambda, spectrum_omega_log


@dataclass
class ConsistencyCheck:
    name: str
    smearing_index: int | None
    ours: float
    reference: float
    abs_diff: float
    rel_diff: float
    tolerance_rel: float
    tolerance_abs: float
    ok: bool
    note: str = ""

    def to_dict(self) -> dict[str, Any]:
        return asdict(self)


def _cmp(name: str, isig: int | None, ours: float, ref: float, tol_rel: float,
         tol_abs: float = 0.0, note: str = "") -> ConsistencyCheck:
    ad = abs(ours - ref)
    rd = ad / abs(ref) if ref != 0 else (0.0 if ad == 0 else math.inf)
    ok = bool(math.isfinite(ours) and math.isfinite(ref) and (ad <= tol_abs or rd <= tol_rel))
    return ConsistencyCheck(name, isig, float(ours), float(ref), float(ad), float(rd),
                            float(tol_rel), float(tol_abs), ok, note)


@dataclass
class Tolerances:
    footer_lambda_rel: float = 2e-4      # parse/units check (E16.6 rounding)
    lambda_file_rel: float = 0.05        # histogram vs direct q-sum
    lambda_file_abs: float = 1e-4        # F12.4 print resolution
    omega_log_rel: float = 0.05
    epw_lambda_rel: float = 1e-3


def check_qe_spectrum(sp: Any, tab: QEA2FDos | None, entry: QELambdaEntry | None,
                      tol: Tolerances | None = None) -> list[ConsistencyCheck]:
    tol = tol or Tolerances()
    isig = sp.provenance.get("smearing_index")
    lam = spectrum_lambda(sp)
    wlog = spectrum_omega_log(sp)
    out: list[ConsistencyCheck] = []
    if tab is not None and tab.lambda_reported is not None:
        # QE footer includes all bins; compare against QE's sum restricted to kept bins too.
        d = tab.delta_Ry or float(np.median(np.diff(tab.E_Ry)))
        pos = tab.E_Ry > 0
        neg_part = float(np.sum(2.0 * tab.a2f_total[~pos] / tab.E_Ry[~pos] * d)) if np.any(~pos) else 0.0
        out.append(_cmp("lambda_vs_a2F_footer", isig, lam, tab.lambda_reported - neg_part,
                        tol.footer_lambda_rel, 1e-6,
                        note=f"QE footer lambda={tab.lambda_reported:.8g}; "
                             f"contribution of E<=0 bins removed: {neg_part:.3g}"))
    if entry is not None:
        out.append(_cmp("lambda_vs_lambda_file", isig, lam, entry.lambda_,
                        tol.lambda_file_rel, tol.lambda_file_abs,
                        note=f"broadening {entry.broadening_Ry} Ry; QE value is a direct q-sum"))
        if entry.omega_ln_K is not None and math.isfinite(wlog):
            out.append(_cmp("omega_log_vs_lambda_file", isig, wlog, entry.omega_ln_K,
                            tol.omega_log_rel, 1e-3,
                            note="QE omega_ln [K] = exp(<ln w>)/K_BOLTZMANN_RY, direct q-sum"))
    return out


def check_epw_spectrum(sp: Any, tol: Tolerances | None = None) -> list[ConsistencyCheck]:
    tol = tol or Tolerances()
    ref = sp.provenance.get("epw_reported_lambda")
    if ref is None:
        return []
    return [_cmp("lambda_vs_epw", sp.provenance.get("smearing_index"), spectrum_lambda(sp), ref,
                 tol.epw_lambda_rel, 1e-5)]


def spread_report(spectra: Sequence[Any]) -> dict[str, Any]:
    """Per-alternative lambda, omega_log, A, omega_max and their min/max/relative spread.

    Pure reporting: nothing is selected (ARCH section 5, 'never silently pick the best').
    """
    rows = []
    for sp in spectra:
        rows.append({"label": sp.label,
                     "smearing_index": sp.provenance.get("smearing_index"),
                     "broadening_Ry": sp.provenance.get("el_ph_broadening_Ry"),
                     "A_K": float(math.fsum(sp.weight_K)),
                     "lambda": spectrum_lambda(sp),
                     "omega_log_K": spectrum_omega_log(sp),
                     "omega_max_K": float(max(sp.omega_K))})
    summ: dict[str, Any] = {"n_alternatives": len(rows), "per_spectrum": rows}
    for key in ("A_K", "lambda", "omega_log_K", "omega_max_K"):
        vals = [r[key] for r in rows if math.isfinite(r[key])]
        if vals:
            lo, hi = min(vals), max(vals)
            mid = 0.5 * (lo + hi)
            summ[key] = {"min": lo, "max": hi, "rel_spread": (hi - lo) / mid if mid else 0.0}
    return summ
