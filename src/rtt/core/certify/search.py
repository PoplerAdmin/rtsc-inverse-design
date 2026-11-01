"""Search drivers: choose cells / Rayleigh temperatures (float steering) and obtain certificates
(every accepted certificate has passed the Arb check).  Nothing here decides a status; the
decision is made in classify.py from evidence re-derived by evidence.py.
"""
from __future__ import annotations

import math
import time

from flint import fmpq

from ..constants import A0_FRACTION
from .arbcore import ArbSpectrum, q_of_float
from .cells import LamCache, certify_cell, head_size
from .floatcore import FloatSpectrum, float_tc
from .rayleigh import check_rayleigh, propose_rayleigh

A0_Q = fmpq(A0_FRACTION.numerator, A0_FRACTION.denominator)


def t_end_for(asp: ArbSpectrum) -> float:
    """Smallest-ish double t_end with A_hi <= a0 * t_end (exact rational check, ERRATA D1)."""
    target = asp.A_hi_q / A0_Q
    t = int(target.p) / int(target.q)          # correctly rounded big-int division
    while not (A0_Q * q_of_float(t) >= asp.A_hi_q):
        t = math.nextafter(t, math.inf)
    return t


def certify_upper(asp: ArbSpectrum, fs: FloatSpectrum, variant: str = "mhat",
                  rel_tol: float = 1e-4, T_star: float | None = None, refine: bool = True,
                  N_cap: int = 2400, time_budget: float = 300.0, log=None) -> dict:
    """Top-down cover of [T_up, t_end] by certified cells (all cells pass check_cell).

    Phase 1 (decision): if T* < t_end, first try to reach just below T* -- a single cell
    [T*(1-2^-30), t_end] is attempted before any float Tc is computed (cheap when Tc << T*).
    Phase 2 (refine=True): continue down to Tc_float(1+e-scaled)*(1+rel_tol) within the time
    budget.  Returns {"t_end","T_up","cells" (ascending),"tc_float_hi","seconds","attempts"}."""
    t0 = time.time()
    t_end = t_end_for(asp)
    e = asp.e_f
    out = {"t_end": t_end, "T_up": t_end, "cells": [], "tc_float_hi": None, "attempts": 0}
    cache = LamCache(asp)
    cells_desc: list = []
    cur = t_end

    def attempt(lo):
        nonlocal cur
        out["attempts"] += 1
        cert, res = certify_cell(asp, fs, lo, cur, variant, cache, N_cap=N_cap)
        if cert is not None:
            cert["check"] = {k: res[k] for k in ("worst_rel", "G_tail", "seconds") if k in res}
            cells_desc.append(cert)
            if log:
                log(f"  cell ok [{lo:.10g}, {cur:.10g}] N={cert['N']} N_e={cert['N_e']} "
                    f"M_v={cert['M_v']} ({res['seconds']:.2f}s)")
            cur = lo
            return True
        if log:
            log(f"  cell FAIL [{lo:.10g}, {cur:.10g}]: {res and res.get('reason')}")
        return False

    # ---- phase 1: one direct decision cell
    if T_star is not None and 0 < T_star < t_end:
        lo1 = T_star * (1 - 2.0 ** -30)
        attempt(lo1)
    if not refine and (T_star is None or cur < T_star):
        out["cells"] = list(reversed(cells_desc))
        out["T_up"] = cur
        out["seconds"] = time.time() - t0
        return out
    # ---- phase 2: float-steered refinement
    tcf = float_tc(fs, factor=1 + e, T_hi=cur, T_min=1e-6 * t_end, rtol=max(1e-12, 1e-3 * rel_tol))
    out["tc_float_hi"] = tcf
    targets = []
    if T_star is not None and tcf < T_star < cur:
        targets.append(T_star * (1 - 2.0 ** -30))       # decision first
    if refine:
        targets.append(tcf * (1 + rel_tol) if tcf > 0 else 1e-3 * t_end)
    for target in targets:
        s = 0.5
        while cur > target and time.time() - t0 < time_budget:
            base = tcf if tcf > 0 else 0.0
            lo = max(target, base + (cur - base) * (1 - s))
            if not lo < cur:
                break
            if attempt(lo):
                s = min(0.9, s * 1.6)
            else:
                s = s / 3
                if s < 1e-3:
                    break
    out["cells"] = list(reversed(cells_desc))
    out["T_up"] = cur
    out["seconds"] = time.time() - t0
    return out


def certify_lower(asp: ArbSpectrum, fs: FloatSpectrum, mu_hi: float, omega_c: float,
                  rel_tol: float = 1e-4, T_floor: float | None = None, T_hi: float | None = None,
                  T_min: float | None = None, time_budget: float = 300.0, log=None) -> dict:
    """Largest T (found) with a certified C-RAYLEIGH certificate.  Returns
    {"cert" or None, "T_low", "tc_float_lo", "seconds"}."""
    t0 = time.time()
    e = asp.e_f
    if T_hi is None:
        T_hi = t_end_for(asp)
    if T_min is None:
        T_min = max(1e-3 * T_hi, 0.05)
    tcl = float_tc(fs, factor=1 - e, mu=mu_hi, omega_c=omega_c, T_hi=T_hi, T_min=T_min,
                   rtol=max(1e-12, 1e-3 * rel_tol))
    out = {"cert": None, "T_low": None, "tc_float_lo": tcl, "attempts": 0}
    if tcl <= 0:
        out["seconds"] = time.time() - t0
        return out
    d = rel_tol
    while d < 0.6 and time.time() - t0 < time_budget:
        T = tcl * (1 - d)
        if T_floor is not None and T < T_floor:
            break
        out["attempts"] += 1
        c = propose_rayleigh(fs, T, mu_hi, omega_c)
        if c is not None:
            res = check_rayleigh(asp, c)
            if res["ok"]:
                c["check"] = {"Q_over_x2": res["Q_over_x2"], "seconds": res["seconds"]}
                out["cert"] = c
                out["T_low"] = T
                if log:
                    log(f"  rayleigh ok at T={T:.10g} (d={d:.1e}) Q/|x|^2={res['Q_over_x2']:.3e}")
                break
            if log:
                log(f"  rayleigh fail at T={T:.10g}: {res.get('reason')}")
        d *= 4
    out["seconds"] = time.time() - t0
    return out
