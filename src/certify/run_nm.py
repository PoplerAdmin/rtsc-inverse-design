"""Certificates for rule NM-G (docs/theory/03_certificates.md, Sec. 2-3).

For each grid value lam_k (or lam = inf) this proves, in Arb ball arithmetic:

  THEOREM NM(lam_k, m0):  every finite positive measure a2F on (0, inf) with
      lambda[a2F] <= lam_k   and   lambda<omega^2>[a2F] <= m0 * T^2
  is stable at temperature T (rho_T < 0, both sectors, every mu* >= 0, every cutoff).

Method: Lagrangian / Jensen reduction to the one-parameter family of Einstein atoms with fixed
moment m0 (MFamily) plus a diagonal shift beta*(lam_k - m0/w^2); certified on dyadic w-cells,
with analytic end ranges. The upper side (Einstein(lam_k, w_hi) unstable) is a Rayleigh
certificate, so the true Einstein frequency Omega_E(lam_k) lies in [Omega_lo, Omega_hi].

usage: python run_nm.py --lams 0.8,1.0,...,inf --delta 1e-4 --out <dir>   (writes <dir>/nm_<lam>.json)
"""
import argparse
import json
import math
import os
import sys
import time
from decimal import Decimal, ROUND_FLOOR

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(HERE, "..", "theory"))

from flint import arb  # noqa: E402
from families import MFamily  # noqa: E402
from cover import cover, dyadic_floor, dyadic_ceil  # noqa: E402
from rayleigh import einstein_unstable  # noqa: E402

TT = 300


def fdown(x):
    """float <= the arb number x (directed rounding: one ulp down after nearest rounding)."""
    return math.nextafter(float(x), -math.inf)


def fup(x):
    return math.nextafter(float(x), math.inf)


# ------------------------------------------------------------------ float proposals
def float_wE(lam):
    from eliashberg import einstein, tc
    return 1.0 / tc(einstein(lam, 1.0), xtol_rel=1e-12)


def float_c2():
    from scipy.optimize import brentq
    from bounds_extra import rho_soft
    Bc = brentq(lambda B: rho_soft(B, 1.0, 0.0, math.inf, 400), 0.1, 5.0, xtol=1e-14)
    return 1 / (2 * math.pi * math.sqrt(Bc))


def float_psi(m0, w):
    from eliashberg import einstein, rho
    return rho(einstein(m0 / w ** 2, w), 1.0)


def propose_beta(lam, m0, ws):
    ws_star = math.sqrt(m0 / lam)
    ws = np.sort(np.concatenate([ws, np.linspace(0.93 * ws_star, 1.07 * ws_star, 4001)]))
    ps = np.array([float_psi(m0, w) for w in ws])
    lp = m0 / ws ** 2
    F = lambda b: float(np.max(ps + b * (lam - lp)))
    lo, hi = 0.0, 0.95 / lam
    for _ in range(80):
        a = lo + (hi - lo) / 3
        b = hi - (hi - lo) / 3
        if F(a) < F(b):
            hi = b
        else:
            lo = a
    beta = 0.5 * (lo + hi)
    return beta, F(beta)


def dec_floor(x, sig=10):
    d = Decimal(repr(float(x)))
    q = Decimal(1).scaleb(d.adjusted() - sig + 1)
    return str(d.quantize(q, rounding=ROUND_FLOOR))


# ------------------------------------------------------------------ analytic ends
def ends_check(fam, eta, ws, W):
    """Lemma E (Sec. 2.3): unit-weight Schur bound x^T M x <= (2L - 1)|x|^2 for any Einstein atom,
    with L <= m0/(4w) (large w) and L <= m0/24 (all w). Returns dict with arb upper bounds."""
    m0, beta = fam.m0, fam.beta
    lam = fam.lam if fam.lam is not None else arb(0)
    out = {}
    large = m0 / (2 * arb(W)) - 1 + beta * lam
    out["large_w"] = {"W": W, "bound": float(large.upper()), "ok": bool(large < -eta)}
    if ws > 0:
        small = m0 / 12 - 1 + beta * lam - beta * m0 / (arb(ws) * arb(ws))
        out["small_w"] = {"w_s": ws, "bound": float(small.upper()), "ok": bool(small < -eta)}
    else:
        out["small_w"] = {"w_s": 0.0, "bound": None, "ok": True, "note": "cells start at w = 0"}
    # moment-0 endpoint (affinity in m, Sec. 2.2): -D + beta*lam <= (beta*lam - 1) I
    zero = beta * lam - 1
    out["m_zero_endpoint"] = {"bound": float(zero.upper()), "ok": bool(zero < -eta)}
    out["ok"] = out["large_w"]["ok"] and out["small_w"]["ok"] and out["m_zero_endpoint"]["ok"]
    return out


def choose_ends(m0, lam, beta, eta):
    W = dyadic_ceil(1.02 * m0 / (2 * (1 - beta * (lam or 0) - 4 * eta)), bits=4)
    if lam is None or beta == 0:
        return 0.0, W
    den = m0 / 12 - 1 + beta * lam + 4 * eta
    ws2 = beta * m0 / den if den > 0 else 1e9
    ws = dyadic_floor(0.98 * math.sqrt(ws2), bits=4)
    return ws, W


# ------------------------------------------------------------------ one lambda
def certify_lambda(lam_s, delta, N, K, eta_s, outdir, up_rel, max_time=None):
    t0 = time.time()
    is_inf = lam_s == "inf"
    if is_inf:
        c2 = float_c2()
        mE = 1 / c2 ** 2
        lam = None
    else:
        lam = float(lam_s)
        wE = float_wE(lam)
        mE = lam * wE ** 2
    m0_s = dec_floor((1 - delta) * mE)
    m0 = float(m0_s)
    eta = arb(eta_s)
    if is_inf:
        beta, Fb = 0.0, float_psi(m0, 1e-3)
    else:
        ws_grid = np.exp(np.linspace(np.log(0.05 * math.sqrt(m0 / lam)), np.log(m0), 1500))
        beta, Fb = propose_beta(lam, m0, ws_grid)
        beta = math.floor(beta * 2 ** 30) / 2 ** 30      # exact dyadic
    fam = MFamily(m0_s, None if is_inf else lam_s, repr(beta))
    ws, W = choose_ends(m0, lam, beta, float(eta_s))
    ends = ends_check(fam, eta, ws, W)
    logpath = os.path.join(outdir, "logs", "nm_%s.log" % lam_s)
    os.makedirs(os.path.dirname(logpath), exist_ok=True)
    with open(logpath, "w") as log:
        log.write("# NM lam=%s m0=%s beta=%r eta=%s N=%d K=%d range=[%r,%r] float F(beta)=%.3e\n"
                  % (lam_s, m0_s, beta, eta_s, N, K, ws, W, Fb))
        ok_c, cells, failed, worst, tc_ = cover(fam, eta, ws, W, N, K, log=log, max_time=max_time)
        log.write("cells ok=%s n=%d failed=%d t=%.0fs\n" % (ok_c, len(cells), len(failed), tc_))
    # upper side: Einstein(lam, w_hi) unstable
    if is_inf:
        m1 = (1 + up_rel) ** 2 * mE
        w_hi = 2.0 ** -6
        lam_up = dec_floor(m1 / w_hi ** 2 * (1 + 1e-9) + 1e-9, sig=12)
        ok_u, r_u, _ = einstein_unstable(lam_up, w_hi, N=300)
        upper = {"spectrum": "Einstein", "lambda": lam_up, "w": w_hi,
                 "m_over_T2": fup((arb(lam_up) * w_hi ** 2).upper()), "Q_over_x2_lower": r_u,
                 "certified_unstable": ok_u,
                 "g_hi_K2": fup((arb(lam_up) * w_hi ** 2 * TT ** 2).upper())}
    else:
        w_hi = dyadic_ceil(wE * (1 + up_rel), bits=30)
        ok_u, r_u, _ = einstein_unstable(lam_s, w_hi, N=300 if lam > 20 else 200)
        upper = {"spectrum": "Einstein", "lambda": lam_s, "w": w_hi, "Omega_hi_K": TT * w_hi,
                 "Q_over_x2_lower": r_u, "certified_unstable": ok_u,
                 "g_hi_K2": fup((arb(lam_s) * (TT * w_hi) ** 2).upper())}
    g_lo = arb(m0_s) * TT ** 2
    res = {
        "rule": "NM-G", "lambda": lam_s, "delta": delta,
        "m0_over_T2": m0_s, "g_lo_K2": fdown(g_lo.lower()), "g_lo_exact": "90000 * " + m0_s,
        "Omega_lo_K": None if is_inf else fdown((TT * (arb(m0_s) / arb(lam_s)).sqrt()).lower()),
        "float_g_K2": mE * TT ** 2,
        "beta": repr(beta), "beta_exact_dyadic": True, "eta": eta_s, "N": N, "K": K,
        "w_range": [ws, W], "ends": ends,
        "cells_ok": ok_c, "n_cells": len(cells), "failed_cells": failed,
        "min_cell_margin_scaled": worst, "cells": cells,
        "upper": upper,
        "verdict": "CERTIFIED" if (ok_c and ends["ok"]) else "NOT CERTIFIED",
        "seconds": round(time.time() - t0, 1),
        "prec_bits": 128,
    }
    with open(os.path.join(outdir, "nm_%s.json" % lam_s), "w") as fh:
        json.dump(res, fh, indent=1)
    print("lam=%s verdict=%s cells=%d failed=%d g_lo=%.6e upper_ok=%s t=%.0fs"
          % (lam_s, res["verdict"], len(cells), len(failed), res["g_lo_K2"], ok_u, res["seconds"]),
          flush=True)
    return res


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--lams", required=True)
    ap.add_argument("--deltas", default="1e-6,1e-5,1e-4,1e-3",
                    help="relative slack in m0; tried in order until certified")
    ap.add_argument("--max_time", type=float, default=1800.0)
    ap.add_argument("--up", type=float, default=1e-4, help="relative w excess for the upper side")
    ap.add_argument("--N", type=int, default=120)
    ap.add_argument("--K", type=int, default=600)
    ap.add_argument("--eta", default="1e-7")
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)
    for lam_s in a.lams.split(","):
        for d in [float(x) for x in a.deltas.split(",")]:
            r = certify_lambda(lam_s, d, a.N, a.K, a.eta, a.out, a.up, a.max_time)
            if r["verdict"] == "CERTIFIED":
                break


if __name__ == "__main__":
    main()
