"""Certificates for N1mu (docs/theory/03_certificates.md, Sec. 4): mu*-dependent A floors.

THEOREM N1mu(mu*, omega_c, A_floor): in the isotropic linearized ME model with the symmetric
Coulomb block of cutoff omega_c (docs/theory/01, Sec. 1), every finite positive measure a2F with
  A = int a2F <= A_floor
is stable at every temperature T >= 300 K.

Proof structure (per segment of constant N_c(T) = j):
  * T in (omega_c/((2j+1)pi), omega_c/((2j-1)pi)] has N_c = j and A/T < a_j := A_floor/T_lo(j),
    T_lo(j) = max(300, omega_c/((2j+1)pi)).  The Einstein family AFamily(a_j, mu*, j) is certified
    stable for all w = Omega/T > 0 (cells + analytic ends); monotonicity in a and the affinity
    reduction (01, Sec. 3.1) extend this to every a2F with A/T <= a_j.
  * T >= A_floor/4.33148145: universal wE1 certificate (A/T <= a0, any mu* >= 0, any cutoff).
usage: python run_n1mu.py --mu 0.10 --wc 20000 --delta 1e-4 --out <dir>
"""
import argparse
import json
import math
import os
import sys
import time
from decimal import Decimal, ROUND_CEILING

import numpy as np
from scipy.optimize import brentq, minimize_scalar

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(HERE, "..", "theory"))
from flint import arb  # noqa: E402
from families import AFamily  # noqa: E402
from cover import cover, dyadic_floor, dyadic_ceil  # noqa: E402
from run_nm import dec_floor  # noqa: E402

A0_UNIV = "4.33148145"   # wE1 universal certificate (R(a0, 1e-9) for all a <= a0, all mu*, N_c)
TT = 300.0


def Nc_of(T, wc):
    return int(math.floor((wc / (math.pi * T) + 1) / 2))


def rho_f(a, w, mu, j):
    from eliashberg import einstein, rho
    return rho(einstein(2 * a / w, w), 1.0, mu, 2 * j * math.pi if j > 0 else 0.5)


def a_crit_float(mu, j):
    """float min over w of the critical a (rho = 0) with N_c = j at T = 1."""
    def ac(w):
        return brentq(lambda a: rho_f(a, w, mu, j), 1.0, 40.0, xtol=1e-12)
    ws = np.linspace(3.0, 14.0, 45)
    vals = [ac(w) for w in ws]
    i = int(np.argmin(vals))
    r = minimize_scalar(ac, bounds=(ws[max(i - 1, 0)], ws[min(i + 1, len(ws) - 1)]),
                        method="bounded", options={"xatol": 1e-9})
    return r.fun, r.x


def dec_ceil(x, sig=10):
    d = Decimal(repr(float(x)))
    q = Decimal(1).scaleb(d.adjusted() - sig + 1)
    return str(d.quantize(q, rounding=ROUND_CEILING))


def ends_A(a, eta):
    """small w: Q <= (a Lang(w/2) - 1)|x|^2 (Lemma 7, Lang increasing);
    large w: weighted Schur test, h_n = H (n < n0), 1 (n >= n0), l(k) <= eps = 2a/W:
       rows n < n0:  R <= eps - 1 + 2 (n0-1) eps + (a + (n0-1) eps)/H
       rows n >= n0: R <= a - (2 n0 + 1) + 2 (H-1) n0 eps          (Sec. 4.2)."""
    af = float(a.upper())
    ws = dyadic_floor(0.95 * 6 / af, bits=4)
    small = a * ((arb(ws) / 2).coth() - 2 / arb(ws)) - 1
    n0 = int(math.ceil(af / 2)) + 1
    H = arb(int(math.ceil(2 * af)))
    W = 2.0 ** 11
    eps = 2 * a / arb(W)
    r_lo = eps - 1 + 2 * (n0 - 1) * eps + (a + (n0 - 1) * eps) / H
    r_hi = a - (2 * n0 + 1) + 2 * (H - 1) * n0 * eps
    ok = bool(small < -eta) and bool(r_lo < -eta) and bool(r_hi < -eta)
    return ws, W, {"w_s": ws, "small_bound": float(small.upper()), "W": W, "n0": n0,
                   "H": float(H.mid()), "row_lt_n0": float(r_lo.upper()),
                   "row_ge_n0": float(r_hi.upper()), "ok": ok}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--mu", required=True)
    ap.add_argument("--wc", type=float, required=True)
    ap.add_argument("--delta", type=float, default=1e-4)
    ap.add_argument("--N", type=int, default=120)
    ap.add_argument("--K", type=int, default=600)
    ap.add_argument("--eta", default="1e-7")
    ap.add_argument("--out", required=True)
    o = ap.parse_args()
    t0 = time.time()
    mu = float(o.mu)
    wc = o.wc
    eta = arb(o.eta)
    j300 = Nc_of(TT, wc)
    # segments j = j300, j300-1, ... with T_lo(j)
    segs = []
    for j in range(j300, 0, -1):
        Tlo = TT if j == j300 else wc / ((2 * j + 1) * math.pi)
        ac, wopt = a_crit_float(mu, j)
        segs.append(dict(j=j, T_lo=Tlo, a_crit_float=ac, w_opt_float=wopt, A_req_float=Tlo * ac))
        # later segments have T_lo * a_crit >= T_lo * a0 > current minimum: they cannot lower the
        # floor and are covered by the universal bound
        if Tlo * 4.33148145 > min(x["A_req_float"] for x in segs):
            break
    A_float = min(s["A_req_float"] for s in segs)
    A_floor_s = dec_floor((1 - o.delta) * A_float, sig=8)
    A_floor = arb(A_floor_s)
    T_u = A_floor / arb(A0_UNIV)
    os.makedirs(os.path.join(o.out, "logs"), exist_ok=True)
    tag = "mu%s_wc%d" % (o.mu, int(wc))
    certs = []
    all_ok = True
    for s in segs:
        j = s["j"]
        # this segment matters only if T_lo(j) < T_u (else the universal bound covers it)
        if not (arb(s["T_lo"]) < T_u):
            s["covered_by"] = "universal a0 bound (T_lo >= A_floor/a0)"
            certs.append(s)
            continue
        if j == j300:
            a_j = dec_ceil(float(A_floor_s) / TT)
        else:   # T_lo exact: omega_c/((2j+1) pi); a_j >= A_floor/T_lo rounded up (checked in arb)
            a_j = dec_ceil(float(A_floor_s) * (2 * j + 1) * math.pi / wc * (1 + 1e-12))
        Tlo_b = arb(TT) if j == j300 else arb(wc) / ((2 * j + 1) * arb.pi())
        assert bool(arb(a_j) * Tlo_b >= A_floor), "a_j rounding"
        fam = AFamily(a_j, o.mu, j)
        ws, W, ends = ends_A(arb(a_j), eta)
        with open(os.path.join(o.out, "logs", "n1mu_%s_j%d.log" % (tag, j)), "w") as log:
            log.write("# N1mu mu=%s wc=%g j=%d a=%s N=%d K=%d range=[%r,%r]\n"
                      % (o.mu, wc, j, a_j, o.N, o.K, ws, W))
            ok_c, cells, failed, worst, tsec = cover(fam, eta, ws, W, o.N, o.K, log=log)
        s.update(a_cert=a_j, ends=ends, cells_ok=ok_c, n_cells=len(cells), failed_cells=failed,
                 min_cell_margin_scaled=worst, cells=cells, seconds=round(tsec, 1))
        all_ok &= ok_c and ends["ok"]
        certs.append(s)
        print("%s j=%d a=%s ok=%s ends=%s cells=%d failed=%d %.0fs" % (tag, j, a_j, ok_c, ends["ok"],
                                                                       len(cells), len(failed), tsec),
              flush=True)
    res = {"rule": "N1mu", "mu_star": o.mu, "omega_c_K": wc, "T_target_K": TT, "delta": o.delta,
           "A_floor_K": A_floor_s, "A_float_K": A_float, "N_c_at_300": j300,
           "universal_a0": A0_UNIV, "T_universal_K": float(T_u.lower()),
           "eta": o.eta, "N": o.N, "K": o.K, "segments": certs,
           "verdict": "CERTIFIED" if all_ok else "NOT CERTIFIED",
           "seconds": round(time.time() - t0, 1), "prec_bits": 128}
    with open(os.path.join(o.out, "n1mu_%s.json" % tag), "w") as fh:
        json.dump(res, fh, indent=1)
    print("%s A_floor=%s verdict=%s" % (tag, A_floor_s, res["verdict"]), flush=True)


if __name__ == "__main__":
    main()
