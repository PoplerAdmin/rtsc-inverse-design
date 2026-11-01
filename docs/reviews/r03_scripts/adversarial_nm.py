"""R03 review: adversarial float checks of the NM-G and N1mu certificates with the INDEPENDENT
solver of docs/reviews/r01_scripts/indep_solver.py (full Matsubara axis, non-symmetric eigenproblem).

For each tested lambda_k, spectra are built with lambda <= lambda_k and m = lambda<omega^2> just
below g_lo(lambda_k) (T = 300 K), and the top eigenvalue of K diag(1/|w~|) is evaluated at
T = 300, 330, 450 K (stable iff < 1). Shapes: Einstein, two-atom (wide frequency ratios, all
weight splits), three-atom, lambda < lambda_k, soft-mode + hard, box spectra; plus a Nelder-Mead
maximisation of the top eigenvalue over two- and three-atom spectra on the constraint surface.
N1mu: Einstein and two-atom spectra with A = 0.9999 A_floor at mu*, sym convention, T on a grid
>= 300 K incl. segment boundaries. Corollary: Tc <= C_cert sqrt(m) for large-lambda spectra.

usage: python adversarial_nm.py > adversarial_nm.log
"""
import json
import math
import os
import sys

import numpy as np
from scipy.optimize import minimize

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "r01_scripts"))
from indep_solver import top_eig, tc_indep  # noqa: E402

CERT = os.path.join(HERE, "..", "..", "theory", "certificates")
TT = 300.0
worst_global = {"ev": -1, "what": None}


def glo(lam_s):
    d = json.load(open(os.path.join(CERT, "nm_%s.json" % lam_s)))
    return float(d["m0_over_T2"]) * TT ** 2


def spec_from(lams, oms):
    """atoms with partial couplings lams at frequencies oms -> (om, w) with w = lam*om/2."""
    oms = np.asarray(oms, float)
    lams = np.asarray(lams, float)
    return oms, lams * oms / 2.0


def m_of(om, w):
    return 2.0 * float(np.sum(om * w))


def lam_of(om, w):
    return 2.0 * float(np.sum(w / om))


def scale_to_m(oms, fr, lam, mtarget):
    """keep couplings fractions fr*lam, rescale all frequencies by s so that m = mtarget."""
    om, w = spec_from(fr * lam, oms)
    s = math.sqrt(mtarget / m_of(om, w))
    return spec_from(fr * lam, oms * s)


def check(om, w, tag, Ts=(300.0, 330.0, 450.0)):
    ev = max(top_eig(om, w, T, Nmult=40.0) for T in Ts)
    if ev > worst_global["ev"]:
        worst_global.update(ev=ev, what=tag)
    flag = "VIOLATION" if ev >= 1.0 else ""
    return ev, flag


def nm_tests(lam_s, frac=0.9999):
    lamk = float(lam_s)
    g = glo(lam_s)
    mt = frac * g
    rows = []
    # 1 Einstein at lambda_k and at smaller lambda
    for lf in (1.0, 0.95, 0.8, 0.5, 0.2):
        lam = lf * lamk
        om, w = spec_from([lam], [math.sqrt(mt / lam)])
        rows.append(("Einstein lam=%.3g" % lam,) + check(om, w, "%s Einstein lf=%g" % (lam_s, lf)))
    # 2 two-atom, wide ratios, lambda = lambda_k
    for r in (1.5, 3, 10, 30, 100, 1000):
        for f in (0.02, 0.1, 0.3, 0.5, 0.7, 0.9, 0.98):
            fr = np.array([f, 1 - f])
            om, w = scale_to_m(np.array([1.0, r]), fr, lamk, mt)
            rows.append(("2atom r=%g f_lo=%g" % (r, f),) + check(om, w, "%s 2atom r=%g f=%g" % (lam_s, r, f)))
    # 3 soft atom with big coupling share + hard atom, lambda < lambda_k
    for lf in (0.6, 0.9, 1.0):
        for r in (20, 200):
            fr = np.array([0.8, 0.2])
            om, w = scale_to_m(np.array([1.0, r]), fr, lf * lamk, mt)
            rows.append(("soft+hard lf=%g r=%g" % (lf, r),) + check(om, w, "%s softhard" % lam_s))
    # 4 box spectrum
    for width in (0.2, 1.0, 1.8):
        oms = np.linspace(1 - width / 2, 1 + width / 2, 41)
        fr = np.ones(41) / 41
        om, w = scale_to_m(oms, fr, lamk, mt)
        rows.append(("box width=%g" % width,) + check(om, w, "%s box" % lam_s))
    # 5 Nelder-Mead over 2- and 3-atom spectra on the constraint surface (lambda=lambda_k, m=mt)
    def obj(p, k):
        p = np.clip(p, -4.5, 4.5)
        lr = np.concatenate([[0.0], p[:k - 1]])
        fr = np.exp(np.concatenate([[0.0], p[k - 1:]]))
        fr = fr / fr.sum()
        om, w = scale_to_m(np.exp(lr), fr, lamk, mt)
        return -top_eig(om, w, TT, Nmult=30.0)
    rng = np.random.default_rng(7)
    for k in (2, 3):
        best = None
        for _ in range(3):
            p0 = np.concatenate([rng.uniform(-3, 3, k - 1), rng.uniform(-2, 2, k - 1)])
            r = minimize(obj, p0, args=(k,), method="Nelder-Mead",
                         options={"maxiter": 150, "xatol": 1e-4, "fatol": 1e-7})
            if best is None or r.fun < best.fun:
                best = r
        ev = -best.fun
        if ev > worst_global["ev"]:
            worst_global.update(ev=ev, what="%s NM-opt k=%d" % (lam_s, k))
        rows.append(("NM-opt %d atoms p=%s" % (k, np.round(best.x, 3)), ev, "VIOLATION" if ev >= 1 else ""))
    print("== lambda_k = %s  g_lo = %.6e K^2, m = %.5f g_lo" % (lam_s, g, frac), flush=True)
    for name, ev, fl in rows:
        print("   %-40s top_ev = %.6f %s" % (name, ev, fl), flush=True)
    return max(r[1] for r in rows)


def sharpness(lam_s):
    """Einstein at lambda_k with m = 1.001 g_lo should be UNSTABLE at 300 K (top_ev > 1)."""
    lamk = float(lam_s)
    mt = 1.001 * glo(lam_s)
    om, w = spec_from([lamk], [math.sqrt(mt / lamk)])
    return top_eig(om, w, TT, Nmult=40.0)


def n1mu_tests(mu, wc):
    d = json.load(open(os.path.join(CERT, "n1mu_mu%s_wc%d.json" % (mu, wc))))
    A = 0.9999 * float(d["A_floor_K"])
    j300 = d["N_c_at_300"]
    Ts = [300.0]
    for j in range(j300, 0, -1):
        b = wc / ((2 * j + 1) * math.pi)
        if b > 300:
            Ts += [b * (1 - 1e-9), b * (1 + 1e-9)]
    Ts += [350.0, 400.0, 450.0, 500.0]
    Ts = sorted(t for t in Ts if 300 <= t <= 520)
    worst = -1
    for Om in (500, 900, 1300, 1700, 2100, 2600, 3500, 5000):
        om, w = np.array([float(Om)]), np.array([A])
        for T in Ts:
            worst = max(worst, top_eig(om, w, T, float(mu), float(wc), "sym", Nmult=30.0))
    for (O1, O2, f) in ((300, 3000, 0.5), (800, 2500, 0.3), (1000, 8000, 0.5), (100, 1500, 0.1)):
        om, w = np.array([float(O1), float(O2)]), np.array([f * A, (1 - f) * A])
        for T in Ts:
            worst = max(worst, top_eig(om, w, T, float(mu), float(wc), "sym", Nmult=30.0))
    # sharpness: Einstein near optimum with A = 1.001 * A_float should be unstable at some T in grid
    sh = -1
    Af = 1.002 * float(d["A_float_K"])
    for Om in np.linspace(1200, 3200, 21):
        for T in Ts:
            sh = max(sh, top_eig(np.array([Om]), np.array([Af]), T, float(mu), float(wc), "sym", Nmult=30.0))
    print("== N1mu mu*=%s wc=%d A=0.9999 A_floor: max top_ev over shapes x T grid = %.6f %s;"
          " sharpness (A=1.002 A_float) max top_ev = %.6f" % (mu, wc, worst, "VIOLATION" if worst >= 1 else "",
                                                               sh), flush=True)
    return worst


def corollary():
    C = 0.18272634
    print("== corollary Tc <= %.8f sqrt(m) (mu*=0)" % C, flush=True)
    cases = [("Einstein lam=10, 1000K", [10.0], [1000.0]), ("Einstein lam=100, 300K", [100.0], [300.0]),
             ("Einstein lam=1000, 100K", [1000.0], [100.0]),
             ("2atom lam 50@50K + 2@2000K", [50.0, 2.0], [50.0, 2000.0]),
             ("2atom lam 500@20K + 0.5@1500K", [500.0, 0.5], [20.0, 1500.0])]
    for name, lams, oms in cases:
        om, w = spec_from(lams, oms)
        m = m_of(om, w)
        tc = tc_indep(om, w, T_hi=0.25 * math.sqrt(m), Nmult=40.0)
        print("   %-34s Tc=%.4f  C*sqrt(m)=%.4f  ratio=%.7f %s" % (name, tc, C * math.sqrt(m), tc / (C * math.sqrt(m)),
                                                              "VIOLATION" if tc > C * math.sqrt(m) else ""), flush=True)


if __name__ == "__main__":
    lams = sys.argv[1].split(",") if len(sys.argv) > 1 else ["0.8", "1.0", "1.5", "2.0", "3.0", "5.0", "10.0", "100.0"]
    for ls in lams:
        nm_tests(ls)
        print("   sharpness: Einstein lam_k, m = 1.001 g_lo: top_ev = %.6f (should be > 1)" % sharpness(ls), flush=True)
    for mu in ("0.10", "0.13"):
        for wc in (5000, 10000, 20000, 50000):
            n1mu_tests(mu, wc)
    corollary()
    print("WORST top_ev over all NM tests: %.6f at %s  -> %s" % (worst_global["ev"], worst_global["what"],
          "VIOLATION FOUND" if worst_global["ev"] >= 1 else "no violation"), flush=True)
