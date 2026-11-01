"""Feasible region for Tc = 300 K in isotropic Migdal-Eliashberg theory.

Outputs (docs/theory/data/):
  t_vs_lambda.csv      per shape, mu*, cutoff multiple c (omega_c = c * omega_max), lambda:
                       Tc/omega_log and the moments required for Tc = 300 K
  min_A.csv            per shape, mu*, c: lambda minimising A/Tc and the required moments there
  einstein_abs.csv     Einstein extremal problem with an ABSOLUTE cutoff omega_c (for the
                       extreme-point bound): min over Omega, lambda of A/Tc and of lam<w2>/Tc^2
  random_spectra.csv   mu* = 0 random atomic spectra: Tc/A and Tc/sqrt(lam<w2>)
  br_derivative.csv    Bergmann-Rainer delta Tc/delta a2F(omega) at Tc = 300 K
Run (WSL): OMP_NUM_THREADS=2 nice -n 19 python scan_targets.py <outdir>
"""
import csv
import math
import os
import sys
import time

import numpy as np
from scipy.optimize import minimize_scalar, brentq

sys.path.insert(0, os.path.dirname(__file__))
from eliashberg import (Spectrum, einstein, gaussians, tc, rho, functional_derivative)  # noqa

OUT = sys.argv[1] if len(sys.argv) > 1 else "."
TT = 300.0
MUS = [0.0, 0.08, 0.10, 0.13, 0.16, 0.20]
CMULTS = [5, 10]
LAMS = [0.75, 1.0, 1.25, 1.5, 1.75, 2.0, 2.5, 3.0, 3.5, 4.0, 5.0, 6.0, 8.0]


def two_delta(r, f_hi):
    """lambda share f_hi at Omega_hi = r * Omega_lo; unit lambda."""
    lo, hi = 1.0, float(r)
    return Spectrum(np.array([lo, hi]), np.array([(1 - f_hi) * lo / 2, f_hi * hi / 2]))


def box(lo_frac=0.1):
    om = np.linspace(lo_frac, 1.0, 400)
    from eliashberg import from_density
    return from_density(om, np.ones_like(om), 1.0)


def debye():
    om = np.linspace(0.005, 1.0, 400)
    from eliashberg import from_density
    return from_density(om, om ** 2, 1.0)


# meV -> K: 11.6045
mK = 11.6045
SHAPES = {
    "einstein": einstein(1.0, 1.0),
    "2delta_r3_f0.5": two_delta(3, 0.5),
    "2delta_r8_f0.5": two_delta(8, 0.5),
    "2delta_r8_f0.8": two_delta(8, 0.8),
    "box_0.1-1": box(0.1),
    "debye_w2": debye(),
    # schematic, not fitted to any calculation (see 01_eliashberg_targets.md, Sec. 5)
    "H3S_like": gaussians([20 * mK, 105 * mK, 175 * mK], [6 * mK, 20 * mK, 25 * mK],
                          [0.12, 0.55, 0.33], 1.0),
    "LaH10_like": gaussians([12 * mK, 60 * mK, 140 * mK], [4 * mK, 20 * mK, 35 * mK],
                            [0.15, 0.45, 0.40], 1.0),
    "MgB2_like": gaussians([35 * mK, 70 * mK], [8 * mK, 5 * mK], [0.35, 0.65], 1.0),
}


def unit_shape(sh: Spectrum, lam: float) -> Spectrum:
    """Shape scaled to the given lambda and to omega_log = 1000 K."""
    s = sh.scaled(lam_factor=lam / sh.lam)
    return s.scaled(freq_factor=1000.0 / s.omega_log)


def tc_shape(sh, lam, mu, cm):
    sp = unit_shape(sh, lam)
    return sp, tc(sp, mu, cm * sp.omega_max, T_lo=0.02 * sp.omega_log)


def row_for(name, sh, mu, cm, lam):
    sp, T = tc_shape(sh, lam, mu, cm)
    if T <= 0:
        return None
    s = TT / T  # stretch factor making Tc = 300 K
    return dict(shape=name, mu=mu, cmult=cm, lam=round(lam, 6), Tc_over_wlog=T / sp.omega_log,
                A_over_Tc=sp.A / T, lamw2_over_Tc2=sp.lam_omega2 / T ** 2,
                wlog_req=sp.omega_log * s, w2_req=sp.omega2 * s, w1_req=sp.omega1 * s,
                wmax_req=sp.omega_max * s, A_req=sp.A * s, lamw2_req=sp.lam_omega2 * s * s,
                sqrt_lamw2_req=math.sqrt(sp.lam_omega2) * s)


def write(path, rows):
    if not rows:
        return
    with open(path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        w.writeheader()
        for r in rows:
            w.writerow({k: (f"{v:.6g}" if isinstance(v, float) else v) for k, v in r.items()})


def main():
    t0 = time.time()
    # ---------------- 1. lambda scans and minimal A per shape
    rows, mins = [], []
    for name, sh in SHAPES.items():
        for mu in MUS:
            for cm in (CMULTS if mu > 0 else [10]):
                for lam in LAMS:
                    r = row_for(name, sh, mu, cm, lam)
                    if r:
                        rows.append(r)

                def aot(lam):
                    sp, T = tc_shape(sh, lam, mu, cm)
                    return sp.A / T if T > 0 else 1e9
                res = minimize_scalar(aot, bounds=(0.5, 25.0), method="bounded",
                                      options=dict(xatol=1e-3))
                r = row_for(name, sh, mu, cm, res.x)
                r["shape"] = name
                mins.append(r)
                print(f"[{time.time()-t0:7.1f}s] {name} mu*={mu} c={cm}: lam_opt={res.x:.3f} "
                      f"A_req={r['A_req']:.1f}K wlog_req={r['wlog_req']:.1f}K", flush=True)
    write(os.path.join(OUT, "t_vs_lambda.csv"), rows)
    write(os.path.join(OUT, "min_A.csv"), mins)

    # ---------------- 2. Einstein extremal problem, absolute cutoff
    # For absolute omega_c, A/Tc for Einstein (lam, Omega) at Tc=300 depends on (lam, Omega).
    # Parametrise by Omega/Tc = w and lam; find lam_c(w) with Tc = 300 K at Omega = 300 w.
    ein = []
    for mu in MUS:
        for wc in ([math.inf] if mu == 0 else [5000.0, 10000.0, 20000.0, 50000.0]):
            def lam_c(w):
                Om = TT * w
                # instability at T = 300 K iff rho >= 0; rho is convex in lam with rho(0) = -1
                f = lambda lam: rho(einstein(lam, Om), TT, mu, wc)
                hi = 2.0
                while f(hi) < 0:
                    hi *= 1.6
                    if hi > 5000:
                        return math.inf
                return brentq(f, 0.05, hi, xtol=1e-7)

            def a_of_w(w):
                lc = lam_c(w)
                return lc * w / 2 if math.isfinite(lc) else 1e9
            rA = minimize_scalar(a_of_w, bounds=(2.0, 20.0), method="bounded",
                                 options=dict(xatol=1e-4))
            lA = lam_c(rA.x)

            def m2_of_w(w):  # lam <w2> / Tc^2 = lam w^2
                lc = lam_c(w)
                return lc * w * w if math.isfinite(lc) else 1e12
            # lam<w2> is minimised as Omega -> 0 (lam -> inf); sample a decreasing w grid
            m2 = [(w, m2_of_w(w)) for w in [8, 4, 2, 1, 0.5, 0.25]]
            ein.append(dict(mu=mu, omega_c=wc, w_opt=rA.x, lam_opt=lA, A_over_Tc=rA.fun,
                            A_req=TT * rA.fun, Omega_opt=TT * rA.x,
                            **{f"lamw2_over_Tc2_at_w{w}": v for w, v in m2}))
            print(f"[{time.time()-t0:7.1f}s] Einstein abs mu*={mu} wc={wc}: A/Tc={rA.fun:.5f} "
                  f"w={rA.x:.3f} lam={lA:.3f}  lamw2/Tc2: {[round(v,2) for _,v in m2]}", flush=True)
    write(os.path.join(OUT, "einstein_abs.csv"), ein)

    # ---------------- 3. random spectra, mu*=0
    rng = np.random.default_rng(12345)
    rnd = []
    for i in range(300):
        k = rng.integers(1, 6)
        om = np.exp(rng.uniform(np.log(50), np.log(3000), k))
        w = rng.dirichlet(np.ones(k)) * rng.uniform(200, 3000)
        sp = Spectrum(om, w)
        T = tc(sp, 0.0, T_lo=0.01 * sp.omega_log)
        if T <= 0:
            continue
        rnd.append(dict(i=i, natoms=int(k), A=sp.A, lam=sp.lam, Tc=T, Tc_over_A=T / sp.A,
                        Tc_over_sqrt_lamw2=T / math.sqrt(sp.lam_omega2),
                        Tc_over_wlog=T / sp.omega_log))
    write(os.path.join(OUT, "random_spectra.csv"), rnd)
    print(f"[{time.time()-t0:7.1f}s] random: max Tc/A={max(r['Tc_over_A'] for r in rnd):.5f} "
          f"max Tc/sqrt(lamw2)={max(r['Tc_over_sqrt_lamw2'] for r in rnd):.5f}", flush=True)

    # ---------------- 4. Bergmann-Rainer derivative at Tc = 300 K
    br = []
    for name, lam, mu in [("einstein", 1.14365, 0.0), ("H3S_like", 2.0, 0.13),
                          ("LaH10_like", 2.5, 0.13), ("box_0.1-1", 2.0, 0.13)]:
        sh = SHAPES[name]
        sp, T = tc_shape(sh, lam, mu, 10)
        sp = sp.scaled(freq_factor=TT / T)
        T = tc(sp, mu, 10 * sp.omega_max)
        oms = np.linspace(0.05, 30, 300) * T
        d = functional_derivative(sp, T, oms, mu, 10 * sp.omega_max)
        ipk = int(np.argmax(d))
        print(f"[{time.time()-t0:7.1f}s] BR {name} lam={lam} mu*={mu}: Tc={T:.2f} "
              f"peak omega/Tc={oms[ipk]/T:.3f}  dTc/da2F(peak)={d[ipk]:.4f}", flush=True)
        for o, v in zip(oms, d):
            br.append(dict(shape=name, lam=lam, mu=mu, Tc=T, omega_over_Tc=o / T, dTc_da2F=v))
    write(os.path.join(OUT, "br_derivative.csv"), br)
    print(f"done in {time.time()-t0:.1f}s", flush=True)


if __name__ == "__main__":
    main()
