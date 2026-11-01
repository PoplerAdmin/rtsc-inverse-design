"""Summarise MatterSim phonon proxies for generated sets and compare each with a reference set.

For each set directory containing relax/relaxed.extxyz and ph/phonons.json (mattersim_screen.py):
  stable fraction (MatterSim, imag tol 15 cm^-1), median omega_max (K), median omega_log of the phonon DOS (K),
  and the calibrated linear proxy log10A_hat (coefficients from calib_proxy.py, 'stable' fit).
Mann-Whitney U (two-sided) against the reference set, on MatterSim-stable structures only.
MatterSim phonons are NOT reliable for hydrides (docs/design/10 §3.4); the stable filter is a ranking aid.

Usage: python -I gen_phonon_proxy.py CALIB.json OUT.json REF_SET SET1 [SET2 ...]
"""
import json
import sys

import numpy as np
from ase.io import read
from pymatgen.core import Composition
from scipy.stats import mannwhitneyu

LIGHT = {"H", "Li", "Be", "B", "C", "N", "O", "F"}
CM1_TO_K = 1.4387769


def load(d, coef):
    rel = read(f"{d}/relax/relaxed.extxyz", ":")
    rows = []
    for r in json.load(open(f"{d}/ph/phonons.json"))["results"]:
        a = rel[r["index"]]
        c = Composition(a.get_chemical_formula())
        n = c.num_atoms
        els = [e.symbol for e in c.elements]
        xl = sum(c[e] for e in els if e in LIGHT) / n
        xh = c["H"] / n if "H" in els else 0.0
        wmax = max(r["max_freq_cm1"], 1.0) * CM1_TO_K
        wlog = r["omega_log_phDOS_K"]
        p = None
        if wlog:
            p = (coef["intercept"] + coef["log10_wmax_K"] * np.log10(wmax) + coef["log10_wlog_phDOS_K"] * np.log10(wlog)
                 + coef["x_light"] * xl + coef["x_H"] * xh)
        rows.append(dict(stable=r["dynamically_stable_ml"], wmax_K=wmax, wlog_K=wlog, proxy_log10A=p,
                         formula=c.reduced_formula))
    return rows


def main():
    calib, out, ref, sets = sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4:]
    coef = json.load(open(calib))["stable"]["coef_full_fit"]
    data = {d: load(d, coef) for d in [ref] + sets}
    res = {}
    for d, rows in data.items():
        st = [r for r in rows if r["stable"] and r["wlog_K"]]
        s = {"n": len(rows), "n_stable": len(st)}
        for k in ("wmax_K", "wlog_K", "proxy_log10A"):
            v = np.array([r[k] for r in st], float)
            s[k + "_median_stable"] = float(np.median(v)) if len(v) else None
            if d != ref and len(v):
                r0 = np.array([r[k] for r in data[ref] if r["stable"] and r["wlog_K"]], float)
                s[k + "_MWU_p"] = float(mannwhitneyu(v, r0, alternative="two-sided").pvalue)
        s["top_proxy"] = sorted(((round(r["proxy_log10A"], 3), r["formula"]) for r in st), reverse=True)[:8]
        res[d] = s
    json.dump(res, open(out, "w"), indent=1)
    print(json.dumps({d.rsplit("/", 1)[-1]: v for d, v in res.items()}, indent=1))


if __name__ == "__main__":
    main()
