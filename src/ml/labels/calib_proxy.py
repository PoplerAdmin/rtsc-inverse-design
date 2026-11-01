"""Calibrate cheap MatterSim-phonon proxies of A on the held-out test split (DFT A known).

Inputs: test.extxyz (info: material_id, log10A), relax_test/relaxed.json, ph_test/phonons.json
(from mattersim_screen.py).  Reports Spearman correlations of log10 A with
  - MatterSim max phonon frequency, MatterSim phonon-DOS omega_log, x_light,
  - and a 5-fold CV linear model log10A ~ log10(omega_max) + log10(omega_log_phDOS) + x_light + x_H,
restricted to MatterSim-stable structures and for all.  The fitted model is saved (coefficients) so
generated structures can be scored with the same proxy.

Usage: python -I calib_proxy.py CALIBDIR OUT.json
"""
import json
import sys

import numpy as np
from ase.io import read
from pymatgen.core import Composition
from scipy.stats import spearmanr

LIGHT = {"H", "Li", "Be", "B", "C", "N", "O", "F"}
CM1_TO_K = 1.4387769


def feats(atoms, ph):
    c = Composition(atoms.get_chemical_formula())
    n = c.num_atoms
    els = [e.symbol for e in c.elements]
    xl = sum(c[e] for e in els if e in LIGHT) / n
    xh = c["H"] / n if "H" in els else 0.0
    wmax = max(ph["max_freq_cm1"], 1.0) * CM1_TO_K
    wlog = ph["omega_log_phDOS_K"] or np.nan
    return [np.log10(wmax), np.log10(wlog), xl, xh]


def fit(X, y):
    X1 = np.c_[np.ones(len(X)), X]
    coef, *_ = np.linalg.lstsq(X1, y, rcond=None)
    return coef


def pred(coef, X):
    return np.c_[np.ones(len(X)), X] @ coef


def main():
    d, out = sys.argv[1], sys.argv[2]
    ats = read(f"{d}/test.extxyz", ":")
    rel = read(f"{d}/relax_test/relaxed.extxyz", ":")
    ph = {r["index"]: r for r in json.load(open(f"{d}/ph_test/phonons.json"))["results"]}
    X, y, stab = [], [], []
    for i, a in enumerate(ats):
        if i not in ph or ph[i]["omega_log_phDOS_K"] is None:
            continue
        X.append(feats(rel[i], ph[i]))
        y.append(a.info["log10A"])
        stab.append(ph[i]["dynamically_stable_ml"])
    X, y, stab = np.array(X), np.array(y), np.array(stab)
    res = {"n": int(len(y)), "n_mattersim_stable": int(stab.sum())}
    names = ["log10_wmax_K", "log10_wlog_phDOS_K", "x_light", "x_H"]
    for subset, m in (("all", np.ones(len(y), bool)), ("stable", stab)):
        r = {"n": int(m.sum())}
        for j, nm in enumerate(names):
            r["spearman_" + nm] = float(spearmanr(X[m, j], y[m])[0])
        idx = np.where(m)[0]
        rng = np.random.default_rng(0)
        rng.shuffle(idx)
        folds = np.array_split(idx, 5)
        p = np.zeros(len(y))
        for k in range(5):
            te = folds[k]
            tr = np.concatenate([folds[j] for j in range(5) if j != k])
            p[te] = pred(fit(X[tr], y[tr]), X[te])
        r["cv_linear_spearman"] = float(spearmanr(p[m], y[m])[0])
        r["cv_linear_mae_log10"] = float(np.abs(p[m] - y[m]).mean())
        r["baseline_mae_log10"] = float(np.abs(y[m] - y[m].mean()).mean())
        r["coef_full_fit"] = dict(zip(["intercept"] + names, map(float, fit(X[m], y[m]))))
        res[subset] = r
    json.dump(res, open(out, "w"), indent=1)
    print(json.dumps(res, indent=1))


if __name__ == "__main__":
    main()
