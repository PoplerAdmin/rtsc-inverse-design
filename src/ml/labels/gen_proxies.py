"""Composition proxies for MatterGen samples conditioned on rtt_log10A, plus a kNN label lookup.

For each generation directory (generated_crystals.extxyz) report per-structure:
  x_H, x_light (H,Li,Be,B,C,N,O,F atom fraction), mean atomic mass per atom, n_elements, nsites,
  exact_formula_in_train (reduced formula present in the training manifest),
  knn_log10A: median log10(A) of the k=5 training entries nearest in element-fraction space (L1).
The kNN lookup is validated on the held-out test split (--validate) to show how much it knows.
Mann-Whitney U tests (two-sided) compare each conditioned set with the reference set.

Usage: python -I gen_proxies.py MANIFEST.tsv OUT.json REF_DIR DIR1 [DIR2 ...]
"""
import json
import sys

import numpy as np
from ase.io import read
from pymatgen.core import Composition, Element
from scipy.stats import mannwhitneyu

LIGHT = {"H", "Li", "Be", "B", "C", "N", "O", "F"}
K = 5


def vec(comp):
    v = np.zeros(100)
    n = comp.num_atoms
    for el, amt in comp.items():
        v[el.Z] = amt / n
    return v


def load_manifest(fn):
    rows = [l.rstrip("\n").split("\t") for l in open(fn)]
    hdr, rows = rows[0], rows[1:]
    out = []
    for r in rows:
        d = dict(zip(hdr, r))
        d["A_K"] = float(d["A_K"])
        d["comp"] = Composition(d["formula"])
        d["vec"] = vec(d["comp"])
        out.append(d)
    return out


def knn(train_X, train_y, v):
    d = np.abs(train_X - v).sum(1)
    idx = np.argsort(d)[:K]
    return float(np.median(train_y[idx])), float(d[idx].mean())


def props(comp):
    n = comp.num_atoms
    els = [e.symbol for e in comp.elements]
    return dict(x_H=comp["H"] / n if "H" in els else 0.0,
                x_light=sum(comp[e] for e in els if e in LIGHT) / n,
                mass=sum(Element(e).atomic_mass * comp[e] for e in els) / n,
                n_el=len(els))


def main():
    man, outfn, ref, dirs = sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4:]
    rows = load_manifest(man)
    tr = [r for r in rows if r["split"] == "train"]
    te = [r for r in rows if r["split"] == "test"]
    X = np.array([r["vec"] for r in tr])
    y = np.log10([r["A_K"] for r in tr])
    trforms = {r["formula"] for r in tr}
    res = {"k": K}
    # validation of the kNN proxy on the held-out test split
    pred = np.array([knn(X, y, r["vec"])[0] for r in te])
    yt = np.log10([r["A_K"] for r in te])
    res["knn_validation_test"] = dict(n=len(te), mae_log10=float(np.abs(pred - yt).mean()),
                                      mae_log10_baseline_mean=float(np.abs(y.mean() - yt).mean()),
                                      spearman=float(__import__("scipy.stats").stats.spearmanr(pred, yt)[0]))
    sets = {}
    for d in [ref] + dirs:
        at = read(f"{d}/generated_crystals.extxyz", ":")
        per = []
        for a in at:
            c = Composition(a.get_chemical_formula())
            p = props(c)
            p["nsites"] = len(a)
            p["formula"] = c.reduced_formula
            p["in_train"] = c.reduced_formula in trforms
            p["knn_log10A"], p["knn_dist"] = knn(X, y, vec(c))
            per.append(p)
        sets[d] = per
    refp = sets[ref]
    summ = {}
    for d, per in sets.items():
        s = {"n": len(per)}
        for key in ("x_H", "x_light", "mass", "n_el", "nsites", "knn_log10A", "knn_dist"):
            v = np.array([p[key] for p in per], float)
            s[key + "_mean"] = float(v.mean())
            s[key + "_median"] = float(np.median(v))
            if d != ref:
                r0 = np.array([p[key] for p in refp], float)
                s[key + "_MWU_p"] = float(mannwhitneyu(v, r0, alternative="two-sided").pvalue)
        s["frac_with_H"] = float(np.mean([p["x_H"] > 0 for p in per]))
        s["frac_formula_in_train"] = float(np.mean([p["in_train"] for p in per]))
        s["knn_A_K_median"] = float(10 ** s["knn_log10A_median"])
        s["top_formulas"] = sorted({p["formula"] for p in per}, key=lambda f: -sum(q["formula"] == f for q in per))[:12]
        summ[d] = s
    res["sets"] = summ
    res["per_structure"] = sets
    json.dump(res, open(outfn, "w"), indent=1)
    print(json.dumps({"knn_validation_test": res["knn_validation_test"],
                      "sets": {d.rsplit("/", 1)[-1]: {k: (round(v, 4) if isinstance(v, float) else v) for k, v in s.items()}
                               for d, s in summ.items()}}, indent=1))


if __name__ == "__main__":
    main()
