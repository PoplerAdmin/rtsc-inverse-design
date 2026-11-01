"""Merge label JSONL files, apply training filters, write summary + MatterGen CSV splits.

Filters for TRAINING rows (all recorded in the summary):
  - pressure <= 1 GPa (all current sources are 0 GPa)
  - nsites <= 20 (MatterGen training domain)
  - finite A > 0
  - dynamical-stability gate: cerq: phonon-dispersion min freq >= WMIN_CM (default -30 cm^-1, i.e.
    tolerating small Gamma/ASR noise); jarvis: 'stability' == 'stable'
  - de-duplication: cerq preferred over jarvis for the same (reduced formula, nsites) is NOT done
    (different functionals / structures); exact duplicates of agm id within cerq are dropped.
Split: grouped by reduced formula (no formula in two splits), 90/5/5, seed 0.
Label column: rtt_log10A = log10(A / K).

Usage: python -I build_dataset.py LABELDIR OUTDIR [--wmin-cm -30]
"""
import argparse
import collections
import csv
import json
import os
import random

import numpy as np


def load(fn):
    return [json.loads(l) for l in open(fn)] if os.path.exists(fn) else []


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("labeldir")
    ap.add_argument("outdir")
    ap.add_argument("--wmin-cm", type=float, default=-30.0)
    ap.add_argument("--no-jarvis", action="store_true")
    a = ap.parse_args()
    os.makedirs(a.outdir, exist_ok=True)
    cerq = load(os.path.join(a.labeldir, "cerq.jsonl"))
    jar = load(os.path.join(a.labeldir, "jarvis.jsonl"))
    summ = collections.OrderedDict()
    rows, seen = [], set()
    cnt = collections.Counter()
    for r in cerq:
        cnt["cerq_total"] += 1
        if "error" in r:
            cnt["cerq_error"] += 1
            continue
        if r["id"] in seen:
            cnt["cerq_dup"] += 1
            continue
        seen.add(r["id"])
        rows.append(r)
    for r in jar:
        cnt["jarvis_total"] += 1
        if "error" in r:
            cnt["jarvis_error"] += 1
            continue
        rows.append(r)

    def stable(r):
        if r["source"] == "cerq":
            return r.get("wmin_cm", -1e9) >= a.wmin_cm
        return r["stability"] == "stable"

    keep = []
    for r in rows:
        src = r["source"]
        if not (np.isfinite(r.get("A_K", np.nan)) and r["A_K"] > 0 and np.isfinite(r.get("lam", np.nan))):
            cnt[f"{src}_drop_nanA"] += 1
            continue
        if r.get("pressure_GPa", 0.0) > 1.0:
            cnt[f"{src}_drop_pressure"] += 1
            continue
        if r["nsites"] > 20:
            cnt[f"{src}_drop_nsites>20"] += 1
            continue
        if not stable(r):
            cnt[f"{src}_drop_unstable"] += 1
            continue
        if a.no_jarvis and src == "jarvis":
            continue
        keep.append(r)
    summ["counts"] = dict(cnt)

    def dist(rs, tag):
        A = np.array([r["A_K"] for r in rs])
        lam = np.array([r["lam"] for r in rs])
        lw2 = np.array([r["lam_w2_K2"] for r in rs])
        H = np.array([r["has_H"] for r in rs])
        d = {"n": len(rs), "n_H": int(H.sum()),
             "A_K_pct_5_25_50_75_95_99": [round(float(x), 1) for x in np.percentile(A, [5, 25, 50, 75, 95, 99])] if len(A) else [],
             "A_K_max": round(float(A.max()), 1) if len(A) else None,
             "n_A_ge": {t: int((A >= t).sum()) for t in (300, 500, 800, 1000, 1300, 2060)},
             "n_A_ge_H": {t: int(((A >= t) & H).sum()) for t in (300, 500, 800, 1000, 1300)},
             "lam_pct_50_95_99": [round(float(x), 3) for x in np.percentile(lam, [50, 95, 99])] if len(lam) else [],
             "n_lamw2_ge_6e6": int((lw2 >= 6e6).sum()),
             "lamw2_max": float(np.nanmax(lw2)) if len(lw2) else None}
        summ[tag] = d

    allrows = [r for r in rows if np.isfinite(r.get("A_K", np.nan)) and r["A_K"] > 0]
    for src in ("cerq", "jarvis"):
        dist([r for r in allrows if r["source"] == src], f"all_{src}")
        dist([r for r in keep if r["source"] == src], f"train_{src}")
    dist(keep, "train_all")
    top = sorted(allrows, key=lambda r: -r["A_K"])[:40]
    summ["top40_A_all"] = [{k: (round(r[k], 3) if isinstance(r.get(k), float) else r.get(k)) for k in
                            ("source", "id", "formula", "nsites", "A_K", "lam", "wlog_K", "lam_w2_K2",
                             "wmin_cm", "stability", "x_H", "x_light", "Tc_AD_mu010", "Tc_El_mu010")}
                           for r in top]

    # grouped split
    forms = sorted({r["formula"] for r in keep})
    random.Random(0).shuffle(forms)
    n = len(forms)
    fval = set(forms[: n // 20])
    ftest = set(forms[n // 20: n // 10])
    split = {"train": [], "val": [], "test": []}
    for r in keep:
        s = "val" if r["formula"] in fval else "test" if r["formula"] in ftest else "train"
        split[s].append(r)
    summ["split_sizes"] = {k: len(v) for k, v in split.items()}
    from pymatgen.core import Structure
    dens = []
    man = open(os.path.join(a.outdir, "manifest.tsv"), "w")
    man.write("split\tsource\tid\tformula\tnsites\tA_K\tlam\twlog_K\tlam_w2_K2\tx_H\tx_light\n")
    os.makedirs(os.path.join(a.outdir, "csv"), exist_ok=True)
    for s, rs in split.items():
        with open(os.path.join(a.outdir, "csv", f"{s}.csv"), "w", newline="") as fo:
            w = csv.writer(fo)
            w.writerow(["material_id", "cif", "rtt_log10A"])
            for r in rs:
                st = Structure.from_dict(r["structure"])
                dens.append(len(st) / st.volume)
                mid = f"{r['source']}:{r['id']}"
                w.writerow([mid, st.to(fmt="cif"), f"{np.log10(r['A_K']):.6f}"])
                man.write(f"{s}\t{r['source']}\t{r['id']}\t{r['formula']}\t{r['nsites']}\t{r['A_K']:.2f}\t"
                          f"{r['lam']:.4f}\t{r['wlog_K']:.2f}\t{r['lam_w2_K2']:.1f}\t{r['x_H']:.3f}\t{r['x_light']:.3f}\n")
    summ["average_density_atoms_per_A3"] = float(np.mean(dens))
    y = np.log10([r["A_K"] for r in split["train"]])
    summ["train_log10A_mean_std"] = [float(y.mean()), float(y.std())]
    json.dump(summ, open(os.path.join(a.outdir, "summary.json"), "w"), indent=1)
    print(json.dumps({k: v for k, v in summ.items() if k != "top40_A_all"}, indent=1))


if __name__ == "__main__":
    main()
