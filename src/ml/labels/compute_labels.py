"""Compute alpha^2F moment labels (A, lambda, omega_log, lambda<omega^2>) per structure.

Definitions (docs/theory/01_eliashberg_targets.md §0), all integrals over omega > 0 only:
  A          = int a2F(w) dw                       [K]   (Tc <= A/4.33148, theorem N1)
  lambda     = 2 int a2F(w)/w dw
  omega_log  = exp( (2/lambda) int a2F(w) ln(w)/w dw )   [K]
  lam_w2     = lambda<w^2> = 2 int w a2F(w) dw     [K^2]
Frequencies converted to Kelvin (1 meV = 11.604518 K, 1 Ry = 13605.693 meV).

Sources:
  cerq  : Materials Cloud 3kbt5-r3n56 (Cerqueira, Sanna, Marques), streamed extract (fetch_cerqueira.py).
          a2F.dosK, K=1..10 are QE double-delta smearings 0.005..0.050 Ry.  The authors' own
          McMillan.dat integrals equal smearing K=6 (0.030 Ry) -- verified on the first entries -- so
          K=6 is the primary label; K=4..8 give the smearing spread.
  jarvis: JARVIS-DFT EPC (figshare 38307921, jarvis_epc_data_figshare_1058.json): 50-bin a2F in meV
          ('a2F_original_x/y').  COMPUTED from a coarse grid.
  alexsc: JARVIS 'alex_supercon' (figshare 55551155): lambda, omega_log only -> only the rigorous
          lower bound A >= lambda*omega_log/2 (Jensen: <w> >= omega_log) is DERIVED.  Not used as a label.

Output: one JSON-lines file per source with structure (pymatgen dict) + labels + flags.
Run with ~/.rtt/envs/py/bin/python -I.
"""
import argparse
import gzip
import json
import os
import re
import sys
import warnings

import numpy as np

warnings.filterwarnings("ignore")
MEV_K = 11.604518
RY_MEV = 13605.693
LIGHT = {"H", "Li", "Be", "B", "C", "N", "O", "F"}


def sfloat(x):
    """float() that maps the archive's '-----' (no Tc / not converged) to None."""
    try:
        return float(x)
    except ValueError:
        return None


def moments(w, f):
    """w in K (array), f = a2F (dimensionless).  Positive frequencies only."""
    p = w > 0
    w, f = w[p], np.clip(f[p], 0, None)
    A = np.trapezoid(f, w)
    lam = 2 * np.trapezoid(f / w, w)
    if lam <= 0:
        return dict(A_K=float(A), lam=float(lam), wlog_K=float("nan"), lam_w2_K2=float("nan"))
    wlog = np.exp(2 / lam * np.trapezoid(f / w * np.log(w), w))
    lw2 = 2 * np.trapezoid(f * w, w)
    return dict(A_K=float(A), lam=float(lam), wlog_K=float(wlog), lam_w2_K2=float(lw2))


def comp_flags(struct):
    from pymatgen.core import Composition
    c = struct.composition
    n = c.num_atoms
    els = [e.symbol for e in c.elements]
    return dict(formula=c.reduced_formula, nsites=len(struct), elements=sorted(els),
                x_H=float(c["H"] / n) if "H" in els else 0.0,
                x_light=float(sum(c[e] for e in els if e in LIGHT) / n),
                has_H="H" in els)


def do_cerq(root, out):
    from pymatgen.core import Structure
    n_ok = n_bad = 0
    with open(out, "w") as fo:
        for batch in sorted(os.listdir(root)):
            bd = os.path.join(root, batch)
            if not os.path.isdir(bd):
                continue
            for mat in sorted(os.listdir(bd)):
                d = os.path.join(bd, mat)
                rec = {"source": "cerq", "id": mat.rsplit("_", 1)[-1], "dir": f"{batch}/{mat}"}
                try:
                    s = Structure.from_file(os.path.join(d, "geo_opt.cif"))
                    labs = {}
                    neg = {}
                    for k in range(1, 11):
                        fn = os.path.join(d, f"a2F.dos{k}.gz")
                        if not os.path.exists(fn):
                            continue
                        a = np.loadtxt(gzip.open(fn, "rt"), comments=["#", " lambda", "lambda"])
                        w = a[:, 0] * RY_MEV * MEV_K
                        f = a[:, 1]
                        labs[k] = moments(w, f)
                        q = w < 0
                        neg[k] = float(np.trapezoid(np.abs(f[q]), w[q])) if q.sum() > 1 else 0.0
                    if 6 not in labs:
                        raise ValueError("no smearing 6")
                    rec.update(labs[6])
                    rec["A_neg_K"] = neg[6]
                    sp = [labs[k]["A_K"] for k in range(4, 9) if k in labs]
                    rec["A_K_smear_min"], rec["A_K_smear_max"] = float(min(sp)), float(max(sp))
                    rec["lam_smear"] = [round(labs[k]["lam"], 4) for k in sorted(labs)]
                    mc = open(os.path.join(d, "McMillan.dat")).read() if os.path.exists(os.path.join(d, "McMillan.dat")) else ""
                    m = re.search(r"integral a2F\s*=\s*(\S+)", mc)
                    rec["A_K_authors"] = float(m.group(1)) * MEV_K if m else None
                    m = re.search(r"#\s*lambda\s*=\s*(\S+)", mc)
                    rec["lam_authors"] = float(m.group(1)) if m else None
                    m = re.search(r"^\s*0\.100\s+(\S+)\s+(\S+)", mc, re.M)
                    rec["Tc_AD_mu010"] = sfloat(m.group(2)) if m else None
                    el = os.path.join(d, "Eliashberg.dat")
                    if os.path.exists(el):
                        m = re.search(r"^\s*0\.100\s+(\S+)", open(el).read(), re.M)
                        rec["Tc_El_mu010"] = sfloat(m.group(1)) if m else None
                    # phonon dispersion min frequency (cm^-1)
                    fg = os.path.join(d, "qe.freq.gp.gz")
                    if os.path.exists(fg):
                        fr = np.loadtxt(gzip.open(fg, "rt"))
                        rec["wmin_cm"] = float(fr[:, 1:].min())
                        rec["wmax_cm"] = float(fr[:, 1:].max())
                    inp = open(os.path.join(d, "geoopt_0.in")).read() if os.path.exists(os.path.join(d, "geoopt_0.in")) else ""
                    m = re.search(r"\bpress\s*=\s*([-\d.eEdD+]+)", inp)
                    rec["pressure_GPa"] = float(m.group(1).replace("d", "e")) / 10 if m else 0.0  # QE press in kbar
                    rec.update(comp_flags(s))
                    rec["structure"] = s.as_dict()
                    n_ok += 1
                except Exception as e:  # noqa
                    rec["error"] = repr(e)[:200]
                    n_bad += 1
                fo.write(json.dumps(rec) + "\n")
    print("cerq ok", n_ok, "bad", n_bad)


def jatoms_to_struct(a):
    from pymatgen.core import Lattice, Structure
    return Structure(Lattice(a["lattice_mat"]), a["elements"], a["coords"],
                     coords_are_cartesian=a.get("cartesian", False))


def do_jarvis(fn, out):
    d = json.load(open(fn))
    with open(out, "w") as fo:
        for e in d:
            rec = {"source": "jarvis", "id": e["jid"], "stability": e["stability"],
                   "pressure_GPa": float(e["press"].split()[0]), "lam_authors": e["lamb"],
                   "wlog_authors_K": e["wlog"], "Tc_authors": e["Tc"]}
            try:
                s = jatoms_to_struct(e["atoms"])
                x = np.array(e["a2F_original_x"]) * MEV_K
                y = np.array(e["a2F_original_y"])
                rec.update(moments(x, y))
                q = x < 0
                rec["A_neg_K"] = float(np.trapezoid(np.abs(y[q]), x[q])) if q.sum() > 1 else 0.0
                rec["wmin_K_a2f"] = float(x.min())
                rec.update(comp_flags(s))
                rec["structure"] = s.as_dict()
            except Exception as ex:  # noqa
                rec["error"] = repr(ex)[:200]
            fo.write(json.dumps(rec) + "\n")
    print("jarvis", len(d))


def do_alexsc(fn, out):
    d = json.load(open(fn))
    with open(out, "w") as fo:
        for e in d:
            rec = {"source": "alexsc", "id": e["id"], "lam_authors": e["la"], "wlog_authors_K": e["wlog"],
                   "Tc_authors": e["Tc"], "dosef": e.get("dosef"), "debye": e.get("debye")}
            rec["A_K_lowerbound"] = 0.5 * e["la"] * e["wlog"]
            try:
                s = jatoms_to_struct(e["atoms"])
                rec.update(comp_flags(s))
                rec["structure"] = s.as_dict()
            except Exception as ex:  # noqa
                rec["error"] = repr(ex)[:200]
            fo.write(json.dumps(rec) + "\n")
    print("alexsc", len(d))


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("which", choices=["cerq", "jarvis", "alexsc"])
    ap.add_argument("inp")
    ap.add_argument("out")
    a = ap.parse_args()
    {"cerq": do_cerq, "jarvis": do_jarvis, "alexsc": do_alexsc}[a.which](a.inp, a.out)
