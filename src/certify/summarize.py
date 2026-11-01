"""Collect the certificate JSONs into docs/theory/certificates/SUMMARY.json and print markdown tables.
usage: python src/certify/summarize.py [docs/theory/certificates]"""
import glob
import json
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
root = sys.argv[1] if len(sys.argv) > 1 else os.path.join(HERE, "..", "..", "docs", "theory", "certificates")
nm, mu = [], []
for p in glob.glob(os.path.join(root, "*.json")):
    d = json.load(open(p))
    if d.get("rule") == "NM-G":
        nm.append(d)
    elif d.get("rule") == "N1mu":
        mu.append(d)
key = lambda d: float("inf") if d["lambda"] == "inf" else float(d["lambda"])
nm.sort(key=key)
mu.sort(key=lambda d: (float(d["mu_star"]), d["omega_c_K"]))
rows = []
print("| lambda_k | g_lo (K^2, certified) | g_hi (K^2, certified) | Omega_lo (K) | Omega_hi (K) | "
      "float g (01 table) | beta | delta | cells | verdict |")
print("|---|---|---|---|---|---|---|---|---|---|")
for d in nm:
    u = d["upper"]
    r = dict(lam=d["lambda"], m0_over_T2=d["m0_over_T2"], g_lo_K2=d["g_lo_K2"], g_hi_K2=u["g_hi_K2"],
             Omega_lo_K=d["Omega_lo_K"], Omega_hi_K=u.get("Omega_hi_K"), float_g_K2=d["float_g_K2"],
             beta=d["beta"], delta=d["delta"], n_cells=d["n_cells"], verdict=d["verdict"],
             upper_certified=u["certified_unstable"])
    rows.append(r)
    f = lambda x, p: "-" if x is None else ("%.*f" % (p, x))
    print("| %s | %.6e | %.6e | %s | %s | %.6e | %.6g | %g | %d | %s%s |" % (
        r["lam"], r["g_lo_K2"], r["g_hi_K2"], f(r["Omega_lo_K"], 3), f(r["Omega_hi_K"], 3),
        r["float_g_K2"], float(r["beta"]), r["delta"], r["n_cells"], r["verdict"],
        "" if r["upper_certified"] else " (upper NOT certified)"))
print()
print("| mu* | omega_c (K) | A_floor (K, certified) | float floor (K) | N_c(300 K) | segments certified | verdict |")
print("|---|---|---|---|---|---|---|")
mrows = []
for d in mu:
    segs = [s["j"] for s in d["segments"] if "a_cert" in s]
    mrows.append(dict(mu=d["mu_star"], omega_c=d["omega_c_K"], A_floor_K=d["A_floor_K"],
                      A_float_K=d["A_float_K"], Nc300=d["N_c_at_300"], segments=segs,
                      verdict=d["verdict"]))
    print("| %s | %g | %s | %.2f | %d | %s | %s |" % (d["mu_star"], d["omega_c_K"], d["A_floor_K"],
                                                    d["A_float_K"], d["N_c_at_300"],
                                                    ",".join(map(str, segs)), d["verdict"]))
json.dump({"NM_G": rows, "N1mu": mrows}, open(os.path.join(root, "SUMMARY.json"), "w"), indent=1)
print()
print("TOTALS: NM-G certificates %d, cells %d; N1mu files %d, segment certificates %d, cells %d"
      % (len(nm), sum(d["n_cells"] for d in nm), len(mu),
         sum(1 for d in mu for s in d["segments"] if "a_cert" in s),
         sum(s["n_cells"] for d in mu for s in d["segments"] if "a_cert" in s)))
