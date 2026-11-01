"""Render the CSVs in docs/theory/data/ as markdown tables (written to data/tables.md)."""
import csv
import os
import sys

D = sys.argv[1] if len(sys.argv) > 1 else "."


def load(name):
    with open(os.path.join(D, name)) as f:
        return list(csv.DictReader(f))


def F(x, nd=0):
    return f"{float(x):.{nd}f}"


out = []
P = out.append

# T1: min A per shape (cutoff 10 omega_max)
rows = load("min_A.csv")
shapes = list(dict.fromkeys(r["shape"] for r in rows))
mus = ["0", "0.08", "0.1", "0.13", "0.16", "0.2"]
P("### Table 1. Minimal A = ∫α²F dω (K) for Tc = 300 K, per shape, ω_c = 10 ω_max "
  "(λ at the optimum in parentheses)\n")
P("| shape | " + " | ".join(f"μ*={m}" for m in mus) + " |")
P("|---|" + "---|" * len(mus))
for s in shapes:
    cells = []
    for m in mus:
        r = [x for x in rows if x["shape"] == s and x["mu"] == m and x["cmult"] == "10"]
        cells.append(f"{F(r[0]['A_req'])} ({F(r[0]['lam'], 2)})" if r else "–")
    P(f"| {s} | " + " | ".join(cells) + " |")
P("\nSame optimum, required ω_log / sqrt⟨ω²⟩ (K):\n")
P("| shape | " + " | ".join(f"μ*={m}" for m in mus) + " |")
P("|---|" + "---|" * len(mus))
for s in shapes:
    cells = []
    for m in mus:
        r = [x for x in rows if x["shape"] == s and x["mu"] == m and x["cmult"] == "10"]
        cells.append(f"{F(r[0]['wlog_req'])} / {F(r[0]['w2_req'])}" if r else "–")
    P(f"| {s} | " + " | ".join(cells) + " |")
P("\nCutoff sensitivity (ω_c = 5 ω_max vs 10 ω_max), relative change of A_min:\n")
P("| shape | " + " | ".join(f"μ*={m}" for m in mus[1:]) + " |")
P("|---|" + "---|" * (len(mus) - 1))
for s in shapes:
    cells = []
    for m in mus[1:]:
        a5 = [x for x in rows if x["shape"] == s and x["mu"] == m and x["cmult"] == "5"]
        a10 = [x for x in rows if x["shape"] == s and x["mu"] == m and x["cmult"] == "10"]
        cells.append(f"{100*(float(a5[0]['A_req'])/float(a10[0]['A_req'])-1):+.1f}%")
    P(f"| {s} | " + " | ".join(cells) + " |")

# T2: Einstein absolute cutoff
rows = load("einstein_abs.csv")
P("\n### Table 2. Extremal (Einstein) minimal A for Tc = 300 K, absolute cutoff ω_c\n")
P("| μ* | ω_c (K) | A_min (K) | A_min/Tc | Ω_opt (K) | Ω_opt/Tc | λ_opt |")
P("|---|---|---|---|---|---|---|")
for r in rows:
    P(f"| {r['mu']} | {r['omega_c']} | {F(r['A_req'],1)} | {F(r['A_over_Tc'],4)} | "
      f"{F(r['Omega_opt'])} | {F(r['w_opt'],3)} | {F(r['lam_opt'],3)} |")

# T3: soft-mode constant
rows = [r for r in load("soft_mode_constant.csv") if r["N"] == "400"]
P("\n### Table 3. Soft-mode constant c2: Tc ≤ c2·sqrt(λ⟨ω²⟩), ω_c = 20000 K at T = 300 K\n")
P("| μ* | B_c | c2 | λ⟨ω²⟩ required (K²) | sqrt(λ⟨ω²⟩) required (K) | ∫ω α²F dω required (K²) |")
P("|---|---|---|---|---|---|")
for r in rows:
    v = float(r["lamw2_req_300"])
    P(f"| {r['mu']} | {F(r['B_c'],5)} | {F(r['c2'],5)} | {v:.3e} | {v**0.5:.0f} | {v/2:.3e} |")

# T4: two-atom fixed lambda
rows = load("two_atom_fixed_lambda.csv")
P("\n### Table 4. Fixed λ (two-atom extremal search, ω_c = 20000 K): minimal sqrt⟨ω²⟩ and "
  "minimal λ⟨ω²⟩ for Tc = 300 K; the minimiser was the Einstein spectrum in every case\n")
P("| μ* | λ | min sqrt⟨ω²⟩ = Ω_E (K) | min λ⟨ω²⟩ (K²) | min A on grid (K) | ω_log of best-A spectrum (K) |")
P("|---|---|---|---|---|---|")
for r in rows:
    P(f"| {r['mu']} | {r['lam']} | {F(r['min_w2'])} | {float(r['min_lamw2']):.3e} | "
      f"{F(r['minA'])} (r={r['minA_r']}, f={r['minA_f']}) | {F(r['minA_wlog'])} |")

# T5: wlog requirement at fixed lambda per shape
rows = load("t_vs_lambda.csv")
for m in ["0.1", "0.13"]:
    P(f"\n### Table 5{'a' if m == '0.1' else 'b'}. Required ω_log (K) for Tc = 300 K at fixed λ, "
      f"μ* = {m}, ω_c = 10 ω_max\n")
    lams = ["1.5", "2", "2.5", "3", "4"]
    P("| shape | " + " | ".join(f"λ={l}" for l in lams) + " |")
    P("|---|" + "---|" * len(lams))
    for s in shapes:
        cells = []
        for l in lams:
            r = [x for x in rows if x["shape"] == s and x["mu"] == m and x["cmult"] == "10"
                 and x["lam"] == l]
            cells.append(F(r[0]["wlog_req"]) if r else "–")
        P(f"| {s} | " + " | ".join(cells) + " |")

# T6: unbounded wlog
rows = load("unbounded_wlog_demo.csv")
P("\n### Table 6. Two atoms, λ = 2 split equally, frequency ratio r, μ* = 0\n")
P("| r | Tc/ω_log | Tc/A | Tc/sqrt(λ⟨ω²⟩) |")
P("|---|---|---|---|")
for r in rows:
    P(f"| {r['r']} | {F(r['Tc_over_wlog'],4)} | {F(r['Tc_over_A'],4)} | "
      f"{F(r['Tc_over_sqrt_lamw2'],4)} |")

# T7: BR peaks
rows = load("br_derivative.csv")
P("\n### Table 7. Bergmann–Rainer derivative δTc/δα²F(ω) at Tc = 300 K (ω_c = 10 ω_max)\n")
P("| shape | λ | μ* | peak ω/Tc | peak ω (K) | peak value | value at ω = 2Tc | value at ω = 20Tc |")
P("|---|---|---|---|---|---|---|---|")
for s in dict.fromkeys(r["shape"] for r in rows):
    rr = [x for x in rows if x["shape"] == s]
    pk = max(rr, key=lambda x: float(x["dTc_da2F"]))
    v2 = min(rr, key=lambda x: abs(float(x["omega_over_Tc"]) - 2))
    v20 = min(rr, key=lambda x: abs(float(x["omega_over_Tc"]) - 20))
    P(f"| {s} | {rr[0]['lam']} | {rr[0]['mu']} | {F(pk['omega_over_Tc'],2)} | "
      f"{300*float(pk['omega_over_Tc']):.0f} | {F(pk['dTc_da2F'],4)} | {F(v2['dTc_da2F'],4)} | "
      f"{F(v20['dTc_da2F'],4)} |")

rnd = load("random_spectra.csv")
P(f"\nRandom μ*=0 spectra (n={len(rnd)}, 1–5 atoms): max Tc/A = "
  f"{max(float(r['Tc_over_A']) for r in rnd):.5f}, max Tc/sqrt(λ⟨ω²⟩) = "
  f"{max(float(r['Tc_over_sqrt_lamw2']) for r in rnd):.5f}, max Tc/ω_log = "
  f"{max(float(r['Tc_over_wlog']) for r in rnd):.4f}")

with open(os.path.join(D, "tables.md"), "w") as f:
    f.write("\n".join(out) + "\n")
print("\n".join(out))
