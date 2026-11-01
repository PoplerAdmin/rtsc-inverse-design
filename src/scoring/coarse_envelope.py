#!/usr/bin/env python3
"""Apply the theory-doc ENVELOPE rule (docs/theory/01_eliashberg_targets.md §0.2(b), §0.4) to a
COARSE ladder and compare the coarse upper envelope with the full reference's nominal value.

Input: the per-sigma CSV written by coarse_calibration.py. For each (material, q_coarse run,
q_fine run) pair of coarse runs:
  X_nom   = X(q_fine) at the smallest sigma admitted in BOTH coarse runs (assumption: the theory
            doc does not define X_nom explicitly; this is the natural reading)
  d_sig   = max_i |X(q_fine, s_i) - X_nom| / X_nom over sigma admitted in q_fine
  d_q     = |X(q_fine, s_nom) - X(q_coarse, s_nom)| / X_nom
  eps_w   = max(|w_log(q_fine) - w_log(q_coarse)| / w_log(q_fine) at s_nom, 0.02), capped at 0.9
  A_hi    = A_nom (1 + d_sig + d_q)(1 + eps_sys) / (1 - eps_w),  eps_sys = 0.10
  m_hi    = m_nom (1 + d_sig + d_q)(1 + eps_sys)
Output columns: ratios A_hi_coarse / A_ref_nom and m_hi_coarse / m_ref_nom (reference nominal =
smallest admitted reference sigma), and A_hi_coarse / A_ref_hi_lb, where
  A_ref_hi_lb = A_ref_nom (1 + d_sig_ref)(1 + eps_sys) / (1 - 0.02)
is a LOWER bound on the reference's own ENVELOPE value (the bench has one q mesh, so its d_q is
unknown and set to 0; eps_w set to its floor 0.02). A_ref_hi is what the full-tier N1-ISO rule
compares with the target, so a coarse screen-out is "safe" for a material when the coarse
statistic times (1+e_c) is >= A_ref_hi.
A second block (q_fine_run = '-') gives, per single coarse run, the doc-11 statistic
max over admitted sigma of A_coarse, divided by A_ref_nom and by A_ref_hi_lb. Informational only.

Usage: python coarse_envelope.py per_sigma.csv MAT:coarse_run:fine_run [...] > out.csv
"""
from __future__ import annotations

import csv
import sys
from collections import defaultdict


def main(argv: list[str]) -> int:
    rows = list(csv.DictReader(open(argv[0])))
    tab: dict[tuple[str, str], dict[int, dict]] = defaultdict(dict)
    for r in rows:
        tab[(r["material"], r["run"])][int(r["isig"])] = r
    out = csv.writer(sys.stdout)
    out.writerow(["material", "q_coarse_run", "q_fine_run", "sigma_nom_Ry", "A_nom_K", "d_sig_A", "d_q_A",
                  "eps_w", "A_hi_K", "A_ref_nom_K", "A_hi_over_ref_nom", "A_nom_over_ref_nom",
                  "m_hi_over_ref_nom", "lambda_nom_over_ref_nom", "A_ref_hi_lb_K", "A_hi_over_ref_hi_lb"])

    def ref_hi_lb(mat: str) -> tuple[float, float]:
        ref = tab[(mat, "ref")]
        adm_r = sorted(i for i, r in ref.items() if r["admitted"] == "True")
        a0 = float(ref[adm_r[0]]["A_K"])
        ds = max(abs(float(ref[i]["A_K"]) - a0) / a0 for i in adm_r)
        return a0, a0 * (1 + ds) * 1.10 / 0.98
    for spec in argv[1:]:
        mat, rc, rf = spec.split(":")
        c, f, ref = tab[(mat, rc)], tab[(mat, rf)], tab[(mat, "ref")]
        adm = lambda t: sorted(i for i, r in t.items() if r["admitted"] == "True")  # noqa: E731
        both = [i for i in adm(f) if i in adm(c)]
        if not both or len(adm(f)) < 3:
            out.writerow([mat, rc, rf, "ladder incomplete (cannot form envelope -> no screen-out)"])
            continue
        s = both[0]
        g = lambda t, i, k: float(t[i][k])  # noqa: E731
        res = {}
        for k in ("A_K", "m_K2", "lambda"):
            xn = g(f, s, k)
            dsig = max(abs(g(f, i, k) - xn) / xn for i in adm(f))
            dq = abs(xn - g(c, s, k)) / xn
            res[k] = (xn, dsig, dq)
        ew = min(max(abs(g(f, s, "wlog_K") - g(c, s, "wlog_K")) / g(f, s, "wlog_K"), 0.02), 0.9)
        A_nom, dsA, dqA = res["A_K"]
        A_hi = A_nom * (1 + dsA + dqA) * 1.10 / (1 - ew)
        m_nom, dsm, dqm = res["m_K2"]
        m_hi = m_nom * (1 + dsm + dqm) * 1.10
        r0 = adm(ref)[0]
        out.writerow([mat, rc, rf, f[s]["sigma_Ry"], f"{A_nom:.4g}", f"{dsA:.3f}", f"{dqA:.3f}", f"{ew:.3f}",
                      f"{A_hi:.4g}", f"{g(ref, r0, 'A_K'):.4g}", f"{A_hi / g(ref, r0, 'A_K'):.3f}",
                      f"{A_nom / g(ref, r0, 'A_K'):.3f}", f"{m_hi / g(ref, r0, 'm_K2'):.3f}",
                      f"{res['lambda'][0] / g(ref, r0, 'lambda'):.3f}",
                      f"{ref_hi_lb(mat)[1]:.4g}", f"{A_hi / ref_hi_lb(mat)[1]:.3f}"])
    # single-run doc-11 statistic
    for (mat, run), t in sorted(tab.items()):
        if run == "ref":
            continue
        adm_c = [i for i, r in t.items() if r["admitted"] == "True"]
        a_nom, a_hi = ref_hi_lb(mat)
        if not adm_c:
            out.writerow([mat, run, "-", "no admitted smearing (cannot screen out)"])
            continue
        amax = max(float(t[i]["A_K"]) for i in adm_c)
        out.writerow([mat, run, "-", "maxadm", f"{amax:.4g}", "", "", "", "", f"{a_nom:.4g}", "",
                      f"{amax / a_nom:.3f}", "", "", f"{a_hi:.4g}", f"{amax / a_hi:.3f}"])
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
