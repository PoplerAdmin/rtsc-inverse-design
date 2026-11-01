"""Re-check every certificate JSON in docs/theory/certificates/ from scratch (Arb, 128 bits).

For each file the exact parameters are re-read (decimal strings / dyadic floats), the analytic end
lemmas are re-evaluated, every stored cell is re-certified without splitting, the cells are
checked to tile the stated w-range exactly, and the Rayleigh (upper-side) certificates are
re-evaluated. Float data in the JSON (float_g, margins, timings) are informational only.

usage (from the repo root, WSL):  python src/certify/verify.py [docs/theory/certificates] [--quick]
--quick re-checks every 10th cell only (smoke test; not a full verification).
Exit code 0 iff everything re-verifies.
"""
import glob
from decimal import Decimal
from fractions import Fraction
import json
import math
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from flint import arb  # noqa: E402
from families import MFamily, AFamily  # noqa: E402
from cover import cover, tiled  # noqa: E402
from rayleigh import einstein_unstable  # noqa: E402
from run_nm import ends_check  # noqa: E402
from run_n1mu import ends_A, A0_UNIV, Nc_of  # noqa: E402


def check_cells(fam, eta, ws, W, cells, N, K, quick):
    cells = [tuple(c) for c in cells]
    if not tiled(cells, ws, W):
        return False, "cells do not tile [%r, %r]" % (ws, W)
    sub = cells[::10] if quick else cells
    ok, done, failed, worst, t = cover(fam, eta, ws, W, N, K, cells=sub)
    if failed:
        return False, "%d cells failed re-check: %s" % (len(failed), failed[:3])
    return True, "%d/%d cells re-certified (min scaled margin %.2e, %.0fs)" % (len(done), len(cells),
                                                                                 worst, t)


def verify_nm(d, quick):
    lam_s = d["lambda"]
    pre = []
    # proof requirements (R03 finding 4): beta >= 0 (needed for "lambda <= lambda_k"),
    # beta = 0 for lambda = inf, upper certificate refers to the same lambda,
    # reported g_lo does not exceed the exact certified 90000*m0.
    if not bool(arb(d["beta"]) >= 0):
        pre.append("beta < 0")
    if lam_s == "inf" and not bool(arb(d["beta"]) == 0):
        pre.append("beta != 0 at lambda = inf")
    if lam_s != "inf" and d["upper"]["lambda"] != lam_s:
        pre.append("upper certificate lambda mismatch")
    if not Fraction(float(d["g_lo_K2"])) <= 90000 * Fraction(Decimal(d["m0_over_T2"])):
        pre.append("g_lo_K2 > 90000*m0")
    if pre:
        return False, False, "NM lam=%s: PRECONDITION FAILED: %s" % (lam_s, "; ".join(pre))
    fam = MFamily(d["m0_over_T2"], None if lam_s == "inf" else lam_s, d["beta"])
    eta = arb(d["eta"])
    ws, W = d["w_range"]
    ends = ends_check(fam, eta, ws, W)
    ok_e = ends["ok"]
    ok_c, msg = check_cells(fam, eta, ws, W, d["cells"], d["N"], d["K"], quick)
    u = d["upper"]
    ok_u, r, _ = einstein_unstable(u["lambda"], u["w"], N=300 if (lam_s == "inf" or float(lam_s) > 20)
                                   else 200)
    ok = ok_e and ok_c
    return ok, ok_u, "NM lam=%s: ends=%s cells: %s; upper(Einstein unstable)=%s" % (lam_s, ok_e, msg, ok_u)


def verify_n1mu(d, quick):
    A_floor = arb(d["A_floor_K"])
    wc = d["omega_c_K"]
    T_u = A_floor / arb(A0_UNIV)
    j300 = Nc_of(300.0, wc)
    assert j300 == d["N_c_at_300"]
    eta = arb(d["eta"])
    ok = True
    msgs = []
    segs = {s["j"]: s for s in d["segments"]}
    # segment j = 0 (T > omega_c/pi, no Coulomb term): must lie in the universal range T >= T_u
    # (R03 finding 6); otherwise a mu* = 0 certificate at a = A_floor*pi/omega_c would be needed.
    if not bool(arb(wc) / arb.pi() >= T_u):
        ok = False
        msgs.append("j=0 segment NOT covered (omega_c/pi < A_floor/a0)")
    else:
        msgs.append("j=0 covered by universal bound")
    for j in range(j300, 0, -1):
        Tlo = arb(300) if j == j300 else arb(wc) / ((2 * j + 1) * arb.pi())
        if bool(Tlo >= T_u):
            msgs.append("j=%d covered by universal bound" % j)
            continue
        s = segs.get(j)
        if s is None or "a_cert" not in s:
            ok = False
            msgs.append("j=%d MISSING" % j)
            continue
        a = arb(s["a_cert"])
        if not bool(a * Tlo >= A_floor):
            ok = False
            msgs.append("j=%d a_cert too small" % j)
            continue
        fam = AFamily(s["a_cert"], d["mu_star"], j)
        ws, W, ends = ends_A(a, eta)
        okc, m = check_cells(fam, eta, ws, W, s["cells"], d["N"], d["K"], quick)
        ok &= ends["ok"] and okc
        msgs.append("j=%d ends=%s %s" % (j, ends["ok"], m))
    return ok, "N1mu mu*=%s wc=%g A_floor=%s: %s" % (d["mu_star"], wc, d["A_floor_K"], "; ".join(msgs))


def main():
    args = [a for a in sys.argv[1:] if not a.startswith("--")]
    quick = "--quick" in sys.argv
    root = args[0] if args else os.path.join(HERE, "..", "..", "docs", "theory", "certificates")
    allok = True
    for p in sorted(glob.glob(os.path.join(root, "*.json"))):
        d = json.load(open(p))
        if d.get("rule") == "NM-G":
            ok, ok_u, msg = verify_nm(d, quick)
            allok &= ok and ok_u
        elif d.get("rule") == "N1mu":
            ok, msg = verify_n1mu(d, quick)
            allok &= ok
        else:
            continue
        print(("PASS " if ok else "FAIL ") + os.path.basename(p) + "  " + msg, flush=True)
    print("ALL CERTIFICATES %s%s" % ("VERIFIED" if allok else "NOT VERIFIED",
                                     " (quick mode: subset of cells)" if quick else ""))
    sys.exit(0 if allok else 1)


if __name__ == "__main__":
    main()
