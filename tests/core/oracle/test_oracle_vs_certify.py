"""Cross-test: rigorous rtt.core.certify vs the independent float oracle (tests/oracle/oracle.py).

All spectra are SYNTHETIC (tests/oracle/cases/*.json, built by make_suite.py).

For each case and each mu_hi in {0, 0.1, 0.15} (omega_c = 10 * max omega, T* = 300 K):
  * containment: T_low <= Tc_oracle(mu_hi) + err   and   Tc_oracle(0) - err <= T_up
    (Tc is non-increasing in mu*, so these are the binding ends of the S5 conjunction over
    [0, mu_hi]; for envelope cases the oracle min/max over sampled corners of the mass box is
    used, since the certified bracket must hold for every spectrum in the box).
  * verdicts at T* = 300 K:  REJECTED_PROVEN => Tc_oracle(mu*=0) < 300 (+ err);
                             CERTIFIED_CANDIDATE => Tc_oracle(mu_hi) >= 300 (- err).
Bracket widths and timings go to tests/oracle/crosstest_table.{json,md}.

Environment knobs: ORACLE_CASES=comma,separated,names   ORACLE_MU=0,0.1   ORACLE_TSTAR=300
The module skips cleanly when rtt.core.certify is not importable.
"""
from __future__ import annotations

import glob
import inspect
import json
import math
import os
import sys
import time

import pytest

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", ".."))
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)

certify = pytest.importorskip("rtt.core.certify", reason="rtt.core.certify not importable yet")
if not hasattr(certify, "classify"):
    pytest.skip("rtt.core.certify has no classify() yet", allow_module_level=True)

from rtt.core.spectrum import Spectrum  # noqa: E402
from rtt.core.status import Status  # noqa: E402

T_STAR = float(os.environ.get("ORACLE_TSTAR", "300"))
MUS = [float(x) for x in os.environ.get("ORACLE_MU", "0,0.1,0.15").split(",")]
_SEL = os.environ.get("ORACLE_CASES")


def _load_cases():
    cs = []
    for f in sorted(glob.glob(os.path.join(HERE, "cases", "*.json"))):
        d = json.load(open(f))
        if _SEL and d["name"] not in _SEL.split(","):
            continue
        cs.append(d)
    return cs


CASES = _load_cases()
ROWS: list[dict] = []


# ---------------------------------------------------------------------------------------------
# oracle side
# ---------------------------------------------------------------------------------------------
def oracle_tc_range(d, mu):
    """(Tc_lo, Tc_hi, err, floor): oracle Tc at mu (min/max over envelope corners if any).
    Tc None => no instability above floor."""
    key = str(float(mu))
    if d.get("envelope"):
        e = d["envelope"][key]
        corners = list(e["corners"].values())
    else:
        corners = [d["oracle"][key]]
    tcs = [c["Tc"] for c in corners]
    errs = [c["err"] for c in corners if c["Tc"] is not None]
    floors = [c["T_floor"] for c in corners if c["Tc"] is None]
    err = max(errs) if errs else 0.0
    lo = None if any(t is None for t in tcs) else min(tcs)
    hi = max((t for t in tcs if t is not None), default=None)
    floor = max(floors) if floors else None
    return lo, hi, err, floor


# ---------------------------------------------------------------------------------------------
# certify side: adapter around the public API (names from ARCH.md section 4 and __init__ docs)
# ---------------------------------------------------------------------------------------------
def _spectrum(d) -> Spectrum:
    return Spectrum.from_dict(d["spectrum"])


TIME_BUDGET = float(os.environ.get("ORACLE_CERT_BUDGET_S", "600"))


def call_classify(spec: Spectrum, mu_hi: float, omega_c: float, T_star: float):
    """certify.classify(spectra, params, material_id) -> Verdict (public API, certify/__init__)."""
    if hasattr(certify, "CertifyParams"):
        P = certify.CertifyParams
        names = set(inspect.signature(P).parameters)
        kw = dict(T_star=T_star, mu_hi=mu_hi, omega_c=omega_c)
        if "time_budget_s" in names:
            kw["time_budget_s"] = TIME_BUDGET
        params = P(**{k: v for k, v in kw.items() if k in names})
        return certify.classify([spec], params, material_id="oracle-" + spec.label[11:40])
    # fallback: match keyword names to the published signature
    sig = inspect.signature(certify.classify)
    names = set(sig.parameters)
    kw = {}

    def put(cands, val):
        for c in cands:
            if c in names:
                kw[c] = val
                return True
        return False

    args = []
    if not put(("spectra", "alternatives"), [spec]):
        if not put(("spectrum", "spec", "s"), spec):
            args.append(spec)
    put(("T_star", "T_star_K", "t_star", "Tstar"), T_star)
    put(("mu_hi", "mu_star_hi", "mu_max"), mu_hi)
    if not put(("omega_c", "omega_c_K"), omega_c):
        put(("omega_c_factor",), 10.0)
    put(("material_id",), "oracle-test")
    return certify.classify(*args, **kw)


def _status_name(v):
    s = getattr(v, "status", v)
    return s.value if hasattr(s, "value") else str(s)


# ---------------------------------------------------------------------------------------------
@pytest.mark.parametrize("mu_hi", MUS)
@pytest.mark.parametrize("case", CASES, ids=[c["name"] for c in CASES])
def test_classify_consistent_with_oracle(case, mu_hi):
    spec = _spectrum(case)
    omega_c = case["omega_c_K"]
    t = time.time()
    v = call_classify(spec, mu_hi, omega_c, T_STAR)
    dt = time.time() - t
    st = _status_name(v)
    T_low, T_up = getattr(v, "T_low_K", None), getattr(v, "T_up_K", None)
    # binding oracle values: upper end at mu* = 0, lower end at mu* = mu_hi
    _, tc_hi0, err0, floor0 = oracle_tc_range(case, 0.0)
    tc_lo_mu, _, err_mu, floor_mu = oracle_tc_range(case, mu_hi)
    row = dict(case=case["name"], mu_hi=mu_hi, status=st, T_low=T_low, T_up=T_up,
               width=(T_up - T_low) if (T_low is not None and T_up is not None) else None,
               Tc_oracle_mu0_max=tc_hi0, Tc_oracle_muhi_min=tc_lo_mu,
               oracle_floor_muhi=floor_mu, seconds=round(dt, 2), failures=[],
               rule_ids=list(getattr(v, "rule_ids", []) or []),
               reasons=[str(x)[:160] for x in (getattr(v, "reasons", []) or [])][:3])
    ROWS.append(row)
    fails = row["failures"]
    # --- containment
    if T_up is not None and tc_hi0 is not None and not (tc_hi0 - err0 <= T_up):
        fails.append(f"CONTAINMENT: Tc_oracle(mu*=0)={tc_hi0!r} (err {err0:.1e}) > T_up={T_up!r}")
    if T_low is not None and T_low > 0:
        if tc_lo_mu is not None:
            if not (T_low <= tc_lo_mu + err_mu):
                fails.append(f"CONTAINMENT: T_low={T_low!r} > Tc_oracle(mu*={mu_hi})="
                             f"{tc_lo_mu!r} (err {err_mu:.1e})")
        elif floor_mu is not None and T_low > floor_mu:
            fails.append(f"CONTAINMENT: T_low={T_low!r} but oracle finds no instability above "
                         f"{floor_mu!r} at mu*={mu_hi}")
    # --- verdict consistency at T*
    if st == "REJECTED_PROVEN" and tc_hi0 is not None and not (tc_hi0 - err0 < T_STAR):
        fails.append(f"VERDICT: REJECTED_PROVEN but Tc_oracle(mu*=0)={tc_hi0!r} >= {T_STAR}")
    if st == "CERTIFIED_CANDIDATE":
        if tc_lo_mu is None or not (tc_lo_mu + err_mu >= T_STAR):
            fails.append(f"VERDICT: CERTIFIED_CANDIDATE but Tc_oracle(mu*={mu_hi})="
                         f"{tc_lo_mu!r} < {T_STAR}")
    assert not fails, f"{case['name']} mu_hi={mu_hi}: " + "; ".join(fails)


def teardown_module(module):
    if not ROWS:
        return
    json.dump(ROWS, open(os.path.join(HERE, "crosstest_table.json"), "w"), indent=1,
              default=str)
    hdr = "| case | mu_hi | status | T_low | T_up | width | Tc_or(mu0) | Tc_or(mu_hi) | s | fail |"
    lines = [hdr, "|" + "---|" * 10]

    def f(x):
        return "-" if x is None else (f"{x:.6g}" if isinstance(x, float) else str(x))
    for r in ROWS:
        lines.append(f"| {r['case']} | {r['mu_hi']} | {r['status']} | {f(r['T_low'])} | "
                     f"{f(r['T_up'])} | {f(r['width'])} | {f(r['Tc_oracle_mu0_max'])} | "
                     f"{f(r['Tc_oracle_muhi_min'])} | {r['seconds']} | "
                     f"{'FAIL' if r['failures'] else ''} |")
    open(os.path.join(HERE, "crosstest_table.md"), "w").write("\n".join(lines) + "\n")
