"""R-MULTI: multiband universal ceiling (ARCH section 4c; wMB Prop. N audited SOUND; opt-in via
CertifyParams.rejection_model = "multiband", not the default).

Model: multiband isotropic-per-band Migdal-Eliashberg (wMB/PROOF.md Theorem 1), band coupling
weights A_ij = int alpha^2F_ij (K), DOS weights N_i > 0, Coulomb matrix with N_i mu*_ij PSD.
Theorem (Lean MultiBand.R_UNIV_MB / R_UNIV_MB_rowsum; hypotheses H(a0,eta) and Prop N):
if some v > 0 has A v <= kappa v, (1+e) kappa <= kappa_hi and kappa_hi < a0 T*, then every
spectrum with A' <= (1+e) A entrywise is stable at every T >= T*, for every PSD Coulomb matrix and
every cutoff.

Checks (all exact rational or Arb):
  route "eig"    : W = N_i A_ij (exact), W-bar = (W + W^T)/2 if W is not exactly symmetric
                   (assumption "W symmetrized"), K_ij = W-bar_ij / (r_i r_j) with r_i = sqrt(N_i)
                   (Arb), u > 0 (doubles), and (1+e) (K u)_i <= kappa_hi u_i for every i  (Lean
                   sub_of_K -> A v <= kappa v with kappa = kappa_hi/(1+e), v = u/r; sub_env).
  route "rowsum" : (1+e) max_i sum_j A_ij <= kappa_hi  (or the stated row-sum upper bounds;
                   Lean R_UNIV_MB_rowsum).
  decision       : kappa_hi < a0 T*   (strict, exact rationals).
  route "cw"     : direct Collatz-Wielandt vector on A (no DOS needed; referee note (a)):
                   v > 0 (doubles) with (1+e) (A v)_i <= kappa_hi v_i in exact rationals (Lean
                   R_UNIV_MB takes any v > 0 with A v <= kappa v).
BandCoupling objects (team DFT, rtt.core.dfpt.band_coupling) are read by duck typing; their float
helpers are never used.

PITFALL (referee note (b), NON-CONSERVATIVE): A_ij must be in the convention of wMB Theorem 1,
lambda_i = sum_j lambda_ij (band i's mass enhancement is the ROW SUM).  Band-resolved alpha^2F_ij
given as "contributions" normalised to sum to the isotropic alpha^2F (i.e. W/N_tot) make kappa
too small, and a superconductor can be wrongly REJECTED.  Data whose provenance carries a
normalisation tag mentioning "contribution" (or "isotropic"-normalised sums) are refused.

ENVELOPE (referee note (c)): rel_envelope e means  true A_ij in [(1-e), (1+e)] x nominal A_ij.
If an error estimate says nominal = (1-delta) x true, the correct e is delta/(1-delta).
"""
from __future__ import annotations

import hashlib
import json
import math

import numpy as np
from flint import arb, fmpq

from .arbcore import q_of_float
from .search import A0_Q

LEAN_MB = ["MultiBand.R_UNIV_MB", "MultiBand.R_UNIV_MB_eig", "MultiBand.R_UNIV_MB_rowsum",
           "MultiBand.sub_of_K", "MultiBand.sub_env", "MultiBand.thm1_even", "MultiBand.thm3_odd"]


def bc_payload(bc) -> dict:
    """The numbers of a BandCoupling that the check uses (exact doubles)."""
    get = lambda name: getattr(bc, name, None) if not isinstance(bc, dict) else bc.get(name)
    A = get("A_K")
    N = get("N")
    R = get("row_sum_upper_K")
    prov = get("provenance") or {}
    return {"kind": get("kind"), "labels": list(get("labels") or []),
            "A_K": [list(map(float, r)) for r in A] if A is not None else None,
            "N": [float(x) for x in N] if N is not None else None,
            "row_sum_upper_K": [float(x) for x in R] if R is not None else None,
            "rel_envelope": float(get("rel_envelope") or 0.0),
            "label": get("label") or "",
            "data_assumptions": _data_assumptions(prov),
            "refusal": normalisation_refusal(prov)}


def bc_digest(payload: dict) -> str:
    d = {k: payload.get(k) for k in ("kind", "A_K", "N", "row_sum_upper_K", "rel_envelope",
                                     "refusal", "data_assumptions")}

    def hexify(o):
        if isinstance(o, float):
            return o.hex()
        if isinstance(o, list):
            return [hexify(x) for x in o]
        return o
    s = json.dumps({k: hexify(v) for k, v in d.items()}, sort_keys=True)
    return hashlib.sha256(s.encode()).hexdigest()


def _data_assumptions(prov) -> list:
    """Assumptions and caveats carried by the data (referee R1): provenance['assumptions'],
    provenance['caveats'], the stated derived bound, and for EPW data the EPW conditions."""
    if not isinstance(prov, dict):
        return []
    out = [str(a) for a in (prov.get("assumptions") or [])]
    out += [f"DATA CAVEAT: {c}" for c in (prov.get("caveats") or [])]
    if prov.get("bound"):
        out.append(f"DATA BOUND: {prov['bound']}")
    if "EPW" in str(prov.get("format", "")):
        out.append("EPW: omega_max bounds all fine-grid phonon frequencies; modes with "
                   "wf <= eps_acustic excluded from lambda_nk")
    return out


_NORM_KEYS = ("normalisation", "normalization", "convention", "a2f_normalisation",
              "a2f_normalization", "pair_normalisation", "pair_normalization")


def normalisation_refusal(payload_prov: dict) -> str | None:
    """Referee note (b): refuse 'contribution'-normalised band data (W/N_tot convention)."""
    if not isinstance(payload_prov, dict):
        return None
    for k, v in payload_prov.items():
        key = str(k).lower()
        val = json.dumps(v).lower() if not isinstance(v, str) else v.lower()
        if any(t in key for t in _NORM_KEYS) and ("contribution" in val or "isotropic" in val
                                                   or "n_tot" in val or "w/n" in val):
            return f"provenance {k}={v!r}: contribution-normalised alpha^2F_ij refused (wMB convention needed)"
    for a in payload_prov.get("assumptions", []) or []:
        if "contribution" in str(a).lower():
            return f"provenance assumption {a!r}: contribution-normalised data refused"
    return None


def _round_up(x: float, ulps: int = 4) -> float:
    for _ in range(ulps):
        x = math.nextafter(x, math.inf)
    return x


def propose_rmulti(payload: dict) -> dict:
    """Float proposal (untrusted): kappa_hi, u (eig route) and the row-sum route value."""
    e = payload["rel_envelope"]
    out = {"routes": []}
    A, N = payload["A_K"], payload["N"]
    if payload["kind"] == "row_sum_upper_bound":
        A = None                                  # referee R3: only the declared R_i count
    if A is not None and N is not None:
        An = np.asarray(A, float)
        Nn = np.asarray(N, float)
        W = Nn[:, None] * An
        sym = bool(np.array_equal(W, W.T))
        Wb = W if sym else 0.5 * (W + W.T)
        r = np.sqrt(Nn)
        K = Wb / r[:, None] / r[None, :]
        k0 = float(np.max(np.linalg.eigvalsh(K)))
        kc = k0 * (1 + 1e-9) + 1e-300
        try:
            u = np.linalg.solve(kc * np.eye(len(Nn)) - K, np.ones(len(Nn)))
        except np.linalg.LinAlgError:
            u = None
        if u is not None and np.all(np.isfinite(u)) and np.all(u > 0):
            out["routes"].append({"route": "eig", "kappa_hi": _round_up((1 + e) * kc * (1 + 1e-12)),
                                  "u": [float(x) for x in u], "symmetrized": not sym})
    if A is not None and N is None:          # cw route only without DOS (else "eig")
        An = np.asarray(A, float)
        B = An.shape[0]
        ev = np.linalg.eigvals(An)
        k0 = float(np.max(np.abs(ev)))
        kc = k0 * (1 + 1e-9) + 1e-300
        try:
            v = np.linalg.solve(kc * np.eye(B) - An, np.ones(B))
        except np.linalg.LinAlgError:
            v = None
        if v is not None and np.all(np.isfinite(v)) and np.all(v > 0):
            out["routes"].append({"route": "cw", "kappa_hi": _round_up((1 + e) * kc * (1 + 1e-12)),
                                  "v": [float(x) for x in v]})
    if A is not None:
        rs = max(math.fsum(row) for row in A)
        out["routes"].append({"route": "rowsum", "source": "A_K",
                              "kappa_hi": _round_up((1 + e) * rs * (1 + 1e-12))})
    if payload["row_sum_upper_K"] is not None:
        rs = max(payload["row_sum_upper_K"])
        out["routes"].append({"route": "rowsum", "source": "row_sum_upper_K",
                              "kappa_hi": _round_up((1 + e) * rs * (1 + 1e-12))})
    return out


def check_route(payload: dict, route: dict) -> dict:
    """Rigorous check that kappa_hi is a valid R-MULTI bound for this data (not the T* test)."""
    dims = [len(x) for x in (payload.get("A_K"), payload.get("N"), payload.get("row_sum_upper_K"))
            if x is not None]
    if not dims or min(dims) == 0:
        return {"ok": False, "reason": "empty band data (B = 0): never fires (G12)"}
    e_q = q_of_float(payload["rel_envelope"])
    if not (0 <= payload["rel_envelope"] < 1):
        return {"ok": False, "reason": "rel_envelope must be in [0,1)"}
    kh = float(route["kappa_hi"])
    if not (math.isfinite(kh) and kh > 0):
        return {"ok": False, "reason": "kappa_hi must be finite and > 0"}
    kh_q = q_of_float(kh)
    A = payload["A_K"]
    if payload["kind"] == "row_sum_upper_bound":
        # referee R3: for this kind only the declared row-sum bounds R_i have a meaning
        if not (route.get("route") == "rowsum" and route.get("source") == "row_sum_upper_K"):
            return {"ok": False, "reason": "kind row_sum_upper_bound: only the declared R_i route"}
        A = None
    if A is not None and any((not math.isfinite(x)) or x < 0 for row in A for x in row):
        return {"ok": False, "reason": "A_ij must be finite and >= 0"}
    if route["route"] == "rowsum":
        if route.get("source") == "row_sum_upper_K":
            vals = payload["row_sum_upper_K"]
            if vals is None or any((not math.isfinite(x)) or x < 0 for x in vals):
                return {"ok": False, "reason": "row_sum_upper_K missing/invalid"}
            rs = [q_of_float(x) for x in vals]
        else:
            if A is None:
                return {"ok": False, "reason": "A_K missing"}
            rs = []
            for row in A:
                s = fmpq(0)
                for x in row:
                    s += q_of_float(x)
                rs.append(s)
        ok = all((1 + e_q) * s <= kh_q for s in rs)
        return {"ok": ok, "reason": None if ok else "(1+e) max row sum > kappa_hi",
                "lean": "MultiBand.R_UNIV_MB_rowsum"}
    if route["route"] == "cw":
        v = route.get("v")
        if A is None or v is None:
            return {"ok": False, "reason": "cw route needs A_K and v"}
        B = len(A)
        if any(len(r) != B for r in A) or len(v) != B:
            return {"ok": False, "reason": "dimension mismatch"}
        if any((not math.isfinite(x)) or x <= 0 for x in v):
            return {"ok": False, "reason": "v must be > 0"}
        vq = [q_of_float(x) for x in v]
        for i in range(B):
            s = fmpq(0)
            for j in range(B):
                s += q_of_float(A[i][j]) * vq[j]
            if not ((1 + e_q) * s <= kh_q * vq[i]):
                return {"ok": False, "reason": f"(1+e)(A v)_{i} <= kappa_hi v_{i} fails"}
        return {"ok": True, "lean": "MultiBand.R_UNIV_MB (v > 0, A v <= kappa v; sub_env)"}
    if route["route"] == "eig":
        N = payload["N"]
        u = route.get("u")
        if A is None or N is None or u is None:
            return {"ok": False, "reason": "eig route needs A_K, N, u"}
        B = len(N)
        if len(A) != B or any(len(r) != B for r in A) or len(u) != B:
            return {"ok": False, "reason": "dimension mismatch"}
        if any((not math.isfinite(x)) or x <= 0 for x in N):
            return {"ok": False, "reason": "N_i must be > 0"}
        if any((not math.isfinite(x)) or x <= 0 for x in u):
            return {"ok": False, "reason": "u must be > 0"}
        Wq = [[q_of_float(N[i]) * q_of_float(A[i][j]) for j in range(B)] for i in range(B)]
        sym = all(Wq[i][j] == Wq[j][i] for i in range(B) for j in range(B))
        Wb = Wq if sym else [[(Wq[i][j] + Wq[j][i]) / 2 for j in range(B)] for i in range(B)]
        r = [arb(q_of_float(x)).sqrt() for x in N]
        e1 = 1 + arb(e_q)
        ka = arb(kh_q)
        for i in range(B):
            s = arb(0)
            for j in range(B):
                s += arb(Wb[i][j]) / (r[i] * r[j]) * arb(u[j])
            if not bool(e1 * s <= ka * arb(u[i])):
                return {"ok": False, "reason": f"(1+e)(K u)_{i} <= kappa_hi u_{i} not certified"}
        return {"ok": True, "symmetrized": not sym, "lean": "MultiBand.R_UNIV_MB (sub_of_K, sub_env)"}
    return {"ok": False, "reason": "unknown route"}


def rmulti_for(bc, T_star: float) -> dict:
    """Propose + check: best certified route for one BandCoupling.  Returns a certificate block
    {payload, digest, route, kappa_hi}.

    PITFALL (non-conservative): A_ij must follow lambda_i = sum_j lambda_ij (wMB Thm 1).  Pair
    spectra normalised as "contributions" summing to the isotropic alpha^2F (W/N_tot) give a
    kappa that is too small -> wrong rejections.  Such data (provenance normalisation tag) are
    refused: the block then has route None and R-MULTI does not fire.
    Envelope: true A in [(1-e),(1+e)] x nominal; nominal = (1-delta) true needs e = delta/(1-delta).
    """
    payload = bc_payload(bc)
    if payload.get("refusal"):
        return {"payload": payload, "digest": bc_digest(payload), "route": None,
                "kappa_hi": None, "refused": payload["refusal"]}
    prop = propose_rmulti(payload)
    best = None
    for rt in prop["routes"]:
        res = check_route(payload, rt)
        if res["ok"] and (best is None or rt["kappa_hi"] < best["kappa_hi"]):
            best = dict(rt)
            best["check"] = res
    return {"payload": payload, "digest": bc_digest(payload), "route": best,
            "kappa_hi": best["kappa_hi"] if best else None}


def rmulti_evidence(block: dict, T_star: float) -> tuple[bool, dict]:
    """Re-derive R-MULTI evidence from a stored block (exact + Arb only)."""
    info = {}
    if not (isinstance(T_star, (int, float)) and math.isfinite(T_star) and T_star > 0):
        return False, {"reason": "T_star must be finite and > 0", "fires": False}   # referee R4
    payload = block["payload"]
    if bc_digest(payload) != block.get("digest"):
        return False, {"reason": "band-coupling digest mismatch"}
    ref = normalisation_refusal({"assumptions": payload.get("data_assumptions", [])}) \
        or payload.get("refusal")
    if ref:
        return False, {"reason": ref}
    rt = block.get("route")
    if rt is None:
        return False, {"reason": "no certified route"}
    res = check_route(payload, rt)
    info["route_check"] = res
    if not res["ok"]:
        return False, info
    kh_q = q_of_float(float(rt["kappa_hi"]))
    fires = bool(kh_q < A0_Q * q_of_float(float(T_star)))       # strict (Lean R_UNIV_MB)
    info["kappa_hi"] = float(rt["kappa_hi"])
    tq = kh_q / A0_Q                                              # exact kappa_hi / a0
    tb = int(tq.p) / int(tq.q)
    while q_of_float(tb) < tq:                                    # round up
        tb = math.nextafter(tb, math.inf)
    info["T_bound"] = tb
    info["fires"] = fires
    assum = ["M-MULTIBAND (isotropic per band; wMB Thm 1; Prop. N audited SOUND)",
             "A_ij in the lambda_i = sum_j lambda_ij convention (not W/N_tot contributions)",
             "Coulomb N_i mu*_ij PSD", f"E-ENV-MB({payload['rel_envelope']}) entrywise A' <= (1+e)A"]
    if rt["route"] == "eig":
        assum.append("detailed balance with the same N_i for the true model; the envelope covers "
                     "A, not N")
    else:
        assum.append("detailed balance with some N > 0 (N_i A_ij = N_j A_ji; wMB Thm 1 hypothesis; "
                     "symmetrisability / Kolmogorov cycle condition NOT checked by the code)")
    if rt["route"] == "eig" and res.get("symmetrized"):
        assum.append("W symmetrized (verdict about the symmetrized model W-bar = (W+W^T)/2)")
    if payload["kind"] != "computed":
        assum.append(f"band data kind = {payload['kind']}: an upper bound only under the data "
                     f"conditions recorded below (DATA BOUND / DATA CAVEAT)")
    assum += list(payload.get("data_assumptions", []))          # e.g. M-BAND-ISO propagated
    info["assumptions"] = assum
    return fires, info
