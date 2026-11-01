#!/usr/bin/env python3
"""Generate lean/TargetLogic/NMTable.lean from docs/theory/certificates/nm_*.json.
Usage (from the repo root): python3 lean/scripts/gen_nm_table.py > lean/TargetLogic/NMTable.lean
Reads the exact decimal `m0_over_T2` of every certificate (never the float g_lo_K2)."""
import glob, hashlib, json, os, sys
from fractions import Fraction

d0 = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "docs", "theory",
                  "certificates")
rows, inf = [], None
for f in sorted(glob.glob(os.path.join(d0, "nm_*.json"))):
    raw = open(f, "rb").read(); d = json.loads(raw)
    if d.get("verdict", "CERTIFIED") != "CERTIFIED":
        sys.exit(f"{f}: not CERTIFIED")
    h = hashlib.sha256(raw).hexdigest(); b = os.path.basename(f)
    if d["lambda"] == "inf":
        inf = (d["m0_over_T2"], h, b)
    else:
        rows.append((Fraction(d["lambda"]), d["lambda"], d["m0_over_T2"], h, b))
rows.sort()
q = lambda s: f"({Fraction(s).numerator} / {Fraction(s).denominator} : K)"
print("/-\nTargetLogic.NMTable: the certified NM-G table (doc 01 §0.3, [THM-C]).\n"
      "GENERATED from docs/theory/certificates/nm_*.json (exact decimal `m0_over_T2`; do not edit).\n"
      "Entry `(λ_k, m0_k)`: certificate Theorem NM(λ_k, m0_k) (CERT-NM); g_lo(λ_k) = 90000·m0_k at 300 K.\n"
      "Source files and sha256:")
for r in rows:
    print(f"  {r[4]}  {r[3]}")
print(f"  {inf[2]}  {inf[1]}")
print("-/\nimport Mathlib.Algebra.Order.Field.Basic\n\nnamespace TargetLogic\n\nsection\n"
      "variable {K : Type*} [Field K]\n\n/-- `(λ_k, m0_k)` from the certificates, sorted by λ_k. -/\n"
      "def nmM0Table : List (K × K) :=\n  ["
      + ",\n   ".join(f"({q(r[1])}, {q(r[2])})" for r in rows) + "]\n\n"
      "/-- `m0` of the λ = ∞ certificate. -/\ndef nmM0Inf : K := " + q(inf[0]) + "\n\nend\n\nend TargetLogic")
