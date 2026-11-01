"""Summarise QE matdyn a2F.dosN files: lambda, omega_log, omega_2, A = int a2F dw, Allen-Dynes Tc.

All numbers here are informational float estimates (tier-2 scoring), never verdicts.
A uses the same definition as the certify core (rtt/core/certify/floatcore.py:
A = sum of alpha^2F point weights in K = int alpha^2F(w) dw), and A/4.33148 is the float value
of the universal upper bound k_B Tc <= A/4.33148 (the certified version needs the Arb path).

Usage: python a2f_summary.py <dir-with-a2F.dosN> [<dir2> ...]
"""
from __future__ import annotations

import glob
import math
import os
import re
import sys

import numpy as np

RY_TO_K = 157887.512  # 1 Ry / k_B in K (13.605693 eV / 8.617333e-5 eV/K)


def read_a2f(path: str) -> tuple[np.ndarray, np.ndarray]:
    w, a = [], []
    with open(path) as fh:
        for line in fh:
            s = line.strip()
            if not s or s.startswith("#") or s.startswith("lambda"):
                continue
            parts = s.split()
            try:
                w.append(float(parts[0]))
                a.append(float(parts[1]))
            except (ValueError, IndexError):
                continue
    return np.array(w) * RY_TO_K, np.array(a)  # omega in K, a2F dimensionless


def allen_dynes(lam: float, wlog: float, w2: float, mu: float) -> float:
    if lam <= mu * (1 + 0.62 * lam):
        return 0.0
    L1 = 2.46 * (1 + 3.8 * mu)
    L2 = 1.82 * (1 + 6.3 * mu) * (w2 / wlog)
    f1 = (1 + (lam / L1) ** 1.5) ** (1 / 3)
    f2 = 1 + (w2 / wlog - 1) * lam ** 2 / (lam ** 2 + L2 ** 2)
    return f1 * f2 * wlog / 1.2 * math.exp(-1.04 * (1 + lam) / (lam - mu * (1 + 0.62 * lam)))


def summarise(w: np.ndarray, a: np.ndarray) -> dict:
    # Integrate on the full (uniform) matdyn grid: dropping alpha^2F <= 0 points first and then
    # taking np.gradient on the uneven remainder inflates every integral across phonon gaps
    # (Mg2IrH6: lambda 3.68 vs QE 2.42; doc 15). Clip negatives to 0, skip w <= 0 only.
    dw = np.gradient(w) if len(w) > 1 else np.ones_like(w)
    keep = w > 0
    w, a, dw = w[keep], np.clip(a[keep], 0.0, None), dw[keep]
    lam = 2 * np.sum(a / w * dw)
    A = np.sum(a * dw)
    wlog = math.exp(2 / lam * np.sum(a / w * np.log(w) * dw))
    w2 = math.sqrt(2 / lam * np.sum(a * w * dw))
    return {"lambda": lam, "wlog_K": wlog, "w2_K": w2, "A_K": A, "A_over_4.33148_K": A / 4.33148145,
            "TcAD_mu0.10_K": allen_dynes(lam, wlog, w2, 0.10), "TcAD_mu0.13_K": allen_dynes(lam, wlog, w2, 0.13)}


def main() -> None:
    for d in sys.argv[1:]:
        files = sorted(glob.glob(os.path.join(d, "a2F.dos*")), key=lambda p: int(re.findall(r"\d+$", p)[0]))
        print(f"== {d}")
        for f in files:
            r = summarise(*read_a2f(f))
            print(f"{os.path.basename(f):10s} " + " ".join(f"{k}={v:.3f}" for k, v in r.items()))


if __name__ == "__main__":
    main()
