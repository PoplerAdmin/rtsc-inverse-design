"""Measured error budget for the S3 ML-hull window (review RC1 F15).

For each H-containing chemical system:

1. ``hull error``: every Alex-MP reference entry is relaxed with MatterSim (the S3 self hull,
   :func:`rtt.mlff.hull.build_self_hull`). For each converged **H-containing** entry,
   err = E_hull(MatterSim-only hull) − E_hull(DFT hull of the same reference set), eV/atom.
2. ``ZPE term``: for H-containing entries with n_atoms <= ``zpe_max_atoms`` whose hull
   decomposition products are also computable, harmonic MatterSim ZPE (positive modes only)
   gives Δ_ZPE = ZPE(entry) − Σ_i x_i ZPE(product_i), eV/atom: the change of the hull distance
   when ZPE is added (doc 02 H7 is ZPE-corrected; the ML hull is not).

Window (eV/atom) = 0.35 (H7) + √2·p95(|err|) + Mg offset (0.10, doc 10 §3.3) + p95(|Δ_ZPE|),
and never below 0.55 (RC1 F15 interim floor). √2 is the reviewer's allowance for the candidate
side of the difference (the reference entries are themselves hull members, so their own error
enters once). Heuristic and MatterSim-based throughout; the numbers are written to
``docs/design/data/21_hull_calibration.json`` and copied into ``tiers.toml [S3]``.

    nice -n 19 ~/.rtt/envs/py/bin/python -m rtt.mlff.calibrate OUT.json H-Ir-Mg H-Li-Mg ...
"""
from __future__ import annotations

import json
import math
import sys
import time
from typing import Any, Sequence

import numpy as np

H7_WINDOW = 0.35
MG_OFFSET = 0.10          # doc 10 §3.3 (meas.): MP-sourced Mg entries +0.10 eV/atom vs MatterSim
FLOOR = 0.55              # RC1 F15 interim minimum


def _p95(x: Sequence[float]) -> float | None:
    return float(np.percentile(np.abs(np.asarray(x, float)), 95)) if len(x) else None


def system_errors(cs: str, *, device: str | None = None, zpe: bool = True, zpe_max_atoms: int = 20,
                  zpe_max_entries: int = 40) -> dict:
    from pymatgen.analysis.phase_diagram import PDEntry
    from pymatgen.core import Composition
    from ..source.external import refquery
    from .hull import build_self_hull
    from .phonons import ml_phonons
    from .relax import relax_batch

    t0 = time.time()
    data = refquery(cs)
    rows = {r["entry_id"]: r for r in data["rows"]}
    hull = build_self_hull(cs, device=device)
    pd = hull.phase_diagram()
    errs, per = [], []
    for e in hull.entries:
        r = rows.get(e["entry_id"])
        if not e["converged"] or r is None or r.get("e_above_hull_dft_eV") is None:
            continue
        if "H" not in r["chemsys"].split("-") or r["chemsys"] == "H":
            continue
        c = Composition(e["formula"])
        eh_ml = pd.get_e_above_hull(PDEntry(c, e["enthalpy_per_atom_eV"] * c.num_atoms), allow_negative=True)
        d = float(eh_ml) - float(r["e_above_hull_dft_eV"])
        errs.append(d)
        per.append({"entry_id": e["entry_id"], "formula": r["formula"], "e_hull_ml": round(float(eh_ml), 4),
                    "e_hull_dft": r["e_above_hull_dft_eV"], "err": round(d, 4)})
    out: dict[str, Any] = {"chemsys": cs, "n_ref": len(hull.entries), "n_H_entries": len(errs), "entries": per,
                           "p95_abs_err": _p95(errs), "mae": float(np.mean(np.abs(errs))) if errs else None}
    if zpe:
        # ZPE for small converged entries (relaxed again: the self hull stores energies only)
        small = sorted([rows[e["entry_id"]] for e in hull.entries
                        if e["converged"] and e["entry_id"] in rows and rows[e["entry_id"]]["n_atoms"] <= zpe_max_atoms],
                       key=lambda r: (r["e_above_hull_dft_eV"] if r["e_above_hull_dft_eV"] is not None else 9, r["n_atoms"]))
        small = small[:zpe_max_entries]
        rel = relax_batch([r["structure"] for r in small], device=device)
        zpe_of: dict[str, float] = {}
        for r, x in zip(small, rel):
            if x.error or not x.converged:
                continue
            try:
                zpe_of[r["entry_id"]] = ml_phonons(x.structure, device=device).zpe_ml_meV_atom / 1000.0
            except Exception:  # noqa: BLE001
                continue
        dz = []
        hfrac = {}
        for r in small:
            eid = r["entry_id"]
            if eid not in zpe_of or "H" not in r["chemsys"].split("-") or r["chemsys"] == "H":
                continue
            e = next(x for x in hull.entries if x["entry_id"] == eid)
            c = Composition(e["formula"])
            decomp = pd.get_decomposition(c)
            if any(k.name not in zpe_of for k in decomp):
                continue
            v = zpe_of[eid] - sum(f * zpe_of[k.name] for k, f in decomp.items())
            dz.append(v)
            hfrac[r["formula"]] = round(v, 4)
        out.update(n_zpe=len(dz), zpe_delta=hfrac, p95_abs_zpe=_p95(dz), zpe_of_meV={k: round(v * 1000, 1) for k, v in zpe_of.items()})
    out["wall_s"] = round(time.time() - t0, 1)
    return out


def window(results: Sequence[dict]) -> dict:
    errs = [e["err"] for r in results for e in r["entries"]]
    dz = [v for r in results for v in (r.get("zpe_delta") or {}).values()]
    p95 = _p95(errs) or 0.0
    pz = _p95(dz) or 0.0
    margin = math.sqrt(2) * p95 + MG_OFFSET + pz
    return {"n_err": len(errs), "p95_abs_hull_err_eV": round(p95, 4), "sqrt2_p95": round(math.sqrt(2) * p95, 4),
            "mg_offset_eV": MG_OFFSET, "n_zpe": len(dz), "p95_abs_zpe_delta_eV": round(pz, 4),
            "margin_eV": round(margin, 4), "window_eV": round(max(FLOOR, H7_WINDOW + margin), 3),
            "floor_applied": H7_WINDOW + margin < FLOOR}


def main(argv: list[str] | None = None) -> int:
    a = list(sys.argv[1:] if argv is None else argv)
    out, systems = a[0], a[1:]
    res = []
    for cs in systems:
        try:
            r = system_errors(cs)
        except Exception as ex:  # noqa: BLE001
            r = {"chemsys": cs, "error": f"{type(ex).__name__}: {ex}", "entries": []}
        res.append(r)
        print(json.dumps({k: r.get(k) for k in ("chemsys", "n_H_entries", "p95_abs_err", "mae", "n_zpe",
                                                 "p95_abs_zpe", "wall_s", "error")}), flush=True)
        with open(out, "w") as fh:
            json.dump({"systems": res, "window": window(res)}, fh, indent=1)
    print(json.dumps(window(res), indent=1))
    return 0


if __name__ == "__main__":
    sys.exit(main())
