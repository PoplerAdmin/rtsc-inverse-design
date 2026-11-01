"""Control lane C (architecture §3.3 table "Control lane C"; contracts R5, R6, R8).

Known structures are injected directly (no generation). Each control runs every stage at its
surrogate target T*_s; a control removed by a heuristic screen is a screen false negative (R8).

Structures (starting geometries; every stage relaxes them itself):

* **Mg2IrH6** Fm-3m: the Alex-MP DFT entry ``msft-2023062809350933252138247157-GGA``,
  a_conv = 6.632576 Å, H 24e x = 0.260649 (read from the reference LMDB, 2026-10-28).
  Dolui et al., PRL 132, 166001 (2024) give PBE a = 6.655 Å.
* **H3S** Im-3m at 200 GPa: the starting geometry ``H3S_Im-3m.vasp`` of an earlier benchmark run
  (a = 3.09 Å conventional, "approximate starting geometry"). Calibration entry only: outside
  the 0-1 GPa scope, pressure set by ``Origin.pressure_override_GPa``.
* **MgB2** P6/mmm a = 3.086, c = 3.524 Å (as in ``src/scoring/alignn_screen.py``; textbook values,
  UNVERIFIED beyond ~1 %). Flag ``BLOCKED_PH_CRASH``: ph.x SIGSEGV (W4) until diagnosed.
* **Nb** bcc a = 3.30, **Pb** fcc a = 4.95, **Al** fcc a = 4.05 Å (same source, same caveat).

T*_s: 80 K (Mg2IrH6), 150 K (H3S), 20 K (MgB2), ½·Tc_exp for the elements (Nb 9.25 K, Pb 7.20 K,
Al 1.18 K; textbook experimental Tc, UNVERIFIED here).

The 6-8-atom light-element anchor control "from lane D" (§3.3) is not chosen yet: it needs an
Alex-MP query plus a person's choice; :func:`anchor_from_lane_d` builds it once a formula/entry id
is given.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any, Callable

from ..records import Origin, RawStructure, StructureData


@dataclass(frozen=True)
class ControlSpec:
    name: str
    T_star_s_K: float
    pressure_GPa: float
    build: Callable[[], Any]                  # -> pymatgen Structure
    why: str
    flags: tuple[str, ...] = ()
    source: str = ""

    def raw(self) -> RawStructure:
        st = self.build()
        sd = StructureData(lattice=[[float(x) for x in r] for r in st.lattice.matrix],
                           species=[site.specie.symbol for site in st],
                           frac_coords=[[float(x) % 1.0 for x in site.frac_coords] for site in st])
        origin = Origin(lane="C", control=self.name, T_star_surrogate_K=float(self.T_star_s_K),
                        generator=None, flags=["KNOWN_COMPOUND", *self.flags],
                        pressure_override_GPa=(self.pressure_GPa if self.pressure_GPa != 0.0 else None))
        return RawStructure(structure=sd, origin=origin, pressure_GPa=float(self.pressure_GPa),
                            source_ref=f"control:{self.name}", meta={"why": self.why, "source": self.source})


def _sg(sg: str, lat, species, coords):
    from pymatgen.core import Lattice, Structure
    return Structure.from_spacegroup(sg, lat, species, coords)


def _mg2irh6():
    from pymatgen.core import Lattice
    return _sg("Fm-3m", Lattice.cubic(6.632576), ["Mg", "Ir", "H"],
               [[0.25, 0.25, 0.25], [0, 0, 0], [0.260649, 0, 0]])


def _h3s():
    from pymatgen.core import Lattice
    return _sg("Im-3m", Lattice.cubic(3.09), ["S", "H"], [[0, 0, 0], [0.5, 0, 0]])


def _mgb2():
    from pymatgen.core import Lattice
    return _sg("P6/mmm", Lattice.hexagonal(3.086, 3.524), ["Mg", "B"], [[0, 0, 0], [1 / 3, 2 / 3, 0.5]])


def _fcc(a, el):
    from pymatgen.core import Lattice
    return lambda: _sg("Fm-3m", Lattice.cubic(a), [el], [[0, 0, 0]])


def _bcc(a, el):
    from pymatgen.core import Lattice
    return lambda: _sg("Im-3m", Lattice.cubic(a), [el], [[0, 0, 0]])


CONTROLS: list[ControlSpec] = [
    ControlSpec("Mg2IrH6", 80.0, 0.0, _mg2irh6,
                "best predicted ambient hydride (Tc 160-175 K, Dolui et al. 2024); MatterSim phonons fail on it; "
                "the critical S3 control and the Gamma diagnostic",
                ("CONTROL_S3_CRITICAL",), "Alex-MP msft-2023062809350933252138247157-GGA"),
    ControlSpec("H3S", 150.0, 200.0, _h3s,
                "only H anchor with full DFPT on disk; calibration entry only",
                ("CALIBRATION_ONLY", "OUT_OF_SCOPE_PRESSURE"), "earlier benchmark H3S_Im-3m.vasp"),
    ControlSpec("MgB2", 20.0, 0.0, _mgb2, "multiband; G-VALID item 5 and N1-MB",
                ("BLOCKED_PH_CRASH", "MULTIBAND"), "textbook lattice"),
    ControlSpec("Nb", 9.25 / 2, 0.0, _bcc(3.30, "Nb"), "calibration coarse vs full (V2/V3)", (), "textbook lattice"),
    ControlSpec("Pb", 7.20 / 2, 0.0, _fcc(4.95, "Pb"), "calibration coarse vs full (V2/V3)", (), "textbook lattice"),
    ControlSpec("Al", 1.18 / 2, 0.0, _fcc(4.05, "Al"), "calibration coarse vs full (V2/V3)", (), "textbook lattice"),
]

CONTROL_NAMES = tuple(c.name for c in CONTROLS)


def control(name: str) -> ControlSpec:
    for c in CONTROLS:
        if c.name == name:
            return c
    raise KeyError(name)


def control_raws(names: tuple[str, ...] | None = None) -> list[RawStructure]:
    return [c.raw() for c in CONTROLS if names is None or c.name in names]


def anchor_from_lane_d(raw_d: RawStructure, T_c_exp_K: float, name: str | None = None) -> ControlSpec:
    """Promote a lane-D entry to the 6-8-atom anchor control (T*_s = ½·Tc_exp)."""
    from ..structure.canonical import to_pymatgen
    st = to_pymatgen(raw_d.structure.to_json())
    nm = name or f"anchor:{raw_d.meta.get('formula', st.composition.reduced_formula)}"
    return ControlSpec(nm, T_c_exp_K / 2, 0.0, lambda: st, "6-8-atom cost anchor (C2/C3)",
                       ("ANCHOR",), str(raw_d.source_ref))


def inject_controls(ledger: Any, names: tuple[str, ...] | None = None, *, targets: Any = None,
                    tiers: Any = None, context: Any = None) -> list[str]:
    """Insert the controls into the ledger and record their S1 outcome; returns candidate ids.

    Controls go through S1 (prefilter) like everything else, so a control that a screen would
    remove shows up as a recorded false negative instead of being skipped. Idempotent: an
    already-present control (same candidate id) is not inserted twice. ``targets``/``tiers``
    default to the loaded configs (or ``context.targets``/``context.tiers``)."""
    from ..config import load_targets, load_tiers
    from .ingest import ingest_raw
    if context is not None:
        targets, tiers = targets or context.targets, tiers or context.tiers
    targets = targets or load_targets()
    tiers = tiers or load_tiers()
    return [ingest_raw(ledger, rs, targets=targets, tiers=tiers) for rs in control_raws(names)]


__all__ = ["ControlSpec", "CONTROLS", "CONTROL_NAMES", "control", "control_raws", "inject_controls",
           "anchor_from_lane_d"]
