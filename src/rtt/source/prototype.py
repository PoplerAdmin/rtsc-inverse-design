"""Deterministic prototype-substitution lane **P** (campaign 2; R30 §3.2; doc 32).

The lane enumerates the K2PtCl6-type (Fm-3m, sg 225) [MH6] complex hydride family around the
Mg2IrH6 positive control (doc 15) by substitution on fixed prototypes. No generator, no
randomness: the same configuration gives the same candidate list, and every candidate records
its prototype and substitution in ``Origin.conditioning`` (provenance, lane ``P``).

Prototypes (all start from an idealised K2PtCl6 geometry; S3 relaxes them with MatterSim):

==========  ==========================  =====  ====================================================
key         formula / cell              atoms  ordering
==========  ==========================  =====  ====================================================
A2MH6       A2 M H6, Fm-3m primitive    9      parent (K2PtCl6 / Mg2FeH6 type)
AA'MH6      A A' M H6, F-43m primitive  9      the A simple-cubic sublattice split into its two fcc
                                               sublattices = rock-salt order of A/A' (8c -> 4c+4d)
A3A'M2H12   A3 A' M2 H12, doubled cell  18     one A of four replaced in the tetragonal
                                               (a/√2, a/√2, a) two-complex cell (1:3 alloy)
A2MM'H12    A4 M M' H12, doubled cell   18     M/M' alternate on (001) planes of the fcc M lattice
                                               (L1_0-type order) in the same two-complex cell
==========  ==========================  =====  ====================================================

Electron-count rule (derived in doc 32 §1). Count H as hydride H⁻ (2-electron σ donor) and the
A cations at their fixed group valence. Per [MH6] complex

    N = g(M) + 6 + q,      q = Σ ox(A) / n(M)       (complex charge [MH6]^q-)

with g(M) the group number (d^g s^0 count of the neutral metal). For an octahedral low-spin
complex the 6 M-H σ bonding orbitals (a1g + t1u + eg) and the non-bonding t2g hold 18 electrons:
N = 18 is the closed shell (Mg2FeH6, K2PtH6, SrPtH6/BaPtH6 in campaign 1: all insulators).
N - 18 > 0 electrons go into the M-H σ* (eg) band (capacity 4 per complex), N - 18 < 0 are
holes in t2g. Mg2IrH6 has N = 9 + 6 + 4 = 19: one electron per complex in the eg σ* band
(band-level check in doc 32 §1.3). Odd total valence count per primitive cell
(``valence_parity``) forces a non-magnetic metal (Kramers pairs cannot fill an odd count).

Targets (``DEFAULT_N_TARGETS``): N ∈ {17, 17.5, 18.5, 19, 19.5, 20}; N = 18 (closed shell) and
N ≤ 16.5 / ≥ 20.5 are not enumerated. Priority = |N - 19| (Mg2IrH6 filling first), then cell
size. The rule is a *ranking/enumeration* rule only (contract R2): it decides what the lane
proposes, it never removes a candidate from another lane.
"""
from __future__ import annotations

import itertools
from dataclasses import dataclass, field
from fractions import Fraction
from typing import Any, Iterable, Sequence

from ..records import Origin, RawStructure, StructureData

LANE = "P"
GENERATOR = "rtt.source.prototype/1"

#: A-site cations and their fixed formal valence (group valence)
A_OX: dict[str, int] = {"Li": 1, "Na": 1, "K": 1, "Rb": 1, "Cs": 1,
                        "Mg": 2, "Ca": 2, "Sr": 2, "Ba": 2,
                        "Sc": 3, "Y": 3, "La": 3}
#: M-site late transition metals: group number. Tc is excluded (radioactive; also outside the
#: MatterGen element domain), see :data:`EXCLUDED_M`.
M_GROUP: dict[str, int] = {"Mn": 7, "Re": 7,
                           "Fe": 8, "Ru": 8, "Os": 8,
                           "Co": 9, "Rh": 9, "Ir": 9,
                           "Ni": 10, "Pd": 10, "Pt": 10}
EXCLUDED_M = {"Tc": "radioactive"}
M_3D = frozenset({"Mn", "Fe", "Co", "Ni"})

#: per-complex electron counts the lane enumerates (doc 32 §1.5)
DEFAULT_N_TARGETS = (Fraction(17), Fraction(35, 2), Fraction(37, 2), Fraction(19), Fraction(39, 2),
                     Fraction(20))
N_REFERENCE = Fraction(19)            # Mg2IrH6
N_CLOSED = 18

#: Starting-geometry radii (Å): Shannon CN8 ionic radii, recalled from memory (UNVERIFIED to the
#: 0.01 Å level; only the start of the MatterSim relax depends on them). The effective H⁻ radius
#: 1.46 Å makes d(Mg-H) = 0.89 + 1.46 = 2.35 Å reproduce a(Mg2IrH6) = 6.65 Å (doc 15).
R_A_CN8 = {"Li": 0.92, "Na": 1.18, "K": 1.51, "Rb": 1.61, "Cs": 1.74, "Mg": 0.89, "Ca": 1.12,
           "Sr": 1.26, "Ba": 1.42, "Sc": 0.87, "Y": 1.019, "La": 1.16}
R_H_EFF = 1.46
D_MH = {"3d": 1.58, "4d5d": 1.70}     # M-H bond (Ir-H 1.72 Å in PBE Mg2IrH6, doc 15)
SQRT_1_8 = (1.0 / 8.0) ** 0.5         # d(A-H)/a at x_H = 1/4

PROTOTYPES = ("A2MH6", "AA'MH6", "A3A'M2H12", "A2MM'H12")


# ----------------------------------------------------------------------------- electron count
def complex_count(a_sites: Sequence[str], m_sites: Sequence[str], h_per_m: int = 6) -> Fraction:
    """N per [MH_n] complex = mean g(M) + n_H + Σ ox(A)/n(M) (doc 32 §1)."""
    n_m = len(m_sites)
    if n_m == 0:
        raise ValueError("no M site")
    g = Fraction(sum(M_GROUP[m] for m in m_sites), n_m)
    q = Fraction(sum(A_OX[a] for a in a_sites), n_m)
    return g + h_per_m + q


def cell_valence(a_sites: Sequence[str], m_sites: Sequence[str], n_h: int) -> int:
    """Total valence electrons of the cell in the same counting (A: group valence, M: group
    number, H: 1). Equal to n(M)·N."""
    return sum(A_OX[a] for a in a_sites) + sum(M_GROUP[m] for m in m_sites) + n_h


def filling(N: Fraction) -> dict[str, Any]:
    """Band-filling labels for a per-complex count N (octahedral low-spin picture)."""
    x = N - N_CLOSED
    d = {"N_per_complex": float(N), "delta_vs_18": float(x), "delta_vs_Mg2IrH6": float(N - N_REFERENCE)}
    if x == 0:
        d["class"] = "closed_shell_18e"
    elif x > 0:
        d["class"] = "eg_sigma_star_electrons"
        d["eg_filling_fraction"] = float(min(x, 4) / 4)
    else:
        d["class"] = "t2g_holes"
        d["t2g_holes"] = float(-x)
    return d


# ----------------------------------------------------------------------------- geometry
def lattice_guess(a_sites: Sequence[str], m_sites: Sequence[str]) -> tuple[float, float]:
    """(a_conv Å, x_H) start values: a from the mean A-H contact at x = 1/4, x from d(M-H)."""
    r = sum(R_A_CN8[a] for a in a_sites) / len(a_sites)
    a = (r + R_H_EFF) / SQRT_1_8
    dmh = sum(D_MH["3d" if m in M_3D else "4d5d"] for m in m_sites) / len(m_sites)
    return round(a, 4), round(dmh / a, 5)


def _conv_k2ptcl6(A: str, A2: str | None, M: str, a: float, x: float):
    """Conventional (36-atom) cubic cell; A2 (if given) on the 4d half of the 8c sites."""
    from pymatgen.core import Lattice, Structure
    fcc = [(0, 0, 0), (0, .5, .5), (.5, 0, .5), (.5, .5, 0)]
    sp, fc = [], []
    for t in fcc:
        sp.append(M)
        fc.append(t)
        for s in (0.25, 0.75):
            sp.append(A if (A2 is None or s == 0.25) else A2)
            fc.append(tuple((ti + s) % 1 for ti in t))
        for i in range(3):
            for sgn in (1, -1):
                v = [0.0, 0.0, 0.0]
                v[i] = sgn * x
                sp.append("H")
                fc.append(tuple((ti + vi) % 1 for ti, vi in zip(t, v)))
    return Structure(Lattice.cubic(a), sp, fc)


def build_a2mh6(A: str, M: str, a: float | None = None, x: float | None = None):
    """Fm-3m A2MH6 primitive cell (9 atoms)."""
    a0, x0 = lattice_guess([A, A], [M])
    st = _conv_k2ptcl6(A, None, M, a or a0, x or x0)
    return st.get_primitive_structure()


def build_aamh6(A: str, A2: str, M: str, a: float | None = None, x: float | None = None):
    """F-43m AA'MH6 primitive cell (9 atoms): rock-salt order of the A simple-cubic sublattice."""
    a0, x0 = lattice_guess([A, A2], [M])
    st = _conv_k2ptcl6(A, A2, M, a or a0, x or x0)
    return st.get_primitive_structure()


#: two-complex (a/√2, a/√2, a) tetragonal cell in units of the fcc primitive vectors
#: p1 = (0,½,½), p2 = (½,0,½), p3 = (½,½,0): t1 = p2 - p1, t2 = p3, t3 = p1 + p2 - p3; det = 2.
DOUBLE_MATRIX = ((-1, 1, 0), (0, 0, 1), (1, 1, -1))


def _doubled(A: str, M: str, a: float, x: float):
    from pymatgen.core import Lattice, Structure
    import numpy as np
    conv = _conv_k2ptcl6(A, None, M, a, x)
    P = np.array([[0, .5, .5], [.5, 0, .5], [.5, .5, 0]]) * a
    T = np.array(DOUBLE_MATRIX, float) @ P
    if np.linalg.det(T) < 0:
        T[0] = -T[0]
    inv = np.linalg.inv(T)
    sp, fc, seen = [], [], set()
    for s in conv:
        for sh in itertools.product((-1, 0, 1), repeat=3):
            cart = s.coords + np.array(sh) * a
            f = cart @ inv
            f = np.round(f % 1.0, 6) % 1.0
            key = (s.specie.symbol, *np.round(f, 4))
            if key in seen:
                continue
            seen.add(key)
            sp.append(s.specie.symbol)
            fc.append(f)
    st = Structure(Lattice(T), sp, fc)
    if len(st) != 18:
        raise RuntimeError(f"doubled cell has {len(st)} atoms, expected 18")
    return st


def build_a3a_m2h12(A: str, A2: str, M: str, a: float | None = None, x: float | None = None):
    """A3A'(MH6)2, 18 atoms: one of the four A sites of the two-complex cell is A'."""
    a0, x0 = lattice_guess([A, A, A, A2], [M, M])
    st = _doubled(A, M, a or a0, x or x0)
    idx = sorted((i for i, s in enumerate(st) if s.specie.symbol == A),
                 key=lambda i: tuple(round(v, 4) for v in st[i].frac_coords))[0]
    st.replace(idx, A2)
    return st


def build_a2mm_h12(A: str, M: str, M2: str, a: float | None = None, x: float | None = None):
    """A4 M M' H12, 18 atoms: M/M' alternate on (001) planes (L1_0-type order of the fcc M lattice)."""
    a0, x0 = lattice_guess([A] * 4, [M, M2])
    st = _doubled(A, M, a or a0, x or x0)
    ms = sorted((i for i, s in enumerate(st) if s.specie.symbol == M),
                key=lambda i: tuple(round(v, 4) for v in st[i].frac_coords))
    st.replace(ms[1], M2)
    return st


# ----------------------------------------------------------------------------- enumeration
@dataclass(frozen=True)
class ProtoCandidate:
    prototype: str
    a_sites: tuple[str, ...]          # A species of the cell (with multiplicity)
    m_sites: tuple[str, ...]
    n_h: int
    substitution: dict = field(default_factory=dict, compare=False, hash=False)

    @property
    def N(self) -> Fraction:
        return complex_count(self.a_sites, self.m_sites)

    @property
    def cell_valence(self) -> int:
        return cell_valence(self.a_sites, self.m_sites, self.n_h)

    @property
    def n_atoms(self) -> int:
        return len(self.a_sites) + len(self.m_sites) + self.n_h

    @property
    def formula(self) -> str:
        """Reduced formula written A-part, M-part, H (e.g. Mg2IrH6, NaMgPtH6, LiMg3Ir2H12)."""
        import math
        a, m = _count(self.a_sites), _count(self.m_sites)
        g = math.gcd(self.n_h, *a.values(), *m.values())
        part = lambda d: "".join(f"{k}{v // g if v // g > 1 else ''}" for k, v in sorted(d.items()))  # noqa: E731
        return f"{part(a)}{part(m)}H{self.n_h // g}"

    @property
    def parent_N(self) -> Fraction:
        """N of the closest ternary parent (majority A on every A site; for mixed M the parent
        nearer to 19): an alloy whose parent is itself a 19-electron A2MH6 is a small doping of
        a Mg2IrH6 analogue, so it is enumerated first."""
        a_major = max(set(self.a_sites), key=lambda a: (self.a_sites.count(a), a))
        ns = [complex_count((a_major, a_major), (m,)) for m in sorted(set(self.m_sites))]
        if len(set(self.a_sites)) == 2 and self.a_sites.count(self.a_sites[0]) == self.a_sites.count(self.a_sites[1]):
            ns = [complex_count((a, a), (m,)) for a in set(self.a_sites) for m in set(self.m_sites)]
        return min(ns, key=lambda n: (abs(n - N_REFERENCE), n))

    @property
    def distance_to_control(self) -> float:
        """Fraction of A sites that are not Mg plus fraction of M sites that are not Ir."""
        return (sum(a != "Mg" for a in self.a_sites) / len(self.a_sites)
                + sum(m != "Ir" for m in self.m_sites) / len(self.m_sites))

    @property
    def priority(self) -> tuple:
        return (abs(self.N - N_REFERENCE), abs(self.parent_N - N_REFERENCE), round(self.distance_to_control, 3),
                self.n_atoms, self.prototype, self.a_sites, self.m_sites)

    def build(self):
        s = self.substitution
        if self.prototype == "A2MH6":
            return build_a2mh6(s["A"], s["M"])
        if self.prototype == "AA'MH6":
            return build_aamh6(s["A"], s["A'"], s["M"])
        if self.prototype == "A3A'M2H12":
            return build_a3a_m2h12(s["A"], s["A'"], s["M"])
        if self.prototype == "A2MM'H12":
            return build_a2mm_h12(s["A"], s["M"], s["M'"])
        raise ValueError(self.prototype)

    def provenance(self) -> dict:
        f = filling(self.N)
        return {"lane": LANE, "generator": GENERATOR, "prototype": self.prototype,
                "parent": "K2PtCl6-type Fm-3m A2MH6 (Mg2FeH6 / Mg2IrH6)",
                "substitution": dict(self.substitution), "formula": self.formula,
                "N_per_complex": f["N_per_complex"], "filling": f["class"],
                "delta_vs_Mg2IrH6": f["delta_vs_Mg2IrH6"], "parent_N": float(self.parent_N),
                "cell_valence_electrons": self.cell_valence,
                "valence_parity": "odd" if self.cell_valence % 2 else "even",
                "nm_metal_by_parity": bool(self.cell_valence % 2),
                "rule": "doc32#1: N = g(M) + 6 + sum ox(A)/n(M)"}


def _count(xs: Iterable[str]) -> dict[str, int]:
    out: dict[str, int] = {}
    for x in xs:
        out[x] = out.get(x, 0) + 1
    return out


def enumerate_lane_p(*, A: Sequence[str] = tuple(A_OX), M: Sequence[str] = tuple(M_GROUP),
                     targets: Sequence[Fraction] = DEFAULT_N_TARGETS,
                     prototypes: Sequence[str] = PROTOTYPES,
                     alloy_targets: Sequence[Fraction] | None = None) -> list[ProtoCandidate]:
    """All prototype substitutions whose per-complex count N is in ``targets``, sorted by
    :attr:`ProtoCandidate.priority`. ``alloy_targets`` (default: N ∈ [18.5, 19.5]) restricts
    the alloy prototypes, whose role is to tune N around the Mg2IrH6 filling (R30 §2.2(1))."""
    tset = set(Fraction(t) for t in targets)
    aset = set(Fraction(t) for t in (alloy_targets if alloy_targets is not None
                                     else (Fraction(37, 2), Fraction(19), Fraction(39, 2))))
    bad = [m for m in M if m in EXCLUDED_M or m not in M_GROUP]
    if bad:
        raise ValueError(f"M not allowed: {bad}")
    out: list[ProtoCandidate] = []
    A, M = list(A), list(M)
    if "A2MH6" in prototypes:
        for a, m in itertools.product(A, M):
            c = ProtoCandidate("A2MH6", (a, a), (m,), 6, {"A": a, "M": m})
            if c.N in tset:
                out.append(c)
    if "AA'MH6" in prototypes:
        for (a1, a2), m in itertools.product(itertools.combinations(A, 2), M):
            c = ProtoCandidate("AA'MH6", (a1, a2), (m,), 6, {"A": a1, "A'": a2, "M": m})
            if c.N in tset and c.N in aset:
                out.append(c)
    if "A3A'M2H12" in prototypes:
        for a1, a2, m in itertools.product(A, A, M):
            if a1 == a2 or A_OX[a1] == A_OX[a2]:
                continue                       # isovalent 1:3 mixing does not change N
            c = ProtoCandidate("A3A'M2H12", (a1, a1, a1, a2), (m, m), 12, {"A": a1, "A'": a2, "M": m})
            if c.N in tset and c.N in aset and c.N.denominator != 1:
                out.append(c)                  # integer-N cases are covered by A2MH6 / AA'MH6
    if "A2MM'H12" in prototypes:
        for a, (m1, m2) in itertools.product(A, itertools.combinations(M, 2)):
            if M_GROUP[m1] == M_GROUP[m2]:
                continue                       # isovalent M mixing: no change of N
            c = ProtoCandidate("A2MM'H12", (a,) * 4, (m1, m2), 12, {"A": a, "M": m1, "M'": m2})
            if c.N in tset and c.N in aset:
                out.append(c)
    out.sort(key=lambda c: c.priority)
    return out


def select(cands: Sequence[ProtoCandidate], *, max_per_prototype: dict[str, int] | None = None,
           max_total: int | None = None) -> list[ProtoCandidate]:
    """Deterministic cap: keep the first ``max_per_prototype[p]`` of each prototype (priority
    order), then the first ``max_total`` overall."""
    caps = dict(max_per_prototype or {})
    seen: dict[str, int] = {}
    out = []
    for c in sorted(cands, key=lambda c: c.priority):
        k = seen.get(c.prototype, 0)
        if c.prototype in caps and k >= caps[c.prototype]:
            continue
        seen[c.prototype] = k + 1
        out.append(c)
    return out[:max_total] if max_total else out


def to_raw(c: ProtoCandidate, batch_id: str = "P-campaign2") -> RawStructure:
    st = c.build()
    sd = StructureData(lattice=[[float(v) for v in r] for r in st.lattice.matrix],
                       species=[s.specie.symbol for s in st],
                       frac_coords=[[float(v) % 1.0 for v in s.frac_coords] for s in st])
    prov = c.provenance()
    flags = ["PROTOTYPE_SUBSTITUTION", f"PROTOTYPE={c.prototype}", f"N_PER_COMPLEX={float(c.N):g}"]
    if prov["nm_metal_by_parity"]:
        flags.append("ODD_VALENCE_PARITY")
    origin = Origin(lane=LANE, generator=GENERATOR, conditioning=prov, batch_id=batch_id, flags=flags)
    sub = ",".join(f"{k}={v}" for k, v in c.substitution.items())
    return RawStructure(structure=sd, origin=origin, pressure_GPa=0.0,
                        source_ref=f"prototype:{c.prototype}:{sub}",
                        meta={"formula": c.formula, **prov})


def lane_p_raws(cands: Sequence[ProtoCandidate] | None = None, **kw: Any) -> list[RawStructure]:
    return [to_raw(c, **kw) for c in (cands if cands is not None else enumerate_lane_p())]


# ----------------------------------------------------------------------------- config / plan
def load_lane_p_config(path: Any = None):
    """``configs/lane_p.toml`` as a LoadedConfig (sources checked, contract R7)."""
    import tomllib
    from ..config import ConfigError, LoadedConfig, _resolve, sha256_of, validate_sources
    p = _resolve(path, "lane_p.toml")
    text = p.read_text(encoding="utf-8")
    raw = tomllib.loads(text)
    if raw.get("schema") != "rtt.lane_p/1":
        raise ConfigError(f"lane_p: schema must be 'rtt.lane_p/1', got {raw.get('schema')!r}")
    validate_sources(raw)
    return LoadedConfig(kind="lane_p", path=str(p), text=text, sha256=sha256_of(text), raw=raw)


def _frac(v: float) -> Fraction:
    return Fraction(v).limit_denominator(4)


def plan_lane_p(cfg: Any = None) -> list[ProtoCandidate]:
    """The lane-P candidate list of a config (default ``configs/lane_p.toml``), capped per
    prototype in priority order. Deterministic."""
    cfg = cfg or load_lane_p_config()
    cands = enumerate_lane_p(A=cfg.get("sites.A"), M=cfg.get("sites.M"),
                             targets=[_frac(x) for x in cfg.get("electron_count.targets")],
                             alloy_targets=[_frac(x) for x in cfg.get("electron_count.alloy_targets")],
                             prototypes=cfg.get("lane.prototypes"))
    caps = {k: int(v) for k, v in cfg.get("caps").items() if k not in ("source", "source_id")}
    return select(cands, max_per_prototype=caps)


def run_lane_p(ledger: Any, *, targets: Any, tiers: Any, cfg: Any = None,
               cands: Sequence[ProtoCandidate] | None = None) -> list[str]:
    """S1 for lane P: build every planned prototype cell and ingest it (S1 outcome per cell)."""
    from .ingest import ingest_many
    cfg = cfg or load_lane_p_config()
    cands = list(cands) if cands is not None else plan_lane_p(cfg)
    return ingest_many(ledger, lane_p_raws(cands, batch_id=cfg.get("lane.batch_id", "P")),
                       targets=targets, tiers=tiers)


__all__ = ["load_lane_p_config", "plan_lane_p", "run_lane_p","LANE", "A_OX", "M_GROUP", "EXCLUDED_M", "PROTOTYPES", "DEFAULT_N_TARGETS", "ProtoCandidate",
           "complex_count", "cell_valence", "filling", "lattice_guess", "build_a2mh6", "build_aamh6",
           "build_a3a_m2h12", "build_a2mm_h12", "enumerate_lane_p", "select", "to_raw", "lane_p_raws"]
