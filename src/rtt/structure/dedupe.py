"""S2 duplicate and novelty checks (architecture §3.1 S2: ``DUPLICATE``; ``novel`` metric).

* Grouping key: (reduced formula, pressure rounded to 0.1 GPa). Only structures inside a group
  are compared.
* Equality: pymatgen ``StructureMatcher`` (ltol 0.2, stol 0.3, angle_tol 5°, primitive cell,
  scale=True), the same matcher doc 10 §3.2 used for its uniqueness/novelty numbers.
* The structure hash (``structure_sha256``) is a fast exact-match path only.

A false "same" would silently lose a candidate, so a duplicate is declared only on a matcher fit
**and** a bonding-geometry guard (review RC1 F16). The matcher is tight (``scale=False``,
``stol=0.1`` of (V/N)^(1/3), ``ltol=0.2``) because S2 sees *unrelaxed* cells; the guard
(:func:`bonding_compatible`) never merges two cells whose H-H motif class differs (H2-molecular
vs dissociated/network/isolated), whose shortest H-H distances differ by more than 0.2 Å, or
whose shortest (non-H)-H distances differ by more than 0.08 Å (doc 10 §3.4's two Fm-3m Mg2IrH6
minima differ by 0.09 Å in Ir-H). A post-relaxation match is recorded at S3 as a metric
(``relaxed_match_of``), never as a removal;
the duplicate keeps its own ledger row (``duplicates`` table) and the outcome ``DUPLICATE`` names
the representative (first seen).

Novelty against Alex-MP (``is_novel``) is recorded as a metric, never a removal: known compounds
are wanted (database lane D and controls are known by construction).
"""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any, Iterable

from .canonical import structure_sha256, to_pymatgen

MATCHER_KW = dict(ltol=0.2, stol=0.1, angle_tol=5.0, primitive_cell=True, scale=False,
                  attempt_supercell=False)


def _matcher():
    from pymatgen.analysis.structure_matcher import StructureMatcher
    return StructureMatcher(**MATCHER_KW)


MAX_DHH_A = 0.2
MAX_DMH_A = 0.08


def _bond_signature(st) -> tuple:
    from ..mlff.features import features
    f = features({"lattice": st.lattice.matrix.tolist(), "species": [x.specie.symbol for x in st],
                  "frac_coords": st.frac_coords.tolist()})
    return f.get("motif"), f.get("min_HH_A"), f.get("min_MH_A")


def bonding_compatible(a: Any, b: Any) -> bool:
    """False when the two cells are different H bonding situations (never merged)."""
    ma, hha, mha = _bond_signature(to_pymatgen(a))
    mb, hhb, mhb = _bond_signature(to_pymatgen(b))
    if ma != mb:
        return False
    if (hha is None) != (hhb is None) or (hha is not None and abs(hha - hhb) > MAX_DHH_A):
        return False
    if (mha is None) != (mhb is None) or (mha is not None and abs(mha - mhb) > MAX_DMH_A):
        return False
    return True


def _key(st, pressure_GPa: float) -> tuple[str, float]:
    return st.composition.reduced_formula, round(float(pressure_GPa or 0.0), 1)


@dataclass
class DedupeIndex:
    """In-memory index of accepted (non-duplicate) structures."""
    groups: dict = field(default_factory=dict)        # key -> list[(cid, Structure)]
    by_sha: dict = field(default_factory=dict)        # (sha, p) -> cid
    _m: Any = None

    @property
    def matcher(self):
        if self._m is None:
            self._m = _matcher()
        return self._m

    def add(self, cid: str, s: Any, pressure_GPa: float = 0.0, sha: str | None = None) -> None:
        st = to_pymatgen(s)
        self.groups.setdefault(_key(st, pressure_GPa), []).append((cid, st))
        self.by_sha[(sha or structure_sha256(st), round(float(pressure_GPa or 0.0), 1))] = cid

    def find(self, s: Any, pressure_GPa: float = 0.0, sha: str | None = None) -> str | None:
        st = to_pymatgen(s)
        hit = self.by_sha.get((sha or structure_sha256(st), round(float(pressure_GPa or 0.0), 1)))
        if hit is not None:
            return hit
        for cid, other in self.groups.get(_key(st, pressure_GPa), []):
            if self.matcher.fit(st, other) and bonding_compatible(st, other):
                return cid
        return None

    def __len__(self) -> int:
        return sum(len(v) for v in self.groups.values())


def index_from(items: Iterable[tuple[str, Any, float]]) -> DedupeIndex:
    idx = DedupeIndex()
    for cid, s, p in items:
        idx.add(cid, s, p)
    return idx


def _ledger_items(ledger: Any) -> Iterable[tuple[str, Any, float]]:
    """(candidate_id, structure, pressure) of every non-duplicate candidate in the ledger.

    Uses ``ledger.iter_structures()`` when the ledger provides it, else the ``candidates`` /
    ``structures`` / ``duplicates`` tables through ``ledger.conn``."""
    if hasattr(ledger, "iter_structures"):
        yield from ledger.iter_structures()
        return
    import json
    conn = ledger.conn if hasattr(ledger, "conn") else ledger
    dups = {r[0] for r in conn.execute("SELECT candidate_id FROM duplicates")}
    for cid, js, p in conn.execute(
            "SELECT c.id, s.json, c.pressure_GPa FROM candidates c JOIN structures s "
            "ON s.sha256 = c.structure_sha256 ORDER BY c.created, c.id"):
        if cid in dups:
            continue
        d = json.loads(js)
        yield cid, d.get("structure", d), p


def find_duplicate(s: Any, ledger: Any, pressure_GPa: float = 0.0, *, exclude: str | None = None,
                   index: DedupeIndex | None = None) -> str | None:
    """Candidate id of an earlier equal structure at the same pressure, or None.

    Pass a prebuilt ``index`` for batches (building it from the ledger costs one pass)."""
    if index is None:
        index = index_from((c, x, p) for c, x, p in _ledger_items(ledger) if c != exclude)
    hit = index.find(s, pressure_GPa)
    return hit if hit != exclude else None


def is_novel(s: Any, reference_structures: Iterable[Any]) -> tuple[bool, str | None]:
    """(novel, matching reference id). ``reference_structures``: iterable of (id, structure)."""
    st = to_pymatgen(s)
    rf = st.composition.reduced_formula
    m = _matcher()
    for rid, ref in reference_structures:
        r = to_pymatgen(ref)
        if r.composition.reduced_formula == rf and m.fit(st, r) and bonding_compatible(st, r):
            return False, rid
    return True, None
