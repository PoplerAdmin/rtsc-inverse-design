"""Python <-> Lean correspondence for decideT (architecture §1.2, Decide.lean v2).

1. ``core_b`` over all Bool rows equals the Lean semantics of ``decideCore`` + ``applyPolicy``
   (transcribed independently below) and, when present, the Lean export
   ``lean/logs/decideB_table.txt``.
2. Drift: every Lean text ``rtt.decide.core`` transcribes still occurs verbatim in Decide.lean.
3. When ``lean/export/decide_cases.json`` exists, ``decide_t`` is run on every exported case
   (``test_decide_lean_cases.py``).
"""
from __future__ import annotations

import itertools
import re
from pathlib import Path

import pytest

from rtt.decide import core as C

ROOT = Path(__file__).resolve().parents[1]
DECIDE_LEAN = ROOT / "lean" / "TargetLogic" / "Decide.lean"
TABLE = ROOT / "lean" / "logs" / "decideB_table.txt"
NAMES = {"rejectedProven", "screenedOut", "undecided", "pass"}


def lean_semantics(gates_ok, first, would, b):
    if gates_ok:
        core = "rejectedProven" if first else ("screenedOut" if would else "pass")
    else:
        core = "undecided"
    rank = {"rejectedProven": 2, "screenedOut": 1, "undecided": 0, "pass": 0}[core]
    return "undecided" if (b and rank >= 1) else core


ROWS = [r for r in itertools.product([False, True], repeat=4) if not (r[1] and not r[2])]   # first => would


@pytest.mark.parametrize("row", ROWS)
def test_core_b_truth_table(row):
    assert C.core_b(*row) == lean_semantics(*row)


def test_table_shape():
    tab = {r: C.core_b(*r) for r in ROWS}
    assert len(tab) == 12 and set(tab.values()) == NAMES
    for (g, first, would, b), o in tab.items():
        if o == "rejectedProven":
            assert g and first and not b
        if o == "screenedOut":
            assert g and would and not first and not b


@pytest.mark.skipif(not TABLE.exists(), reason="lean/logs/decideB_table.txt not exported")
def test_core_b_matches_lean_export():
    n = 0
    for line in TABLE.read_text(encoding="utf-8").splitlines():
        bools = re.findall(r"\b(true|false)\b", line)
        outs = [w for w in re.findall(r"[A-Za-z]+", line) if w in NAMES]
        if len(bools) == 4 and len(outs) == 1:
            assert C.core_b(*[b == "true" for b in bools]) == outs[0], line
            n += 1
    assert n > 0


def _norm(s: str) -> str:
    return re.sub(r"\s+", " ", s).strip()


@pytest.mark.parametrize("i", range(len(C.LEAN_TEXTS)))
def test_no_drift_against_decide_lean(i):
    lean = _norm(DECIDE_LEAN.read_text(encoding="utf-8"))
    assert _norm(C.LEAN_TEXTS[i]) in lean, "Decide.lean changed: re-check rtt.decide.core against it"
