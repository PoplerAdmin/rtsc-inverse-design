"""Evidence statement per finalist (architecture §3.6, S9; R12 #10).

Every item of the §3.6 template is rendered as either

    checked (<value>; <record id>)        or        NOT CHECKED (<what is missing>)

so that an absent computation can never read as a passed check. :func:`gather` assembles an
:class:`EvidenceInputs` from the ledger (``rtt.ledger.store.Ledger``); :func:`render` is pure
and works on hand-built inputs too (tests, other stores). :func:`evidence_statement` does both.

Metric keys read from stage records (``tier_runs.metrics_json``), by stage:

* S2: ``novel`` (bool), ``dup_of``; S4: ``gap_eV``, ``moment_muB``;
* S3/S5: ``gamma_min_freq_cm1`` (Gamma dynamical-matrix diagnostic, optional);
* S6: ``min_freq_cm1`` (full q grid of the final structure), ``A_mid_K``, ``A_lo_K``;
* S7: ``e_hull_zpe_eV``, ``hull_ref_profile``, ``synthesis_route``, ``anharmonic_method``,
  ``anharmonic_shift``, ``alpha2F_recomputed``, ``md_method``, ``md_persistent``, ``md_T_K``.

G-VALID items and ENVELOPE inputs come from the ``decide_inputs`` row of the latest decision.
"""
from __future__ import annotations

import json
from dataclasses import dataclass, field
from typing import Any

NC = "NOT CHECKED"
TEMPLATE_ITEMS = (
    "structure", "origin lane", "novelty", "H1 metallic (S4)", "H2 non-magnetic (S4)",
    "harmonic stability, full q grid of the final structure (S6)", "Gamma diagnostic",
    "G-VALID 1 non-adiabatic", "G-VALID 2 stable/anharmonic", "G-VALID 3 Coulomb sign",
    "G-VALID 4 convergence data", "G-VALID 5 multiband data", "ENVELOPE inputs",
    "decide outcome and firing rule", "core verdict", "ZPE-corrected E_hull (H7)",
    "anharmonic check", "MD/PIMD persistence (H6)", "mu* sensitivity",
    "multiband / k-resolved status", "provisional margins used upstream",
    "audit and control recall")


@dataclass
class EvidenceInputs:
    candidate_id: str
    formula: str | None = None
    structure: dict | None = None
    structure_sha256: str | None = None
    origin: dict | None = None
    pressure_GPa: float | None = None
    stages: dict = field(default_factory=dict)        # stage -> {outcome, tier_run_id, metrics, criterion, provisional, tags, threshold, source_id}
    decisions: list = field(default_factory=list)     # [{id, stage, tier, outcome, firing_rule, decide_input_id, assumptions}]
    decide_inputs: dict | None = None                 # DecideInputs.to_json() of the latest decision
    verdicts: list = field(default_factory=list)      # ledger verdict rows (dicts)
    audits: list = field(default_factory=list)        # [{audit, inclusion_prob}]
    control_recall: dict | None = None                # from rtt.report.summary.control_recall
    ranking: dict | None = None                       # {decide_outcome, T_low_K, A_mid_K, A_lo_K}


def _fmt(v: Any) -> str:
    if v is None:
        return "None"
    if isinstance(v, float):
        return f"{v:.4g}"
    return str(v)


def _ok(value: str, rec: str | None) -> str:
    return f"checked ({value}{'; ' + rec if rec else ''})"


def _nc(why: str = "") -> str:
    return f"{NC}" + (f" ({why})" if why else "")


def _stage(ev: EvidenceInputs, s: str) -> dict:
    return ev.stages.get(s) or {}


def _rid(st: dict, stage: str) -> str | None:
    if not st:
        return None
    tr = st.get("tier_run_id")
    oid = st.get("outcome_id")
    parts = [f"{stage} outcome {st.get('outcome')}"]
    if oid is not None:
        parts.append(f"outcome #{oid}")
    if tr is not None:
        parts.append(f"tier_run #{tr}")
    return ", ".join(parts)


def _metric(ev: EvidenceInputs, stage: str, key: str) -> Any:
    return (_stage(ev, stage).get("metrics") or {}).get(key)


def _fields(d: dict | None, names: list[str]) -> tuple[dict, list[str]]:
    d = d or {}
    got = {n: d.get(n) for n in names}
    return got, [n for n, v in got.items() if v is None]


def items(ev: EvidenceInputs) -> list[tuple[str, str]]:
    """(item, text) for every template item, in template order."""
    out: list[tuple[str, str]] = []
    st = ev.structure or {}
    if ev.structure_sha256:
        out.append(("structure", _ok(f"{ev.formula}, spacegroup {st.get('spacegroup')}, "
                                     f"{st.get('n_atoms')} atoms, P = {_fmt(ev.pressure_GPa)} GPa",
                                     f"structure sha256 {ev.structure_sha256[:16]}")))
    else:
        out.append(("structure", _nc("no structure record")))
    o = ev.origin or {}
    out.append(("origin lane", _ok(f"lane {o.get('lane')}" + (f", control {o['control']}"
                                                               if o.get("control") else ""),
                                   None) if o.get("lane") else _nc()))
    s2 = _stage(ev, "S2")
    if s2:
        nov = (s2.get("metrics") or {}).get("novel")
        out.append(("novelty", _ok(f"S2 {s2.get('outcome')}, novel={_fmt(nov)}", _rid(s2, "S2"))
                    if nov is not None else _nc(f"S2 {s2.get('outcome')} but no novelty metric")))
    else:
        out.append(("novelty", _nc("S2 not reached")))
    s4 = _stage(ev, "S4")
    m4 = s4.get("metrics") or {}
    tags4 = s4.get("tags") or []
    crit4 = s4.get("criterion") or ""
    for keys, item, crit in ((("N_EF_per_eV_atom", "gap_eV"), "H1 metallic (S4)", "H1"),
                             (("moment_per_atom_muB", "abs_moment_muB", "moment_muB",
                               "magnetic_energy_meV_atom"), "H2 non-magnetic (S4)", "H2")):
        vals = {k: m4.get(k) for k in keys if m4.get(k) is not None}
        if not s4:
            out.append((item, _nc("S4 not run")))
        elif f"{crit}_NOT_CHECKED" in tags4 or not vals:
            out.append((item, _nc("threshold or data missing, or filter disabled")))
        elif crit4.startswith(f"{crit}-"):          # DEFERRED(H1-marginal), DEFERRED(H2-...)
            out.append((item, _nc(f"S4 {s4.get('outcome')}({crit4}): "
                                  + ", ".join(f"{k} = {_fmt(v)}" for k, v in vals.items()))))
        else:
            out.append((item, _ok(", ".join(f"{k} = {_fmt(v)}" for k, v in vals.items()),
                                  _rid(s4, "S4"))))
    s6 = _stage(ev, "S6")
    mf = (s6.get("metrics") or {}).get("min_freq_cm1")
    out.append(("harmonic stability, full q grid of the final structure (S6)",
                _ok(f"min frequency {_fmt(mf)} cm^-1", _rid(s6, "S6")) if mf is not None
                else _nc("no S6 phonons" if not s6 else "S6 min_freq_cm1 missing")))
    g = None
    for s in ("S6", "S5", "S3"):
        g = _metric(ev, s, "gamma_min_freq_cm1")
        if g is not None:
            out.append(("Gamma diagnostic", _ok(f"min Gamma frequency {_fmt(g)} cm^-1",
                                                _rid(_stage(ev, s), s))))
            break
    if g is None:
        out.append(("Gamma diagnostic", _nc()))
    di = ev.decide_inputs or {}
    val = di.get("validity")
    env = di.get("envelope")
    dsrc = f"decide_input #{ev.decisions[-1].get('decide_input_id')}" if ev.decisions else None
    gv = [
        ("G-VALID 1 non-adiabatic", ["EF", "omegaMax"], env and ["lamNom"]),
        ("G-VALID 2 stable/anharmonic", ["maxImagModeK", "containsHLiBe", "anharmonicCheck"], None),
        ("G-VALID 3 Coulomb sign", ["muStarNegativeSuspected"], None),
        ("G-VALID 4 convergence data", ["admittedSmearings", "qMeshes"], None),
        ("G-VALID 5 multiband data", ["nSheets", "bandResolvedAvailable", "dosFractionsKnown"], None),
    ]
    for name, need, _ in gv:
        if val is None:
            out.append((name, _nc("no decide inputs")))
            continue
        got, miss = _fields(val, need)
        if name.startswith("G-VALID 2"):
            got["anharmonicShift"] = val.get("anharmonicShift")
            got["alpha2FRecomputedAnharmonic"] = val.get("alpha2FRecomputedAnharmonic")
        if name.startswith("G-VALID 4") and env:
            for k in ("dA", "dLam", "dM"):
                got[k] = env.get(k)
                if env.get(k) is None:
                    miss.append(k)
        if name.startswith("G-VALID 1") and env:
            got["lamNom"] = env.get("lamNom")
        txt = ", ".join(f"{k}={_fmt(v)}" for k, v in got.items())
        out.append((name, _nc("missing: " + ", ".join(miss) + "; G-VALID fails")
                    if miss else _ok(txt, dsrc)))
    if env:
        got, miss = _fields(env, ["Anom", "lamNom", "mNom", "dA", "dLam", "dM", "epsSys", "epsOmega"])
        txt = ", ".join(f"{k}={_fmt(v)}" for k, v in got.items()) + f", tier={env.get('tier')}"
        out.append(("ENVELOPE inputs", _nc("missing: " + ", ".join(miss)) if miss else _ok(txt, dsrc)))
    else:
        out.append(("ENVELOPE inputs", _nc("no decide inputs")))
    if ev.decisions:
        d = ev.decisions[-1]
        out.append(("decide outcome and firing rule",
                    _ok(f"{d.get('outcome')} at tier {d.get('tier')} (stage {d.get('stage')}), "
                        f"rule {d.get('firing_rule') or 'none'}"
                        + (f", assumptions {d.get('assumptions')}" if d.get("assumptions") else ""),
                        f"decision #{d.get('id')}")))
    else:
        out.append(("decide outcome and firing rule", _nc("no decision")))
    vs = [v for v in ev.verdicts if v.get("spectrum_tier") in ("full", "converged")]
    if vs:
        v = vs[-1]
        out.append(("core verdict",
                    _ok(f"{v.get('status')} [{_fmt(v.get('T_low_K'))}, {_fmt(v.get('T_up_K'))}] K "
                        f"at mu*_hi={_fmt(v.get('mu_hi'))}, tier {v.get('spectrum_tier')}, rules "
                        f"{_j(v.get('rule_ids'))}, assumptions {_j(v.get('assumptions'))}, flags "
                        f"{_j(v.get('flags'))}, recheck {'ok' if v.get('recheck_ok') else 'NOT ok'}",
                        f"verdict #{v.get('id')}")))
    else:
        coarse = [v for v in ev.verdicts if v.get("spectrum_tier") == "coarse"]
        out.append(("core verdict", _nc("no full/converged verdict"
                                            + (f"; {len(coarse)} coarse verdict(s) stored as "
                                               "evidence only" if coarse else ""))))
    eh = _metric(ev, "S7", "e_hull_zpe_eV")
    if eh is not None:
        route = _metric(ev, "S7", "synthesis_route")
        txt = f"E_hull+ZPE = {_fmt(eh)} eV/atom vs own-DFT refs ({_metric(ev, 'S7', 'hull_ref_profile')})"
        if eh > 0.1:
            txt += f"; synthesis route: {route if route else 'NOT GIVEN (required above 0.1 eV/atom)'}"
        out.append(("ZPE-corrected E_hull (H7)", _ok(txt, _rid(_stage(ev, "S7"), "S7"))))
    else:
        out.append(("ZPE-corrected E_hull (H7)", _nc("S7a not run")))
    am = _metric(ev, "S7", "anharmonic_method")
    out.append(("anharmonic check",
                _ok(f"method {am}, shift {_fmt(_metric(ev, 'S7', 'anharmonic_shift'))}, alpha2F "
                    f"recomputed: {_fmt(_metric(ev, 'S7', 'alpha2F_recomputed'))}",
                    _rid(_stage(ev, "S7"), "S7")) if am else _nc()))
    md = _metric(ev, "S7", "md_method")
    out.append(("MD/PIMD persistence (H6)",
                _ok(f"{md} at {_fmt(_metric(ev, 'S7', 'md_T_K'))} K, persistent: "
                    f"{_fmt(_metric(ev, 'S7', 'md_persistent'))}", _rid(_stage(ev, "S7"), "S7"))
                if md else _nc()))
    mus = {}
    for v in vs:
        mus.setdefault(round(float(v.get("mu_hi") or 0), 3), v)
    want = (0.10, 0.13, 0.15)
    have = [m for m in want if m in mus]
    if have:
        txt = "; ".join(f"mu*_hi={m}: {mus[m].get('status')} [{_fmt(mus[m].get('T_low_K'))}, "
                        f"{_fmt(mus[m].get('T_up_K'))}] K (verdict #{mus[m].get('id')})" for m in have)
        missing = [m for m in want if m not in mus]
        out.append(("mu* sensitivity", _ok(txt + (f"; NOT CHECKED at {missing}" if missing else ""),
                                           None)))
    else:
        out.append(("mu* sensitivity", _nc("no full/converged verdicts")))
    if val is not None and val.get("nSheets") is not None:
        out.append(("multiband / k-resolved status",
                    _ok(f"n_sheets={val.get('nSheets')}, band-resolved data "
                        f"{_fmt(val.get('bandResolvedAvailable'))}, k-resolved check "
                        f"{_fmt(val.get('kResolvedCheckDone'))}", dsrc)))
    else:
        ns = _metric(ev, "S4", "n_sheets")
        out.append(("multiband / k-resolved status",
                    _ok(f"n_sheets={ns} (S4 only; no band-/k-resolved check)", _rid(s4, "S4"))
                    if ns is not None else _nc()))
    prov = []
    for s, st2 in sorted(ev.stages.items()):
        if st2.get("provisional") or "PROVISIONAL_MARGIN" in (st2.get("tags") or []):
            prov.append(f"{s}: {st2.get('outcome')} "
                        f"({st2.get('criterion') or ','.join(t for t in st2.get('tags') or [] if t.startswith('PROVISIONAL'))})")
    out.append(("provisional margins used upstream",
                _ok("; ".join(prov), None) if prov else "checked (none used)"))
    au = "; ".join(f"{a.get('audit')} p={_fmt(a.get('inclusion_prob'))}" for a in ev.audits)
    cr = ev.control_recall
    if cr is None and not ev.audits:
        out.append(("audit and control recall", _nc("no audit or control data")))
    else:
        parts = []
        if au:
            parts.append(f"audits of this candidate: {au}")
        if cr is not None:
            parts.append("control recall: " + (", ".join(
                f"{k} {v['passed']}/{v['total']}" for k, v in sorted(cr.items())) or "no controls"))
        out.append(("audit and control recall", _ok("; ".join(parts), None)))
    return out


def _j(x: Any) -> Any:
    if isinstance(x, str):
        try:
            return json.loads(x)
        except ValueError:
            return x
    return x


def render(ev: EvidenceInputs) -> str:
    lines = [f"# Evidence statement: {ev.formula or '?'} ({ev.candidate_id})", ""]
    if ev.ranking:
        r = ev.ranking
        lines += [f"Ranking keys (§3.6): decide outcome **{r.get('decide_outcome')}**, certified "
                  f"T_low = {_fmt(r.get('T_low_K'))} K, A_mid = {_fmt(r.get('A_mid_K'))} K "
                  f"(A_lo = {_fmt(r.get('A_lo_K'))} K). T_up and A_hi are ceilings only.", ""]
    s5 = ev.stages.get("S5")
    if s5:      # campaign 1: the coarse number is what most candidates end with; show it (not a ranking key)
        m5 = s5.get("metrics") or {}
        lines += [f"Coarse DFPT (S5, q 3³ + 4³, not a ranking key): A_coarse,mid = {_fmt(m5.get('A_mid_K'))} K, "
                  f"A_coarse,lo = {_fmt(m5.get('A_lo_K'))} K; S5 outcome {s5.get('outcome')}"
                  + (f" ({s5.get('criterion')})" if s5.get("criterion") else "")
                  + f"; tier_run #{s5.get('tier_run_id')}.", ""]
    its = items(ev)
    n_nc = sum(1 for _, t in its if t.startswith(NC))
    for name, text in its:
        lines.append(f"- **{name}**: {text}")
    lines += ["", f"{len(its) - n_nc} of {len(its)} items checked; {n_nc} NOT CHECKED."]
    return "\n".join(lines) + "\n"


def gather(ledger: Any, cid: str) -> EvidenceInputs:
    """Assemble the inputs from the ledger (read-only)."""
    rec = ledger.candidate_record(cid)
    ev = EvidenceInputs(candidate_id=cid, formula=rec.formula,
                        structure=rec.structure.to_json(), structure_sha256=rec.structure_sha256,
                        origin=rec.origin.to_json(), pressure_GPa=rec.pressure_GPa)
    for stage, o in ledger.current_outcomes(cid).items():
        e = rec.stages.get(stage)
        met = dict((e.metrics if e else {}) or {})
        if stage == "S2" and isinstance(o.get("value"), dict):
            # campaign 1: S2 stores novelty in the PASS outcome's value (doc 21 §2 item 11),
            # not in tier_run metrics; without this the item always read NOT CHECKED
            met = {**o["value"], **met}
        ev.stages[stage] = {"outcome": o["outcome"], "outcome_id": o.get("id"),
                            "tier_run_id": o.get("tier_run_id"),
                            "metrics": met,
                            "criterion": o.get("criterion"), "provisional": bool(o.get("provisional")),
                            "tags": o.get("tags") or [], "threshold": o.get("threshold"),
                            "source_id": o.get("source_id")}
    rows = ledger.qa("SELECT id, stage, tier, decide_input_id, outcome, firing_rule, assumptions "
                     "FROM decisions WHERE candidate_id=? ORDER BY id", (cid,))
    for r in rows:
        ev.decisions.append({"id": r[0], "stage": r[1], "tier": r[2], "decide_input_id": r[3],
                             "outcome": r[4], "firing_rule": r[5],
                             "assumptions": json.loads(r[6]) if r[6] else []})
    if ev.decisions:
        ev.decide_inputs = ledger.decide_input(ev.decisions[-1]["decide_input_id"]).to_json()
    ev.verdicts = ledger.verdicts_for(cid)
    ev.audits = [{"audit": a, "inclusion_prob": p}
                 for a, p in ledger.qa("SELECT audit, inclusion_prob FROM audits WHERE candidate_id=?",
                                       (cid,))]
    from .summary import control_recall
    ev.control_recall = control_recall(ledger)
    # §3.6: the latest full/converged verdict at mu*_hi = 0.13; its T_low is already the
    # combined bracket's lower end (min over alternatives); missing counts as 0.
    cand = [v for v in ev.verdicts if v.get("spectrum_tier") in ("full", "converged")
            and abs(float(v.get("mu_hi") or 0) - 0.13) < 1e-9]
    t_low = float(cand[-1]["T_low_K"]) if cand and cand[-1].get("T_low_K") is not None else 0.0
    ev.ranking = {"decide_outcome": ev.decisions[-1]["outcome"] if ev.decisions else None,
                  "T_low_K": t_low, "A_mid_K": _metric(ev, "S6", "A_mid_K"),
                  "A_lo_K": _metric(ev, "S6", "A_lo_K")}
    return ev


def evidence_statement(cid: str, ledger: Any) -> str:
    """Markdown evidence statement for one candidate (architecture §5.1 ``evidence_statement``)."""
    return render(gather(ledger, cid))
