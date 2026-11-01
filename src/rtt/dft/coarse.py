"""S5 coarse-DFPT metrics and the coarse margin rule (architecture §3.1 S5, §3.4; doc 11 §4).

Metrics (from :func:`rtt.dft.collect.collect_run`): per smearing A, lambda, omega_log, omega_2
(float, informational), which smearings were admitted (mode-positive), the sigma plateau
(admitted sigma >= ``plateau_sigma_min_Ry``), the lowest phonon frequency on the matdyn mesh and
per irreducible q, and the harmonic ZPE from ``phonon.dos``.

Rule (:func:`decide_coarse`), the **only** removal S5 may make on its own, is doc 13 §5
T2-SCREEN on the q 3^3 + q 4^3 ladders as revised by doc 16: SCREENED_OUT(coarse-A envelope
(all-sigma)) iff (1 + e_c) * A_hi,coarse < A_target, A_hi,coarse being the ENVELOPE built from
EVERY sigma row of both meshes, admitted or not, with a sigma -> 0 extrapolation
(:func:`coarse_envelope_all`), e_c = 0.5 provisional. Doc 13's admitted-only statistic
(:func:`coarse_envelope`) is reported, not used. With one mesh only, the fallback
M_c * A_hi,coarse(1 mesh) < A_target (M_c = 2) applies at q >= 4^3 only. Preconditions (else
DEFERRED): every mesh at or above the per-element cutoff floor (H 80 Ry, doc 16 §3), admission
on every mesh, >= 3 admitted sigma on q 4^3, non-empty plateau, S4 n_sheets == 1, no
imaginary modes or dropped omega <= 0 weight, no q 2^3. Never on lambda, omega_log, m or
imaginary modes. e_c / M_c are never passed to the core engine. Comparisons are exact
(``fractions.Fraction``). The metrics also carry ``min_freq_K`` (S6 full-grid stability).
"""
from __future__ import annotations

import re
from dataclasses import asdict, dataclass, field
from fractions import Fraction
from pathlib import Path
from typing import Any

from .profiles import th_parts

CM1_TO_K = 1.438776877   # hc/k_B in K cm (CODATA 2018)
CM1_TO_MEV = 0.12398419843320026
IMAG_TOL_K = 30.0        # doc 01 G-VALID item 2 (numerical-noise tolerance), used for the tag only

_FREQ = re.compile(r"freq\s*\(\s*\d+\)\s*=\s*(\S+)\s*\[THz\]\s*=\s*(\S+)\s*\[cm-1\]")


def dyn_min_freqs_cm1(workdir: Path) -> dict[str, float]:
    """Lowest frequency (cm^-1, negative = imaginary) per dyn file (irreducible q)."""
    out = {}
    for f in sorted(Path(workdir).glob("dyn[0-9]*"), key=lambda p: int(p.name[3:])):
        if f.name == "dyn0":
            continue
        vals = [float(m.group(2)) for m in _FREQ.finditer(f.read_text(errors="replace"))]
        if vals:
            out[f.name] = min(vals)
    return out


def zpe_meV_per_atom(phonon_dos: Path, nat: int) -> float | None:
    """Harmonic ZPE = (1/2) sum_modes hbar omega, from matdyn's ``fldos`` (col 1 = cm^-1,
    col 2 = total DOS), normalised to 3*nat modes; imaginary (negative) frequencies dropped."""
    import numpy as np
    p = Path(phonon_dos)
    if not p.is_file():
        return None
    rows = []
    for ln in p.read_text(errors="replace").splitlines():
        s = ln.split()
        if not s or s[0].startswith("#"):
            continue
        try:
            rows.append((float(s[0]), float(s[1])))
        except (ValueError, IndexError):
            continue
    if len(rows) < 2:
        return None
    a = np.array(rows)
    w, g = a[:, 0], a[:, 1]
    keep = w > 0
    if not keep.any():
        return None
    dw = np.gradient(w)
    norm = np.sum(g[keep] * dw[keep])
    if norm <= 0:
        return None
    mean_w = np.sum(w[keep] * g[keep] * dw[keep]) / norm
    return float(0.5 * 3 * nat * mean_w * CM1_TO_MEV / nat)


@dataclass
class SigmaRow:
    index: int
    sigma_Ry: float | None
    admitted: bool
    A_K: float | None = None
    lambda_: float | None = None
    omega_log_K: float | None = None
    omega2_K: float | None = None
    on_plateau: bool = False


@dataclass
class CoarseMetrics:
    rows: list[SigmaRow]
    plateau_sigma_min_Ry: float | None
    min_freq_matdyn_cm1: float | None
    min_freq_dyn_cm1: dict = field(default_factory=dict)
    zpe_meV_atom: float | None = None
    q: list | None = None
    flags: list[str] = field(default_factory=list)
    scf_pressure_kbar: float | None = None
    omega_max_K: float | None = None          # max phonon frequency carrying alpha^2F weight
    admission_ok: bool = True                 # mode-positive admission succeeded on this mesh
    ecutwfc_Ry: float | None = None           # cutoff of this run (doc 16 floor check)
    elements: list[str] = field(default_factory=list)

    @classmethod
    def from_metrics(cls, d: dict[str, Any], *, elements: Any = None,
                     ecutwfc_Ry: float | None = None) -> "CoarseMetrics":
        """Rebuild from a stored :meth:`to_metrics` dict (ledger ``tier_runs.metrics_json``
        or one ``by_mesh`` entry), for ``rtt redecide --stage S5``. ``elements`` / ``ecutwfc_Ry``
        fill in records written before doc 16 (they did not store them)."""
        keep = set(SigmaRow.__dataclass_fields__)
        rows = [SigmaRow(**{k: v for k, v in r.items() if k in keep}) for r in d.get("rows", [])]
        mfd = d.get("min_freq_dyn_cm1") or {}
        return cls(rows=rows, plateau_sigma_min_Ry=None,
                   min_freq_matdyn_cm1=d.get("min_freq_matdyn_cm1"),
                   min_freq_dyn_cm1=dict(mfd), zpe_meV_atom=d.get("zpe_meV_atom"), q=d.get("q"),
                   flags=list(d.get("flags") or []), scf_pressure_kbar=d.get("scf_pressure_kbar"),
                   omega_max_K=d.get("omega_max_K"),
                   admission_ok=bool(d.get("admission_ok", bool(d.get("admitted_sigmas")))),
                   ecutwfc_Ry=d.get("ecutwfc_Ry", ecutwfc_Ry),
                   elements=list(d.get("elements") or elements or []))

    @property
    def plateau(self) -> list[SigmaRow]:
        return [r for r in self.rows if r.on_plateau]

    @property
    def min_freq_cm1(self) -> float | None:
        vals = [v for v in [self.min_freq_matdyn_cm1, *self.min_freq_dyn_cm1.values()] if v is not None]
        return min(vals) if vals else None

    @property
    def min_freq_K(self) -> float | None:
        """Lowest harmonic frequency in K, signed (negative = imaginary, QE convention), over
        every irreducible q of the DFPT mesh (dyn files: the dynamical matrices of the run's own,
        final relaxed structure) and the matdyn interpolation mesh (uniform k_matdyn grid,
        ``matdyn_factor`` x q). None if neither is available. This is the S6 full-grid
        stability input (``Validity.maxImagModeK`` = max(0, -min_freq_K))."""
        m = self.min_freq_cm1
        return None if m is None else m * CM1_TO_K

    @property
    def needs_phonon_check(self) -> bool:
        m = self.min_freq_cm1
        return m is not None and -m * CM1_TO_K > IMAG_TOL_K

    def to_metrics(self) -> dict[str, Any]:
        """S5 metrics with the §6.1 key names (plus per-sigma detail)."""
        pl = self.plateau
        return {"A_K_by_sigma": [r.A_K for r in pl],
                "plateau_sigmas": [r.index for r in pl],
                "plateau_sigma_Ry": [r.sigma_Ry for r in pl],
                "lambda_by_sigma": [r.lambda_ for r in pl],
                "omega_log_K_by_sigma": [r.omega_log_K for r in pl],
                "omega2_K_by_sigma": [r.omega2_K for r in pl],
                "A_mid_K": _median([r.A_K for r in pl if r.A_K is not None]),
                "A_lo_K": min([r.A_K for r in self.rows if r.admitted and r.A_K is not None],
                              default=None),
                "omega_max_K": self.omega_max_K,
                "min_freq_cm1": self.min_freq_cm1,
                "min_freq_K": self.min_freq_K,
                "max_imag_mode_K": (None if self.min_freq_K is None
                                    else max(0.0, -self.min_freq_K)),
                "min_freq_sources": {"dyn_irreducible_q": len(self.min_freq_dyn_cm1),
                                     "matdyn_mesh": self.min_freq_matdyn_cm1 is not None},
                "min_freq_matdyn_cm1": self.min_freq_matdyn_cm1,
                "min_freq_dyn_cm1": self.min_freq_dyn_cm1,
                "zpe_meV_atom": self.zpe_meV_atom, "q": self.q,
                "admitted_sigmas": [r.index for r in self.rows if r.admitted],
                "admission_ok": self.admission_ok, "ecutwfc_Ry": self.ecutwfc_Ry,
                "elements": list(self.elements),
                "needs_phonon_check": self.needs_phonon_check,
                "scf_pressure_kbar": self.scf_pressure_kbar,
                "rows": [asdict(r) for r in self.rows], "flags": list(self.flags)}


def _median(xs: list[float]) -> float | None:
    if not xs:
        return None
    s = sorted(xs)
    n = len(s)
    return s[n // 2] if n % 2 else 0.5 * (s[n // 2 - 1] + s[n // 2])


def coarse_metrics(collected: Any, run: Any, plateau_sigma_min_Ry: float | None = None
                   ) -> CoarseMetrics:
    """Metrics from a :class:`rtt.dft.collect.Collected` and its :class:`DFPTRun`."""
    import numpy as np

    from rtt.core.dfpt.parsers import parse_qe_a2f_dos, parse_qe_lambda_file
    from rtt.core.dfpt.units import to_kelvin_factor
    w = Path(run.workdir)
    if plateau_sigma_min_Ry is None:
        plateau_sigma_min_Ry = run.profile.get("plateau_sigma_min_Ry")
    nsig = int((run.settings or {}).get("el_ph_nsigma", 10))
    lam = parse_qe_lambda_file(w / "lambda") if (w / "lambda").exists() else []
    admitted_idx = {sp.provenance.get("smearing_index"): sp for sp in collected.spectra}
    rows = []
    f = to_kelvin_factor("Ry")
    for i in range(1, nsig + 1):
        sig = lam[i - 1].broadening_Ry if len(lam) >= i else None
        sp = admitted_idx.get(i)
        r = SigmaRow(index=i, sigma_Ry=sig, admitted=sp is not None)
        if sp is not None:
            r.A_K, r.lambda_ = sp.A(), sp.lambda_()
            r.omega_log_K, r.omega2_K = sp.omega_log(), sp.omega2()
        elif (w / f"a2F.dos{i}").exists():      # excluded smearing: kept as upper evidence (doc 16)
            try:
                t = parse_qe_a2f_dos(w / f"a2F.dos{i}")
                pos = t.E_Ry > 0
                a, E = t.a2f_total[pos], t.E_Ry[pos]
                r.A_K = float(np.sum(a) * t.delta_Ry * f)
                r.lambda_ = float(2.0 * np.sum(a / E) * t.delta_Ry)
                if r.lambda_ > 0:
                    r.omega_log_K = float(np.exp(2.0 / r.lambda_ * np.sum(a * np.log(E * f) / E)
                                                 * t.delta_Ry))
            except Exception:
                pass
        r.on_plateau = bool(r.admitted and sig is not None and plateau_sigma_min_Ry is not None
                            and sig >= float(plateau_sigma_min_Ry) - 1e-12)
        rows.append(r)
    fmin = None
    for sp in collected.spectra:
        fmin = sp.provenance.get("min_freq_matdyn_grid_cm1", fmin)
    if fmin is None and (w / "freq.dos").is_file():
        # read directly so the stability check survives a failed a2F collection
        try:
            from rtt.core.dfpt.workflow import min_frequency_cm1
            fmin = min_frequency_cm1(w / "freq.dos")
        except Exception:
            fmin = None
    nat = len((run.cell0 or {}).get("species", [])) or 1
    m = CoarseMetrics(rows=rows, plateau_sigma_min_Ry=plateau_sigma_min_Ry,
                      min_freq_matdyn_cm1=fmin, min_freq_dyn_cm1=dyn_min_freqs_cm1(w),
                      zpe_meV_atom=zpe_meV_per_atom(w / "phonon.dos", nat), q=run.q,
                      flags=list(collected.flags), admission_ok=bool(collected.ok),
                      ecutwfc_Ry=(run.settings or {}).get("ecutwfc_Ry"),
                      elements=sorted(set((run.cell0 or {}).get("species", []))))
    if not m.plateau:
        m.flags.append("NO_PLATEAU")
    if collected.spectra:
        m.omega_max_K = max(sp.omega_max() for sp in collected.spectra)
    for sp in collected.spectra:
        ass = sp.provenance.get("assumptions", []) or []
        if any(str(a).startswith("A-DROP-NONPOS") for a in ass) and "A_DROP_NONPOS" not in m.flags:
            m.flags.append("A_DROP_NONPOS")      # weight at omega <= 0 dropped: lowers A (RC1 F9)
        if any(str(a).startswith("A-NEG-CLIP") for a in ass) and "A_NEG_CLIP" not in m.flags:
            m.flags.append("A_NEG_CLIP")         # raises A: not a removal hazard, recorded
    sc = w / "scf.out"
    if sc.is_file():
        from .electronic import PULAY_TOL_KBAR, parse_pressure_kbar
        p = parse_pressure_kbar(sc.read_text(errors="replace"))
        m.scf_pressure_kbar = p
        target = 10.0 * float((run.settings or {}).get("pressure_GPa", 0.0))
        if (p is not None and (run.settings or {}).get("do_relax", False)
                and abs(p - target) > max(PULAY_TOL_KBAR, 0.05 * abs(target))):
            m.flags.append(f"RELAX_PRESSURE_MISMATCH(P={p:.1f}kbar,target={target:.1f})")
    return m


_th = th_parts          # backwards-compatible name

#: flags meaning that part of the coupling sits at omega <= 0 or the a2F is unreliable
PHONON_HAZARD_PREFIXES = ("IMAGINARY_MODES", "NEGATIVE_FREQ_A2F_WEIGHT", "A_DROP_NONPOS",
                          "NEEDS_PHONON_CHECK")
EPS_OMEGA_FLOOR, EPS_OMEGA_CAP = 0.02, 0.9     # doc 01 eps_omega_min; doc 13 §5 cap
MIN_ADMITTED_FINE = 3                         # doc 13 §5: >= 3 admitted sigma on q 4^3
FALLBACK_MIN_Q = 4                            # doc 13 §5: single-mesh fallback only at q >= 4^3


def _qmin(m: CoarseMetrics) -> int | None:
    return min(m.q) if m.q else None


def coarse_envelope(fine: CoarseMetrics, coarse: CoarseMetrics, eps_sys: float) -> dict[str, Any]:
    """doc 01 §0.4 ENVELOPE(A) on the coarse ladders (doc 13 §5; port of
    ``src/scoring/coarse_envelope.py``)::

        s_nom  = smallest sigma admitted on BOTH meshes
        A_nom  = A(fine, s_nom);  d_sig = max_i |A(fine, s_i) - A_nom| / A_nom, admitted fine s_i
        d_q    = |A(fine, s_nom) - A(coarse, s_nom)| / A_nom
        eps_w  = |w_log(fine) - w_log(coarse)| / w_log(fine) at s_nom, floor 0.02, cap 0.9
        A_hi   = A_nom (1 + d_sig + d_q)(1 + eps_sys) / (1 - eps_w)

    A and omega_log are integrated from the alpha^2F atoms (core Spectrum methods), never
    taken from QE's ``lambda`` footer (wrong with imaginary modes: it sums E <= 0 bins).
    Returns {"ok": False, "why": ...} when the ladder cannot form an envelope."""
    adm_f = {r.index: r for r in fine.rows if r.admitted and r.A_K is not None}
    adm_c = {r.index: r for r in coarse.rows if r.admitted and r.A_K is not None}
    both = sorted(set(adm_f) & set(adm_c))
    if not both:
        return {"ok": False, "why": "no sigma admitted on both meshes"}
    s = both[0]
    a_nom = float(adm_f[s].A_K)
    if a_nom <= 0:
        return {"ok": False, "why": "A_nom <= 0"}
    d_sig = max(abs(float(r.A_K) - a_nom) / a_nom for r in adm_f.values())
    d_q = abs(a_nom - float(adm_c[s].A_K)) / a_nom
    wf, wc = adm_f[s].omega_log_K, adm_c[s].omega_log_K
    if wf is None or wc is None or not wf > 0:
        return {"ok": False, "why": "omega_log missing at s_nom"}
    eps_w = min(max(abs(wf - wc) / wf, EPS_OMEGA_FLOOR), EPS_OMEGA_CAP)
    a_hi = a_nom * (1 + d_sig + d_q) * (1 + eps_sys) / (1 - eps_w)
    return {"ok": True, "sigma_nom_index": s, "sigma_nom_Ry": adm_f[s].sigma_Ry, "A_nom_K": a_nom,
            "d_sig_A": d_sig, "d_q_A": d_q, "eps_omega": eps_w, "eps_sys": eps_sys, "A_hi_K": a_hi}


def _row_A(m: CoarseMetrics) -> dict[int, SigmaRow]:
    """Every sigma row with a finite A > 0 and a known sigma, admitted or not (doc 16 §4)."""
    return {r.index: r for r in m.rows
            if r.A_K is not None and r.A_K == r.A_K and r.A_K > 0 and r.sigma_Ry is not None}


def a_top(m: CoarseMetrics) -> float | None:
    """Upper A estimate of one mesh (doc 16 §4): max over ALL sigma rows (admitted or not) of A,
    and the linear extrapolation to sigma = 0 from the two narrowest smearings,
    2 A(s_1) - A(s_2). Mode-positive admission removes exactly the narrow-sigma rows where A of a
    sharp-N(E_F) hydride is largest (Mg2IrH6: A 905 K at 0.01 Ry vs 304 K at 0.10 Ry), so the
    excluded rows are kept here as upper evidence. Never used as a proof input."""
    rows = sorted(_row_A(m).values(), key=lambda r: r.sigma_Ry)
    if not rows:
        return None
    top = max(float(r.A_K) for r in rows)
    if len(rows) >= 2:
        top = max(top, 2.0 * float(rows[0].A_K) - float(rows[1].A_K))
    return top


def coarse_envelope_all(fine: CoarseMetrics, coarse: CoarseMetrics | None,
                        eps_sys: float) -> dict[str, Any]:
    """Doc 16 T2-SCREEN statistic (replaces :func:`coarse_envelope` as the removal statistic)::

        A_nom  = max over meshes of a_top(mesh)         (all rows; sigma -> 0 extrapolation)
        d_sig  = (A_nom - min_i A(fine, s_i)) / A_nom    (all fine rows)
        s*     = smallest sigma with A and omega_log on BOTH meshes (all rows)
        d_q    = |A(fine, s*) - A(coarse, s*)| / A(fine, s*)
        eps_w  = |w_log(fine, s*) - w_log(coarse, s*)| / w_log(fine, s*), floor 0.02, cap 0.9
        A_hi   = A_nom (1 + d_sig + d_q)(1 + eps_sys) / (1 - eps_w)

    With one mesh (``coarse`` None): d_q = 0, eps_w = 0.02 floor. Returns {"ok": False, ...}
    when no row has A. The doc 13 statistic (:func:`coarse_envelope`) is reported alongside."""
    rf = _row_A(fine)
    if not rf:
        return {"ok": False, "why": "no sigma row with A > 0 on the fine mesh"}
    tops = [t for t in (a_top(fine), a_top(coarse) if coarse is not None else None) if t is not None]
    a_nom = max(tops)
    d_sig = (a_nom - min(float(r.A_K) for r in rf.values())) / a_nom
    d_q, eps_w, s_star = 0.0, EPS_OMEGA_FLOOR, None
    if coarse is not None:
        rc = _row_A(coarse)
        both = sorted(i for i in set(rf) & set(rc)
                      if rf[i].omega_log_K and rc[i].omega_log_K and rf[i].omega_log_K > 0)
        if not both:
            return {"ok": False, "why": "no sigma with A and omega_log on both meshes"}
        s_star = both[0]
        af, ac = float(rf[s_star].A_K), float(rc[s_star].A_K)
        d_q = abs(af - ac) / af
        wf, wc = float(rf[s_star].omega_log_K), float(rc[s_star].omega_log_K)
        eps_w = min(max(abs(wf - wc) / wf, EPS_OMEGA_FLOOR), EPS_OMEGA_CAP)
    a_hi = a_nom * (1 + d_sig + d_q) * (1 + eps_sys) / (1 - eps_w)
    return {"ok": True, "A_nom_K": a_nom, "A_top_fine_K": a_top(fine),
            "A_top_coarse_K": a_top(coarse) if coarse is not None else None,
            "sigma_star_index": s_star, "d_sig_A": d_sig, "d_q_A": d_q, "eps_omega": eps_w,
            "eps_sys": eps_sys, "A_hi_K": a_hi}


def ecut_floor_blockers(meshes: list[CoarseMetrics], floors: Any) -> list[tuple[str, str]]:
    """Doc 16 §3: a removal needs every mesh at or above the per-element cutoff floor
    (tiers ``S5.ecut_floor_Ry``). Unknown cutoff or elements with a floor table -> blocked."""
    from .profiles import ecut_floor_Ry
    if not floors:
        return []
    out = []
    for m in meshes:
        if not m.elements:
            out.append(("ecut_unknown", f"q = {m.q}: elements unknown, cutoff floor not checkable"))
            continue
        need = ecut_floor_Ry(m.elements, floors)
        if need is None:
            continue
        if m.ecutwfc_Ry is None or float(m.ecutwfc_Ry) < need - 1e-9:
            out.append(("ecut_below_floor", f"q = {m.q}: ecutwfc {m.ecutwfc_Ry} Ry < floor "
                                            f"{need:g} Ry for {sorted(m.elements)}"))
    return out


def decide_coarse(metrics: CoarseMetrics | None, *, floor_K: Any, eps_sys: Any,
                  margin_M_c: Any = None, e_c: Any = None,
                  metrics_coarser: CoarseMetrics | None = None,
                  n_sheets: int | None = None, w_min: float | None = None,
                  tier_run_id: int | None = None, error: str | None = None,
                  ecut_floors: Any = None) -> Any:
    """S5 coarse screen -> :class:`rtt.records.Outcome`: PASS, SCREENED_OUT, DEFERRED or ERROR
    (never a rejection). ``metrics`` is the finest coarse mesh (q 4^3), ``metrics_coarser`` the
    q 3^3 mesh; ``floor_K`` = A_target = a0*T* (``Targets.floor``; the surrogate for controls).

    Rule T2-SCREEN (doc 13 §5 as revised by doc 16 §4), both meshes present::

        SCREENED_OUT(coarse-A envelope (all-sigma))  iff  (1 + e_c) * A_hi,coarse < A_target

    with A_hi,coarse from :func:`coarse_envelope_all` (every sigma row, admitted or not, and the
    sigma -> 0 extrapolation; doc 13's admitted-only :func:`coarse_envelope` under-covered
    Mg2IrH6 by 0.60x) and e_c = tiers ``S5.e_c`` (provisional 0.5). Fallback with one mesh,
    only at q >= 4^3::

        SCREENED_OUT(coarse-A max (all-sigma, single mesh))  iff  M_c * A_hi,coarse(1 mesh) < A_target

    Doc 16 §3 adds the blocker ``ecut_below_floor``: every mesh must run at or above the
    per-element cutoff floor ``ecut_floors`` (tiers ``S5.ecut_floor_Ry``; H 80 Ry).

    A removal fires only if ALL preconditions hold, otherwise it becomes DEFERRED(<first
    blocker>) with every blocker tagged: no mesh below q 3^3 (q 2^3 is never used); admission
    succeeded on every mesh; >= 3 admitted sigma on the fine mesh (two-mesh rule); a non-empty
    sigma plateau; S4 ``n_sheets == 1`` (the isotropic A of a multiband metal is not a valid
    input, doc 01 §5.1; ``w_min`` is recorded only); no imaginary modes / dropped omega <= 0
    weight / unreliable a2F on any mesh (R2). eps_sys is required (ValueError). Never screens on
    lambda, omega_log, m or imaginary modes. Outcomes made with a provisional e_c or M_c carry
    ``PROVISIONAL_MARGIN``."""
    from rtt.records import DEFERRED, ERROR, PASS, SCREENED_OUT, Outcome
    if error or metrics is None:
        return Outcome(stage="S5", outcome=ERROR, reason=error or "no metrics",
                       tier_run_id=tier_run_id)
    es, es_src, es_prov = th_parts(eps_sys)
    if es is None or float(es) < 0:
        raise ValueError("eps_sys is required (>= 0) for the coarse screen")
    two = metrics_coarser is not None
    meshes = [metrics] + ([metrics_coarser] if two else [])
    if two:
        ec, src, prov = th_parts(e_c)
        if ec is None or float(ec) < 0:
            raise ValueError("e_c is required for the two-mesh T2-SCREEN rule")
        factor = 1 + Fraction(str(float(ec)))
    else:
        mc, src, prov = th_parts(margin_M_c)
        if mc is None or float(mc) <= 0:
            raise ValueError("margin M_c must be > 0 for the single-mesh fallback")
        factor = Fraction(str(float(mc)))
    tags = ["PROVISIONAL_MARGIN"] if prov else []
    if es_prov:
        tags.append("PROVISIONAL_EPS_SYS")
    target = Fraction(floor_K) if not isinstance(floor_K, Fraction) else floor_K

    blockers: list[tuple[str, str]] = []
    # doc 16 §4: the removal statistic uses EVERY sigma row (admitted or not) plus the sigma -> 0
    # extrapolation; the doc 13 statistic (admitted rows only) is reported in the reason only.
    env = coarse_envelope_all(metrics, metrics_coarser if two else None, float(es))
    if not env["ok"]:
        return Outcome(stage="S5", outcome=PASS, tier_run_id=tier_run_id,
                       provisional=bool(prov), tags=tags + ["ENVELOPE_INCOMPLETE"],
                       reason=f"T2-SCREEN cannot form the envelope ({env['why']}): no removal")
    stat = Fraction(env["A_hi_K"])
    parts = (f"A_nom {env['A_nom_K']:.1f} K (all-sigma, sigma->0), d_sig {env['d_sig_A']:.3f}, "
             f"d_q {env['d_q_A']:.3f}, eps_w {env['eps_omega']:.3f}")
    if two:
        crit = "coarse-A envelope (all-sigma)"
        old = coarse_envelope(metrics, metrics_coarser, float(es))
        doc13 = f"; doc13 A_hi {old['A_hi_K']:.1f} K" if old.get("ok") else ""
        detail = (f"(1+e_c)*A_hi,coarse = {float(factor * stat):.1f} K vs A_target "
                  f"{float(target):.1f} K ({parts}{doc13})")
    else:
        crit = "coarse-A max (all-sigma, single mesh)"
        detail = (f"M_c*A_hi,coarse = {float(factor * stat):.1f} K vs A_target "
                  f"{float(target):.1f} K (single-mesh fallback; {parts})")
        q = _qmin(metrics)
        if q is None or q < FALLBACK_MIN_Q:
            blockers.append(("single_mesh_q_below_4",
                             f"single-mesh fallback needs q >= 4^3 (q = {metrics.q})"))
    blockers += ecut_floor_blockers(meshes, ecut_floors)
    would_remove = factor * stat < target

    for m in meshes:
        q = _qmin(m)
        if q is not None and q < 3:
            blockers.append(("q2_not_allowed", f"mesh q = {m.q} below 3^3 is never used"))
        if not m.admission_ok:
            blockers.append(("admission_failed", f"smearing admission failed on q = {m.q}"))
        if m.needs_phonon_check or any(f.startswith(PHONON_HAZARD_PREFIXES) for f in m.flags):
            blockers.append(("imaginary_modes",
                             f"imaginary modes / dropped omega<=0 weight on q = {m.q}"))
    if two and sum(1 for r in metrics.rows if r.admitted) < MIN_ADMITTED_FINE:
        blockers.append(("ladder_short", f"< {MIN_ADMITTED_FINE} admitted sigma on the fine mesh"))
    if not metrics.plateau:
        blockers.append(("no_plateau", "no admitted sigma on the plateau"))
    if n_sheets != 1:
        what = "unknown" if n_sheets is None else str(n_sheets)
        blockers.append(("multiband", f"n_sheets = {what}: the screen needs a single-sheet metal "
                                      f"(doc 01 §5.1; w_min = {w_min})"))
    if any(b[0] == "imaginary_modes" for b in blockers):
        tags.append("NEEDS_PHONON_CHECK")

    if not would_remove:
        return Outcome(stage="S5", outcome=PASS, tier_run_id=tier_run_id, provisional=bool(prov),
                       tags=tags, reason=detail + ": not below target")
    if blockers:
        return Outcome(stage="S5", outcome=DEFERRED, criterion=blockers[0][0],
                       value=float(factor * stat), threshold=float(target),
                       source_id=src or "doc16#4", provisional=bool(prov),
                       tier_run_id=tier_run_id, tags=tags + [f"BLOCKED:{b[0]}" for b in blockers],
                       reason=detail + "; removal blocked: " + "; ".join(b[1] for b in blockers))
    return Outcome(stage="S5", outcome=SCREENED_OUT, criterion=crit, value=float(factor * stat),
                   threshold=float(target), source_id=src or "doc16#4", provisional=bool(prov),
                   tier_run_id=tier_run_id, tags=tags, reason=detail)
