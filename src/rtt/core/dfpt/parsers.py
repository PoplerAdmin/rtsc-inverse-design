"""Parsers for alpha^2F outputs and conversion to the atomic measure ``rtt.core.spectrum.Spectrum``.

FORMAT SOURCES (QE source fetched from github.com/QEF/q-e at tags ``qe-6.7.0`` and ``develop``;
the a2Fdos routine is identical between them apart from ASR arguments)
-------------------------------------------------------------------------------------------
QE ``a2F.dos<isig>`` -- written by SUBROUTINE a2Fdos in ``PHonon/PH/matdyn.f90``::

    write(ifn,*)
    write(ifn,*) '# Eliashberg function a2F (per both spin)'
    write(ifn,*) '#  frequencies in Rydberg  '
    write(ifn,*) '# DOS normalized to E in Rydberg: a2F_total, a2F(mode) '
    write(ifn,*)
    do n = 1, ndos
       E = Emin + (n-1)*DeltaE + 0.5d0*DeltaE           ! bin centre, Rydberg
       dos_a2F(j) = dos_gam(...)/dos_ee(isig)/2/pi      ! tetrahedron density at E
       lambda = lambda + 2.d0 * dos_tot/E * DeltaE
       write (ifn, '(3X,1000E16.6)') E, dos_tot, dos_a2F(1:nmodes)
    enddo
    write(ifn,*) " lambda =",lambda,'   Delta = ',DeltaE

  => column 1 = phonon energy in Ry (bin centres of a uniform grid of width Delta, in Ry);
     column 2 = alpha^2F_total (dimensionless density per unit energy-in-Ry, i.e.
     int a2F dE(Ry) is the Ry-measure); columns 3.. = per-mode a2F. QE's own lambda is the
     midpoint sum sum_n 2 a2F(E_n)/E_n * Delta. Hence the atomic measure that reproduces QE's
     lambda exactly is: atom at E_n, mass a2F(E_n)*Delta, both converted Ry -> K
     (factor RYDBERG_SI/K_BOLTZMANN_SI = 157887.512 K, see units.py).
     Emin = min(0, min freq) (matdyn.f90 main program), so with imaginary modes the grid
     extends to E < 0; matdyn zeroes gamma for freq <= 1e-7 Ry ("correction for small
     frequencies") but tetrahedra straddling 0 can leave weight at E <= 0. Those atoms are
     dropped from the Spectrum and the dropped weight is recorded as a WARNING.
     Delta is converted from the input DeltaE in cm^-1 by /RY_TO_CMM1 before the call.

QE ``lambda`` file -- same routine::

    write(400,'(" Broadening ",F8.4," lambda ",F12.4," dos(Ef)",F8.4," omega_ln [K]",F12.4)') &
         deg(isig), lambda, dos_ee(isig), E
  with lambda = sum_q wq sum_nu gamma/omega^2 / (N(Ef) pi) (direct q sum, NOT the histogram),
  E = exp(sum log(omega) gamma/omega^2 / sum gamma/omega^2) / K_BOLTZMANN_RY (kelvin),
  modes with freq <= 1e-5 Ry excluded; deg(isig) = isig*el_ph_sigma in Ry (elphon.f90:
  ``deg(isig) = isig * el_ph_sigma``); dos(Ef) in states/spin/Ry/cell (elphon.f90 FORMAT 9005).

EPW ``<prefix>.a2f`` -- SUBROUTINE evaluate_a2f_lambda, ``EPW/src/supercond.f90`` (qe-6.7.0)::

    WRITE(iua2ffil, '(" w[meV] a2f and integrated 2*a2f/w for ", i4, " smearing values")')
    l_a2f_tmp = l_a2f_tmp + 2.0d0 * (a2f(iwph, ismear) / wsph(iwph)) * dwsph
    WRITE(iua2ffil, '(f12.7, 20f12.7)') wsph(iwph) * 1000.d0, a2f(iwph, :), l_a2f(iwph, :)
    "Integrated el-ph coupling" / "  #   <lambda per smearing>" / "Phonon smearing (meV)" / ...
  with wsph(iwph) = iwph * dwsph, dwsph = wsphmax/nqstep (io_eliashberg.f90) -> frequency in
  meV on a uniform grid starting at dw; atom at w_i with mass a2F_i*dw reproduces EPW's lambda.
  ``<prefix>.a2f_iso``: header "w[meV] a2f" then (w[meV], a2f[, mode-projected...]) rows.
"""
from __future__ import annotations

import math
import re
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Sequence

import numpy as np

from .units import canonical_unit, to_kelvin_factor

try:  # shared type owned by another team (ARCH section 2)
    from rtt.core.spectrum import Spectrum
except Exception:  # pragma: no cover - only while the shared module is missing
    Spectrum = None  # type: ignore[assignment]


# ------------------------------------------------------------------------------------------
# number parsing (Fortran E/D formats, including 3-digit exponents written without 'E')
# ------------------------------------------------------------------------------------------
_FORT_NOE = re.compile(r"^([-+]?\d*\.\d+)([-+]\d{3})$")


def ffloat(tok: str) -> float:
    t = tok.strip().replace("D", "E").replace("d", "e")
    try:
        return float(t)
    except ValueError:
        m = _FORT_NOE.match(t)
        if m:
            return float(m.group(1) + "E" + m.group(2))
        raise


def _numeric_row(line: str) -> list[float] | None:
    toks = line.split()
    if not toks:
        return None
    try:
        return [ffloat(t) for t in toks]
    except ValueError:
        return None


# ------------------------------------------------------------------------------------------
# QE matdyn a2F.dos / lambda
# ------------------------------------------------------------------------------------------
@dataclass
class QEA2FDos:
    path: str
    isig: int | None
    E_Ry: np.ndarray
    a2f_total: np.ndarray
    a2f_modes: np.ndarray            # (n, nmodes)
    lambda_reported: float | None    # footer " lambda = ... Delta = ..."
    delta_Ry: float | None
    header: list[str] = field(default_factory=list)
    complete: bool = True            # False only via allow_incomplete=True (diagnostics)

    @property
    def nmodes(self) -> int:
        return int(self.a2f_modes.shape[1]) if self.a2f_modes.ndim == 2 else 0

    def lambda_all_bins(self) -> float:
        """QE's histogram formula over ALL bins (including E <= 0 bins), for format checks."""
        d = self.delta_Ry if self.delta_Ry else float(np.median(np.diff(self.E_Ry)))
        return float(np.sum(2.0 * self.a2f_total / self.E_Ry * d))


_FOOTER = re.compile(r"lambda\s*=\s*(\S+)\s+Delta\s*=\s*(\S+)", re.I)


class IncompleteFileError(ValueError):
    """The file lacks the end-of-data marker the code always writes (truncated/incomplete).
    A truncated alpha^2F silently lowers A and lambda, which is anti-conservative for
    REJECTED_PROVEN, so this is always an error unless ``allow_incomplete=True``."""


def parse_qe_a2f_dos(path: str | Path, allow_incomplete: bool = False) -> QEA2FDos:
    """Parse a matdyn a2F.dos file. The footer ``lambda = ... Delta = ...`` (written after the
    last row by matdyn.f90) is REQUIRED: without it the file may be truncated. With
    ``allow_incomplete=True`` (diagnostics only) a footer-less file is returned with
    ``complete=False``; spectra built from it carry the assumption 'A-FILE-INCOMPLETE' and the
    flag 'INCOMPLETE_INPUT_FILE'."""
    p = Path(path)
    m = re.search(r"a2F\.dos(\d+)$", p.name)
    isig = int(m.group(1)) if m else None
    rows: list[list[float]] = []
    header: list[str] = []
    lam = delta = None
    for line in p.read_text(errors="replace").splitlines():
        s = line.strip()
        if not s:
            continue
        mf = _FOOTER.search(s)
        if mf:
            if lam is not None:
                raise ValueError(f"{p}: more than one 'lambda = ... Delta =' footer")
            lam, delta = ffloat(mf.group(1)), ffloat(mf.group(2))
            continue
        if lam is not None:
            raise ValueError(f"{p}: content after the 'lambda = ... Delta =' footer: {s[:60]!r}")
        if s.startswith("#"):
            header.append(s)
            continue
        r = _numeric_row(s)
        if r is None:
            header.append(s)
            continue
        rows.append(r)
    if not rows:
        raise ValueError(f"{p}: no data rows")
    ncol = len(rows[0])
    if any(len(r) != ncol for r in rows) or ncol < 2:
        raise ValueError(f"{p}: ragged or too-narrow data table")
    a = np.asarray(rows, dtype=float)
    if not any("Rydberg" in h for h in header):
        raise ValueError(f"{p}: header does not state 'frequencies in Rydberg' -- not a QE "
                         "matdyn a2F.dos file (or an unsupported version)")
    complete = True
    if lam is None:
        if not allow_incomplete:
            raise IncompleteFileError(
                f"{p}: no 'lambda = ... Delta =' footer -- the file is truncated or incomplete "
                "(matdyn writes the footer after the last row); refusing to build a spectrum")
        complete = False
    elif not (math.isfinite(delta) and delta > 0):
        raise ValueError(f"{p}: invalid footer Delta={delta}")
    return QEA2FDos(str(p), isig, a[:, 0], a[:, 1], a[:, 2:], lam, delta, header, complete)


@dataclass
class QELambdaEntry:
    broadening_Ry: float
    lambda_: float
    dos_ef: float | None              # states/spin/Ry/cell (None if Fortran overflowed '****')
    omega_ln_K: float | None


_LAMBDA_LINE = re.compile(
    r"Broadening\s+(\S+)\s+lambda\s+(\S+)\s+dos\(Ef\)\s*(\S+)\s+omega_ln\s*\[K\]\s*(\S+)")


def parse_qe_lambda_file(path: str | Path) -> list[QELambdaEntry]:
    out = []
    for line in Path(path).read_text(errors="replace").splitlines():
        m = _LAMBDA_LINE.search(line)
        if not m:
            continue

        def g(i: int) -> float | None:
            t = m.group(i)
            return None if "*" in t else ffloat(t)

        b, lam = g(1), g(2)
        if b is None or lam is None:
            raise ValueError(f"{path}: unreadable line {line!r}")
        out.append(QELambdaEntry(b, lam, g(3), g(4)))
    return out


def parse_qe_version(text: str) -> str | None:
    m = re.search(r"Program\s+(\S+)\s+v\.(\S+)", text)
    return f"{m.group(1)} {m.group(2)}" if m else None


# ------------------------------------------------------------------------------------------
# EPW
# ------------------------------------------------------------------------------------------
@dataclass
class EPWA2F:
    path: str
    omega_meV: np.ndarray
    a2f: np.ndarray                  # (n, nsmear)
    lambda_cumulative: np.ndarray    # (n, nsmear)
    lambda_reported: list[float]
    phonon_smearing_meV: list[float]
    electron_smearing_eV: float | None
    fermi_window_eV: float | None
    complete: bool = True


def parse_epw_a2f(path: str | Path, allow_incomplete: bool = False) -> EPWA2F:
    """Parse EPW <prefix>.a2f. The 'Integrated el-ph coupling' block (one value per smearing,
    written after the last row) is REQUIRED; ``allow_incomplete=True`` is for diagnostics and
    marks the result ``complete=False``."""
    p = Path(path)
    lines = p.read_text(errors="replace").splitlines()
    nsm = None
    rows = []
    lam_rep: list[float] = []
    smear: list[float] = []
    esm = fw = None
    mode = None
    for line in lines:
        s = line.strip()
        if not s:
            continue
        mh = re.search(r"for\s+(\d+)\s+smearing values", s)
        if mh:
            nsm = int(mh.group(1))
            continue
        if s.startswith("Integrated el-ph coupling"):
            mode = "lam"
            continue
        if s.startswith("Phonon smearing"):
            mode = "smear"
            continue
        if s.startswith("Electron smearing"):
            esm = ffloat(s.split()[-1])
            continue
        if s.startswith("Fermi window"):
            fw = ffloat(s.split()[-1])
            continue
        if s.startswith("Summed el-ph coupling"):
            continue
        if s.startswith("#"):
            vals = [ffloat(t) for t in s[1:].split()]
            if mode == "lam":
                lam_rep = vals
            elif mode == "smear":
                smear = vals
            continue
        r = _numeric_row(s)
        if r is not None and mode is None:
            rows.append(r)
    if nsm is None:
        raise ValueError(f"{p}: missing EPW header 'w[meV] a2f and integrated 2*a2f/w for N "
                         "smearing values'")
    a = np.asarray(rows, dtype=float)
    if a.ndim != 2 or a.shape[1] != 1 + 2 * nsm:
        raise ValueError(f"{p}: expected {1 + 2 * nsm} columns, got {a.shape}")
    complete = len(lam_rep) == nsm
    if not complete and not allow_incomplete:
        raise IncompleteFileError(
            f"{p}: missing or short 'Integrated el-ph coupling' block ({len(lam_rep)} of {nsm} "
            "values) -- the file is truncated or incomplete; refusing to build spectra")
    if complete:
        # the block repeats the last cumulative-lambda row: a mismatch means damaged rows
        last = a[-1, 1 + nsm:]
        for j in range(nsm):
            if abs(last[j] - lam_rep[j]) > 2e-7 + 1e-9 * abs(lam_rep[j]):
                raise ValueError(f"{p}: last cumulative lambda {last[j]} != integrated "
                                 f"coupling {lam_rep[j]} (smearing {j + 1})")
    return EPWA2F(str(p), a[:, 0], a[:, 1:1 + nsm], a[:, 1 + nsm:], lam_rep, smear, esm, fw,
                  complete)


# ------------------------------------------------------------------------------------------
# generic readers
# ------------------------------------------------------------------------------------------
def read_columns(path: str | Path, cols: Sequence[int] = (0, 1),
                 comments: str = "#!%") -> np.ndarray:
    """Numeric rows of a whitespace/comma table; non-numeric lines are skipped."""
    rows = []
    for line in Path(path).read_text(errors="replace").splitlines():
        s = line.strip()
        if not s or s[0] in comments:
            continue
        r = _numeric_row(s.replace(",", " "))
        if r is None or len(r) <= max(cols):
            continue
        rows.append([r[c] for c in cols])
    if not rows:
        raise ValueError(f"{path}: no numeric rows")
    return np.asarray(rows, dtype=float)


# ------------------------------------------------------------------------------------------
# sampled density -> atomic measure
# ------------------------------------------------------------------------------------------
RULES = ("bins", "trapezoid", "midpoint")


def grid_to_atoms(omega: Sequence[float], a2f: Sequence[float], rule: str = "bins",
                  delta: float | None = None) -> tuple[np.ndarray, np.ndarray]:
    """Atoms (positions, masses) of a sampled density, in the units of ``omega``.

    rule="bins"      : atom at omega_i, mass a2F_i * width_i; width = ``delta`` if given
                       (uniform histogram, e.g. QE a2F.dos / EPW), else the Voronoi cell of
                       omega_i (ends extended by half a spacing). Reproduces QE/EPW lambda.
    rule="trapezoid" : atom at omega_i, trapezoid quadrature weight (half cells at the ends).
    rule="midpoint"  : one atom per interval at its midpoint with mass
                       (a2F_i + a2F_{i+1})/2 * (omega_{i+1} - omega_i).
    No point is dropped here (see ``atoms_to_spectrum``).
    """
    x = np.asarray(omega, dtype=float)
    y = np.asarray(a2f, dtype=float)
    if x.shape != y.shape or x.ndim != 1 or len(x) < 2:
        raise ValueError("omega and a2f must be 1-D arrays of equal length >= 2")
    order = np.argsort(x, kind="stable")
    x, y = x[order], y[order]
    if np.any(np.diff(x) <= 0):
        raise ValueError("duplicate frequency points")
    if rule == "bins":
        if delta is not None:
            w = np.full_like(x, float(delta))
        else:
            edges = np.empty(len(x) + 1)
            edges[1:-1] = 0.5 * (x[1:] + x[:-1])
            edges[0] = x[0] - 0.5 * (x[1] - x[0])
            edges[-1] = x[-1] + 0.5 * (x[-1] - x[-2])
            w = np.diff(edges)
        return x, y * w
    if rule == "trapezoid":
        h = np.diff(x)
        w = np.zeros_like(x)
        w[:-1] += 0.5 * h
        w[1:] += 0.5 * h
        return x, y * w
    if rule == "midpoint":
        h = np.diff(x)
        return 0.5 * (x[1:] + x[:-1]), 0.5 * (y[1:] + y[:-1]) * h
    raise ValueError(f"unknown rule {rule!r}; choose from {RULES}")


def atoms_to_spectrum(pos: np.ndarray, mass: np.ndarray, *, unit: str, label: str = "",
                      provenance: dict | None = None, rel_envelope: float = 0.0,
                      negative_policy: str = "error", negative_tol_rel: float = 1e-6,
                      drop_zero: bool = True) -> "Spectrum":
    """Convert atoms in ``unit`` to a kelvin ``Spectrum``; drop omega <= 0 atoms (recorded).

    Negative masses at omega > 0 (numerical noise in a2F): if their total magnitude is at most
    ``negative_tol_rel`` * A_kept they are clipped to 0 (recorded); otherwise
    ``negative_policy="error"`` raises and ``"clip"`` clips with a WARNING. Clipping raises A
    and lambda, so it is recorded in provenance for downstream review.
    """
    if Spectrum is None:  # pragma: no cover
        raise ImportError("rtt.core.spectrum.Spectrum is not available")
    f = to_kelvin_factor(unit)
    pos_K = np.asarray(pos, dtype=float) * f
    mass_K = np.asarray(mass, dtype=float) * f
    prov = dict(provenance or {})
    bad = ~(np.isfinite(pos_K) & np.isfinite(mass_K))
    if np.any(bad):
        # never drop NaN/inf silently: e.g. QE matdyn dos_gam gives 0/0 = NaN when a mode has
        # freq exactly 0 at a tetrahedron vertex (gamma is zeroed there, Y = gamma/freq)
        raise ValueError(f"{int(np.count_nonzero(bad))} non-finite (NaN/inf) alpha^2F values in "
                         f"{prov.get('source_file', 'input')}: spectrum unusable")
    warnings: list[str] = list(prov.pop("warnings", []))
    assumptions: list[str] = list(prov.pop("assumptions", []))
    nonpos = pos_K <= 0
    dropped_signed = float(np.sum(mass_K[nonpos]))
    dropped_abs = float(np.sum(np.abs(mass_K[nonpos])))
    keep_pos, keep_mass = pos_K[~nonpos], mass_K[~nonpos]
    A_kept_pos = float(np.sum(keep_mass[keep_mass > 0]))
    neg = keep_mass < 0
    neg_mass = float(-np.sum(keep_mass[neg])) if np.any(neg) else 0.0
    if neg_mass > 0:
        if neg_mass > negative_tol_rel * max(A_kept_pos, 1e-300):
            msg = (f"negative alpha^2F mass {neg_mass:.6g} K at omega > 0 "
                   f"({neg_mass / max(A_kept_pos, 1e-300):.3g} of A)")
            if negative_policy == "error":
                raise ValueError(msg + " -- refusing to clip (negative_policy='error')")
            warnings.append("WARNING: " + msg + "; clipped to 0 (raises A and lambda)")
        keep_mass = np.where(neg, 0.0, keep_mass)
        # clipping raises A: anti-conservative for CERTIFIED_CANDIDATE -> assumption-level
        assumptions.append(f"A-NEG-CLIP({neg_mass:.6g}K): negative alpha^2F mass at omega>0 "
                           "clipped to 0; raises A and lambda (anti-conservative for CERTIFIED)")
    if dropped_abs > 0:
        frac = dropped_abs / max(A_kept_pos, 1e-300)
        warnings.append(f"WARNING: dropped {int(np.count_nonzero(nonpos & (mass_K != 0)))} atoms "
                        f"at omega <= 0 carrying |mass| {dropped_abs:.6g} K "
                        f"(signed {dropped_signed:.6g} K, {frac:.3g} of kept A); "
                        "nonzero weight at omega <= 0 indicates imaginary/unstable modes or "
                        "interpolation artefacts")
        if dropped_signed > 0:
            # dropping positive weight lowers A: anti-conservative for REJECTED_PROVEN
            assumptions.append(f"A-DROP-NONPOS({dropped_signed:.6g}K): alpha^2F weight at "
                               "omega<=0 (imaginary/unstable modes) excluded from the model "
                               "spectrum; lowers A (anti-conservative for REJECTED)")
    nz = keep_mass > 0 if drop_zero else np.ones_like(keep_mass, dtype=bool)
    om, w = keep_pos[nz], keep_mass[nz]
    if len(om) == 0:
        raise ValueError(f"no positive weight at omega > 0 in {prov.get('source_file', 'input')} "
                         "(alpha^2F identically zero after dropping omega <= 0): cannot form a "
                         "Spectrum")
    prov.update({
        "input_unit": canonical_unit(unit),
        "to_kelvin_factor": f,
        "n_points_in": int(len(pos_K)),
        "n_atoms": int(len(om)),
        "dropped_nonpositive_omega": {"count": int(np.count_nonzero(nonpos)),
                                      "abs_mass_K": dropped_abs, "signed_mass_K": dropped_signed},
        "clipped_negative_mass_K": neg_mass,
        "warnings": warnings,
        "assumptions": assumptions,
    })
    return Spectrum(tuple(float(v) for v in om), tuple(float(v) for v in w),
                    rel_envelope, label, prov)


def mode_negative_fractions(tab: QEA2FDos) -> list[float]:
    """Negative alpha^2F mass of each phonon mode at omega > 0, as a fraction of the total positive
    mass (uniform matdyn bins, so Delta cancels).

    alpha^2F_nu(omega) = (1/N_F) sum_q gamma_qnu/(pi omega_qnu) delta(omega - omega_qnu) with
    gamma_qnu = sum_k |g|^2 x (non-negative double-delta weights) is non-negative mode by mode, so
    any negative value is numerical: a double-delta smearing the k grid does not resolve, made
    visible by the Fourier interpolation of gamma near Gamma. Checking modes, not only the total,
    catches ringing hidden by cancellation between branches.
    """
    pos = tab.E_Ry > 0
    A = float(np.sum(np.clip(tab.a2f_total[pos], 0.0, None)))
    modes = tab.a2f_modes[pos]
    return [float(-np.sum(np.clip(modes[:, k], None, 0.0))) / max(A, 1e-300)
            for k in range(modes.shape[1])]


def qe_a2f_dos_to_spectrum(tab: QEA2FDos, *, label: str = "", broadening_Ry: float | None = None,
                           extra_provenance: dict | None = None, **kw: Any) -> "Spectrum":
    """QE a2F.dos -> Spectrum with QE's own bin rule (atom at E_n, mass a2F_n * Delta)."""
    delta = tab.delta_Ry
    extra_assumptions: list[str] = []
    if delta is None:
        if tab.complete:
            raise ValueError(f"{tab.path}: no Delta from the footer")
        # diagnostics only (allow_incomplete=True): never silently
        delta = float(np.median(np.diff(tab.E_Ry)))
        extra_assumptions.append("A-FILE-INCOMPLETE: QE a2F.dos without footer (possibly "
                                 "truncated); Delta from median spacing -- NOT a valid input "
                                 "for any proof")
    d = np.diff(tab.E_Ry)
    # E16.6 keeps 6 mantissa digits: absolute print error per value <= 5e-6 |E|
    if len(d) and np.max(np.abs(d - delta)) > 1e-4 * abs(delta) + 1.2e-5 * np.max(np.abs(tab.E_Ry)):
        raise ValueError(f"{tab.path}: frequency grid is not uniform with Delta={delta}")
    pos, mass = grid_to_atoms(tab.E_Ry, tab.a2f_total, "bins", delta=delta)
    prov = {"source_file": tab.path, "format": "QE matdyn a2F.dos (PHonon/PH/matdyn.f90 a2Fdos)",
            "smearing_index": tab.isig, "el_ph_broadening_Ry": broadening_Ry,
            "rule": "bins: atom at bin centre E_n, mass a2F(E_n)*Delta (QE histogram rule)",
            "delta_Ry": delta, "qe_reported_lambda_histogram": tab.lambda_reported,
            "complete": tab.complete}
    prov.update(extra_provenance or {})
    if extra_assumptions:
        prov["assumptions"] = list(prov.get("assumptions", [])) + extra_assumptions
        prov["flags"] = list(prov.get("flags", [])) + ["INCOMPLETE_INPUT_FILE"]
    return atoms_to_spectrum(pos, mass, unit="Ry", label=label, provenance=prov, **kw)


def epw_a2f_to_spectra(tab: EPWA2F, *, label_prefix: str = "",
                       extra_provenance: dict | None = None, **kw: Any) -> list["Spectrum"]:
    """EPW <prefix>.a2f -> one Spectrum per phonon smearing (atom at w_i, mass a2F_i*dw)."""
    w = tab.omega_meV
    dw = float(w[0]) if w[0] > 0 else float(np.median(np.diff(w)))
    if np.max(np.abs(np.diff(w) - dw)) > 2e-6 * max(1.0, float(np.max(w))) + 1e-6:
        dw_opt: float | None = None   # non-uniform: Voronoi cells
    else:
        dw_opt = dw
    out = []
    for j in range(tab.a2f.shape[1]):
        pos, mass = grid_to_atoms(w, tab.a2f[:, j], "bins", delta=dw_opt)
        prov = {"source_file": tab.path, "format": "EPW <prefix>.a2f (EPW/src/supercond.f90)",
                "smearing_index": j + 1,
                "phonon_smearing_meV": tab.phonon_smearing_meV[j] if j < len(tab.phonon_smearing_meV) else None,
                "epw_reported_lambda": tab.lambda_reported[j] if j < len(tab.lambda_reported) else None,
                "rule": "bins: atom at w_i, mass a2F_i*dw (EPW lambda rule)", "dw_meV": dw_opt,
                "complete": tab.complete}
        prov.update(extra_provenance or {})
        if not tab.complete:
            prov["assumptions"] = list(prov.get("assumptions", [])) + [
                "A-FILE-INCOMPLETE: EPW .a2f without the integrated-coupling block (possibly "
                "truncated) -- NOT a valid input for any proof"]
            prov["flags"] = list(prov.get("flags", [])) + ["INCOMPLETE_INPUT_FILE"]
        out.append(atoms_to_spectrum(pos, mass, unit="meV",
                                     label=f"{label_prefix}epw_smear{j + 1}", provenance=prov, **kw))
    return out


def read_two_column(path: str | Path, unit: str, *, rule: str = "bins",
                    cols: Sequence[int] = (0, 1), label: str = "",
                    extra_provenance: dict | None = None, **kw: Any) -> "Spectrum":
    """Generic (omega, a2F) table -> Spectrum. ``unit`` (meV, THz, cm-1, Ry, K, eV, ...) is
    REQUIRED and refers to the omega column; a2F is dimensionless (per unit energy of the same
    axis). ``rule`` as in ``grid_to_atoms``.
    """
    canonical_unit(unit)
    a = read_columns(path, cols)
    pos, mass = grid_to_atoms(a[:, 0], a[:, 1], rule)
    prov = {"source_file": str(path), "format": "two-column (omega, a2F)", "rule": rule,
            "columns": list(cols)}
    prov.update(extra_provenance or {})
    return atoms_to_spectrum(pos, mass, unit=unit, label=label, provenance=prov, **kw)


def spectrum_lambda(sp: "Spectrum") -> float:
    o = np.asarray(sp.omega_K)
    w = np.asarray(sp.weight_K)
    return float(2.0 * math.fsum(w / o))


def spectrum_omega_log(sp: "Spectrum") -> float:
    o = np.asarray(sp.omega_K)
    w = np.asarray(sp.weight_K)
    lam = 2.0 * math.fsum(w / o)
    if lam <= 0:
        return float("nan")
    return float(math.exp(2.0 * math.fsum(w / o * np.log(o)) / lam))
