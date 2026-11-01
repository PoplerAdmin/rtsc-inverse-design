"""Band-resolved coupling data for the multiband rejection rule (wMB; Lean library MultiBand).

Convention (wMB Theorem 1 / section 7): band i's mass enhancement is lambda_i = sum_j lambda_ij with
lambda_ij = 2 int alpha^2F_ij(w)/w dw, detailed balance N_i alpha^2F_ij = N_j alpha^2F_ji, and

    A_ij := int_0^inf alpha^2F_ij(w) dw     (kelvin),   kappa = Perron root of A.

Rejection (team CERT, rule R-MULTI) needs an UPPER bound on kappa, so this module carries data of
three kinds (``kind``):

* ``"computed"``               A_ij computed from band-resolved spectra (best; any code).
* ``"entrywise_upper_bound"``  every A_ij is >= the true value, e.g. from lambda_ij and omega_max
                               via A_ij <= (omega_max/2) lambda_ij (wMB section 7 "Fallback").
* ``"row_sum_upper_bound"``    only upper bounds on the row sums A_i = sum_j A_ij are known
                               (``row_sum_upper_K``); ``A_K`` must then be None (forbidden)
                               (kappa <= max_i A_i, wMB section 4.7); e.g. from EPW's per-state
                               lambda_nk: the anisotropic corollary kappa <= max_k A_k (wMB section 6)
                               with A_k <= (omega_max/2) lambda_k.

The float helpers below (``kappa_float`` etc.) are informational; R-MULTI must re-evaluate with
rigorous arithmetic. ``rel_envelope`` e is an assumption: true A_ij <= (1+e) * stated A_ij.

EPW SOURCES CHECKED (github QEF/q-e, tag qe-6.7.0 = EPW 5.3, and develop):
EPW/src/supercond.f90 SUBROUTINE evaluate_a2f_lambda::

    lambda_eph = sum_{imode: wf > eps_acustic} g2(ik,iq,ibnd,jbnd,imode) / wf(imode,iq0)
    lambda_k(ik,ibnd) += wqf(iq) * w0g(jbnd, k+q) * lambda_eph      ; then lambda_k *= 2
    WRITE '(3f12.6, i8, f12.6, f24.15)' x1,x2,x3, ibnd, ekfs-ef0, lambda_k   -> <prefix>.lambda_FS
    WRITE '(2f12.7, 2i7)' ekfs-ef0, lambda_k, ik, ibnd                       -> <prefix>.lambda_aniso
          (iverbosity = 2 only; ik = irreducible index)
  and EPW/src/io_eliashberg.f90 (6.7.0; io/io_supercond.f90 in develop)::

    w0g(n, ik) = w0gauss((ekfs(n,ik) - ef0)/degaussw, 0)/degaussw  ! Gaussian FS weight
    n counts the bands INSIDE the Fermi window at that k (so 'ibnd' is NOT a global band label)
    wsphmax = 1.1d0 * MAXVAL(wf(:, :))   (read_frequencies; last point of <prefix>.a2f grid)

So EPW exposes only per-state lambda_nk (= sum over partner states, i.e. the row sum in the wMB
convention). It writes no per-state alpha^2F, no A_nk and no pair matrix: ``.lambda_pairs`` and
``.lambda_k_pairs`` are HISTOGRAMS of values, not matrices. The usable bound is therefore
kappa <= max_nk A_nk <= (omega_max/2) max_nk lambda_nk with omega_max >= every phonon frequency.
Caveat recorded in provenance: EPW drops modes with wf <= eps_acustic from lambda_nk.
"""
from __future__ import annotations

import json
import math
import re
import struct
from fractions import Fraction
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Callable, Sequence

import numpy as np

from .parsers import grid_to_atoms, read_columns, atoms_to_spectrum, ffloat
from .units import canonical_unit, to_kelvin_factor

try:
    from rtt.core.spectrum import Spectrum
except Exception:  # pragma: no cover
    Spectrum = None  # type: ignore[assignment]

KINDS = ("computed", "entrywise_upper_bound", "row_sum_upper_bound")


def _tuple2(m: Any) -> tuple[tuple[float, ...], ...]:
    return tuple(tuple(float(x) for x in row) for row in m)


@dataclass(frozen=True)
class BandCoupling:
    labels: tuple[str, ...]
    kind: str = "computed"
    A_K: tuple[tuple[float, ...], ...] | None = None          # A_ij in kelvin
    N: tuple[float, ...] | None = None                        # DOS weights N_i > 0 (ratios matter)
    row_sum_upper_K: tuple[float, ...] | None = None          # for kind == row_sum_upper_bound
    a2f_pairs: dict = field(default_factory=dict)             # "i,j" (0-based) -> Spectrum
    rel_envelope: float = 0.0
    label: str = ""
    provenance: dict = field(default_factory=dict, compare=False, hash=False)

    def __post_init__(self) -> None:
        labels = tuple(str(x) for x in self.labels)
        object.__setattr__(self, "labels", labels)
        B = len(labels)
        if B < 1:
            raise ValueError("BandCoupling needs at least one band")
        if self.kind not in KINDS:
            raise ValueError(f"kind must be one of {KINDS}, got {self.kind!r}")
        e = float(self.rel_envelope)
        if not (math.isfinite(e) and 0.0 <= e < 1.0):
            raise ValueError("rel_envelope must lie in [0, 1)")
        object.__setattr__(self, "rel_envelope", e)
        if self.A_K is not None:
            A = _tuple2(self.A_K)
            if len(A) != B or any(len(r) != B for r in A):
                raise ValueError(f"A_K must be {B}x{B}")
            for r in A:
                for x in r:
                    if not (math.isfinite(x) and x >= 0):
                        raise ValueError(f"A_ij must be finite and >= 0 (alpha^2F_ij are positive "
                                         f"measures, wMB Thm 1), got {x}")
            object.__setattr__(self, "A_K", A)
        if self.N is not None:
            N = tuple(float(x) for x in self.N)
            if len(N) != B or any(not (math.isfinite(x) and x > 0) for x in N):
                raise ValueError(f"N must be {B} finite positive weights")
            object.__setattr__(self, "N", N)
        if self.row_sum_upper_K is not None:
            R = tuple(float(x) for x in self.row_sum_upper_K)
            if len(R) != B or any(not (math.isfinite(x) and x >= 0) for x in R):
                raise ValueError(f"row_sum_upper_K must be {B} finite values >= 0")
            object.__setattr__(self, "row_sum_upper_K", R)
        if self.kind == "row_sum_upper_bound":
            if self.row_sum_upper_K is None:
                raise ValueError("kind 'row_sum_upper_bound' needs row_sum_upper_K")
            if self.A_K is not None:
                raise ValueError("kind 'row_sum_upper_bound' must not carry A_K (its meaning "
                                 "would be ambiguous); use kind 'computed' or "
                                 "'entrywise_upper_bound' with row_sum_upper_K instead")
        else:
            if self.A_K is None:
                raise ValueError(f"kind {self.kind!r} needs A_K")
        if self.A_K is not None and self.row_sum_upper_K is not None:
            for i, (r, u) in enumerate(zip(self.A_K, self.row_sum_upper_K)):
                # exact rational comparison of the stated doubles (no tolerance)
                if sum((Fraction(x) for x in r), Fraction(0)) > Fraction(u):
                    raise ValueError(f"row {i}: sum A_ij = {math.fsum(r)} exceeds row_sum_upper_K {u}")
        pairs = {}
        for key, sp in dict(self.a2f_pairs).items():
            i, j = (int(x) for x in str(key).split(","))
            if not (0 <= i < B and 0 <= j < B):
                raise ValueError(f"a2f_pairs key {key!r} out of range")
            pairs[f"{i},{j}"] = sp
            if self.A_K is not None and self.kind == "computed":
                a = math.fsum(sp.weight_K)
                if abs(a - self.A_K[i][j]) > 1e-9 * max(1.0, abs(a)):
                    raise ValueError(f"A_K[{i}][{j}] = {self.A_K[i][j]} but its alpha^2F has "
                                     f"total weight {a}")
        object.__setattr__(self, "a2f_pairs", pairs)

    # ------------------------------------------------------------- informational float helpers
    @property
    def n_bands(self) -> int:
        return len(self.labels)

    def row_sums_K(self) -> tuple[float, ...]:
        """A_i = sum_j A_ij (or the stated upper bounds when only those are known)."""
        if self.A_K is not None:
            return tuple(math.fsum(r) for r in self.A_K)
        return tuple(self.row_sum_upper_K or ())

    def max_row_sum_K(self) -> float:
        return max(self.row_sums_K())

    def detailed_balance_asymmetry(self) -> float | None:
        """||W - W^T||_F / ||W||_F with W_ij = N_i A_ij (wMB section 7 data-quality check)."""
        if self.A_K is None or self.N is None:
            return None
        W = np.asarray(self.N)[:, None] * np.asarray(self.A_K)
        n = np.linalg.norm(W)
        return float(np.linalg.norm(W - W.T) / n) if n > 0 else 0.0

    def kappa_float(self) -> float:
        """Float estimate of the bound used by R-MULTI (wMB section 7 steps 1-3):
        lambda_max(N^-1/2 sym(N A) N^-1/2) if A and N are known (>= rho(A)), else max row sum.
        Informational only; not a rigorous evaluation."""
        if self.A_K is not None and self.N is not None:
            N = np.asarray(self.N)
            W = N[:, None] * np.asarray(self.A_K)
            Wb = 0.5 * (W + W.T)
            s = 1.0 / np.sqrt(N)
            K = s[:, None] * Wb * s[None, :]
            return float(np.max(np.linalg.eigvalsh(K)))
        if self.A_K is not None:
            return float(max(abs(np.linalg.eigvals(np.asarray(self.A_K)))))
        return self.max_row_sum_K()

    # ------------------------------------------------------------- serialisation
    def to_dict(self) -> dict[str, Any]:
        return {"labels": list(self.labels), "kind": self.kind,
                "A_K": [list(r) for r in self.A_K] if self.A_K is not None else None,
                "N": list(self.N) if self.N is not None else None,
                "row_sum_upper_K": list(self.row_sum_upper_K) if self.row_sum_upper_K else None,
                "a2f_pairs": {k: v.to_dict() for k, v in self.a2f_pairs.items()},
                "rel_envelope": self.rel_envelope, "label": self.label,
                "provenance": dict(self.provenance), "unit": "K"}

    def to_json(self, **kw: Any) -> str:
        return json.dumps(self.to_dict(), **kw)

    @classmethod
    def from_dict(cls, d: dict[str, Any], base_dir: str | Path | None = None) -> "BandCoupling":
        """Read the user format (JSON or TOML table). Keys::

            labels = ["sigma", "pi"]                # required
            unit = "meV"                            # energy unit of A / row_sum_upper (default K)
            kind = "computed"                       # see KINDS (default computed)
            A = [[a11, a12], [a21, a22]]            # or A_K (already kelvin)
            N = [n1, n2]                            # optional DOS weights
            row_sum_upper = [r1, r2]                # or row_sum_upper_K
            rel_envelope = 0.1
            [[a2f_pairs]]                           # optional, per pair (0-based i, j)
            i = 0; j = 1; unit = "meV"; rule = "bins"
            omega = [...]; a2f = [...]              # or file = "pair01.dat" (two columns)
            # alternative upper-bound input:  lambda = [[...]] with omega_max = <in unit>
        """
        d = dict(d)
        unit = canonical_unit(d.get("unit", "K"))
        f = to_kelvin_factor(unit)
        labels = d["labels"]
        kind = d.get("kind")
        prov = dict(d.get("provenance", {}))
        prov.setdefault("input_unit", unit)
        A = d.get("A_K")
        if A is None and d.get("A") is not None:
            A = [[float(x) * f for x in r] for r in d["A"]]
        R = d.get("row_sum_upper_K")
        if R is None and d.get("row_sum_upper") is not None:
            R = [float(x) * f for x in d["row_sum_upper"]]
        if A is None and d.get("lambda") is not None:
            if d.get("omega_max") is None:
                raise ValueError("'lambda' input needs 'omega_max' (same unit)")
            wmax = float(d["omega_max"]) * f
            A = [[0.5 * wmax * float(x) for x in r] for r in d["lambda"]]
            kind = kind or "entrywise_upper_bound"
            if kind != "entrywise_upper_bound":
                raise ValueError("A from lambda and omega_max is an entrywise upper bound")
            prov["derived"] = "A_ij <= (omega_max/2) lambda_ij (wMB section 7 fallback)"
            prov["omega_max_K"] = wmax
        pairs: dict[str, Any] = {}
        for p in d.get("a2f_pairs", []) or []:
            pu = p.get("unit", unit)
            rule = p.get("rule", "bins")
            if "file" in p:
                path = Path(p["file"])
                if base_dir is not None and not path.is_absolute():
                    path = Path(base_dir) / path
                arr = read_columns(path, tuple(p.get("columns", (0, 1))))
                om, a2 = arr[:, 0], arr[:, 1]
                src = str(path)
            else:
                om, a2 = np.asarray(p["omega"], float), np.asarray(p["a2f"], float)
                src = "inline"
            pos, mass = grid_to_atoms(om, a2, rule)
            sp = atoms_to_spectrum(pos, mass, unit=pu, label=f"a2F_{p['i']}{p['j']}",
                                   provenance={"source_file": src, "rule": rule},
                                   negative_policy=p.get("negative_policy", "error"))
            pairs[f"{int(p['i'])},{int(p['j'])}"] = sp
        if A is None and pairs and len(pairs) == len(labels) ** 2:
            B = len(labels)
            A = [[math.fsum(pairs[f"{i},{j}"].weight_K) for j in range(B)] for i in range(B)]
        if kind is None:
            kind = "computed" if A is not None else "row_sum_upper_bound"
        return cls(tuple(labels), kind, A, d.get("N"), R, pairs, float(d.get("rel_envelope", 0.0)),
                   str(d.get("label", "")), prov)

    @classmethod
    def from_json(cls, s: str) -> "BandCoupling":
        return cls.from_dict(json.loads(s))

    @classmethod
    def from_file(cls, path: str | Path, table: str | None = None) -> "BandCoupling":
        """Read ``.json`` or ``.toml`` (optionally a sub-table, e.g. ``table='band_coupling'``)."""
        p = Path(path)
        if p.suffix.lower() == ".toml":
            import tomllib
            with open(p, "rb") as fh:
                d = tomllib.load(fh)
        else:
            d = json.loads(p.read_text())
        if table:
            for k in table.split("."):
                d = d[k]
        d = dict(d)
        d.setdefault("provenance", {})
        d["provenance"] = dict(d["provenance"], source_file=str(p))
        return cls.from_dict(d, base_dir=p.parent)


# ------------------------------------------------------------------------------------------
# EPW anisotropic outputs
# ------------------------------------------------------------------------------------------
@dataclass
class EPWLambdaStates:
    path: str
    k: np.ndarray            # (n, 3) cartesian, units 2pi/alat (lambda_FS) or empty (lambda_aniso)
    band: np.ndarray         # (n,) EPW index of the band INSIDE the Fermi window at that k
    e_eV: np.ndarray         # (n,) E_nk - E_F (eV)
    lam: np.ndarray          # (n,) lambda_nk
    kindex: np.ndarray | None = None   # lambda_aniso only: irreducible k index
    full_grid: bool = True   # lambda_FS lists the full nkf grid (uniform k weights)

    def fs_weights(self, degaussw_eV: float) -> np.ndarray:
        """EPW's Fermi-surface weight w0g = exp(-(E/s)^2)/(sqrt(pi) s), s = degaussw (eV)."""
        x = self.e_eV / degaussw_eV
        return np.exp(-x * x) / (math.sqrt(math.pi) * degaussw_eV)


def parse_epw_lambda_fs(path: str | Path) -> EPWLambdaStates:
    """<prefix>.lambda_FS (or <prefix>.lambda): rows (kx, ky, kz, band, E-Ef[eV], lambda_nk)."""
    rows = []
    for line in Path(path).read_text(errors="replace").splitlines():
        s = line.strip()
        if not s or s.startswith("#"):
            continue
        t = s.split()
        if len(t) != 6:
            raise ValueError(f"{path}: expected 6 columns, got {line!r}")
        rows.append([ffloat(x) for x in t])
    if not rows:
        raise ValueError(f"{path}: no states")
    a = np.asarray(rows)
    return EPWLambdaStates(str(path), a[:, :3], a[:, 3].astype(int), a[:, 4], a[:, 5])


def parse_epw_lambda_aniso(path: str | Path) -> EPWLambdaStates:
    """<prefix>.lambda_aniso (iverbosity=2): rows (E-Ef[eV], lambda_nk, ik_irr, band).
    Irreducible k only and no multiplicities: usable for max_nk, not for DOS averages."""
    rows = []
    for line in Path(path).read_text(errors="replace").splitlines():
        s = line.strip()
        if not s or s.startswith("#"):
            continue
        rows.append([ffloat(x) for x in s.split()[:4]])
    a = np.asarray(rows)
    return EPWLambdaStates(str(path), np.zeros((len(a), 0)), a[:, 3].astype(int), a[:, 0],
                           a[:, 1], kindex=a[:, 2].astype(int), full_grid=False)


def epw_omega_max_meV(a2f_path: str | Path) -> float:
    """Last frequency of an EPW <prefix>.a2f / .a2f_iso grid (meV). When EPW computed the grid
    itself this is wsphmax = 1.1 * max phonon frequency on the fine q grid (io_eliashberg.f90,
    read_frequencies), so it bounds every frequency entering lambda_nk."""
    w = None
    for line in Path(a2f_path).read_text(errors="replace").splitlines():
        t = line.split()
        if not t:
            continue
        try:
            v = [float(x) for x in t]
        except ValueError:
            if w is not None and not line.strip()[0].isdigit():
                break
            continue
        w = v[0]
    if w is None:
        raise ValueError(f"{a2f_path}: no frequency rows")
    return float(w)


def band_coupling_from_epw_lambda(states: EPWLambdaStates, *, omega_max: float,
                                  omega_unit: str = "meV",
                                  groups: Callable[[np.ndarray, int, float], str] | None = None,
                                  reduction: str = "max", degaussw_eV: float | None = None,
                                  eps_acustic_cm1: float | None = None,
                                  freq_file: str | Path | None = None,
                                  strict_modes: bool = True,
                                  label: str = "") -> BandCoupling:
    """Row-sum upper bounds A_G <= (omega_max/2) * lambda_G from EPW per-state lambda_nk.

    reduction="max" (default): lambda_G = max over states in group G. Valid for the discretised
      anisotropic model (wMB section 6: kappa <= max_k A_k); the group split is cosmetic.
    reduction="dos_average": lambda_G = sum_k N_k lambda_k / sum_k N_k with EPW's Gaussian FS
      weights (needs degaussw_eV and the full-grid lambda_FS). This is the row sum of the
      BAND-ISOTROPIC model built from these states -- a different model; recorded as assumption
      'M-BAND-ISO'.
    ``groups(k_cart, band, e_eV) -> name`` assigns states to sheets (default: one group). EPW's
    band index counts bands inside the Fermi window at each k, so it is not a global label.

    SOUNDNESS CONDITIONS (recorded in provenance["assumptions"]; verified when possible):
      EPW-OMEGA-MAX     omega_max >= every phonon frequency of the fine q grid entering lambda_nk;
      EPW-NO-IMAGINARY  no imaginary modes (EPW stores them as wf < 0 and drops them);
      EPW-NO-EXCLUDED   no modes with 0 < wf <= eps_acustic besides the q=Gamma acoustic modes
                        (EPW drops all wf <= eps_acustic from lambda_nk).
    ``freq_file`` = <outdir>/<prefix>.ephmat/freq (written by EPW with ephwrite=.true.) lets us
    check all three. With ``strict_modes=True`` (default) any violation raises ValueError; the
    verified conditions are then marked '(verified: <file>)'. Without ``freq_file`` they stay
    unverified assumptions.
    """
    wmax_K = float(omega_max) * to_kelvin_factor(omega_unit)
    if not (wmax_K > 0 and math.isfinite(wmax_K)):
        raise ValueError("omega_max must be > 0")
    lam = states.lam
    if not np.all(np.isfinite(lam)) or np.any(lam < 0):
        raise ValueError(f"{states.path}: lambda_nk must be finite and >= 0")
    n = len(lam)
    if groups is None:
        names = np.array(["all_states"] * n, dtype=object)
    else:
        kk = states.k if states.k.shape[1] == 3 else np.full((n, 3), np.nan)
        names = np.array([str(groups(kk[i], int(states.band[i]), float(states.e_eV[i])))
                          for i in range(n)], dtype=object)
    order = list(dict.fromkeys(names.tolist()))
    prov: dict[str, Any] = {
        "source_file": states.path, "format": "EPW per-state lambda_nk (supercond.f90)",
        "omega_max_K": wmax_K, "omega_max_input": [omega_max, omega_unit],
        "bound": "A_G <= (omega_max/2) * lambda_G (wMB section 7 fallback, row sums)",
        "n_states": n, "reduction": reduction,
        "caveats": ["EPW excludes modes with wf <= eps_acustic from lambda_nk"
                    + (f" (eps_acustic = {eps_acustic_cm1} cm^-1)" if eps_acustic_cm1 else ""),
                    "omega_max must bound every phonon frequency on the fine q grid"],
        "assumptions": [],
    }
    conds = {
        "EPW-OMEGA-MAX": f"EPW-OMEGA-MAX({wmax_K:.6g}K): omega_max bounds every fine-grid phonon "
                         "frequency entering lambda_nk",
        "EPW-NO-IMAGINARY": "EPW-NO-IMAGINARY: no imaginary modes on the fine q grid (EPW drops "
                            "wf < 0 from lambda_nk; the bound covers only real modes)",
        "EPW-NO-EXCLUDED": "EPW-NO-EXCLUDED: no modes with 0 < wf <= eps_acustic other than the "
                           "q=Gamma acoustic modes (EPW drops them from lambda_nk)",
    }
    if freq_file is not None:
        if eps_acustic_cm1 is None:
            raise ValueError("freq_file check needs eps_acustic_cm1 (EPW input eps_acustic)")
        rep = check_epw_frequencies(read_epw_freq(freq_file), omega_max_K=wmax_K,
                                    eps_acustic_cm1=eps_acustic_cm1)
        prov["mode_check"] = rep
        problems = []
        if rep["max_freq_K"] > wmax_K:
            problems.append(f"max fine-grid frequency {rep['max_freq_K']:.6g} K > omega_max")
        else:
            conds["EPW-OMEGA-MAX"] += f" (verified: {freq_file})"
        if rep["n_imaginary"]:
            problems.append(f"{rep['n_imaginary']} imaginary modes (min {rep['min_freq_K']:.4g} K)")
        else:
            conds["EPW-NO-IMAGINARY"] += f" (verified: {freq_file})"
        if rep["n_excluded_nongamma"]:
            problems.append(f"{rep['n_excluded_nongamma']} modes with 0 < wf <= eps_acustic "
                            "away from Gamma")
        else:
            conds["EPW-NO-EXCLUDED"] += f" (verified: {freq_file})"
        if problems and strict_modes:
            raise ValueError(f"{freq_file}: EPW data refused for the lambda_nk bound: "
                             + "; ".join(problems))
        prov["mode_problems"] = problems
    prov["assumptions"].extend(conds.values())
    lamG, NG = [], []
    if reduction == "max":
        for g in order:
            lamG.append(float(np.max(lam[names == g])))
        N = None
    elif reduction == "dos_average":
        if degaussw_eV is None or not states.full_grid:
            raise ValueError("dos_average needs degaussw_eV and a full-grid lambda_FS file")
        w = states.fs_weights(degaussw_eV)
        for g in order:
            m = names == g
            NG.append(float(np.sum(w[m])))
            lamG.append(float(np.sum(w[m] * lam[m]) / np.sum(w[m])))
        N = tuple(x / sum(NG) for x in NG)
        prov["assumptions"].append("M-BAND-ISO: isotropic within each group (band-averaged model)")
        prov["degaussw_eV"] = degaussw_eV
    else:
        raise ValueError("reduction must be 'max' or 'dos_average'")
    prov["lambda_G"] = dict(zip(order, lamG))
    R = tuple(0.5 * wmax_K * x for x in lamG)
    return BandCoupling(tuple(order), "row_sum_upper_bound", None, N, R, {}, 0.0, label, prov)


def read_epw_freq(path: str | Path) -> dict[str, Any]:
    """Read EPW's unformatted <prefix>.ephmat/freq (io_eliashberg.f90 write_ephmat, qe-6.7.0):
    record 1 = (nqtotf, nqf1, nqf2, nqf3, nmodes) int32; then per selected q: record xqf(3)
    float64 and nmodes records wf (float64, Ry; imaginary modes negative, ephwann_shuffle.f90).
    Sequential gfortran records with 4-byte little-endian length markers (checked)."""
    data = Path(path).read_bytes()
    pos = 0
    recs: list[bytes] = []
    while pos < len(data):
        if pos + 4 > len(data):
            raise ValueError(f"{path}: truncated record marker")
        (n,) = struct.unpack_from("<i", data, pos)
        if n < 0 or pos + 8 + n > len(data):
            raise ValueError(f"{path}: bad record length {n} (not a gfortran sequential file?)")
        body = data[pos + 4:pos + 4 + n]
        (n2,) = struct.unpack_from("<i", data, pos + 4 + n)
        if n2 != n:
            raise ValueError(f"{path}: record markers disagree ({n} vs {n2})")
        recs.append(body)
        pos += 8 + n
    if not recs or len(recs[0]) != 20:
        raise ValueError(f"{path}: unexpected header record")
    nqtotf, nq1, nq2, nq3, nmodes = struct.unpack("<5i", recs[0])
    body = recs[1:]
    if len(body) % (1 + nmodes):
        raise ValueError(f"{path}: {len(body)} records is not a multiple of 1+nmodes")
    xq, wf = [], []
    for i in range(0, len(body), 1 + nmodes):
        if len(body[i]) != 24 or any(len(r) != 8 for r in body[i + 1:i + 1 + nmodes]):
            raise ValueError(f"{path}: unexpected record sizes")
        xq.append(struct.unpack("<3d", body[i]))
        wf.append([struct.unpack("<d", r)[0] for r in body[i + 1:i + 1 + nmodes]])
    return {"nqtotf": nqtotf, "nqf": (nq1, nq2, nq3), "nmodes": nmodes,
            "xq": np.asarray(xq), "wf_Ry": np.asarray(wf)}


def check_epw_frequencies(freq: dict[str, Any], *, omega_max_K: float,
                          eps_acustic_cm1: float, gamma_tol: float = 1e-8) -> dict[str, Any]:
    """Summary of the fine-grid modes behind lambda_nk: max/min frequency (K), imaginary modes
    (wf < 0, excluding |wf| < 1e-12 Ry round-off), and modes 0 <= wf <= eps_acustic away from
    q = Gamma (EPW drops them)."""
    from .units import RY_TO_K, CMM1_TO_K
    wf = freq["wf_Ry"]
    eps_K = eps_acustic_cm1 * CMM1_TO_K
    wK = wf * RY_TO_K
    gamma = np.all(np.abs(freq["xq"]) < gamma_tol, axis=1)
    imag = wf < -1e-12
    small = (wK <= eps_K) & ~imag
    return {"n_q_selected": int(wf.shape[0]), "nmodes": int(wf.shape[1]),
            "max_freq_K": float(wK.max()), "min_freq_K": float(wK.min()),
            "omega_max_K": float(omega_max_K), "eps_acustic_K": eps_K,
            "n_imaginary": int(np.count_nonzero(imag)),
            "n_excluded_gamma": int(np.count_nonzero(small[gamma])),
            "n_excluded_nongamma": int(np.count_nonzero(small[~gamma]))}


def epw_lambda_dos_average(states: EPWLambdaStates, degaussw_eV: float) -> float:
    """sum_k w0g_k lambda_k / sum_k w0g_k over a full-grid lambda_FS: compare with EPW's
    isotropic lambda ('Summed el-ph coupling' in <prefix>.a2f) as a parse/normalisation check."""
    w = states.fs_weights(degaussw_eV)
    return float(np.sum(w * states.lam) / np.sum(w))


def parse_epw_summed_lambda(a2f_path: str | Path) -> float | None:
    m = re.search(r"Summed el-ph coupling\s+(\S+)", Path(a2f_path).read_text(errors="replace"))
    return ffloat(m.group(1)) if m else None


def band_coupling_from_spectra(labels: Sequence[str], pairs: dict[tuple[int, int], Any],
                               N: Sequence[float] | None = None, label: str = "",
                               provenance: dict | None = None,
                               rel_envelope: float = 0.0) -> BandCoupling:
    """Build a 'computed' BandCoupling from per-pair Spectra alpha^2F_ij (kelvin); missing pairs
    are zero couplings."""
    B = len(labels)
    A = [[0.0] * B for _ in range(B)]
    sp = {}
    for (i, j), s in pairs.items():
        A[i][j] = math.fsum(s.weight_K)
        sp[f"{i},{j}"] = s
    return BandCoupling(tuple(labels), "computed", A, tuple(N) if N is not None else None, None,
                        sp, rel_envelope, label, dict(provenance or {}))
