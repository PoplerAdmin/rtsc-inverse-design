"""DFT/DFPT settings profiles (architecture §7.2, §3.1 S4-S7).

A :class:`Profile` fixes every numerical setting of one tier: plane-wave cutoff, k densities,
q mesh(es), the e-ph smearing ladder and its admission rule, relax settings and the scratch
pruning policy. Four profiles exist:

* ``s4``       pw.x vc-relax -> scf -> ELF -> nscf (denser k) -> ILDOS -> spin-polarised scf
               (S4 dft-electronic, :mod:`rtt.dft.electronic`).
* ``coarse``   interpolated e-ph DFPT at q 3^3-4^3 (S5; doc 11 §3.1 measured settings).
* ``full``     interpolated e-ph DFPT on two q meshes (S6).
* ``ref_hull`` vc-relax + scf at the converged cutoff for own-DFT reference hulls (S7a).

Pseudopotentials: **SG15 ONCV PBE-1.2 norm-conserving files only** (`<El>_ONCV_PBE-1.2.upf`,
scalar-relativistic). QE 7.6 segfaults in e-ph with ultrasoft/PAW (architecture W4), so
:func:`sg15_pseudos` refuses any file whose UPF header is not ``pseudo_type="NC"`` and any
fully-relativistic (``_FR``) file. ecutrho = 4 x ecutwfc (norm-conserving).

Numbers and their sources:

* coarse: ecut 60 Ry (~2/3 of the converged 90 Ry of the H3S reference), k_scf spacing
  0.25 1/A and k_fit spacing 0.125 1/A (H3S at 200 GPa: k 12^3 / k_fit 24^3, i.e. the
  "k_fit >= 24^3-equivalent density" of §3.1), q from tiers.toml (default 4^3), el_ph_sigma
  0.01 Ry x 10 (doc 11 §3.1 ladder 0.01-0.10 Ry), plateau sigma >= 0.03 Ry (tiers.toml
  ``S5.plateau_sigma_min_Ry``). Source: doc 11 §3.1 (measured on H3S only).
* full: ecut 90 Ry, k spacing 0.17 1/A (H3S reference k 18^3), fit factor 2 (k_fit 36^3),
  el_ph_sigma 0.005 Ry x 10 (rtt.core.dfpt default), q meshes from tiers.toml ``S6.q_meshes``.
* s4: ecut 60 Ry, k spacing 0.25 1/A (scf), 0.15 1/A (nscf, for gap / n_sheets / N(E_F)),
  degauss 0.02 Ry Marzari-Vanderbilt (rtt.core.dfpt default).
* ref_hull: ecut 90 Ry, k spacing 0.20 1/A.

Every number may be overridden from ``tiers.toml`` (table ``[profiles.<name>]``, keys equal to
the :class:`Profile` field names) through :func:`profiles_from_tiers`; overrides are recorded in
``Profile.overrides`` and so in ``tier_runs.params_json``.
"""
from __future__ import annotations

import hashlib
import json
import os
import re
from dataclasses import asdict, dataclass, field, replace
from pathlib import Path
from typing import Any

QE_BIN_DEFAULT = "~/micromamba/envs/qe/bin"
PSEUDO_DIR_DEFAULT = "~/pseudo/sg15"
SG15_PATTERN = "{el}_ONCV_PBE-1.2.upf"


def qe_bin() -> Path:
    return Path(os.environ.get("RTT_QE_BIN", QE_BIN_DEFAULT)).expanduser()


def pseudo_dir() -> Path:
    return Path(os.environ.get("RTT_PSEUDO_DIR", PSEUDO_DIR_DEFAULT)).expanduser()


class PseudoError(ValueError):
    """An element has no admissible (SG15 NC scalar-relativistic) pseudopotential."""


@dataclass(frozen=True)
class Profile:
    name: str
    kind: str                                  # "electronic" | "dfpt" | "hull"
    spectrum_tier: str | None = None           # "coarse" | "full" for DFPT profiles
    # plane waves / electrons
    ecutwfc_Ry: float = 60.0
    ecutrho_factor: float = 4.0                # NC pseudopotentials
    smearing: str = "mv"
    degauss_Ry: float = 0.02
    conv_thr: float = 1e-10
    mixing_beta: float = 0.5
    # relax
    do_relax: bool = True
    press_conv_thr_kbar: float = 0.5
    etot_conv_thr: float = 1e-5
    forc_conv_thr: float = 1e-4
    relax_nstep: int = 100
    cell_dofree: str = "all"
    # k grids (spacing in 1/A including 2 pi; rtt.core.dfpt.grids convention)
    k_spacing: float = 0.25
    k_fit_spacing: float | None = None
    fit_factor: int = 2
    nscf_k_spacing: float | None = None        # s4 only
    matdyn_factor: int = 4                     # k_matdyn = matdyn_factor * q (rtt.core.dfpt default)
    # DFPT
    q: tuple[int, int, int] | None = None
    q_meshes: tuple[tuple[int, int, int], ...] = ()
    el_ph_sigma_Ry: float = 0.005
    el_ph_nsigma: int = 10
    smearing_admission: str = "mode_positive"
    tr2_ph: float = 1e-14
    plateau_sigma_min_Ry: float | None = None
    # electronic (s4)
    starting_magnetization: float = 0.5        # ferromagnetic start for the spin check
    dos_sigma_eV: float = 0.136                # Gaussian width for N(E_F) (0.01 Ry)
    ildos_window_eV: float = 0.136             # E_F +- window for the ILDOS (Belli H_DOS)
    # scratch policy
    prune_per_q: bool = True                   # hash-then-delete wfc of finished q dirs
    prune_dvscf_per_q: bool = False            # also dvscf of finished q (breaks EPW reuse)
    prune_final: str = "wfc"                   # "none" | "wfc" | "tmp" after success
    timeouts_s: dict = field(default_factory=dict)
    # pseudopotentials
    pseudo_pattern: str = SG15_PATTERN
    source_id: str = ""
    overrides: tuple[tuple[str, Any], ...] = ()

    # ------------------------------------------------------------------ helpers
    @property
    def ecutrho_Ry(self) -> float:
        return float(self.ecutwfc_Ry) * float(self.ecutrho_factor)

    def to_dict(self) -> dict[str, Any]:
        d = asdict(self)
        d["overrides"] = [list(x) for x in self.overrides]
        return d

    def sha256(self) -> str:
        return hashlib.sha256(json.dumps(self.to_dict(), sort_keys=True, default=list)
                              .encode()).hexdigest()

    def with_overrides(self, **kw: Any) -> "Profile":
        known = set(self.__dataclass_fields__) - {"name", "kind", "overrides"}
        bad = set(kw) - known
        if bad:
            raise ValueError(f"profile {self.name}: unknown override keys {sorted(bad)}")
        norm: dict[str, Any] = {}
        for k, v in kw.items():
            if k == "q" and v is not None:
                v = _triple(v)
            elif k == "q_meshes":
                v = tuple(_triple(x) for x in v)
            elif k == "timeouts_s":
                v = dict(v)
            norm[k] = v
        ov = dict(self.overrides)
        ov.update(norm)
        return replace(self, **norm, overrides=tuple(sorted(ov.items(), key=lambda t: t[0])))

    def eph_settings(self, *, q: tuple[int, int, int] | None = None, pressure_GPa: float = 0.0,
                     do_relax: bool | None = None) -> Any:
        """Core :class:`EPhSettings` for one DFPT run of this profile."""
        from rtt.core.dfpt.workflow import EPhSettings
        if self.kind != "dfpt":
            raise ValueError(f"profile {self.name} is not a DFPT profile")
        qq = _triple(q) if q is not None else self.q
        if qq is None:
            raise ValueError(f"profile {self.name}: no q mesh given")
        return EPhSettings(
            pressure_GPa=float(pressure_GPa),
            do_relax=self.do_relax if do_relax is None else bool(do_relax),
            ecutwfc_Ry=float(self.ecutwfc_Ry), ecutrho_Ry=self.ecutrho_Ry,
            smearing=self.smearing, degauss_Ry=float(self.degauss_Ry), conv_thr=self.conv_thr,
            mixing_beta=self.mixing_beta, press_conv_thr_kbar=self.press_conv_thr_kbar,
            etot_conv_thr=self.etot_conv_thr, forc_conv_thr=self.forc_conv_thr,
            relax_nstep=self.relax_nstep, cell_dofree=self.cell_dofree,
            k_spacing=self.k_spacing, fit_factor=self.fit_factor,
            k_fit_spacing=self.k_fit_spacing, q=qq,
            k_matdyn=tuple(self.matdyn_factor * x for x in qq),
            el_ph_sigma_Ry=self.el_ph_sigma_Ry, el_ph_nsigma=self.el_ph_nsigma,
            smearing_admission=self.smearing_admission, tr2_ph=self.tr2_ph,
            timeouts_s=dict(self.timeouts_s))


_ELEMENT = re.compile(r"^[A-Z][a-z]?$")


def ecut_floor_Ry(elements: Any, floors: Any) -> float | None:
    """Largest per-element wavefunction-cutoff floor (Ry) among ``elements``, or None.

    ``floors`` is the tiers.toml table ``S5.ecut_floor_Ry`` (``{H = 80.0, source = ...}``) or a
    plain ``{element: Ry}`` dict; non-element keys (source, provisional) are ignored. Doc 16 §3:
    SG15 H is not converged below 80 Ry (doc 15 §1), so any H-containing coarse run uses >= 80."""
    if not floors:
        return None
    els = {str(e) for e in (elements or ())}
    vals = [float(v) for k, v in dict(floors).items()
            if _ELEMENT.match(str(k)) and k in els and isinstance(v, (int, float))
            and not isinstance(v, bool)]
    return max(vals) if vals else None


def apply_ecut_floor(profile: "Profile", elements: Any, floors: Any) -> "Profile":
    """``profile`` with ecutwfc raised to :func:`ecut_floor_Ry` when below it (recorded in
    ``Profile.overrides``, hence in the run's params); unchanged otherwise. Never lowers."""
    f = ecut_floor_Ry(elements, floors)
    if f is None or float(profile.ecutwfc_Ry) >= f:
        return profile
    return profile.with_overrides(ecutwfc_Ry=f)


def _triple(v: Any) -> tuple[int, int, int]:
    if isinstance(v, int):
        return (v, v, v)
    t = tuple(int(x) for x in v)
    if len(t) != 3 or min(t) < 1:
        raise ValueError(f"bad grid {v!r}")
    return t  # type: ignore[return-value]


S4 = Profile(name="s4", kind="electronic", ecutwfc_Ry=60.0, k_spacing=0.25,
             nscf_k_spacing=0.15, do_relax=True, prune_per_q=False, prune_final="wfc",
             source_id="doc11#3.3 (Belli chain, H3S 60 Ry); arch#3.1 S4")
COARSE = Profile(name="coarse", kind="dfpt", spectrum_tier="coarse", ecutwfc_Ry=60.0,
                 k_spacing=0.25, k_fit_spacing=0.125, q=(4, 4, 4), q_meshes=((3, 3, 3), (4, 4, 4)), el_ph_sigma_Ry=0.01,
                 el_ph_nsigma=10, tr2_ph=1e-14, plateau_sigma_min_Ry=0.03, do_relax=True,
                 prune_per_q=True, prune_dvscf_per_q=True, prune_final="wfc",
                 source_id="doc11#3.1 (H3S c2/c3 settings); tiers.toml [S5]")
FULL = Profile(name="full", kind="dfpt", spectrum_tier="full", ecutwfc_Ry=90.0,
               k_spacing=0.17, fit_factor=2, q=(4, 4, 4), q_meshes=((4, 4, 4), (6, 6, 6)),
               el_ph_sigma_Ry=0.005, el_ph_nsigma=10, tr2_ph=1e-16, plateau_sigma_min_Ry=None,
               do_relax=True, prune_per_q=True, prune_dvscf_per_q=False, prune_final="wfc",
               source_id="rtt.core.dfpt defaults; H3S reference (90 Ry, k 18/36); tiers.toml [S6]")
REF_HULL = Profile(name="ref_hull", kind="hull", ecutwfc_Ry=90.0, k_spacing=0.20, do_relax=True,
                   prune_per_q=False, prune_final="wfc",
                   source_id="arch#3.1 S7a; same cutoff as full")

PROFILES: dict[str, Profile] = {p.name: p for p in (S4, COARSE, FULL, REF_HULL)}


def profiles_from_tiers(tiers: Any | None) -> dict[str, Profile]:
    """Apply tiers.toml to the defaults: ``[S5] q``, ``[S5] plateau_sigma_min_Ry``,
    ``[S6] q_meshes`` and any ``[profiles.<name>]`` keys. ``tiers`` is a
    :class:`rtt.config.LoadedConfig` or a plain dict."""
    out = dict(PROFILES)
    if tiers is None:
        return out
    raw = getattr(tiers, "raw", tiers)

    def val(x: Any) -> Any:
        return x["value"] if isinstance(x, dict) and "value" in x else x

    s5 = raw.get("S5", {}) or {}
    ov: dict[str, Any] = {}
    if "q" in s5:
        ov["q"] = val(s5["q"])
    if "q_meshes" in s5:
        ov["q_meshes"] = val(s5["q_meshes"])
        ov.setdefault("q", ov["q_meshes"][-1])
    if "plateau_sigma_min_Ry" in s5:
        ov["plateau_sigma_min_Ry"] = float(val(s5["plateau_sigma_min_Ry"]))
    if ov:
        out["coarse"] = out["coarse"].with_overrides(**ov)
    s6 = raw.get("S6", {}) or {}
    if "q_meshes" in s6:
        qm = val(s6["q_meshes"])
        out["full"] = out["full"].with_overrides(q_meshes=qm, q=qm[-1])
    for name, kv in (raw.get("profiles", {}) or {}).items():
        if name not in out:
            raise ValueError(f"tiers.toml [profiles.{name}]: unknown profile")
        kv = {k: val(v) for k, v in kv.items() if k not in ("source", "source_id")}
        out[name] = out[name].with_overrides(**kv)
    return out


# ---------------------------------------------------------------------- pseudopotentials
_PSEUDO_TYPE = re.compile(r'pseudo_type\s*=\s*"\s*([A-Za-z0-9-]+)\s*"')
_REL = re.compile(r'relativistic\s*=\s*"\s*([A-Za-z-]+)\s*"')


def check_nc_upf(path: Path) -> dict[str, Any]:
    """Header check: norm-conserving and not fully relativistic. Returns header facts."""
    head = Path(path).read_text(errors="replace")[:20000]
    m = _PSEUDO_TYPE.search(head)
    ptype = m.group(1).upper() if m else None
    r = _REL.search(head)
    rel = r.group(1).lower() if r else None
    if ptype not in ("NC", "SL"):
        raise PseudoError(f"{Path(path).name}: pseudo_type {ptype!r} is not norm-conserving "
                          "(QE 7.6 e-ph needs NC; architecture W4)")
    if rel == "full" or "_FR" in Path(path).name:
        raise PseudoError(f"{Path(path).name}: fully relativistic pseudopotential not allowed")
    return {"pseudo_type": ptype, "relativistic": rel}


def available_elements(directory: str | os.PathLike | None = None,
                       pattern: str = SG15_PATTERN) -> set[str]:
    """Elements with an SG15 file present (the S1 whitelist, tiers.toml S1.pseudo_dir)."""
    d = Path(directory).expanduser() if directory else pseudo_dir()
    rx = re.compile("^" + re.escape(pattern).replace(r"\{el\}", "([A-Z][a-z]?)") + "$")
    out = set()
    if d.is_dir():
        for p in d.iterdir():
            m = rx.match(p.name)
            if m:
                out.add(m.group(1))
    return out


def sg15_pseudos(elements: list[str] | tuple[str, ...], directory: str | os.PathLike | None = None,
                 pattern: str = SG15_PATTERN) -> Any:
    """Core :class:`PseudoTable` mapping each element to its SG15 NC file, after checking
    that the file exists and is norm-conserving. Raises :class:`PseudoError` otherwise."""
    from rtt.core.dfpt.pseudo import PseudoTable
    d = Path(directory).expanduser() if directory else pseudo_dir()
    files = {}
    for el in elements:
        f = pattern.format(el=el)
        p = d / f
        if not p.is_file():
            raise PseudoError(f"no SG15 pseudopotential for {el}: {f} not in the pseudo dir")
        check_nc_upf(p)
        files[el] = f
    return PseudoTable(str(d), files)


def pseudo_md5s(table: Any, elements: list[str]) -> dict[str, str]:
    out = {}
    for el in elements:
        h = hashlib.md5()
        with open(table.path(el), "rb") as fh:
            for chunk in iter(lambda: fh.read(1 << 20), b""):
                h.update(chunk)
        out[el] = h.hexdigest()
    return out


# ---------------------------------------------------------------------- thresholds
#: keys accepted in place of ``value`` in inline TOML tables (targets.toml [hard] uses them)
NAMED_VALUE_KEYS = ("gap_eV_max", "moment_muB_per_atom_tol", "moment_muB_tol", "window_eV",
                    "value_K", "value_eV")


def th_parts(x: Any) -> tuple[Any, str | None, bool]:
    """(value, source_id, provisional) from an :class:`rtt.config.Threshold`, an inline TOML
    table (``{value=..., source=..., provisional=...}`` or one named key such as
    ``{gap_eV_max=0.0, source=...}``), or a bare number (source None).

    Raises ``ValueError`` for a table that has neither ``value`` nor exactly one named key, so a
    threshold can never silently become None (review RC1 F1)."""
    if x is None:
        return None, None, False
    if hasattr(x, "value") and hasattr(x, "source_id"):
        return x.value, x.source_id, bool(getattr(x, "provisional", False))
    if isinstance(x, dict):
        src = x.get("source_id") or x.get("source")
        prov = bool(x.get("provisional", False))
        if "value" in x:
            return x["value"], src, prov
        named = [k for k in NAMED_VALUE_KEYS if k in x]
        if len(named) != 1:
            raise ValueError(f"threshold table without 'value' or a single named key: {sorted(x)}")
        return x[named[0]], src, prov
    if isinstance(x, (int, float)) and not isinstance(x, bool):
        return x, None, False
    raise ValueError(f"unusable threshold {x!r}")
