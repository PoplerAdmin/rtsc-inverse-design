"""Out-of-model risk flags (ARCH section 6, stage 8; contract S6).

Verdicts are statements about the MODEL: isotropic Migdal-Eliashberg theory, linearised gap
equation, block-form Coulomb cutoff, with the given alpha^2F.  The flags below mark physics the
model leaves out.  They NEVER change a status or a bound (``apply_flags`` asserts this): a
REJECTED_PROVEN material with an ANHARMONIC_H flag is still "rejected within the model", and the
flag says the model may not describe that material well.

Codes (detail strings and the numbers used are kept in the record's ``validity`` section):
  SYNTHETIC_DATA         spectra are invented demo data (never a statement about a real material)
  NONADIABATIC           omega_max / E_F > nonadiabatic_ratio (Migdal's theorem needs << 1)
  ANHARMONIC_H           hydrogen fraction >= anharmonic_h_fraction (light H in shallow wells;
                         harmonic DFPT phonons and alpha^2F may be off), extra detail if pressure
                         >= anharmonic_pressure_GPa or soft modes are present
  SOFT_MODES             lowest real phonon frequency below soft_mode_THz (near an instability)
  MAGNETIC_PROXIMITY     magnetic elements present, or a small non-zero moment, or Stoner-like
                         N(E_F) above stoner_dos (spin fluctuations suppress pairing)
  MULTIBAND_ANISOTROPY   lattice anisotropy max(a,b,c)/min(a,b,c) >= anisotropy_ratio, or
                         n_fermi_sheets >= 2 (isotropic theory averages gap anisotropy away)
  MU_STAR_UNCERTAIN      mu_hi below mu_star_reference while d/f elements are present
                         (REJECTED_PROVEN is unaffected: it holds for every mu* >= 0)
  SPECTRUM_SPREAD        alternative spectra (smearings / q-grids) differ in A or lambda by more
                         than spectrum_spread_rel (convergence risk)
  STRONG_COUPLING        lambda > lambda_large (lattice-instability / polaronic regime)
  PARTIAL_OCCUPANCY      the input structure is disordered (computed on an ordered approximant)
  HYPOTHETICAL_STRUCTURE unrelaxed prototype substitution
"""
from __future__ import annotations

import dataclasses
import math
from dataclasses import dataclass, field
from typing import Any, Sequence

from ..sources.base import Candidate
from ..sources.composition import LANTHANIDES_ACTINIDES, MAGNETIC_ELEMENTS, TRANSITION_METALS, atomic_fractions
from ..spectrum import Spectrum, unit_factor
from ..status import Verdict

__all__ = ["ValidityConfig", "Flag", "compute_flags", "apply_flags"]

_K_PER_EV = unit_factor("eV")


@dataclass
class ValidityConfig:
    enabled: bool = True
    nonadiabatic_ratio: float = 0.1
    anharmonic_h_fraction: float = 0.25
    anharmonic_pressure_GPa: float = 50.0
    soft_mode_THz: float = 1.0
    magnetic_moment_warn_muB: float = 0.005
    stoner_dos: float = 1.5            # states/eV/atom/spin-ish; crude Stoner proxy
    anisotropy_ratio: float = 2.5
    mu_star_reference: float = 0.2
    spectrum_spread_rel: float = 0.1
    lambda_large: float = 3.0

    @classmethod
    def from_config(cls, cfg: dict[str, Any] | None) -> "ValidityConfig":
        cfg = dict(cfg or {})
        bad = set(cfg) - set(cls.__dataclass_fields__)
        if bad:
            raise ValueError(f"unknown [validity] keys: {sorted(bad)}")
        return cls(**cfg)


@dataclass
class Flag:
    code: str
    detail: str
    values: dict[str, Any] = field(default_factory=dict)

    def to_dict(self) -> dict[str, Any]:
        return dataclasses.asdict(self)


def _rel_spread(xs: list[float]) -> float:
    xs = [x for x in xs if math.isfinite(x)]
    if len(xs) < 2:
        return 0.0
    m = max(abs(x) for x in xs)
    return (max(xs) - min(xs)) / m if m > 0 else 0.0


def compute_flags(c: Candidate, spectra: Sequence[Spectrum], cfg: ValidityConfig | None = None,
                  *, mu_hi: float | None = None, phonon_min_THz: float | None = None
                  ) -> tuple[list[Flag], list[str]]:
    """Return (flags, unchecked) where ``unchecked`` lists risks with no data to assess."""
    cfg = cfg or ValidityConfig()
    flags: list[Flag] = []
    unchecked: list[str] = []
    if not cfg.enabled:
        return flags, ["validity flags disabled"]
    props = c.properties
    fr = atomic_fractions(c.composition)

    if c.synthetic or any(s.provenance.get("synthetic") for s in spectra):
        flags.append(Flag("SYNTHETIC_DATA", "invented demo spectra; not data for any material"))

    # --- non-adiabatic (Migdal) ratio ---------------------------------------------------------
    w_max = max((s.omega_max() for s in spectra), default=None)
    ef = props.get("fermi_energy_eV", props.get("electronic_bandwidth_eV"))
    if w_max is not None and ef:
        ratio = w_max / (float(ef) * _K_PER_EV)
        if ratio > cfg.nonadiabatic_ratio:
            flags.append(Flag("NONADIABATIC", f"omega_max/E_F = {ratio:.3g} > {cfg.nonadiabatic_ratio}",
                              {"omega_max_K": w_max, "E_F_eV": float(ef), "ratio": ratio}))
    else:
        unchecked.append("NONADIABATIC (no electronic energy scale fermi_energy_eV)")

    # --- phonon softness ---------------------------------------------------------------------
    pmin = phonon_min_THz
    if pmin is None and isinstance(props.get("phonon"), dict) and "min_frequency" in props["phonon"]:
        ph = props["phonon"]
        pmin = float(ph["min_frequency"]) * unit_factor(ph.get("unit", "THz")) / unit_factor("THz")
    soft = pmin is not None and pmin < cfg.soft_mode_THz
    if soft:
        flags.append(Flag("SOFT_MODES", f"lowest phonon frequency {pmin:.3g} THz < {cfg.soft_mode_THz} THz",
                          {"min_frequency_THz": pmin}))

    # --- anharmonic hydrogen -----------------------------------------------------------------
    xh = fr.get("H", 0.0)
    if xh >= cfg.anharmonic_h_fraction:
        p = float(props.get("target_pressure_GPa", c.pressure_GPa or 0.0))
        extra = []
        if p >= cfg.anharmonic_pressure_GPa:
            extra.append(f"pressure {p} GPa >= {cfg.anharmonic_pressure_GPa}")
        if soft:
            extra.append("soft modes present")
        flags.append(Flag("ANHARMONIC_H", f"H fraction {xh:.2f} >= {cfg.anharmonic_h_fraction}"
                          + (f"; {', '.join(extra)}" if extra else ""),
                          {"h_fraction": xh, "pressure_GPa": p}))

    # --- magnetism proximity -----------------------------------------------------------------
    mag_el = sorted(e for e in fr if e in MAGNETIC_ELEMENTS)
    m_pa = props.get("magnetization_per_atom")
    if m_pa is None and props.get("total_magnetization") is not None and c.n_sites:
        m_pa = float(props["total_magnetization"]) / float(props.get("nsites") or c.n_sites)
    dos = props.get("dos_at_fermi")
    reasons = []
    if mag_el:
        reasons.append(f"magnetic elements {mag_el}")
    if m_pa is not None and abs(float(m_pa)) > cfg.magnetic_moment_warn_muB:
        reasons.append(f"|m| = {abs(float(m_pa)):.3g} mu_B/atom")
    if dos is not None and float(dos) >= cfg.stoner_dos:
        reasons.append(f"N(E_F) = {float(dos):.3g} >= {cfg.stoner_dos} (Stoner proxy)")
    if reasons:
        flags.append(Flag("MAGNETIC_PROXIMITY", "; ".join(reasons),
                          {"magnetic_elements": mag_el, "m_per_atom": m_pa, "dos_at_fermi": dos}))

    # --- multiband / anisotropy --------------------------------------------------------------
    abc = c.lattice_lengths()
    sheets = props.get("n_fermi_sheets")
    reasons = []
    if abc and min(abc) > 0 and max(abc) / min(abc) >= cfg.anisotropy_ratio:
        reasons.append(f"lattice anisotropy {max(abc) / min(abc):.2f} >= {cfg.anisotropy_ratio}")
    if sheets is not None and int(sheets) >= 2:
        reasons.append(f"{int(sheets)} Fermi-surface sheets")
    if reasons:
        flags.append(Flag("MULTIBAND_ANISOTROPY", "; ".join(reasons), {"abc": abc, "n_fermi_sheets": sheets}))
    elif abc is None and sheets is None:
        unchecked.append("MULTIBAND_ANISOTROPY (no structure / Fermi-surface data)")

    # --- mu* range ---------------------------------------------------------------------------
    df = sorted(e for e in fr if e in TRANSITION_METALS or e in LANTHANIDES_ACTINIDES)
    if mu_hi is not None and df and mu_hi < cfg.mu_star_reference:
        flags.append(Flag("MU_STAR_UNCERTAIN",
                          f"mu_hi = {mu_hi} < {cfg.mu_star_reference} with d/f elements {df}",
                          {"mu_hi": mu_hi, "df_elements": df}))

    # --- spread across alternative spectra and coupling strength -----------------------------
    if len(spectra) >= 2:
        sa = _rel_spread([s.A() for s in spectra])
        sl = _rel_spread([s.lambda_() for s in spectra])
        if max(sa, sl) > cfg.spectrum_spread_rel:
            flags.append(Flag("SPECTRUM_SPREAD", f"relative spread A {sa:.3g}, lambda {sl:.3g} "
                              f"> {cfg.spectrum_spread_rel}", {"spread_A": sa, "spread_lambda": sl}))
    lam = max((s.lambda_() for s in spectra), default=None)
    if lam is not None and lam > cfg.lambda_large:
        flags.append(Flag("STRONG_COUPLING", f"lambda = {lam:.3g} > {cfg.lambda_large}", {"lambda": lam}))

    if props.get("partial_occupancy"):
        flags.append(Flag("PARTIAL_OCCUPANCY", "input structure is disordered"))
    if c.source == "prototype_substitution":
        flags.append(Flag("HYPOTHETICAL_STRUCTURE", "unrelaxed prototype substitution"))
    return flags, unchecked


def apply_flags(verdict: Verdict, flags: Sequence[Flag]) -> Verdict:
    """Return a copy of ``verdict`` with flag codes appended.  Status and bounds are untouched (S6)."""
    codes = list(verdict.flags)
    for f in flags:
        if f.code not in codes:
            codes.append(f.code)
    out = dataclasses.replace(verdict, flags=codes)
    assert (out.status, out.T_low_K, out.T_up_K, out.rule_ids) == \
        (verdict.status, verdict.T_low_K, verdict.T_up_K, verdict.rule_ids), "S6 violated"
    return out
