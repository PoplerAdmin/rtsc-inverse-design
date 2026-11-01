"""Single entry point for user-attached alpha^2F files: parse + completeness + consistency checks.

``load_attached_spectra(path, kind, unit=None, ...) -> list[Spectrum]`` raises on any problem
that could make the spectrum too small or wrongly scaled (truncated file, failed lambda check,
non-finite values, negative weight without an explicit policy). Every returned Spectrum carries
``provenance["assumptions"]`` (possibly empty) and ``provenance["checks"]``.

kinds
-----
``"qe_a2f_dos"``  matdyn a2F.dos<i>. Footer required; lambda recomputed from the Spectrum must
                  match QE's footer (check ``lambda_vs_a2F_footer``). Optionally ``lambda_file``
                  (matdyn's ``lambda``) adds the lambda/omega_log checks for that smearing.
``"epw"``         EPW <prefix>.a2f, one Spectrum per smearing; the 'Integrated el-ph coupling'
                  block is required and must match the Spectrum (``lambda_vs_epw``).
``"two_column"``  generic (omega, a2F) table with an explicit ``unit``. Such a file has no
                  end-of-data marker, so completeness cannot be verified from the file. Pass
                  ``expected_lambda`` and/or ``expected_A_K`` (from the producing code) to verify;
                  otherwise the Spectrum carries the assumption 'A-FILE-COMPLETE(unverified)'.
"""
from __future__ import annotations

import math
from pathlib import Path
from typing import Any

from .checks import ConsistencyCheck, Tolerances, _cmp, check_epw_spectrum, check_qe_spectrum
from .parsers import (epw_a2f_to_spectra, parse_epw_a2f, parse_qe_a2f_dos, parse_qe_lambda_file,
                      qe_a2f_dos_to_spectrum, read_two_column, spectrum_lambda)

KINDS = ("qe_a2f_dos", "epw", "two_column")


class SpectrumCheckError(ValueError):
    """A consistency check between the parsed Spectrum and the producing code failed."""


def _attach(sp: Any, checks: list[ConsistencyCheck], extra_assumptions: list[str]) -> Any:
    sp.provenance["checks"] = [c.to_dict() for c in checks]
    if extra_assumptions:
        sp.provenance["assumptions"] = list(sp.provenance.get("assumptions", [])) + extra_assumptions
    return sp


def _raise_failed(path: Any, checks: list[ConsistencyCheck]) -> None:
    bad = [c for c in checks if not c.ok]
    if bad:
        raise SpectrumCheckError(
            f"{path}: consistency check(s) failed: "
            + "; ".join(f"{c.name}[sigma{c.smearing_index}] ours={c.ours:.8g} "
                        f"ref={c.reference:.8g} rel={c.rel_diff:.3g}" for c in bad))


def load_attached_spectra(path: str | Path, kind: str, unit: str | None = None, *,
                          label: str | None = None, rule: str = "bins",
                          negative_policy: str = "error", lambda_file: str | Path | None = None,
                          expected_lambda: float | None = None,
                          expected_A_K: float | None = None, expected_rel_tol: float = 1e-3,
                          tol: Tolerances | None = None, columns: tuple[int, int] = (0, 1),
                          extra_provenance: dict | None = None,
                          allow_incomplete: bool = False) -> list:
    """Parse ``path`` and verify it; return the alternative spectra it contains.

    ``allow_incomplete=True`` is for diagnostics only: footer-less files are accepted, the
    footer check is skipped, and each Spectrum gets assumption 'A-FILE-INCOMPLETE' and flag
    'INCOMPLETE_INPUT_FILE' in provenance (it must never feed a proof).
    """
    if kind not in KINDS:
        raise ValueError(f"unknown kind {kind!r}; choose from {KINDS}")
    p = Path(path)
    lab = label if label is not None else p.name
    prov = dict(extra_provenance or {})
    prov.setdefault("attached_file", str(p))
    prov.setdefault("format_requested", kind)
    if kind == "qe_a2f_dos":
        if unit not in (None, "Ry", "ry"):
            raise ValueError("QE a2F.dos frequencies are in Rydberg by format; do not pass unit")
        tab = parse_qe_a2f_dos(p, allow_incomplete=allow_incomplete)
        entry = None
        if lambda_file is not None:
            entries = parse_qe_lambda_file(lambda_file)
            if tab.isig is None or tab.isig > len(entries):
                raise ValueError(f"{lambda_file}: no entry for smearing index {tab.isig}")
            entry = entries[tab.isig - 1]
        sp = qe_a2f_dos_to_spectrum(tab, label=lab, extra_provenance=prov,
                                    broadening_Ry=entry.broadening_Ry if entry else None,
                                    negative_policy=negative_policy)
        checks = check_qe_spectrum(sp, tab, entry, tol)
        if tab.complete and not any(c.name == "lambda_vs_a2F_footer" for c in checks):
            raise SpectrumCheckError(f"{p}: footer check could not be run")
        _raise_failed(p, checks)
        return [_attach(sp, checks, [])]
    if kind == "epw":
        if unit not in (None, "meV", "mev"):
            raise ValueError("EPW .a2f frequencies are in meV by format; do not pass unit")
        tab = parse_epw_a2f(p, allow_incomplete=allow_incomplete)
        out = []
        for sp in epw_a2f_to_spectra(tab, label_prefix=lab + ":", extra_provenance=prov,
                                     negative_policy=negative_policy):
            checks = check_epw_spectrum(sp, tol)
            if tab.complete and not checks:
                raise SpectrumCheckError(f"{p}: EPW lambda check could not be run")
            _raise_failed(p, checks)
            out.append(_attach(sp, checks, []))
        return out
    # two_column
    if not unit:
        raise ValueError("two_column files need an explicit unit (meV, THz, cm-1, Ry, K, eV)")
    sp = read_two_column(p, unit, rule=rule, cols=columns, label=lab, extra_provenance=prov,
                         negative_policy=negative_policy)
    checks: list[ConsistencyCheck] = []
    if expected_lambda is not None:
        checks.append(_cmp("lambda_vs_expected", None, spectrum_lambda(sp), float(expected_lambda),
                           expected_rel_tol, 0.0, note="user-supplied lambda of the full spectrum"))
    if expected_A_K is not None:
        checks.append(_cmp("A_vs_expected", None, math.fsum(sp.weight_K), float(expected_A_K),
                           expected_rel_tol, 0.0, note="user-supplied A (kelvin)"))
    _raise_failed(p, checks)
    extra = [] if checks else [
        "A-FILE-COMPLETE(unverified): two-column file has no end-of-data marker; a truncated "
        "file would lower A (anti-conservative for REJECTED). Pass expected_lambda/expected_A_K "
        "to verify"]
    return [_attach(sp, checks, extra)]
