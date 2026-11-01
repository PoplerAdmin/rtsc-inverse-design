"""Pseudopotential table: one directory + per-element file mapping (SSSP-style).

Config forms accepted::

    [dfpt.pseudo]
    dir = "/data/SSSP_1.3_PBE_efficiency"
    sssp_json = "/data/SSSP_1.3_PBE_efficiency.json"   # optional: filenames + cutoffs
    [dfpt.pseudo.files]                                # optional: overrides / no-JSON setups
    Al = "Al.pbe-n-kjpaw_psl.1.0.0.UPF"

SSSP JSON layout (as distributed by Materials Cloud): ``{"Al": {"filename": ..., "md5": ...,
"pseudopotential": ..., "cutoff_wfc": <Ry>, "cutoff_rho": <Ry>}, ...}``.
UPF v2 headers carry ``element``, ``functional``, ``wfc_cutoff``/``rho_cutoff`` (Ry) attributes;
we read them for validation and as cutoff fallbacks.
"""
from __future__ import annotations

import hashlib
import json
import os
import re
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any


@dataclass
class UPFInfo:
    element: str | None
    functional: str | None
    wfc_cutoff_Ry: float | None
    rho_cutoff_Ry: float | None
    upf_version: str | None


def read_upf_header(path: str | os.PathLike) -> UPFInfo:
    with open(path, "r", errors="replace") as fh:
        head = fh.read(20000)
    ver = None
    mv = re.search(r'<UPF\s+version\s*=\s*"([^"]+)"', head)
    if mv:
        ver = mv.group(1)
    mh = re.search(r"<PP_HEADER(.*?)(/>|>)", head, re.S)

    def attr(name: str) -> str | None:
        if not mh:
            return None
        m = re.search(name + r'\s*=\s*"([^"]*)"', mh.group(1))
        return m.group(1).strip() if m else None

    el = attr("element")
    fn = attr("functional")
    wc, rc = attr("wfc_cutoff"), attr("rho_cutoff")

    def num(x: str | None) -> float | None:
        try:
            v = float(x.replace("D", "E").replace("d", "e")) if x else None
        except ValueError:
            return None
        return v if v and v > 0 else None

    if el is None:  # UPF v1: "  Al   Element" style line inside <PP_HEADER>
        m1 = re.search(r"<PP_HEADER>(.*?)</PP_HEADER>", head, re.S)
        if m1:
            for line in m1.group(1).splitlines():
                if "element" in line.lower():
                    el = line.split()[0]
                    break
                    # (v1 headers do not reliably carry suggested cutoffs)
    return UPFInfo(el, fn, num(wc), num(rc), ver)


def _md5(path: Path) -> str:
    h = hashlib.md5()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


@dataclass
class PseudoTable:
    directory: str
    files: dict[str, str] = field(default_factory=dict)            # element -> filename
    cutoffs_Ry: dict[str, tuple[float, float]] = field(default_factory=dict)  # el -> (wfc, rho)

    @classmethod
    def from_sssp_json(cls, directory: str, json_path: str,
                       overrides: dict[str, str] | None = None) -> "PseudoTable":
        with open(json_path) as fh:
            data = json.load(fh)
        files, cut = {}, {}
        for el, rec in data.items():
            if not isinstance(rec, dict) or "filename" not in rec:
                continue
            files[el] = rec["filename"]
            if "cutoff_wfc" in rec and "cutoff_rho" in rec:
                cut[el] = (float(rec["cutoff_wfc"]), float(rec["cutoff_rho"]))
        files.update(overrides or {})
        return cls(directory, files, cut)

    @classmethod
    def from_config(cls, cfg: dict[str, Any]) -> "PseudoTable":
        d = cfg["dir"]
        if cfg.get("sssp_json"):
            t = cls.from_sssp_json(d, cfg["sssp_json"], cfg.get("files"))
        else:
            t = cls(d, dict(cfg.get("files", {})))
        for el, c in (cfg.get("cutoffs_Ry") or {}).items():
            t.cutoffs_Ry[el] = (float(c[0]), float(c[1]))
        return t

    def path(self, element: str) -> Path:
        if element not in self.files:
            raise KeyError(f"no pseudopotential mapped for element {element!r} "
                           f"(table has {sorted(self.files)})")
        return Path(self.directory) / self.files[element]

    def validate(self, elements: list[str]) -> dict[str, Any]:
        """Check files exist and headers match; return provenance (file, md5, functional)."""
        prov: dict[str, Any] = {}
        functionals = set()
        problems = []
        for el in elements:
            p = self.path(el)
            if not p.is_file():
                problems.append(f"{el}: missing pseudopotential file {p}")
                continue
            info = read_upf_header(p)
            if info.element and info.element.strip().capitalize() != el:
                problems.append(f"{el}: UPF header says element {info.element!r} ({p.name})")
            if info.functional:
                functionals.add(" ".join(info.functional.upper().split()))
            prov[el] = {"file": str(p), "md5": _md5(p), "functional": info.functional,
                        "upf_version": info.upf_version}
        if problems:
            raise ValueError("pseudopotential problems: " + "; ".join(problems))
        prov["_warnings"] = ([f"mixed exchange-correlation functionals in pseudopotentials: "
                              f"{sorted(functionals)}"] if len(functionals) > 1 else [])
        return prov

    def cutoffs(self, elements: list[str]) -> tuple[float | None, float | None]:
        """max over elements of (ecutwfc, ecutrho) in Ry from SSSP json, else UPF header."""
        wfc, rho = [], []
        for el in elements:
            if el in self.cutoffs_Ry:
                w, r = self.cutoffs_Ry[el]
            else:
                p = self.path(el)
                info = read_upf_header(p) if p.is_file() else UPFInfo(None, None, None, None, None)
                w, r = info.wfc_cutoff_Ry, info.rho_cutoff_Ry
            if w:
                wfc.append(w)
            if r:
                rho.append(r)
        return (max(wfc) if len(wfc) == len(elements) and wfc else None,
                max(rho) if len(rho) == len(elements) and rho else None)
