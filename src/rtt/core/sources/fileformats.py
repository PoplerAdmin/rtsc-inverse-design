"""Dependency-free readers for the structure formats the sources need.

These are fallbacks: when pymatgen or ASE is installed, ``sources.local`` prefers them.  Each
reader returns a list of ``RawStructure`` (lattice rows in angstrom, species, fractional coords,
plus any scalar metadata found in the file).  They are deliberately strict: anything they cannot
interpret raises ``ValueError`` (the pipeline turns that into an ERROR verdict / skipped file
with a logged reason; nothing is silently invented).
"""
from __future__ import annotations

import math
import re
import shlex
from dataclasses import dataclass, field
from typing import Any


@dataclass
class RawStructure:
    lattice: list[list[float]]
    species: list[str]
    frac_coords: list[list[float]]
    info: dict[str, Any] = field(default_factory=dict)


def cellpar_to_matrix(a, b, c, alpha, beta, gamma) -> list[list[float]]:
    """Lattice rows from cell parameters (a along x, b in the xy plane)."""
    al, be, ga = (math.radians(x) for x in (alpha, beta, gamma))
    ca, cb, cg, sg = math.cos(al), math.cos(be), math.cos(ga), math.sin(ga)
    cx = c * cb
    cy = c * (ca - cb * cg) / sg
    cz2 = c * c - cx * cx - cy * cy
    if cz2 <= 0:
        raise ValueError("invalid cell parameters")
    return [[a, 0.0, 0.0], [b * cg, b * sg, 0.0], [cx, cy, math.sqrt(cz2)]]


def _inv3(m):
    a, b, c = m
    det = (a[0] * (b[1] * c[2] - b[2] * c[1]) - a[1] * (b[0] * c[2] - b[2] * c[0])
           + a[2] * (b[0] * c[1] - b[1] * c[0]))
    if abs(det) < 1e-12:
        raise ValueError("singular lattice")
    inv = [[(b[1] * c[2] - b[2] * c[1]) / det, (a[2] * c[1] - a[1] * c[2]) / det,
            (a[1] * b[2] - a[2] * b[1]) / det],
           [(b[2] * c[0] - b[0] * c[2]) / det, (a[0] * c[2] - a[2] * c[0]) / det,
            (a[2] * b[0] - a[0] * b[2]) / det],
           [(b[0] * c[1] - b[1] * c[0]) / det, (a[1] * c[0] - a[0] * c[1]) / det,
            (a[0] * b[1] - a[1] * b[0]) / det]]
    return inv


def cart_to_frac(lattice, cart):
    inv = _inv3(lattice)   # rows: x_frac = cart @ inv
    return [[sum(p[k] * inv[k][j] for k in range(3)) for j in range(3)] for p in cart]


# --------------------------------------------------------------------------- POSCAR / CONTCAR
def read_poscar(text: str, species_hint: list[str] | None = None) -> RawStructure:
    lines = [ln.rstrip() for ln in text.splitlines()]
    if len(lines) < 8:
        raise ValueError("POSCAR too short")
    title = lines[0].strip()
    scale = [float(x) for x in lines[1].split()[:3]]
    lat = [[float(x) for x in lines[i].split()[:3]] for i in (2, 3, 4)]
    idx = 5
    tokens = lines[idx].split()
    if tokens and not _is_number(tokens[0]):
        symbols = [t.split("/")[0].split("_")[0] for t in tokens]
        idx += 1
        counts = [int(x) for x in lines[idx].split()]
    else:
        counts = [int(x) for x in tokens]
        symbols = species_hint or [t for t in title.split() if re.fullmatch(r"[A-Z][a-z]?", t)]
        if len(symbols) != len(counts):
            raise ValueError("VASP4 POSCAR without species line and no usable title / hint")
    idx += 1
    mode = lines[idx].strip()
    if mode[:1] in "sS":
        idx += 1
        mode = lines[idx].strip()
    if len(scale) == 1:
        s = scale[0]
        if s < 0:  # negative scale = target volume
            vol = abs(_det(lat))
            s = (abs(s) / vol) ** (1.0 / 3.0)
        svec = [s, s, s]
    else:
        svec = scale
    lat = [[x * svec[j] for j, x in enumerate(r)] for r in lat]
    n = sum(counts)
    coords = [[float(x) for x in lines[idx + 1 + i].split()[:3]] for i in range(n)]
    if mode[:1] in "cCkK":   # cartesian coordinates are scaled like the lattice
        coords = cart_to_frac(lat, [[x * svec[j] for j, x in enumerate(p)] for p in coords])
    species = [s for s, c in zip(symbols, counts) for _ in range(c)]
    return RawStructure(lat, species, coords, {"title": title})


def _det(m):
    a, b, c = m
    return (a[0] * (b[1] * c[2] - b[2] * c[1]) - a[1] * (b[0] * c[2] - b[2] * c[0])
            + a[2] * (b[0] * c[1] - b[1] * c[0]))


def _is_number(s: str) -> bool:
    try:
        float(s)
        return True
    except ValueError:
        return False


def split_concatenated_poscars(text: str) -> list[str]:
    """Split a file of concatenated POSCARs (USPEX gatheredPOSCARS, CALYPSO) into blocks."""
    lines = text.splitlines()
    blocks, i = [], 0
    while i < len(lines):
        if not lines[i].strip():
            i += 1
            continue
        start = i
        try:
            j = i + 5
            toks = lines[j].split()
            if toks and not _is_number(toks[0]):
                j += 1
            n = sum(int(x) for x in lines[j].split())
            j += 1
            if lines[j].strip()[:1] in "sS":
                j += 1
            end = j + 1 + n
        except (IndexError, ValueError) as exc:
            raise ValueError(f"cannot split POSCAR block at line {i + 1}: {exc}") from None
        blocks.append("\n".join(lines[start:end]))
        i = end
    return blocks


# --------------------------------------------------------------------------- extended XYZ
_KV = re.compile(r'(\w+)=("[^"]*"|\S+)')


def read_extxyz(text: str) -> list[RawStructure]:
    lines = text.splitlines()
    out, i = [], 0
    while i < len(lines):
        if not lines[i].strip():
            i += 1
            continue
        n = int(lines[i].split()[0])
        comment = lines[i + 1]
        info: dict[str, Any] = {}
        for k, v in _KV.findall(comment):
            v = v.strip('"')
            info[k] = v
        if "Lattice" not in info:
            raise ValueError("extxyz frame without Lattice= (non-periodic structures unsupported)")
        lv = [float(x) for x in info.pop("Lattice").split()]
        lat = [lv[0:3], lv[3:6], lv[6:9]]
        props = info.pop("Properties", "species:S:1:pos:R:3")
        cols = props.split(":")
        layout, off = {}, 0
        for k in range(0, len(cols), 3):
            name, cnt = cols[k], int(cols[k + 2])
            layout[name] = (off, cnt)
            off += cnt
        so, _ = layout.get("species", (0, 1))
        po, _ = layout.get("pos", (1, 3))
        species, cart = [], []
        for ln in lines[i + 2:i + 2 + n]:
            t = ln.split()
            species.append(t[so])
            cart.append([float(x) for x in t[po:po + 3]])
        meta: dict[str, Any] = {}
        for k, v in info.items():
            meta[k] = float(v) if _is_number(v) else v
        out.append(RawStructure(lat, species, cart_to_frac(lat, cart), meta))
        i += 2 + n
    if not out:
        raise ValueError("no frames in extxyz")
    return out


# --------------------------------------------------------------------------- CIF (minimal)
def _cif_num(s: str) -> float:
    return float(re.sub(r"\(.*?\)", "", s))


def _cif_tokens(text: str) -> list[str]:
    toks: list[str] = []
    in_semi = False
    buf: list[str] = []
    for line in text.splitlines():
        if in_semi:
            if line.startswith(";"):
                toks.append("\n".join(buf))
                in_semi, buf = False, []
            else:
                buf.append(line)
            continue
        if line.startswith(";"):
            in_semi, buf = True, [line[1:]]
            continue
        s = line.split("#", 1)[0] if not line.lstrip().startswith("'") else line
        try:
            toks.extend(shlex.split(s, posix=True))
        except ValueError:
            toks.extend(s.split())
    return toks


def _parse_symop(op: str):
    terms = op.lower().replace(" ", "").split(",")
    if len(terms) != 3:
        raise ValueError(f"bad symmetry operation {op!r}")
    rows = []
    for t in terms:
        rot = [0.0, 0.0, 0.0]
        trans = 0.0
        for m in re.finditer(r"([+-]?)([^+-]+)", t):
            sign = -1.0 if m.group(1) == "-" else 1.0
            body = m.group(2)
            if body in "xyz":
                rot["xyz".index(body)] += sign
            elif re.fullmatch(r"\d*\.?\d+(/\d*\.?\d+)?", body):
                if "/" in body:
                    p, q = body.split("/")
                    trans += sign * float(p) / float(q)
                else:
                    trans += sign * float(body)
            elif "*" in body:
                f, v = body.split("*")
                rot["xyz".index(v)] += sign * float(f)
            else:
                raise ValueError(f"bad symmetry term {body!r} in {op!r}")
        rows.append((rot, trans))
    return rows


def read_cif(text: str, symprec: float = 1e-3) -> list[RawStructure]:
    """Minimal CIF reader: cell, symmetry operations, atom sites (with occupancy).

    Applies the listed symmetry operations to the asymmetric unit and removes duplicates.
    Partial occupancies are kept as fractional composition in ``info['partial_occupancy']``.
    """
    toks = _cif_tokens(text)
    blocks: list[dict[str, Any]] = []
    cur: dict[str, Any] | None = None
    i = 0
    while i < len(toks):
        t = toks[i]
        tl = t.lower()
        if tl.startswith("data_"):
            cur = {"_name": t[5:], "_loops": []}
            blocks.append(cur)
            i += 1
        elif tl == "loop_":
            i += 1
            names = []
            while i < len(toks) and toks[i].startswith("_"):
                names.append(toks[i].lower())
                i += 1
            vals = []
            while i < len(toks) and not toks[i].startswith("_") and toks[i].lower() != "loop_" \
                    and not toks[i].lower().startswith("data_"):
                vals.append(toks[i])
                i += 1
            if cur is not None and names:
                rows = [vals[k:k + len(names)] for k in range(0, len(vals) - len(names) + 1, len(names))]
                cur["_loops"].append((names, rows))
        elif t.startswith("_"):
            if cur is not None and i + 1 < len(toks):
                cur[tl] = toks[i + 1]
            i += 2
        else:
            i += 1
    out = []
    for b in blocks:
        if "_cell_length_a" not in b:
            continue
        lat = cellpar_to_matrix(*(_cif_num(b[k]) for k in (
            "_cell_length_a", "_cell_length_b", "_cell_length_c",
            "_cell_angle_alpha", "_cell_angle_beta", "_cell_angle_gamma")))
        ops = None
        sites = None
        for names, rows in b["_loops"]:
            for key in ("_symmetry_equiv_pos_as_xyz", "_space_group_symop_operation_xyz"):
                if key in names:
                    k = names.index(key)
                    ops = [_parse_symop(r[k]) for r in rows]
            if "_atom_site_fract_x" in names:
                sites = (names, rows)
        if sites is None:
            raise ValueError(f"CIF block {b['_name']}: no _atom_site_fract_x loop")
        if ops is None:
            ops = [_parse_symop("x,y,z")]
        names, rows = sites
        ix, iy, iz = (names.index(f"_atom_site_fract_{c}") for c in "xyz")
        isym = names.index("_atom_site_type_symbol") if "_atom_site_type_symbol" in names else None
        ilab = names.index("_atom_site_label") if "_atom_site_label" in names else None
        iocc = names.index("_atom_site_occupancy") if "_atom_site_occupancy" in names else None
        species, coords, partial = [], [], False
        for r in rows:
            raw = r[isym] if isym is not None else r[ilab]
            m = re.match(r"([A-Z][a-z]?)", raw)
            if not m:
                raise ValueError(f"cannot read element from {raw!r}")
            el = m.group(1)
            occ = _cif_num(r[iocc]) if iocc is not None and r[iocc] not in ("?", ".") else 1.0
            if abs(occ - 1.0) > 1e-6:
                partial = True
            p = [_cif_num(r[ix]), _cif_num(r[iy]), _cif_num(r[iz])]
            for op in ops:
                q = [(sum(rot[k] * p[k] for k in range(3)) + tr) % 1.0 for rot, tr in op]
                dup = False
                for s2, c2 in zip(species, coords):
                    if s2 == el:
                        d = [((q[k] - c2[k] + 0.5) % 1.0) - 0.5 for k in range(3)]
                        if max(abs(x) for x in d) < symprec:
                            dup = True
                            break
                if not dup:
                    species.append(el)
                    coords.append(q)
        info: dict[str, Any] = {"name": b["_name"]}
        for k in ("_symmetry_space_group_name_h-m", "_space_group_name_h-m_alt"):
            if k in b:
                info["spacegroup"] = b[k]
        for k in ("_symmetry_int_tables_number", "_space_group_it_number"):
            if k in b:
                info["spacegroup_number"] = int(_cif_num(b[k]))
        if "_chemical_formula_sum" in b:
            info["formula_sum"] = b["_chemical_formula_sum"]
        if partial:
            info["partial_occupancy"] = True
        out.append(RawStructure(lat, species, coords, info))
    if not out:
        raise ValueError("no crystal data block in CIF")
    return out


# --------------------------------------------------------------------------- AIRSS .res
def read_res(text: str) -> RawStructure:
    """AIRSS/SHELX .res: ``TITL seed pressure volume enthalpy spin |spin| nions (sg) n - times``,
    ``CELL lambda a b c alpha beta gamma``, ``SFAC``, atom lines ``El idx x y z occ``, ``END``.
    Pressure in GPa, volume in A^3, enthalpy in eV for the whole cell (AIRSS convention)."""
    info: dict[str, Any] = {}
    lat = None
    species, coords = [], []
    in_atoms = False
    for line in text.splitlines():
        t = line.split()
        if not t:
            continue
        key = t[0].upper()
        if key == "TITL":
            info["seed"] = t[1] if len(t) > 1 else ""
            nums = []
            for x in t[2:]:
                if _is_number(x):
                    nums.append(float(x))
                else:
                    break
            for name, v in zip(("pressure_GPa", "volume_A3", "enthalpy_eV", "spin", "abs_spin",
                                "n_ions"), nums):
                info[name] = v
            m = re.search(r"\(([^)]*)\)", line)
            if m:
                info["spacegroup"] = m.group(1)
            m = re.search(r"n\s*-\s*(\d+)", line)
            if m:
                info["times_found"] = int(m.group(1))
        elif key == "CELL":
            lat = cellpar_to_matrix(*(float(x) for x in t[2:8]))
        elif key == "SFAC":
            in_atoms = True
        elif key == "END":
            break
        elif in_atoms and len(t) >= 5 and re.fullmatch(r"[A-Z][a-z]?", t[0]):
            species.append(t[0])
            coords.append([float(x) for x in t[2:5]])
    if lat is None or not species:
        raise ValueError("incomplete .res file (no CELL or no atoms)")
    return RawStructure(lat, species, coords, info)
