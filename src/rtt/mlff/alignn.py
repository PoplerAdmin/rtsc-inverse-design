"""S3(c): public ALIGNN superconductivity numbers -- **recorded only** (doc 11 §2.3; contract R2).

Ported from ``src/scoring/alignn_screen.py`` (CPU only, 2 threads). ALIGNN is not installed in
the main venv; it is used only when ``RTT_ALIGNN_MODELS`` points at a directory of extracted
ALIGNN 2.0 checkpoints (``<name>/config.json`` + ``best_model.pt``) and ``alignn`` /
``jarvis-tools`` are importable. Otherwise the result is ``{"available": False, ...}``.

The weight is 0 for any structure containing H or another light element (H, Li, Be, B, C, N):
doc 11 measured that these models fail on hydrides. Even at weight > 0 they never remove a
candidate.
"""
from __future__ import annotations

import os
from typing import Any

LIGHT = ("H", "Li", "Be", "B", "C", "N")
_MODELS: dict | None = None


def weight_for(species: list[str], light_weight: float = 0.0, other_weight: float = 0.0) -> float:
    return light_weight if any(s in LIGHT for s in species) else other_weight


def _load(models_dir: str) -> dict:
    global _MODELS
    if _MODELS is None:
        from alignn.deprecated.pretrained import load_pure_torch_model  # type: ignore
        out = {}
        for name in sorted(os.listdir(models_dir)):
            p = os.path.join(models_dir, name)
            if os.path.isdir(p) and os.path.exists(os.path.join(p, "best_model.pt")):
                m, c = load_pure_torch_model(p, device="cpu")
                m.eval()
                out[name] = (m, c)
        _MODELS = out
    return _MODELS


def alignn_numbers(struct: Any, models_dir: str | None = None) -> dict:
    from ..structure.canonical import to_pymatgen
    d = struct.to_json() if hasattr(struct, "to_json") else struct
    st = to_pymatgen(d)
    species = [s.specie.symbol for s in st]
    res: dict[str, Any] = {"weight": weight_for(species), "source_id": "doc11#2.3"}
    models_dir = models_dir or os.environ.get("RTT_ALIGNN_MODELS")
    if not models_dir or not os.path.isdir(models_dir):
        res.update(available=False, reason="RTT_ALIGNN_MODELS not set")
        return res
    try:
        import torch
        torch.set_num_threads(2)
        from alignn.deprecated.pretrained import _atom_features_from_config  # type: ignore
        from alignn.torch_graph_builder import build_pure_torch_graph  # type: ignore
        from pymatgen.io.jarvis import JarvisAtomsAdaptor
        atoms = JarvisAtomsAdaptor.get_atoms(st)
        for name, (m, c) in _load(models_dir).items():
            g, lg = build_pure_torch_graph(atoms=atoms, two_body_cutoff=float(c.get("cutoff", 5.0)),
                                           max_neighbors=int(c.get("max_neighbors", 12)),
                                           atom_features=_atom_features_from_config(c),
                                           compute_line_graph=True, device="cpu")
            lat = torch.tensor(atoms.lattice_mat).type(torch.get_default_dtype())
            with torch.no_grad():
                o = m([g, lg, lat])
            o = o["out"] if isinstance(o, dict) else o
            res[name] = round(float(o.detach().cpu().numpy().flatten()[0]), 3)
        res["available"] = True
    except Exception as ex:  # noqa: BLE001
        res.update(available=False, reason=f"{type(ex).__name__}: {ex}")
    return res


__all__ = ["alignn_numbers", "weight_for"]
