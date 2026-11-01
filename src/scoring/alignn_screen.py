"""Tier-1 ML scoring probe: run ALIGNN 2.0 superconductivity checkpoints on CPU.

Purpose: measure (not trust) what the public ALIGNN Tc / lambda / omega_log models say
about reference structures with known DFPT answers, to see how they behave out of
distribution (high-pressure hydrides with Tc ~ 200 K). Output is a ranking hint only;
it never issues a verdict (see docs/design/11_scoring_tiers.md).

Usage (CPU only, never the GPU):
  CUDA_VISIBLE_DEVICES= OMP_NUM_THREADS=2 nice -n 19 \
    PYTHONPATH=<dir with alignn + jarvis-tools installed via pip --target> \
    ~/.rtt/envs/py/bin/python alignn_screen.py --models-dir <dir> [--poscar f1 f2 ...]

Model checkpoints: ALIGNN 2.0 Figshare registry (alignn/pretrained.py, ALIGNN2_MODELS),
downloaded as zips and extracted into <models-dir>/<short-name>/ (config.json + best_model.pt):
  tc_hyd    = tc_supercon_hydride_radius  (https://ndownloader.figshare.com/files/67163621)
  tc_sc     = tc_supercon_radius          (https://ndownloader.figshare.com/files/67163615)
  alex_tc   = alex_supercon_Tc            (https://ndownloader.figshare.com/files/67163582)
  alex_la   = alex_supercon_la            (https://ndownloader.figshare.com/files/67163591)
  alex_wlog = alex_supercon_wlog          (https://ndownloader.figshare.com/files/67163597)
"""
from __future__ import annotations

import argparse
import json
import os
import time

os.environ.setdefault("CUDA_VISIBLE_DEVICES", "")

import torch  # noqa: E402

torch.set_num_threads(2)

from pymatgen.core import Lattice, Structure  # noqa: E402
from pymatgen.io.jarvis import JarvisAtomsAdaptor  # noqa: E402


def reference_structures() -> dict[str, tuple[Structure, str]]:
    """Reference set: structure + literature/bench DFT Tc note (for comparison only).

    Lattice parameters are approximate textbook/DFT values (UNVERIFIED to better than ~1-2 %).
    """
    s = {}
    s["Al_fcc_0GPa"] = (Structure.from_spacegroup("Fm-3m", Lattice.cubic(4.05), ["Al"], [[0, 0, 0]]),
                        "exp Tc 1.2 K")
    s["Pb_fcc_0GPa"] = (Structure.from_spacegroup("Fm-3m", Lattice.cubic(4.95), ["Pb"], [[0, 0, 0]]),
                        "exp Tc 7.2 K")
    s["Nb_bcc_0GPa"] = (Structure.from_spacegroup("Im-3m", Lattice.cubic(3.30), ["Nb"], [[0, 0, 0]]),
                        "exp Tc 9.2 K")
    s["MgB2_0GPa"] = (Structure.from_spacegroup("P6/mmm", Lattice.hexagonal(3.086, 3.524), ["Mg", "B"],
                                                [[0, 0, 0], [1 / 3, 2 / 3, 0.5]]),
                      "exp Tc 39 K (isotropic ME ~20-25 K)")
    # H3S Im-3m: cell taken from an earlier benchmark relax at 200 GPa (primitive bcc, a_conv = 2*1.494124 A)
    s["H3S_Im-3m_200GPa"] = (Structure.from_spacegroup("Im-3m", Lattice.cubic(2 * 1.494124044), ["S", "H"],
                                                       [[0, 0, 0], [0, 0.5, 0.5]]),
                             "DFT Tc ~190-200 K; bench DFPT lambda 2.40, wlog 1244 K (sigma 0.01 Ry)")
    s["CaH6_Im-3m_150GPa"] = (Structure.from_spacegroup("Im-3m", Lattice.cubic(3.50), ["Ca", "H"],
                                                        [[0, 0, 0], [0.25, 0, 0.5]]),
                              "DFT Tc ~215-235 K (Wang 2012 PNAS); a UNVERIFIED")
    s["LaH10_Fm-3m_200GPa"] = (Structure.from_spacegroup("Fm-3m", Lattice.cubic(5.10), ["La", "H", "H"],
                                                         [[0, 0, 0], [0.25, 0.25, 0.25],
                                                          [0.12, 0.12, 0.12]]),
                               "DFT Tc ~250-280 K; a and x(H 32f) UNVERIFIED")
    return s


def load_models(models_dir: str) -> dict:
    from alignn.deprecated.pretrained import load_pure_torch_model  # type: ignore
    out = {}
    for name in sorted(os.listdir(models_dir)):
        p = os.path.join(models_dir, name)
        if os.path.isdir(p) and os.path.exists(os.path.join(p, "best_model.pt")):
            model, config = load_pure_torch_model(p, device="cpu")
            model.eval()
            out[name] = (model, config)
    return out


def predict(model, config, atoms) -> float:
    from alignn.deprecated.pretrained import _atom_features_from_config  # type: ignore
    from alignn.torch_graph_builder import build_pure_torch_graph  # type: ignore
    g, lg = build_pure_torch_graph(atoms=atoms, two_body_cutoff=float(config.get("cutoff", 5.0)),
                                   max_neighbors=int(config.get("max_neighbors", 12)),
                                   atom_features=_atom_features_from_config(config),
                                   compute_line_graph=True, device="cpu")
    lat = torch.tensor(atoms.lattice_mat).type(torch.get_default_dtype())
    with torch.no_grad():
        out = model([g, lg, lat])
    out = out["out"] if isinstance(out, dict) else out
    return float(out.detach().cpu().numpy().flatten()[0])


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--models-dir", required=True)
    ap.add_argument("--poscar", nargs="*", default=[])
    ap.add_argument("--out", default="alignn_predictions.json")
    a = ap.parse_args()
    structs = reference_structures()
    for f in a.poscar:
        structs[os.path.basename(f)] = (Structure.from_file(f), "user")
    models = load_models(a.models_dir)
    rows = []
    for sname, (st, note) in structs.items():
        atoms = JarvisAtomsAdaptor.get_atoms(st)
        row = {"structure": sname, "natoms": len(st), "note": note}
        t0 = time.time()
        for mname, (m, c) in models.items():
            row[mname] = round(predict(m, c, atoms), 3)
        row["seconds_all_models"] = round(time.time() - t0, 3)
        rows.append(row)
        print(json.dumps(row))
    with open(a.out, "w") as fh:
        json.dump(rows, fh, indent=1)


if __name__ == "__main__":
    main()
