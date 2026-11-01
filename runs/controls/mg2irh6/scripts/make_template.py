"""Write an ibrav=0 pw.x scf template (CELL_PARAMETERS + positions in angstrom) for Fm-3m Mg2IrH6,
for src/scoring/coarse_dfpt.py --pw-template.  Usage: make_template.py a_A x_H ecut out"""
import sys
a, x, ecut, out = float(sys.argv[1]), float(sys.argv[2]), float(sys.argv[3]), sys.argv[4]
h = a / 2
cell = [(-h, 0, h), (0, h, h), (-h, h, 0)]
pos = [("Ir", 0, 0, 0), ("Mg", .25, .25, .25), ("Mg", -.25, -.25, -.25),
       ("H", x, 0, 0), ("H", -x, 0, 0), ("H", 0, x, 0), ("H", 0, -x, 0), ("H", 0, 0, x), ("H", 0, 0, -x)]
txt = f""" &control
    calculation = 'scf'
 /
 &system
    ibrav = 0
    nat = 9
    ntyp = 3
    ecutwfc = {ecut}
    ecutrho = {4 * ecut}
 /
 &electrons
 /
ATOMIC_SPECIES
 Mg  24.30500 Mg_ONCV_PBE-1.2.upf
 Ir  192.21700 Ir_ONCV_PBE-1.2.upf
 H   1.00794 H_ONCV_PBE-1.2.upf
CELL_PARAMETERS angstrom
""" + "\n".join("  " + " ".join(f"{v:18.12f}" for v in r) for r in cell) + "\nATOMIC_POSITIONS angstrom\n" + \
    "\n".join(f" {s:<3s} " + " ".join(f"{v * a:18.12f}" for v in p) for s, *p in pos) + "\nK_POINTS automatic\n 8 8 8 0 0 0\n"
open(out, "w").write(txt)
