"""Human-readable S4 rows (formula, sg, n, outcome, gap, N(E_F), sheets, H_DOS, E_form) for doc 30."""
import collections, json, sqlite3, sys
c = sqlite3.connect(sys.argv[1])
for line in open(sys.argv[2]):
    r = json.loads(line)
    sha = c.execute("select structure_sha256 from candidates where id=?", (r["cid"],)).fetchone()[0]
    s = json.loads(c.execute("select json from structures where sha256=?", (sha,)).fetchone()[0])
    f = "".join(f"{e}{n if n > 1 else ''}" for e, n in sorted(collections.Counter(s["species"]).items()))
    rd = lambda v, k=3: None if v is None else round(v, k)  # noqa: E731
    print("|", r["cid"], "|", r["lane"], "|", f, "|", s.get("spacegroup"), "|", len(s["species"]), "|",
          r["S4"] + (f"({r['crit']})" if r["crit"] else ""), "|", rd(r["gap_eV"], 2), "|", rd(r["N_EF_per_eV_atom"]),
          "|", r["n_sheets"], "|", rd(r["H_DOS"], 2), "|", rd(r["e_form_eV"]), "|")
