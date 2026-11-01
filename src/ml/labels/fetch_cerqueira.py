"""Stream one batch of the Cerqueira-Sanna-Marques EPC archive and keep only label-relevant files.

Source: Materials Cloud Archive record 3kbt5-r3n56, DOI 10.24435/materialscloud:qv-bq (CC BY 4.0).
Each batch-X.tar.bz2 (~3 GB) holds ~500 QE electron-phonon calculations.  The full archive is ~50 GB,
so we stream it (curl | this script) and keep only small files per material:
  Tc.dat, lambda, McMillan.dat, Eliashberg.dat, geo_opt.cif, start.cif, scf_fine.in, geoopt_0.in,
  qe.freq.gp (phonon dispersion, for imaginary-mode checks), matdyn.dos, a2F.dos1..10 (gzipped).
The md5 (to compare with the archive's published checksum) and sha256 of the whole stream are
written to <out>/<batch>.stream_hash.json.

Usage:  curl -sL URL | python -I fetch_cerqueira.py OUTDIR BATCHNAME
"""
import gzip
import hashlib
import json
import os
import sys
import tarfile

KEEP = {"Tc.dat", "lambda", "McMillan.dat", "Eliashberg.dat", "geo_opt.cif", "start.cif",
        "scf_fine.in", "geoopt_0.in", "qe.freq.gp", "matdyn.dos"}


class HashReader:
    def __init__(self, f):
        self.f, self.md5, self.sha, self.n = f, hashlib.md5(), hashlib.sha256(), 0

    def read(self, n=-1):
        b = self.f.read(n)
        self.md5.update(b)
        self.sha.update(b)
        self.n += len(b)
        return b


def main():
    out, batch = sys.argv[1], sys.argv[2]
    os.makedirs(out, exist_ok=True)
    hr = HashReader(sys.stdin.buffer)
    nmat, nfile = set(), 0
    with tarfile.open(fileobj=hr, mode="r|bz2") as t:
        for m in t:
            if not m.isfile():
                continue
            parts = m.name.split("/")
            if len(parts) != 3 or ".." in parts:
                continue
            _, mat, fn = parts
            keep_a2f = fn.startswith("a2F.dos")
            if fn not in KEEP and not keep_a2f:
                continue
            d = os.path.join(out, batch, mat)
            os.makedirs(d, exist_ok=True)
            data = t.extractfile(m).read()
            if keep_a2f or fn in ("qe.freq.gp", "matdyn.dos"):
                with gzip.open(os.path.join(d, fn + ".gz"), "wb") as g:
                    g.write(data)
            else:
                with open(os.path.join(d, fn), "wb") as g:
                    g.write(data)
            nmat.add(mat)
            nfile += 1
    # drain any trailing bytes so the hash covers the whole file
    while hr.read(1 << 20):
        pass
    rec = {"batch": batch, "bytes": hr.n, "md5": hr.md5.hexdigest(), "sha256": hr.sha.hexdigest(),
           "materials": len(nmat), "files_kept": nfile}
    json.dump(rec, open(os.path.join(out, f"{batch}.stream_hash.json"), "w"), indent=1)
    print(json.dumps(rec), flush=True)


if __name__ == "__main__":
    main()
