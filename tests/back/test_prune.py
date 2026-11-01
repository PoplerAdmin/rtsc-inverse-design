"""rtt.dft.prune: hash-then-delete, per-q selection, manifests, verification, disk guard."""
from __future__ import annotations

import os
import shutil
from pathlib import Path

import pytest

from rtt.dft import prune as P


def _ph_tree(w: Path, nq: int = 4, done_through: int = 3) -> None:
    """Fake ldisp ph.x scratch: dyn1..dyn<done_through> complete, q dirs 2..nq with wfc+dvscf."""
    (w / "tmp" / "pwscf.save").mkdir(parents=True)
    (w / "tmp" / "pwscf.save" / "wfc1.dat").write_bytes(b"keep-save")
    (w / "tmp" / "pwscf.wfc1").write_bytes(b"scf-wfc")
    for n in range(2, nq + 1):
        d = w / "tmp" / "_ph0" / f"pwscf.q_{n}"
        (d / "_ph0").mkdir(parents=True)
        (d / f"pwscf.wfc1").write_bytes(b"w" * (100 * n))
        (d / "_ph0" / "pwscf.dvscf1").write_bytes(b"d" * 50)
        (d / "data-file-schema.xml").write_text("<xml/>")
    for n in range(1, done_through + 1):
        (w / f"dyn{n}").write_text(f"dyn {n}\n")
    (w / "dyn0").write_text("4 4 4\n4\n")


def test_finished_q_numbers(tmp_path):
    _ph_tree(tmp_path, nq=4, done_through=3)
    # dyn3 exists but dyn4 not yet -> q3 still running; q1, q2 finished
    assert P.finished_q_numbers(tmp_path) == {1, 2}


def test_prune_q_scratch_hash_then_delete(tmp_path, monkeypatch):
    monkeypatch.setenv("RTT_MACHINE_LABEL", "test-machine")
    _ph_tree(tmp_path, nq=4, done_through=4)          # dyn4 exists -> q2, q3 finished
    central = tmp_path / "central"
    keep = tmp_path / "copy"
    shutil.copytree(tmp_path / "tmp", keep)
    freed = P.prune_q_scratch(tmp_path, central_dir=central)
    assert freed == 200 + 300
    assert not (tmp_path / "tmp/_ph0/pwscf.q_2/pwscf.wfc1").exists()
    assert not (tmp_path / "tmp/_ph0/pwscf.q_3/pwscf.wfc1").exists()
    assert (tmp_path / "tmp/_ph0/pwscf.q_4/pwscf.wfc1").exists()          # running q untouched
    assert (tmp_path / "tmp/_ph0/pwscf.q_2/_ph0/pwscf.dvscf1").exists()   # dvscf kept by default
    assert (tmp_path / "tmp/_ph0/pwscf.q_2/data-file-schema.xml").exists()
    man = tmp_path / "PRUNED_MANIFEST.tsv"
    rows = man.read_text().splitlines()
    assert rows[0] == P.HEADER.strip()
    assert len(rows) == 3 and all("test-machine" in r for r in rows[1:])
    cent = list(central.glob("*.tsv"))
    assert len(cent) == 1 and cent[0].read_text() == man.read_text()
    # outside verification: restore the deleted files from the copy and re-hash
    for n in (2, 3):
        dst = tmp_path / f"tmp/_ph0/pwscf.q_{n}/pwscf.wfc1"
        shutil.copy(keep / f"_ph0/pwscf.q_{n}/pwscf.wfc1", dst)
    base = Path.home()
    counts = P.verify_manifest(man, base)
    assert counts == {"MATCH": 2, "MISMATCH": 0, "MISSING": 0}
    (tmp_path / "tmp/_ph0/pwscf.q_2/pwscf.wfc1").write_bytes(b"tampered")
    assert P.verify_manifest(man, base)["MISMATCH"] == 1


def test_prune_dvscf_and_final(tmp_path):
    _ph_tree(tmp_path, nq=3, done_through=3)
    P.prune_q_scratch(tmp_path, with_dvscf=True)
    assert not (tmp_path / "tmp/_ph0/pwscf.q_2/_ph0/pwscf.dvscf1").exists()
    P.prune_final(tmp_path, "wfc")
    assert not (tmp_path / "tmp/pwscf.wfc1").exists()
    assert (tmp_path / "tmp/pwscf.save/wfc1.dat").exists()      # .save never touched by wfc mode
    P.prune_final(tmp_path, "tmp")
    assert not any(p.is_file() for p in (tmp_path / "tmp").rglob("*")) if (tmp_path / "tmp").exists() else True
    n_rows = len((tmp_path / "PRUNED_MANIFEST.tsv").read_text().splitlines()) - 1
    assert n_rows >= 5


def test_manifest_written_before_delete(tmp_path, monkeypatch):
    f = tmp_path / "a.wfc1"
    f.write_bytes(b"x")
    calls = []
    real_unlink = Path.unlink

    def spy(self, *a, **k):
        calls.append(("unlink", (tmp_path / "PRUNED_MANIFEST.tsv").read_text().count("\n")))
        return real_unlink(self, *a, **k)

    monkeypatch.setattr(Path, "unlink", spy)
    P.prune_files([f], tmp_path, "test")
    assert calls == [("unlink", 2)]          # header + record already on disk at unlink time


def test_prune_requires_reason(tmp_path):
    with pytest.raises(ValueError):
        P.prune_files([], tmp_path, "")


def test_disp_hides_home():
    assert P.disp(Path.home() / "x" / "y").startswith("~/")


def test_disk_guard(tmp_path):
    P.disk_guard(tmp_path, 0.0)
    with pytest.raises(P.DiskSpaceLow):
        P.disk_guard(tmp_path, 1e12)


def test_inventory(tmp_path):
    (tmp_path / "a").write_text("1")
    (tmp_path / "sub").mkdir()
    (tmp_path / "sub" / "b").write_text("2")
    out = P.inventory(tmp_path)
    rows = out.read_text().splitlines()
    assert len(rows) == 3 and rows[0].endswith("status")


@pytest.mark.skipif(not Path("/proc").is_dir(), reason="needs /proc")
def test_descendants_and_watcher(tmp_path):
    import subprocess
    p = subprocess.Popen(["sleep", "5"])
    try:
        assert p.pid in P.descendants(os.getpid())
    finally:
        p.kill()
        p.wait()
    _ph_tree(tmp_path, nq=3, done_through=3)
    with P.ScratchWatcher(tmp_path, interval_s=0.05) as w:
        pass
    assert w.pruned_bytes == 200
