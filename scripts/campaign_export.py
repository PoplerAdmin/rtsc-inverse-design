"""Export a campaign ledger to JSONL for the repo, privacy-sanitized (docs/CONTEXT.md):
the home directory becomes ``~`` and the machine's host name ``host``; the manifest is
recomputed over the sanitized files and the substitution is recorded in SANITIZED.txt.

    python scripts/campaign_export.py LEDGER OUTDIR
"""
import hashlib, os, re, socket, sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from rtt.ledger.store import Ledger

ledger, out = Path(sys.argv[1]).expanduser(), Path(sys.argv[2])
tmp = out.with_name(out.name + ".raw")
Ledger(ledger).export_jsonl(tmp)
home, host, user = str(Path.home()), socket.gethostname(), os.path.basename(str(Path.home()))
out.mkdir(parents=True, exist_ok=True)
lines = []
for p in sorted(tmp.glob("*.jsonl")):
    t = p.read_text(encoding="utf-8").replace(home, "~")
    if host:
        t = t.replace(host, "host")
    t = re.sub(re.escape(user), "user", t, flags=re.IGNORECASE) if len(user) >= 4 else t
    q = out / p.name
    q.write_text(t, encoding="utf-8", newline="\n")
    lines.append(f"{hashlib.sha256(q.read_bytes()).hexdigest()}  {q.name}")
(out / "manifest.sha256").write_text("\n".join(lines) + "\n", encoding="utf-8")
(out / "SANITIZED.txt").write_text(
    "JSONL exported by rtt.ledger.store.Ledger.export_jsonl, then: home directory -> '~', host name -> 'host',\n"
    "user name -> 'user' (docs/CONTEXT.md privacy rule). manifest.sha256 covers the sanitized files.\n",
    encoding="utf-8")
for p in tmp.iterdir():
    p.unlink()
tmp.rmdir()
print(out / "manifest.sha256")
