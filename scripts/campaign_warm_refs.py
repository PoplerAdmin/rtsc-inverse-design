"""Pre-compute S4-profile elemental references (cache shared with S4Electronic.formation_energy)."""
import json, sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from rtt.config import load_machine, load_tiers
from rtt.dft.element_refs import element_refs
from rtt.dft.profiles import profiles_from_tiers
from rtt.dft.workflow import default_launcher
work = Path(sys.argv[1]).expanduser()
els = sys.argv[3:]
prof = profiles_from_tiers(load_tiers())["s4"]
for el in els:
    mu, info = element_refs([el], prof, work / "element_refs", launcher=default_launcher(int(sys.argv[2])),
                            cache_dir=work / "element_refs" / "cache", machine=load_machine())
    print(el, json.dumps(info), flush=True)
