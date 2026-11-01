"""DFT / DFPT stages S4-S7 (architecture §5.1 ``rtt.dft``). Builder notes:
docs/design/22_build_back.md."""
from .profiles import PROFILES, Profile, PseudoError, profiles_from_tiers, sg15_pseudos

__all__ = ["PROFILES", "Profile", "PseudoError", "profiles_from_tiers", "sg15_pseudos"]
