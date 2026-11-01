"""S2 structure stage helpers: canonical form + hash, validity, duplicates, novelty."""
from .canonical import canonicalize, structure_sha256, to_pymatgen, to_sdict  # noqa: F401
from .dedupe import DedupeIndex, find_duplicate, is_novel  # noqa: F401
from .validate import validate  # noqa: F401
