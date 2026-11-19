"""Shared loader for packaged files under mhctools/data/.

Every reader of curated data (motif rule coefficients, reference catalogs,
lineage inventories, benchmark panels, licenses) goes through this module so a
missing or corrupt packaged file fails the same clear way everywhere.
"""

import hashlib
from importlib.resources import files
import json


def load_resource_bytes(filename):
    """Read raw bytes of one file packaged under ``mhctools/data/<filename>``."""
    return files("mhctools.data").joinpath(filename).read_bytes()


def load_json_resource(filename, *, expected_sha256=None, error_message=None):
    """Parse one JSON file packaged under ``mhctools/data/<filename>``.

    When ``expected_sha256`` is supplied, verify the raw file bytes against
    that SHA-256 hex digest before parsing.
    """
    raw = load_resource_bytes(filename)
    if expected_sha256 is not None and hashlib.sha256(raw).hexdigest() != expected_sha256:
        raise ValueError(error_message or ("Bundled resource checksum mismatch: " + filename))
    return json.loads(raw.decode("utf-8"))
