"""Shared loader for packaged JSON files under mhctools/data/.

Every reader of curated data (motif rule coefficients, reference catalogs,
lineage inventories, benchmark panels) goes through this one function so a
missing or corrupt packaged file fails the same clear way everywhere.
"""

from importlib.resources import files
import json


def load_json_resource(filename):
    """Parse one JSON file packaged under ``mhctools/data/<filename>``."""
    return json.loads(files("mhctools").joinpath("data/" + filename).read_text(encoding="utf-8"))
