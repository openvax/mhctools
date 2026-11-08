"""Small on-disk snapshot layouts for discovery and wrapper plumbing tests."""

import json
from pathlib import Path

from mhctools.artifacts import _SNAPSHOTS


def populate_snapshot(root, name, managed=False):
    """Create required asset placeholders without downloading model weights."""
    root = Path(root)
    snapshot = _SNAPSHOTS[name]
    for relative in snapshot.required_paths:
        path = root / relative
        if path.suffix or path.name == "LICENSE":
            path.parent.mkdir(parents=True, exist_ok=True)
            path.touch()
        else:
            path.mkdir(parents=True, exist_ok=True)
    for pattern in snapshot.required_globs:
        path = root / pattern.replace("*", "0")
        path.parent.mkdir(parents=True, exist_ok=True)
        path.touch()
    if managed:
        (root / ".mhctools-artifact.json").write_text(json.dumps({
            "name": name, "repository": snapshot.repository,
            "revision": snapshot.revision,
        }))
    return root
