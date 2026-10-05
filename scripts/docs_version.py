"""Render the documentation header with the version in the checked-out source."""

import ast
from pathlib import Path


def on_config(config):
    """Read the package version without importing optional predictor runtimes."""
    source = Path(config.config_file_path).parent / "mhctools" / "__init__.py"
    tree = ast.parse(source.read_text(encoding="utf-8"))
    for node in tree.body:
        if isinstance(node, ast.Assign) and any(
                isinstance(target, ast.Name) and target.id == "__version__"
                for target in node.targets):
            config.extra["package_version"] = ast.literal_eval(node.value)
            return config
    raise ValueError("No mhctools package version found in %s" % source)
