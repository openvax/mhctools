"""Documentation metadata follows the checked-out package without importing it."""

from pathlib import Path
import runpy
from types import SimpleNamespace


def test_header_version_tracks_source_changes_without_importing(tmp_path):
    on_config = runpy.run_path(
        str(Path(__file__).resolve().parents[1] / 'scripts/docs_version.py'))['on_config']
    package = tmp_path / 'mhctools'
    package.mkdir()
    source = package / '__init__.py'
    config = SimpleNamespace(config_file_path=str(tmp_path / 'mkdocs.yml'), extra={})
    for version in ['3.46.9', '3.47.0']:
        source.write_text(
            'raise RuntimeError("the docs hook must not import predictors")\n'
            '__version__ = "%s"\n' % version)
        assert on_config(config) is config
        assert config.extra['package_version'] == version
