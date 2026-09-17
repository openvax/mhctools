import importlib.util
from io import BytesIO
import json
from pathlib import Path
from urllib.error import HTTPError

import pytest


ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "scripts" / "check_pypi_release.py"
SPEC = importlib.util.spec_from_file_location("check_pypi_release", SCRIPT)
PREFLIGHT = importlib.util.module_from_spec(SPEC)
assert SPEC.loader is not None
SPEC.loader.exec_module(PREFLIGHT)


def response(payload, status=200):
    result = BytesIO(json.dumps(payload).encode())
    result.status = status
    return result


def test_published_release_is_detected(monkeypatch):
    requested = []

    def get(request, timeout):
        requested.append(request)
        assert timeout == 3
        return response({"info": {"name": "MHC.Tools", "version": "3.44.42"}})

    monkeypatch.setattr(PREFLIGHT, "urlopen", get)
    assert PREFLIGHT.pypi_release_exists("mhc_tools", "3.44.42", timeout=3)
    assert requested[0].full_url.endswith("/mhc-tools/3.44.42/json")
    assert requested[0].get_header("Cache-control") == "no-cache"


def test_only_404_means_release_is_available(monkeypatch):
    def missing(request, timeout):
        raise HTTPError(request.full_url, 404, "Not Found", {}, None)

    monkeypatch.setattr(PREFLIGHT, "urlopen", missing)
    assert not PREFLIGHT.pypi_release_exists("mhctools", "3.44.43")


@pytest.mark.parametrize(
    "error",
    [
        TimeoutError("lookup timed out"),
        HTTPError("https://pypi.org", 503, "Service Unavailable", {}, None),
    ],
)
def test_lookup_failure_is_not_treated_as_available(monkeypatch, error):
    def fail(*args, **kwargs):
        raise error

    monkeypatch.setattr(PREFLIGHT, "urlopen", fail)
    with pytest.raises(RuntimeError, match="could not verify mhctools 3.44.43"):
        PREFLIGHT.pypi_release_exists("mhctools", "3.44.43")


@pytest.mark.parametrize(
    "payload",
    [
        {},
        {"info": None},
        {"info": {"name": "other", "version": "3.44.43"}},
        {"info": {"name": "mhctools", "version": "3.44.430"}},
    ],
)
def test_malformed_or_mismatched_metadata_fails_closed(monkeypatch, payload):
    monkeypatch.setattr(PREFLIGHT, "urlopen", lambda *args, **kwargs: response(payload))
    with pytest.raises(RuntimeError, match="could not verify"):
        PREFLIGHT.pypi_release_exists("mhctools", "3.44.43")


@pytest.mark.parametrize("version", ["", "v3.44.43", "3.44", "03.44.43"])
def test_invalid_version_fails_before_network(monkeypatch, version):
    def unexpected_request(*args, **kwargs):
        pytest.fail("invalid version reached the network")

    monkeypatch.setattr(PREFLIGHT, "urlopen", unexpected_request)
    with pytest.raises(ValueError):
        PREFLIGHT.pypi_release_exists("mhctools", version)
