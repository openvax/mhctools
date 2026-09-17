#!/usr/bin/env python
"""Fail closed unless PyPI confirms that an exact release is absent."""

import argparse
from http.client import HTTPException
import json
import re
import sys
from urllib.error import HTTPError
from urllib.parse import quote
from urllib.request import Request, urlopen


PROJECT_PATTERN = re.compile(r"^[A-Za-z0-9]+(?:[-_.][A-Za-z0-9]+)*$")
VERSION_PATTERN = re.compile(
    r"^(?:0|[1-9][0-9]*)[.](?:0|[1-9][0-9]*)[.](?:0|[1-9][0-9]*)$"
)


def canonical_project_name(value: str) -> str:
    """Return a normalized distribution name or reject unsafe input."""
    if not PROJECT_PATTERN.fullmatch(value):
        raise ValueError(f"invalid project name: {value!r}")
    return re.sub(r"[-_.]+", "-", value).lower()


def validate_release_version(value: str) -> str:
    """Validate the X.Y.Z version form accepted by deploy.sh."""
    if not VERSION_PATTERN.fullmatch(value):
        raise ValueError(f"version must be X.Y.Z without leading zeroes: {value!r}")
    return value


def pypi_release_exists(project: str, version: str, *, timeout: float = 10) -> bool:
    """Return false only when PyPI reports HTTP 404 for the exact release."""
    canonical_project = canonical_project_name(project)
    release_version = validate_release_version(version)
    url = (
        "https://pypi.org/pypi/"
        f"{quote(canonical_project, safe='')}/{quote(release_version, safe='')}/json"
    )
    request = Request(
        url,
        headers={"Accept": "application/json", "Cache-Control": "no-cache"},
    )
    try:
        with urlopen(request, timeout=timeout) as response:
            if response.status != 200:
                raise ValueError(f"unexpected HTTP status {response.status}")
            payload = json.load(response)
        if not isinstance(payload, dict) or not isinstance(payload.get("info"), dict):
            raise ValueError("response has no release metadata")
        info = payload["info"]
        if (
            not isinstance(info.get("name"), str)
            or canonical_project_name(info["name"]) != canonical_project
        ):
            raise ValueError("response names a different project")
        if info.get("version") != release_version:
            raise ValueError("response names a different version")
        return True
    except HTTPError as error:
        error.close()
        if error.code == 404:
            return False
        raise RuntimeError(
            f"could not verify {canonical_project} {release_version} on PyPI: {error}"
        ) from error
    except (OSError, HTTPException, UnicodeError, ValueError) as error:
        raise RuntimeError(
            f"could not verify {canonical_project} {release_version} on PyPI: {error}"
        ) from error


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("project")
    parser.add_argument("version")
    args = parser.parse_args(argv)
    try:
        exists = pypi_release_exists(args.project, args.version)
    except (ValueError, RuntimeError) as error:
        print(f"deploy.sh: {error}", file=sys.stderr)
        return 1
    if exists:
        print(
            f"deploy.sh: {args.project} {args.version} already exists on PyPI",
            file=sys.stderr,
        )
        return 1
    print(f"Confirmed {args.project} {args.version} is not published on PyPI")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
