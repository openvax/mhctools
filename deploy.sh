#!/bin/bash
set -euo pipefail

PYTHON="${PYTHON:-python}"

usage() {
    echo "Usage: $0 [version]" >&2
    echo "Example: $0 3.13.5" >&2
}

die() {
    echo "deploy.sh: $*" >&2
    exit 1
}

current_version() {
    "$PYTHON" - <<'PY'
import re
from pathlib import Path

text = Path("mhctools/__init__.py").read_text()
match = re.search(r'^__version__ = "([^"]+)"$', text, re.MULTILINE)
if not match:
    raise SystemExit("Unable to find mhctools.__version__")
print(match.group(1))
PY
}

bump_version() {
    local version="$1"
    "$PYTHON" - "$version" <<'PY'
import re
import sys
from pathlib import Path

version = sys.argv[1]
if not re.fullmatch(r"[0-9]+[.][0-9]+[.][0-9]+", version):
    raise SystemExit(
        "Version must be X.Y.Z without a leading 'v', got %r" % version)

path = Path("mhctools/__init__.py")
text = path.read_text()
new_text, count = re.subn(
    r'^__version__ = "[^"]+"$',
    '__version__ = "%s"' % version,
    text,
    count=1,
    flags=re.MULTILINE)
if count != 1:
    raise SystemExit("Unable to update mhctools.__version__")
path.write_text(new_text)
PY
}

require_clean_tree() {
    git update-index -q --refresh
    git diff --quiet || die "working tree has unstaged changes"
    git diff --cached --quiet || die "working tree has staged changes"
}

require_release_branch() {
    local branch
    branch="$(git rev-parse --abbrev-ref HEAD)"
    if [[ "$branch" != "main" && "$branch" != "master" ]]; then
        die "must deploy from main or master, currently on $branch"
    fi
    echo "$branch"
}

require_tooling() {
    "$PYTHON" -m build --version >/dev/null ||
        die "python module 'build' is unavailable; install the dev extra"
    "$PYTHON" -m twine --version >/dev/null ||
        die "python module 'twine' is unavailable; install the dev extra"
}

if [[ "$#" -gt 1 ]]; then
    usage
    exit 2
fi

if [[ "${1:-}" == "-h" || "${1:-}" == "--help" ]]; then
    usage
    exit 0
fi

branch="$(require_release_branch)"
require_clean_tree
require_tooling

requested_version="${1:-}"
if [[ -n "$requested_version" ]]; then
    requested_version="${requested_version#v}"
    old_version="$(current_version)"
    if [[ "$requested_version" != "$old_version" ]]; then
        bump_version "$requested_version"
        ./lint.sh
        ./test.sh
        git add mhctools/__init__.py
        git commit -m "Bump to $requested_version"
        git push origin "$branch"
    else
        ./lint.sh
        ./test.sh
    fi
else
    ./lint.sh
    ./test.sh
fi

release_version="$(current_version)"
tag="v$release_version"

git fetch --tags origin
if git rev-parse -q --verify "refs/tags/$tag" >/dev/null; then
    die "tag $tag already exists"
fi

rm -rf dist
"$PYTHON" -m build
"$PYTHON" -m twine upload dist/*

git tag -a "$tag" -m "Release $tag"
git push origin "$tag"
