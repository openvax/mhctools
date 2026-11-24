#!/bin/bash
set -euo pipefail

PYTHON="${PYTHON:-python}"
DEPLOY_SSH_CONNECT_TIMEOUT="${DEPLOY_SSH_CONNECT_TIMEOUT:-15}"
DEPLOY_SSH_SERVER_ALIVE_INTERVAL="${DEPLOY_SSH_SERVER_ALIVE_INTERVAL:-15}"
DEPLOY_SSH_SERVER_ALIVE_COUNT_MAX="${DEPLOY_SSH_SERVER_ALIVE_COUNT_MAX:-3}"
DEPLOY_GIT_REMOTE_TIMEOUT="${DEPLOY_GIT_REMOTE_TIMEOUT:-120}"

usage() {
    echo "Usage: $0 [version]" >&2
    echo "Example: $0 3.13.5" >&2
}

die() {
    echo "deploy.sh: $*" >&2
    exit 1
}

require_positive_int() {
    local name="$1"
    local value="$2"
    if [[ ! "$value" =~ ^[1-9][0-9]*$ ]]; then
        die "$name must be a positive integer, got '$value'"
    fi
}

require_positive_number() {
    local name="$1"
    local value="$2"
    if [[ ! "$value" =~ ^([0-9]+([.][0-9]+)?|[.][0-9]+)$ ]] || \
            ! "$PYTHON" - "$value" <<'PY'
import math
import sys

value = float(sys.argv[1])
raise SystemExit(0 if math.isfinite(value) and value > 0 else 1)
PY
    then
        die "$name must be a positive number, got '$value'"
    fi
}

build_git_ssh_command() {
    local base="${GIT_SSH_COMMAND:-ssh}"
    # OpenSSH uses the first occurrence of each option and still reads
    # ~/.ssh/config when -F is not overridden, so appending these flags
    # respects user SSH config and any explicit GIT_SSH_COMMAND overrides
    # while bounding connect and liveness stalls by default.
    printf '%s -o BatchMode=yes -o ConnectTimeout=%s -o ServerAliveInterval=%s -o ServerAliveCountMax=%s' \
        "$base" \
        "$DEPLOY_SSH_CONNECT_TIMEOUT" \
        "$DEPLOY_SSH_SERVER_ALIVE_INTERVAL" \
        "$DEPLOY_SSH_SERVER_ALIVE_COUNT_MAX"
}

run_git_remote() {
    local step="$1"
    shift
    GIT_SSH_COMMAND="$deploy_git_ssh_command" \
    GIT_TERMINAL_PROMPT=0 \
    "$PYTHON" - "$DEPLOY_GIT_REMOTE_TIMEOUT" "$step" "$@" <<'PY'
import os
import shlex
import signal
import subprocess
import sys
import time

timeout = float(sys.argv[1])
step = sys.argv[2]
command = sys.argv[3:]

process = subprocess.Popen(
    command,
    stdin=subprocess.DEVNULL,
    start_new_session=True,
)


def terminate_group():
    try:
        os.killpg(process.pid, signal.SIGTERM)
    except ProcessLookupError:
        return
    deadline = time.monotonic() + 2.0
    while time.monotonic() < deadline:
        if process.poll() is not None:
            try:
                os.killpg(process.pid, 0)
            except ProcessLookupError:
                return
        time.sleep(0.05)
    try:
        os.killpg(process.pid, signal.SIGKILL)
    except ProcessLookupError:
        pass
    try:
        process.wait(timeout=5.0)
    except subprocess.TimeoutExpired:
        pass


def handle_signal(signum, _frame):
    terminate_group()
    raise SystemExit(128 + signum)


signal.signal(signal.SIGTERM, handle_signal)
signal.signal(signal.SIGINT, handle_signal)

try:
    returncode = process.wait(timeout=timeout)
except subprocess.TimeoutExpired:
    terminate_group()
    print(
        f"deploy.sh: timed out after {timeout:g}s during {step} "
        f"({shlex.join(command)})",
        file=sys.stderr,
        flush=True,
    )
    raise SystemExit(124)

raise SystemExit(returncode)
PY
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

require_unpublished() {
    local phase="${1:-preflight}"
    run_git_remote "fetch tags from origin ($phase)" git fetch --tags origin
    if git rev-parse -q --verify "refs/tags/$tag" >/dev/null; then
        die "tag $tag already exists"
    fi
    "$PYTHON" scripts/check_pypi_release.py mhctools "$release_version"
}

release_lock() {
    local status=$?
    trap - EXIT
    if ! run_git_remote "release lock ($lock_ref)" \
            git push --no-follow-tags --force-with-lease="$lock_ref:$lock_commit" \
            origin ":$lock_ref"; then
        echo "deploy.sh: could not release lock $lock_commit; inspect $lock_ref on origin before recovery (docs/testing.md#release-ownership)" >&2
        if [[ "$status" -eq 0 ]]; then status=1; fi
    fi
    exit "$status"
}

if [[ "$#" -gt 1 ]]; then
    usage
    exit 2
fi

if [[ "${1:-}" == "-h" || "${1:-}" == "--help" ]]; then
    usage
    exit 0
fi

require_positive_int DEPLOY_SSH_CONNECT_TIMEOUT "$DEPLOY_SSH_CONNECT_TIMEOUT"
require_positive_int DEPLOY_SSH_SERVER_ALIVE_INTERVAL "$DEPLOY_SSH_SERVER_ALIVE_INTERVAL"
require_positive_int DEPLOY_SSH_SERVER_ALIVE_COUNT_MAX "$DEPLOY_SSH_SERVER_ALIVE_COUNT_MAX"
require_positive_number DEPLOY_GIT_REMOTE_TIMEOUT "$DEPLOY_GIT_REMOTE_TIMEOUT"
deploy_git_ssh_command="$(build_git_ssh_command)"

branch="$(require_release_branch)"
require_clean_tree

requested_version="${1:-}"
release_version="$(current_version)"
if [[ -n "$requested_version" ]]; then
    requested_version="${requested_version#v}"
    release_version="$requested_version"
fi
tag="v$release_version"

# Stop on an existing tag, published release, or unverifiable PyPI lookup
# before tooling checks, expensive gates, build, upload, or version mutation.
require_unpublished "initial preflight"
require_tooling

# A unique commit makes the remote ref an owner token, even when two deploys
# start from the same HEAD. An empty lease permits creation only, never takeover.
lock_ref="refs/heads/release-lock"
lock_token="$("$PYTHON" -c 'import uuid; print(uuid.uuid4())')"
lock_commit="$(git commit-tree 'HEAD^{tree}' -p HEAD \
    -m "Deploy mhctools $release_version" -m "Owner: $lock_token")"
if ! run_git_remote "acquire release lock ($lock_ref)" \
        git push --no-follow-tags --force-with-lease="$lock_ref:" \
        origin "$lock_commit:$lock_ref"; then
    run_git_remote "inspect release lock ($lock_ref)" \
        git ls-remote origin "$lock_ref" >&2 || true
    die "cannot acquire release lock on origin; inspect $lock_ref (docs/testing.md#release-ownership)"
fi
trap release_lock EXIT
trap 'exit 130' INT
trap 'exit 143' TERM
echo "deploy.sh: acquired release lock $lock_commit for $release_version"

if [[ -n "$requested_version" ]]; then
    old_version="$(current_version)"
    if [[ "$requested_version" != "$old_version" ]]; then
        bump_version "$requested_version"
        ./lint.sh
        ./test.sh
        git add mhctools/__init__.py
        git commit -m "Bump to $requested_version"
        run_git_remote "push version bump to origin/$branch" \
            git push origin "$branch"
    else
        ./lint.sh
        ./test.sh
    fi
else
    ./lint.sh
    ./test.sh
fi

if [[ "$(current_version)" != "$release_version" ]]; then
    die "release version changed unexpectedly during preflight"
fi

rm -rf dist
"$PYTHON" -m build

# Recheck after the long gates/build, including uploads outside this protocol.
remote_lock="$(run_git_remote "verify release lock ownership ($lock_ref)" \
    git ls-remote origin "$lock_ref")"
if [[ "$remote_lock" != "$lock_commit"$'\t'"$lock_ref" ]]; then
    die "release lock ownership changed; refusing to upload"
fi
require_unpublished "pre-upload recheck"
"$PYTHON" -m twine upload dist/*

git tag -a "$tag" -m "Release $tag"
run_git_remote "push release tag $tag to origin" git push origin "$tag"
