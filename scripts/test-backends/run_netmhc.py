#!/usr/bin/env python3
"""Run installed legacy NetMHC tools in the integration-test Linux runtime."""

import os
from pathlib import Path
import signal
import subprocess
import sys
from uuid import uuid4


def run_container(arguments):
    """Own the container independently of the attached Docker client."""
    name = "mhctools-netmhc-" + uuid4().hex

    def interrupted(signum, frame):
        raise SystemExit(128 + signum)

    previous = signal.signal(signal.SIGTERM, interrupted)
    returncode = 1
    try:
        # Create before starting: interruption during creation cannot leave a
        # running container whose identity the launcher has not recorded yet.
        returncode = subprocess.run(
            ["docker", "create", "--name", name] + arguments,
            stdout=subprocess.DEVNULL).returncode
        if returncode == 0:
            returncode = subprocess.call(["docker", "start", "--attach", name])
    finally:
        signal.signal(signal.SIGTERM, signal.SIG_IGN)
        try:
            removed = subprocess.run(
                ["docker", "rm", "--force", name], stdout=subprocess.DEVNULL,
                # Finish within the batch runner's five-second signal grace.
                stderr=subprocess.PIPE, text=True, timeout=3)
            if removed.returncode and "No such container" not in removed.stderr:
                print("Could not remove %s: %s" % (name, removed.stderr[-500:].strip()),
                      file=sys.stderr)
                returncode = returncode or 1
        except (OSError, subprocess.TimeoutExpired) as error:
            print("Could not remove %s: %s" % (name, error), file=sys.stderr)
            returncode = returncode or 1
        finally:
            signal.signal(signal.SIGTERM, previous)
    return returncode


def main():
    if len(sys.argv) < 2 or sys.argv[1] not in ("netMHC-3.4", "netMHCcons"):
        raise SystemExit("Usage: run_netmhc.py {netMHC-3.4,netMHCcons} [arguments]")
    configured = os.environ.get("NETMHC_BUNDLE_HOME")
    if not configured:
        raise SystemExit("Set NETMHC_BUNDLE_HOME to your installed licensed bundle")
    bundle = Path(configured).expanduser().resolve()
    backend = sys.argv[1]
    if not (bundle / "bin" / backend).is_file():
        raise SystemExit("Missing bundle launcher: %s" % (bundle / "bin" / backend))
    command = [
        "--pull", "never",
        "--platform", "linux/amd64", "--network", "none", "--read-only",
        "--cap-drop", "ALL", "--security-opt", "no-new-privileges",
        "--tmpfs", "/tmp:rw,nosuid,size=256m",
        "--volume", "%s:/netmhc-bundle:ro" % bundle,
    ]
    arguments = []
    # Predictors pass input files by absolute path. Mount just those files;
    # all upstream scratch data stays inside the disposable container.
    original_arguments = sys.argv[2:]
    for index, argument in enumerate(original_arguments):
        path = Path(argument)
        if index and original_arguments[index - 1] == "-tdir":
            arguments.append("/tmp/netmhc")
        elif path.is_file():
            destination = "/inputs/%d%s" % (index, path.suffix)
            command.extend(("--volume", "%s:%s:ro" % (path.resolve(), destination)))
            arguments.append(destination)
        else:
            arguments.append(argument)
    command.extend((
        os.environ.get("MHCTOOLS_LEGACY_NETMHC_IMAGE", "mhctools-test-netmhc-legacy:1"),
        "/bin/sh", "-c", 'mkdir -p /tmp/netmhc; exec "$@"', "sh",
        "/netmhc-bundle/bin/%s" % backend,
    ))
    return run_container(command + arguments)


if __name__ == "__main__":
    sys.exit(main())
