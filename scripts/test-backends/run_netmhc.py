#!/usr/bin/env python3
"""Run installed legacy NetMHC tools in the integration-test Linux runtime."""

import os
from pathlib import Path
import subprocess
import sys


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
        "docker", "run", "--rm", "--pull", "never",
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
    return subprocess.call(command + arguments)


if __name__ == "__main__":
    sys.exit(main())
