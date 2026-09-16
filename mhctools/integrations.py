# Copyright (c) 2026 Mount Sinai School of Medicine
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0

"""Bounded capability checks for optional predictor integrations.

Artifact discovery, process launch, and reproduced inference are deliberately
separate claims. In particular, a checkout or executable being present never
implies that its dependencies work or that it can reproduce a prediction.
"""

from dataclasses import asdict, dataclass
import math
import subprocess
import sys
from typing import Optional

from .artifacts import artifact_status, list_artifacts
from .optional_backend import probe_executable


CAPABILITY_LEVELS = ("located", "runnable", "reproduced")


@dataclass(frozen=True)
class IntegrationStatus:
    """Observed capability of one predictor integration.

    ``None`` means a capability was not established; it is distinct from a
    probe that ran and returned ``False``.
    """

    name: str
    located: bool
    runnable: Optional[bool]
    reproduced: Optional[bool]
    capability: str
    path: str
    detail: str

    def to_dict(self):
        """Return a JSON-serializable representation."""
        return asdict(self)

    def meets(self, level):
        """Return whether this record establishes *level*."""
        if level not in CAPABILITY_LEVELS:
            raise ValueError("Unknown capability level %r" % level)
        return getattr(self, level) is True


@dataclass(frozen=True)
class _ProbeResult:
    runnable: bool
    reproduced: Optional[bool]
    detail: str


_EXECUTABLE_ARGS = {
    "mixmhc2pred": ("-h",),
    "mixmhcpred": ("--help",),
    "netmhc": ("-h",),
    "netmhccons": ("-h",),
    "netmhciipan": ("-h",),
    "netmhcpan": ("-h",),
    "netmhcstabpan": ("-h",),
    "prime": ("-h",),
}

_LAUNCH_FAILURE_PATTERNS = (
    "no binaries found",
    "cannot execute binary file",
    "Exec format error",
    "Bad CPU type",
)


def _capability(located, runnable, reproduced):
    if reproduced is True:
        return "reproduced"
    if runnable is True:
        return "runnable"
    if located:
        return "located"
    return "missing"


def _probe_command(path, args, timeout):
    result = probe_executable(
        path,
        args=args,
        timeout=timeout,
        failure_patterns=_LAUNCH_FAILURE_PATTERNS,
    )
    return _ProbeResult(
        runnable=result.runnable,
        reproduced=None,
        detail=result.reason,
    )


def _probe_python_import(module_name, timeout):
    try:
        completed = subprocess.run(
            [sys.executable, "-c", "import %s" % module_name],
            capture_output=True,
            text=True,
            timeout=timeout,
        )
    except (OSError, subprocess.TimeoutExpired) as error:
        return _ProbeResult(False, None, "import probe failed: %s" % error)
    if completed.returncode:
        diagnostic = (completed.stderr or completed.stdout).strip()[-500:]
        return _ProbeResult(
            False, None, "import probe exited with code %d: %s" % (
                completed.returncode, diagnostic or "no diagnostic output"))
    return _ProbeResult(True, None, "package import probe succeeded")


def _probe_calis():
    try:
        from .calis import Calis
        score = Calis().predict(["GILGFVFTL"])[0].preds[0].score
    except Exception as error:  # report optional integration failures verbatim
        return _ProbeResult(False, False, "reference inference failed: %s" % error)
    expected = 0.30484
    reproduced = math.isclose(score, expected, abs_tol=1e-5)
    detail = (
        "reproduced canonical IEDB score %.5f" % expected
        if reproduced else
        "inference ran but score %.8g did not match reference %.5f"
        % (score, expected)
    )
    return _ProbeResult(True, reproduced, detail)


def _probe_netchop():
    # IEDB's public NetChop example and its recorded 3.1 Cterm score. Running
    # this tiny fixture is the only reliable launch probe because netChop has
    # no side-effect-free help/version mode.
    sequence = (
        "MSLLTEVETPIRNEWGCRCNDSSDPLVVAASIIGIVHLILWIIDRLFSKSIYRIFKHGL"
        "KRGPSTEGVPESMREEYREEQQNAVDADDGHFVSIELE"
    )
    try:
        from .netchop import NetChop
        scores = NetChop(model_variant=0).cleavage_probs(sequence)
    except Exception as error:  # report optional integration failures verbatim
        return _ProbeResult(False, False, "reference inference failed: %s" % error)
    expected = 0.976629
    reproduced = (
        len(scores) == len(sequence)
        and math.isclose(scores[95], expected, abs_tol=5e-7)
    )
    detail = (
        "reproduced NetChop 3.1 Cterm reference score %.6f" % expected
        if reproduced else
        "inference ran but did not reproduce the NetChop 3.1 Cterm reference"
    )
    return _ProbeResult(True, reproduced, detail)


def _run_probe(status, timeout):
    if status.name == "calis":
        return _probe_calis()
    if status.name == "netchop":
        return _probe_netchop()
    if status.name in _EXECUTABLE_ARGS:
        return _probe_command(
            status.path, _EXECUTABLE_ARGS[status.name], timeout)
    if status.name in ("mhcflurry", "mhcflurry-affinity"):
        return _probe_python_import("mhcflurry", timeout)
    if status.name == "pepsickle":
        return _probe_python_import("pepsickle", timeout)
    return None


def integration_status(name, check="runnable", data_dir=None, timeout=10):
    """Report located, runnable, and reproduced state for one integration.

    Parameters
    ----------
    name : str
        Artifact/integration name accepted by :func:`artifact_status`.
    check : {"located", "runnable", "reproduced"}
        Highest capability to attempt. Reproduction probes are intentionally
        sparse and reference-backed; unsupported checks remain ``None``.
    data_dir : path-like, optional
        Override the mhctools-managed artifact directory.
    timeout : float
        Bound for command launch and import probes.
    """
    if check not in CAPABILITY_LEVELS:
        raise ValueError(
            "check must be one of %s, got %r" % (CAPABILITY_LEVELS, check))
    if timeout <= 0:
        raise ValueError("timeout must be greater than zero")
    status = artifact_status(name, data_dir=data_dir)
    located = status.status == "ready"
    runnable = None
    reproduced = None
    detail = status.detail
    if not located:
        if check != "located":
            runnable = False
        if check == "reproduced":
            reproduced = False
    elif check != "located":
        result = _run_probe(status, timeout)
        if result is None:
            detail = (
                "No bounded runtime probe is registered; files are located "
                "but execution has not been established")
        else:
            runnable = result.runnable
            detail = result.detail
            if check == "reproduced" or result.reproduced is not None:
                reproduced = result.reproduced
    return IntegrationStatus(
        name=status.name,
        located=located,
        runnable=runnable,
        reproduced=reproduced,
        capability=_capability(located, runnable, reproduced),
        path=status.path if located else "",
        detail=detail,
    )


def list_integrations(
        names=None, check="runnable", data_dir=None, timeout=10):
    """Return capability reports for known or selected integrations."""
    statuses = list_artifacts(names, data_dir=data_dir)
    return [
        integration_status(
            status.name, check=check, data_dir=data_dir, timeout=timeout)
        for status in statuses
    ]
