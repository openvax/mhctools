"""Per-process interpreter selection for the official MixMHCpred launcher."""

from contextlib import contextmanager
from dataclasses import asdict, dataclass
import json
import os
from pathlib import Path
import shlex
import shutil
import subprocess
from tempfile import TemporaryDirectory


_SETUP = ("Install numpy 'pandas>=2,<3' scipy logomaker matplotlib in a separate venv; "
          "set MIXMHCPRED_PYTHON to its Python executable or pass python_executable.")
_PROBE = """
import json, platform, sys
import numpy, pandas, scipy, logomaker, matplotlib
print(json.dumps(dict(
    reported_executable=sys.executable, python_version=platform.python_version(),
    pandas_version=pandas.__version__, numpy_version=numpy.__version__,
    scipy_version=scipy.__version__, logomaker_version=logomaker.__version__,
    matplotlib_version=matplotlib.__version__)))
"""


def _isolated_environment():
    environment = dict(os.environ)
    # Packages belong to the selected interpreter, not the caller's import path.
    environment.pop("PYTHONPATH", None)
    environment.pop("PYTHONHOME", None)
    environment.update(PYTHONNOUSERSITE="1", MPLBACKEND="Agg")
    return environment


@dataclass(frozen=True)
class MixMHCpredRuntime:
    """Verified interpreter and dependency versions for isolated inference."""

    python_executable: str
    reported_executable: str
    python_version: str
    pandas_version: str
    numpy_version: str
    scipy_version: str
    logomaker_version: str
    matplotlib_version: str

    def to_dict(self):
        """Return runtime diagnostics without conflating host and backend."""
        return asdict(self)

    @contextmanager
    def environment(self):
        """Bind upstream's python3 calls to this interpreter for one run.

        Keep the shim alive through nested PRIME execution. Do not change
        the parent PATH or assume the chosen executable is named python3.
        """
        with TemporaryDirectory(prefix="mhctools-mixmhcpred-python-") as directory:
            launcher = Path(directory) / "python3"
            launcher.write_text("#!/bin/sh\nexec " + shlex.quote(self.python_executable) + ' "$@"\n')
            launcher.chmod(0o755)
            environment = _isolated_environment()
            environment["PATH"] = directory + os.pathsep + environment.get("PATH", os.defpath)
            yield environment


def resolve_runtime(python_executable=None, timeout=30):
    """Resolve an explicit interpreter or MIXMHCPRED_PYTHON; otherwise inherit.

    Existing upstream/user launchers retain control when no interpreter is
    configured. An explicit selection is probed before any prediction runs.
    """
    candidate = python_executable if python_executable is not None else os.environ.get("MIXMHCPRED_PYTHON")
    if not candidate:
        return None
    executable = shutil.which(str(Path(candidate).expanduser()))
    if executable is None:
        raise FileNotFoundError("MixMHCpred backend Python was not found or is not executable: %s. %s"
                                % (candidate, _SETUP))
    # Resolving this symlink would escape a venv into its base interpreter.
    executable = str(Path(executable).absolute())
    try:
        completed = subprocess.run([executable, "-c", _PROBE],
                                   env=_isolated_environment(), capture_output=True,
                                   text=True, timeout=timeout)
    except (OSError, subprocess.TimeoutExpired) as error:
        raise RuntimeError("Could not probe MixMHCpred backend Python %s: %s. %s"
                           % (executable, error, _SETUP)) from error
    if completed.returncode:
        detail = (completed.stderr or completed.stdout).strip()[-1200:]
        raise RuntimeError("MixMHCpred backend dependencies failed at %s: %s. %s"
                           % (executable, detail, _SETUP))
    try:
        info = json.loads(completed.stdout)
        runtime = MixMHCpredRuntime(python_executable=executable, **info)
        pandas_major = int(runtime.pandas_version.split(".")[0])
    except (ValueError, TypeError, AttributeError) as error:
        raise RuntimeError("Invalid MixMHCpred runtime diagnostics from %s. %s"
                           % (executable, _SETUP)) from error
    if pandas_major >= 3:
        raise RuntimeError("MixMHCpred 3.0 requires backend pandas<3; selected %s has pandas %s. %s"
                           % (executable, runtime.pandas_version, _SETUP))
    return runtime


@contextmanager
def runtime_environment(runtime):
    """Enter an isolated runtime, or retain a user-managed launcher's behavior."""
    if runtime is None:
        yield None
    else:
        with runtime.environment() as environment:
            yield environment
