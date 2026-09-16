# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

import logging
import os
from pathlib import Path
import platform
import shutil
import subprocess
import tempfile

from .proteasome_predictor import ProteasomePredictor

logger = logging.getLogger(__name__)

# Timeout in seconds for a single netChop invocation.
NETCHOP_TIMEOUT_SECONDS = 120

# NetChop 3.1 is distributed as 32-bit x86 Linux binaries.  A digest-pinned
# userspace image lets Docker Desktop/qemu run those user-supplied binaries on
# modern Macs and ARM Linux without redistributing NetChop itself.
NETCHOP_CONTAINER_IMAGE = (
    "i386/debian@sha256:"
    "75efd55b326373cf69989912388c0d50c5390638af7378d2fedc3aeb9d100e46"
)


def _netchop_installation_candidates(program_name, netchop_dir=None):
    """Yield explicit and conventional NetChop 3.1 installation roots."""
    if netchop_dir:
        yield Path(netchop_dir).expanduser()
    configured = os.environ.get("NETCHOP_HOME")
    if configured:
        yield Path(configured).expanduser()
    bundle_home = os.environ.get("NETMHC_BUNDLE_HOME")
    if bundle_home:
        yield Path(bundle_home).expanduser() / "netchop-3.1"
    resolved = shutil.which(str(program_name))
    if resolved:
        executable = Path(resolved).resolve()
        if executable.parent.name == "bin":
            yield executable.parent.parent


def resolve_netchop_dir(program_name="netChop", netchop_dir=None):
    """Resolve a licensed NetChop tree containing ``bin/netChop``."""
    for candidate in _netchop_installation_candidates(
            program_name, netchop_dir=netchop_dir):
        candidate = candidate.resolve()
        if (candidate / "bin" / "netChop").is_file():
            return candidate
    return None


class NetChop(ProteasomePredictor):
    """
    Wrapper around the netChop command-line tool.

    On native x86 Linux, this runs ``netChop`` directly. On other
    architectures, it automatically runs a locally installed, licensed
    NetChop tree in a network-disabled, digest-pinned Docker container.

    Parameters
    ----------
    default_peptide_lengths : list of int, optional
        Peptide lengths used when scanning proteins. Default ``[9]``.

    scoring : callable, optional
        See :class:`ProcessingPredictor`.  Default:
        ``score_cterm_anti_max_internal``.

    program_name : str
        Name or path of the netChop executable (default ``"netChop"``).

    netchop_dir : path-like, optional
        Root of the licensed NetChop installation (the directory containing
        ``bin/netChop``). It can also be configured with ``NETCHOP_HOME`` or
        ``NETMHC_BUNDLE_HOME``.

    execution : {"auto", "native", "container"}
        Execution backend. ``auto`` uses the container on non-x86 Linux and
        macOS, where the upstream 32-bit Linux executable cannot run natively.

    model_variant : {0, 1}, optional
        NetChop's ``-v`` model: 0 is C-terminal epitope-trained and 1 is the
        in-vitro 20S model. Omit to retain the upstream default.
    """

    def __init__(
            self,
            default_peptide_lengths=None,
            scoring=None,
            program_name="netChop",
            netchop_dir=None,
            execution="auto",
            model_variant=None,
            container_image=NETCHOP_CONTAINER_IMAGE):
        ProteasomePredictor.__init__(
            self,
            default_peptide_lengths=default_peptide_lengths,
            scoring=scoring,
        )
        if execution not in ("auto", "native", "container"):
            raise ValueError(
                "execution must be 'auto', 'native', or 'container', got %r"
                % execution)
        if model_variant not in (None, 0, 1):
            raise ValueError("model_variant must be 0 (Cterm), 1 (20S), or None")
        self.program_name = program_name
        self.netchop_dir = resolve_netchop_dir(
            program_name=program_name, netchop_dir=netchop_dir)
        self.model_variant = model_variant
        self.container_image = container_image
        if execution == "auto":
            native_x86_linux = (
                platform.system() == "Linux" and
                platform.machine().lower() in ("i386", "i686", "x86_64", "amd64")
            )
            execution = "native" if native_x86_linux else "container"
        self.execution = execution
        if self.execution == "native":
            if self.program_name == "netChop" and self.netchop_dir:
                resolved = str(self.netchop_dir / "bin" / "netChop")
            else:
                resolved = shutil.which(str(self.program_name))
            if not resolved:
                raise FileNotFoundError(
                    "Could not find '%s' on PATH. Is NetChop installed? "
                    "Pass a full path via program_name=." % self.program_name)
            self.program_name = resolved
        else:
            if not self.netchop_dir:
                raise FileNotFoundError(
                    "Could not find a licensed NetChop installation containing "
                    "bin/netChop. Pass netchop_dir= or set NETCHOP_HOME.")
            if not shutil.which("docker"):
                raise FileNotFoundError(
                    "Docker is required to run the 32-bit NetChop binaries on "
                    "this platform, but 'docker' was not found on PATH")

    def __str__(self):
        return "%s(program_name=%r, execution=%r, scoring=%s)" % (
            self.__class__.__name__,
            self.program_name,
            self.execution,
            getattr(self.scoring, "__name__", repr(self.scoring)))

    def _predictor_name(self):
        return "netchop"

    def cleavage_probs(self, sequence):
        """
        Run netChop on a single sequence.

        Returns
        -------
        list of float
            Per-position cleavage probabilities.
        """
        return self.cleavage_probs_many([sequence])[sequence]

    def cleavage_probs_many(self, sequences):
        """Run NetChop once for a sequence collection.

        Parameters
        ----------
        sequences : iterable of str
            Protein or peptide sequences, retained in input order.

        Returns
        -------
        dict of str to list of float
            Per-position score vectors keyed by input sequence.
        """
        sequences = list(dict.fromkeys(sequences))
        if not sequences:
            return {}
        with tempfile.TemporaryDirectory(prefix="mhctools_netchop_") as tmp:
            fasta_path = Path(tmp) / "sequence.fsa"
            fasta_path.write_text("".join(
                ">seq_%d\n%s\n" % (index, sequence)
                for index, sequence in enumerate(sequences)), encoding="ascii")
            if self.execution == "container":
                command = self._container_command(Path(tmp), fasta_path.name)
                environment = None
            else:
                command = [self.program_name]
                if self.model_variant is not None:
                    command.extend(("-v", str(self.model_variant)))
                command.append(str(fasta_path))
                environment = os.environ.copy()
                if self.netchop_dir:
                    environment["NETCHOP"] = str(self.netchop_dir)
                environment["TMPDIR"] = tmp
            try:
                result = subprocess.run(
                    command,
                    capture_output=True,
                    timeout=NETCHOP_TIMEOUT_SECONDS,
                    env=environment,
                )
            except subprocess.TimeoutExpired:
                raise RuntimeError(
                    "%s timed out after %d seconds on a sequence of "
                    "collection containing %d residues"
                    % (self.execution, NETCHOP_TIMEOUT_SECONDS,
                       sum(map(len, sequences))))
            except FileNotFoundError:
                raise FileNotFoundError(
                    "Could not execute NetChop using the %s backend"
                    % self.execution)
        stderr_text = result.stderr.decode("utf-8", errors="replace").strip()
        if stderr_text:
            logger.warning("%s stderr:\n%s", self.program_name, stderr_text)
        if result.returncode != 0:
            image_hint = ""
            if (self.execution == "container" and
                    b"No such image" in result.stderr):
                image_hint = (
                    "\nPreload the pinned compatibility image with: docker "
                    "pull --platform linux/386 %s" % self.container_image)
            raise RuntimeError(
                "%s exited with code %d.\nstdout: %s\nstderr: %s%s"
                % (self.program_name, result.returncode,
                   result.stdout.decode("utf-8", errors="replace").strip(),
                   stderr_text, image_hint))
        parsed = self.parse_netchop(result.stdout)
        if len(parsed) != len(sequences):
            raise ValueError(
                "Expected %d result sequences from %s, got %d. "
                "stdout: %s\nstderr: %s"
                % (len(sequences), self.program_name, len(parsed),
                   result.stdout.decode("utf-8", errors="replace").strip(),
                   stderr_text))
        for index, (sequence, scores) in enumerate(zip(sequences, parsed)):
            if len(scores) != len(sequence):
                raise ValueError(
                    "Expected %d per-position scores from %s for sequence %d, "
                    "got %d. This usually means netChop's internal temp files "
                    "are broken — try reinstalling NetChop or use pepsickle "
                    "as an alternative (--mhc-predictor pepsickle).\n"
                    "stdout: %s\nstderr: %s"
                    % (len(sequence), self.program_name, index, len(scores),
                       result.stdout.decode(
                           "utf-8", errors="replace").strip(), stderr_text))
        return dict(zip(sequences, parsed))

    def _container_command(self, work_dir, fasta_name):
        """Build the offline container invocation for licensed local files."""
        command = [
            "docker", "run", "--rm",
            "--pull", "never",
            "--platform", "linux/386",
            "--network", "none",
            "--read-only",
            "--cap-drop", "ALL",
            "--security-opt", "no-new-privileges",
            "--tmpfs", "/tmp:rw,nosuid,size=64m",
            "-e", "NETCHOP=/netchop",
            "-e", "TMPDIR=/tmp",
            "-v", "%s:/netchop:ro" % self.netchop_dir,
            "-v", "%s:/work:ro" % work_dir.resolve(),
            self.container_image,
            "/netchop/bin/netChop",
        ]
        if self.model_variant is not None:
            command.extend(("-v", str(self.model_variant)))
        command.append("/work/%s" % fasta_name)
        return command

    @staticmethod
    def parse_netchop(netchop_output):
        """
        Parse netChop stdout.

        Returns
        -------
        list of list of float
            One inner list per input sequence, with per-position
            cleavage scores.
        """
        text = netchop_output.decode("utf-8", errors="replace")
        lines = text.split("\n")
        line_iterator = iter(lines)
        scores = []
        for line in line_iterator:
            if "pos" in line and 'AA' in line and 'score' in line:
                scores.append([])
                dashes_line = next(line_iterator, "")
                if "----" not in dashes_line:
                    raise ValueError(
                        "Expected dashes after netChop header, got: %r"
                        % dashes_line)
                line = next(line_iterator, "-------")
                while '-------' not in line:
                    parts = line.split()
                    if len(parts) < 4:
                        raise ValueError(
                            "Unexpected netChop output line: %r" % line)
                    try:
                        score = float(parts[3])
                    except ValueError:
                        raise ValueError(
                            "Could not parse score from netChop "
                            "output line: %r" % line)
                    scores[-1].append(score)
                    line = next(line_iterator, "-------")
        return scores
