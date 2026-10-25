"""Batch cancellation must reach containers, not just their host clients."""

import os
from pathlib import Path
import re
import shutil
import subprocess
import sys

import pytest

from mhctools.process_helpers import run_multiple_commands_redirect_stdout


LAUNCHER = Path(__file__).resolve().parents[1] / "scripts/test-backends/run_netmhc.py"


@pytest.fixture(scope="module")
def legacy_image():
    if not shutil.which("docker"):
        pytest.skip("Docker is not installed")
    image = os.environ.get("MHCTOOLS_LEGACY_NETMHC_IMAGE", "mhctools-test-netmhc-legacy:1")
    probe = subprocess.run(
        ["docker", "run", "--rm", "--pull", "never", "--platform", "linux/amd64",
         "--network", "none", image, "true"], capture_output=True, text=True, timeout=30)
    if "No such image" in probe.stderr or "Unable to find image" in probe.stderr:
        pytest.skip("Provision the legacy Docker runtime with setup_test_backends.py")
    assert probe.returncode == 0, probe.stderr
    return image


def _bundle(tmp_path, source):
    bundle = tmp_path / "bundle"
    launcher = bundle / "bin/netMHCcons"
    launcher.parent.mkdir(parents=True)
    launcher.write_text("#!/bin/sh\n" + source)
    launcher.chmod(0o755)
    return bundle


@pytest.mark.requires_external_tool
def test_failed_batch_removes_stubborn_container(tmp_path, legacy_image):
    bundle = _bundle(
        tmp_path, "trap '' TERM INT\nprintf 'container-ready:%s\\n' \"$(hostname)\"\n"
        "while :; do sleep 1; done\n")
    environment = dict(os.environ, NETMHC_BUNDLE_HOME=str(bundle),
                       MHCTOOLS_LEGACY_NETMHC_IMAGE=legacy_image)
    output_path = tmp_path / "container-output"
    failure = tmp_path / "fail.py"
    failure.write_text(
        "from pathlib import Path\nimport sys, time\n"
        "deadline = time.monotonic() + 30\n"
        "while 'container-ready:' not in Path(sys.argv[1]).read_text():\n"
        "    if time.monotonic() > deadline: sys.exit(8)\n"
        "    time.sleep(.05)\n"
        "sys.exit(7)\n")
    try:
        with output_path.open("w") as output, (tmp_path / "failure-output").open("w") as failed:
            with pytest.raises(subprocess.CalledProcessError) as error:
                run_multiple_commands_redirect_stdout(
                    {output: [sys.executable, str(LAUNCHER), "netMHCcons"],
                     failed: [sys.executable, str(failure), str(output_path)]},
                    env=environment, process_limit=2, polling_freq=0.05)
            assert error.value.returncode == 7, output_path.read_text()
        container = re.search(r"container-ready:([0-9a-f]+)", output_path.read_text())[1]
        remaining = subprocess.check_output(
            ["docker", "ps", "--all", "--quiet", "--filter", "id=" + container], text=True)
        assert not remaining.strip(), "Batch returned with its container still present"
    finally:
        match = re.search(r"container-ready:([0-9a-f]+)", output_path.read_text())
        if match:
            # A failed regression may remove only the container it just created.
            subprocess.run(["docker", "rm", "--force", match[1]],
                           stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, timeout=15)


@pytest.mark.parametrize("exit_code", [0, 13])
@pytest.mark.requires_external_tool
def test_container_launcher_preserves_output_and_exit_status(tmp_path, legacy_image, exit_code):
    bundle = _bundle(tmp_path, f"echo 'predictor output'\nexit {exit_code}\n")
    completed = subprocess.run(
        [sys.executable, str(LAUNCHER), "netMHCcons"],
        env=dict(os.environ, NETMHC_BUNDLE_HOME=str(bundle),
                 MHCTOOLS_LEGACY_NETMHC_IMAGE=legacy_image),
        capture_output=True, text=True, timeout=30)
    assert completed.returncode == exit_code, completed.stderr
    assert completed.stdout == "predictor output\n"


@pytest.mark.parametrize("create_status,start_status,cleanup_status,expected", [
    (23, 0, 0, 23), (0, 13, 55, 13), (0, 0, 55, 1),
])
def test_launcher_preserves_errors_and_reports_failed_cleanup(
        tmp_path, create_status, start_status, cleanup_status, expected):
    bundle = _bundle(tmp_path, "exit 0\n")
    commands = tmp_path / "docker-commands"
    fake_bin = tmp_path / "bin"
    fake_bin.mkdir()
    docker = fake_bin / "docker"
    docker.write_text(
        f"#!{sys.executable}\nimport sys\nfrom pathlib import Path\n"
        f"with Path({str(commands)!r}).open('a') as output: output.write(sys.argv[1] + '\\n')\n"
        f"statuses = {{'create': {create_status}, 'start': {start_status}, 'rm': {cleanup_status}}}\n"
        "status = statuses[sys.argv[1]]\n"
        "if status: print(sys.argv[1] + ' failed', file=sys.stderr)\n"
        "sys.exit(status)\n")
    docker.chmod(0o755)
    completed = subprocess.run(
        [sys.executable, str(LAUNCHER), "netMHCcons"],
        env=dict(os.environ, NETMHC_BUNDLE_HOME=str(bundle),
                 PATH=str(fake_bin) + os.pathsep + os.environ.get("PATH", os.defpath)),
        capture_output=True, text=True, timeout=10)
    assert completed.returncode == expected, completed.stderr
    assert commands.read_text().splitlines() == (
        ["create", "rm"] if create_status else ["create", "start", "rm"])
    if cleanup_status:
        assert "Could not remove mhctools-netmhc-" in completed.stderr
        assert "rm failed" in completed.stderr
    else:
        assert "create failed" in completed.stderr
