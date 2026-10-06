#!/usr/bin/env python3
"""Provision optional integration tests; see docs/testing.md for prerequisites."""

import argparse
import hashlib
import json
import os
from pathlib import Path
import platform
import shlex
import shutil
import subprocess
import sys
import tarfile
import urllib.request
import zipfile


ROOT = Path(__file__).resolve().parents[1]


def run(*args, **kwargs):
    subprocess.run([str(arg) for arg in args], check=True, **kwargs)


def capture(*args):
    return subprocess.check_output([str(arg) for arg in args], text=True).strip()


def checkout(destination, repository, revision, sparse=None):
    if not (destination / ".git").exists():
        run("git", "init", destination)
        run("git", "-C", destination, "remote", "add", "origin", repository)
    if capture("git", "-C", destination, "status", "--porcelain", "--untracked-files=no"):
        raise SystemExit("Refusing to overwrite modified checkout: %s" % destination)
    run("git", "-C", destination, "fetch", "--depth=1", "origin", revision)
    if sparse:
        run("git", "-C", destination, "sparse-checkout", "set", sparse)
    run("git", "-C", destination, "checkout", "--detach", revision)


def make_env(root, name, python, requirements, torch=False):
    destination = root / name
    interpreter = destination / "bin/python"
    if not interpreter.exists():
        run(python, "-m", "venv", destination)
    # uv-created environments may not contain pip yet.
    run(interpreter, "-m", "ensurepip")
    if torch and platform.system() == "Linux":
        run(interpreter, "-m", "pip", "install", "torch", "torchvision",
            "--index-url", "https://download.pytorch.org/whl/cpu")
    run(interpreter, "-m", "pip", "install", *requirements)
    return interpreter


def launcher(path, command, extra=""):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("#!/bin/sh\n" + extra + "exec " + shlex.join(map(str, command)) + ' "$@"\n')
    path.chmod(0o755)


def fetch(backend, *options):
    return json.loads(capture("mhctools", "fetch", backend, "--json", *options))["path"]


def half_life(root, python, config):
    from mhctools.peptiverse import (
        ESM2_REVISION, UPSTREAM_REVISION, _ESM2_ARTIFACTS, _PEPTIVERSE_ARTIFACTS,
    )

    runtime = make_env(root, "peptiverse-env", python, [
        "torch>=2.1", "transformers==4.46.0", "lightning==2.5.5", "xgboost",
        "scikit-learn", "joblib", "mapie", "pandas", "SmilesPE", "rdkit",
    ], torch=True)
    source = root / "PeptiVerse"
    # Run the download client in its isolated environment as well.
    code = (
        "from huggingface_hub import snapshot_download; import json,sys; "
        "print(snapshot_download(**json.loads(sys.argv[1])))")
    capture(runtime, "-c", code, json.dumps(dict(
        repo_id="ChatterjeeLab/PeptiVerse", revision=UPSTREAM_REVISION,
        local_dir=str(source),
        allow_patterns=list(_PEPTIVERSE_ARTIFACTS) + ["tokenizer/*.py"])))
    esm = capture(runtime, "-c", code, json.dumps(dict(
        repo_id="facebook/esm2_t33_650M_UR50D", revision=ESM2_REVISION,
        allow_patterns=list(_ESM2_ARTIFACTS) + ["model.safetensors"])))
    config.update(PEPTIVERSE_HOME=str(source), PEPTIVERSE_PYTHON=str(runtime),
                  PEPTIVERSE_ESM_HOME=esm)

    runtime = make_env(root, "plifepred2", python, [
        "plifepred2==1.0", "numpy<2", "pandas<3", "tqdm"])
    package = capture(runtime, "-c",
                      "import plifepred2; print(next(iter(plifepred2.__path__)))")
    features = root / "Pfeature"
    checkout(features, "https://github.com/raghavagps/Pfeature.git",
             "93636eb95bed9df2893b7a0c56b1215e648ecdbf", "Standalone")
    config.update(PLIFEPRED2_HOME=package, PLIFEPRED2_PYTHON=str(runtime),
                  PFEATURE_HOME=str(features / "Standalone"))


def recognition(root, python, config):
    # CapHLA executes in-process; install mhctools[caphla] in the caller's env.
    fetch("caphla")
    config["TULIP_HOME"] = fetch("tulip")
    runtime = make_env(root, "tulip", python, [
        "torch", "transformers==4.32.1", "scikit-learn", "pandas", "numpy"], torch=True)
    config["TULIP_PYTHON"] = str(runtime)
    fetch("mixtcrpred", "--accept-license")
    runtime = make_env(root, "mixtcrpred", python, [
        "torch", "torchvision", "pytorch-lightning", "scipy", "scikit-learn", "pandas",
    ], torch=True)
    config["MIXTCRPRED_PYTHON"] = str(runtime)


def keras(root, python, config):
    # DeepImmuno and TLimmuno2 ship Keras 2 era weights. Modern TensorFlow
    # reaches that API through the tf-keras shim, so one runtime serves both;
    # pin the pair, because a mismatched tensorflow/tf-keras imports
    # tensorflow fine and then raises on tensorflow.keras.
    config["DEEPIMMUNO_HOME"] = fetch("deepimmuno")
    config["TLIMMUNO2_HOME"] = fetch("tlimmuno2", "--accept-license")
    config["NETCLEAVE_DIR"] = fetch("netcleave", "--accept-license")
    runtime = make_env(root, "keras-env", python, [
        "tensorflow==2.17.0", "tf-keras==2.17.0", "numpy<2", "pandas<3", "pyarrow",
        "scikit-learn", "biopython", "matplotlib",
    ])
    # NetCleave uses modern Keras in the same environment; only the two legacy
    # wrappers set TF_USE_LEGACY_KERAS for their individual subprocesses.
    config.update(DEEPIMMUNO_PYTHON=str(runtime), TLIMMUNO2_PYTHON=str(runtime),
                  NETCLEAVE_PYTHON=str(runtime))


def cleavenet(root, python, config):
    source = Path(fetch("cleavenet"))
    runtime = make_env(root, "cleavenet-env", python,
                       ["-r", str(source / "requirements.txt")])
    config.update(CLEAVENET_HOME=str(source), CLEAVENET_PYTHON=str(runtime))


def nettcr(root, python, config):
    # NetTCR runs in-process. Install mhctools[nettcr] in the caller's env;
    # its LiteRT interpreter needs no TensorFlow or Keras environment.
    config["NETTCR_DIR"] = fetch("nettcr", "--accept-license")


def phagescout(root, python, config):
    # Data only: preserve the host interpreter and dependency versions.
    config["PHAGESCOUT_HOME"] = fetch("phagescout")


def gfeller(root, python, config):
    runtime = make_env(root, "gfeller-env", python, [
        "numpy", "pandas<3", "scipy", "logomaker", "matplotlib"])
    mix = root / "MixMHCpred"
    checkout(mix, "https://github.com/GfellerLab/MixMHCpred.git",
             "0a7f9b9e20d1cf02236f4a0a90d16735be879b38")
    # Exercise the supported runtime configuration, with no custom launcher.
    config.update(MIXMHCPRED_PATH=str(mix / "MixMHCpred"),
                  MIXMHCPRED_V3_PATH=str(mix / "MixMHCpred"),
                  MIXMHCPRED_PYTHON=str(runtime))
    prime = root / "PRIME"
    checkout(prime, "https://github.com/GfellerLab/PRIME.git",
             "7b18d4e11042141e7102f7c69be2b0e03d138dab")
    if platform.system() == "Linux":
        # Upstream ships a macOS binary. Build in a separate runtime tree so
        # the pinned source checkout stays clean and setup remains repeatable.
        runtime_tree = root / "PRIME-linux"
        library = runtime_tree / "lib"
        library.mkdir(parents=True, exist_ok=True)
        (runtime_tree / "temp").mkdir(exist_ok=True)
        shutil.copy2(prime / "PRIME", runtime_tree / "PRIME")
        for asset in (prime / "lib").iterdir():
            destination = library / asset.name
            if asset.name != "PRIME.x" and not destination.exists():
                destination.symlink_to(asset)
        run("g++", "-O3", prime / "lib/PRIME.cc", "-o", library / "PRIME.x")
        prime = runtime_tree
    config["PRIME_EXECUTABLE"] = str(prime / "PRIME")

    archive = root / "MixMHC2pred-2.1-beta1.zip"
    if not archive.exists():
        urllib.request.urlretrieve(
            "https://github.com/GfellerLab/MixMHC2pred/releases/download/"
            "v2.1.beta1.2/MixMHC2pred-2.1-beta1.zip", archive)
    digest = hashlib.sha256(archive.read_bytes()).hexdigest()
    if digest != "fe86e0390c96ca7b4e7b8b68d563c717ce43f4091fe460fc283a33ccd716be74":
        raise SystemExit("Unexpected MixMHC2pred release checksum: %s" % digest)
    mix2 = root / "MixMHC2pred"
    if not (mix2 / "bin/Makefile").exists():
        with zipfile.ZipFile(archive) as release:
            release.extractall(mix2)
    binary = "MixMHC2pred_unix" if platform.system() == "Linux" else "MixMHC2pred"
    (mix2 / binary).chmod(0o755)
    config["MIXMHC2PRED_EXECUTABLE"] = str(mix2 / binary)


# A verbatim subset of the official IEDB MHC-I 3.1.7 bundle, hosted on this
# repo's releases. The upstream 341 MB archive became unreachable from GitHub
# runners on 2026-09-28 (curl exit 28, connection timeout, four attempts) with
# the Actions cache evicted, which blocked unrelated merges. The release
# unpacks to 1031 MB, mostly bundled DTU executables SMM never runs, against
# 9.4 MB it needs. scripts/build_iedb_smm_subset.py
# derives this file reproducibly and documents exactly what it contains; the
# upstream Non-Profit OSL 3.0 license ships inside it.
SMM_SUBSET_URL = (
    "https://github.com/openvax/mhctools/releases/download/"
    "iedb-smm-subset-3.1.7/IEDB_MHC_I-3.1.7-smm-subset.tar.gz")
SMM_SUBSET_SHA256 = (
    "eef720a71991a76ceee64ce32fbb52e19684b91d9d44eb078c18c96f395e192e")


def smm(root, python, config):
    """Install the official Python-only SMM runtime, without DTU executables."""
    archive = root / "IEDB_MHC_I-3.1.7-smm-subset.tar.gz"
    candidate = archive
    if not archive.exists():
        candidate = archive.with_name(archive.name + ".part")
        # curl retries transient transfer failures (including connection
        # timeouts), but not authorization errors such as HTTP 403.
        run("curl", "--fail", "--location", "--retry", "3", "--retry-delay", "2",
            "--retry-max-time", "300", "--connect-timeout", "20", "--max-time", "300",
            "--output", candidate, SMM_SUBSET_URL)
    digest = hashlib.sha256(candidate.read_bytes()).hexdigest()
    if digest != SMM_SUBSET_SHA256:
        raise SystemExit("Unexpected IEDB SMM subset checksum: %s" % digest)
    if candidate != archive:
        candidate.replace(archive)
    installation = root / "iedb-3.1.7"
    # Machines provisioned from the full release still hold ~190 MB of method
    # data this subset does not install. Leaving it in place wastes the disk
    # and, worse, makes record_smm_fixtures.py hash a hybrid tree that no
    # fresh install can reproduce, so the recorded provenance would be
    # machine-specific. Start from an empty directory instead.
    if installation.exists():
        shutil.rmtree(installation)
    # Copy regular files only; compatible with Python 3.9 and no symlink or
    # archive-path traversal. The subset holds only the allowlisted paths and
    # is pinned above by SHA-256; the prefix check is a second constraint so
    # that changing the checksum alone cannot change what lands on disk.
    prefixes = ("mhc_i/src/", "mhc_i/method/allele-info/",
                "mhc_i/method/iedbtools-utilities/",
                "mhc_i/data/MHCI_mhcibinding20130222/")
    permitted_files = ("mhc_i/LIAI_license.txt", "mhc_i/Copenhagen_license.txt",
                       "mhc_i/README")
    with tarfile.open(archive) as release:
        for member in release:
            if not member.isfile():
                continue
            if not (member.name.startswith(prefixes)
                    or member.name in permitted_files):
                raise SystemExit(
                    "Unexpected path in the SMM subset: %s" % member.name)
            destination = installation / member.name
            if not destination.resolve().is_relative_to(installation.resolve()):
                raise ValueError("Unsafe IEDB archive path: %s" % member.name)
            destination.parent.mkdir(parents=True, exist_ok=True)
            with release.extractfile(member) as source, destination.open("wb") as target:
                shutil.copyfileobj(source, target)
    source = installation / "mhc_i"
    # The upstream configure script only substitutes this installation path.
    # Its other steps configure licensed DTU binaries that SMM does not use.
    template = (source / "src/setupinfo.template").read_text()
    (source / "src/setupinfo.py").write_text(template % str(source))
    entry = root / "bin/iedb-mhci"
    launcher(entry, [sys.executable, source / "src/predict_binding.py"])
    config["IEDB_MHCI_EXECUTABLE"] = str(entry)


def legacy(root, python, config):
    bundle = os.environ.get("NETMHC_BUNDLE_HOME")
    if not bundle or not (Path(bundle) / "bin/netMHCcons").is_file():
        raise SystemExit("Set NETMHC_BUNDLE_HOME to your installed licensed bundle")
    config["NETMHC_BUNDLE_HOME"] = str(Path(bundle).resolve())
    run("docker", "build", "--platform", "linux/amd64", "-t",
        "mhctools-test-netmhc-legacy:1", "-f",
        ROOT / "scripts/test-backends/Dockerfile.netmhc-legacy",
        ROOT / "scripts/test-backends")
    for name in ("netMHC-3.4", "netMHCcons"):
        launcher(root / "bin" / name,
                 [sys.executable, ROOT / "scripts/test-backends/run_netmhc.py", name])


def pepsickle(root, python, config):
    image = "mhctools-pepsickle-legacy:0.23.2"
    run("docker", "build", "--platform", "linux/amd64", "-t", image, "-f",
        ROOT / "scripts/test-backends/Dockerfile.pepsickle-legacy",
        ROOT / "scripts/test-backends")
    # Bind the launcher to the built image, not a mutable tag. Inference needs
    # no host mounts, network, or installation of mhctools inside the container.
    image_id = capture("docker", "image", "inspect", image, "--format", "{{.Id}}")
    path = root / "bin" / "pepsickle-gb-python"
    launcher(path, ["docker", "run", "--rm", "--platform", "linux/amd64",
                    "--network", "none", "--read-only", "--tmpfs", "/tmp",
                    "-i", image_id, "python"])
    config["PEPSICKLE_GB_PYTHON"] = str(path)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "groups",
        nargs="+",
        choices=[
            "half-life", "recognition", "gfeller", "keras", "nettcr",
            "legacy", "smm", "pepsickle", "cleavenet", "phagescout",
        ],
    )
    parser.add_argument("--python", default="python3.11", help="Python 3.11 interpreter for isolated runtimes")
    parser.add_argument("--root", type=Path, default=ROOT / "env/test-backends")
    parser.add_argument("--accept-license", action="store_true",
                        help="Accept upstream license terms (recognition/gfeller/keras/nettcr/smm)")
    args = parser.parse_args()
    if set(args.groups) & {"recognition", "gfeller", "keras", "nettcr", "smm"} and not args.accept_license:
        parser.error(
            "recognition/gfeller/keras/nettcr/smm require --accept-license; see docs/testing.md")
    root = args.root.expanduser().resolve()
    root.mkdir(parents=True, exist_ok=True)
    state = root / "config.json"
    config = json.loads(state.read_text()) if state.exists() else {}
    functions = {"half-life": half_life, "recognition": recognition,
                 "gfeller": gfeller, "keras": keras, "nettcr": nettcr,
                 "legacy": legacy, "smm": smm, "pepsickle": pepsickle, "cleavenet": cleavenet,
                 "phagescout": phagescout}
    for group in args.groups:
        functions[group](root, args.python, config)
        state.write_text(json.dumps(config, indent=2) + "\n")
    activate = root / "activate.sh"
    activate.write_text(
        "# Generated by scripts/setup_test_backends.py; do not commit.\n"
        + "".join("export %s=%s\n" % (key, shlex.quote(value))
                  for key, value in sorted(config.items()))
        + 'export PATH=%s:"$PATH"\n' % shlex.quote(str(root / "bin")))
    print("Test environment: %s" % activate)


if __name__ == "__main__":
    main()
