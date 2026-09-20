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


def gfeller(root, python, config):
    runtime = make_env(root, "gfeller-env", python, [
        "numpy", "pandas<3", "scipy", "logomaker", "matplotlib"])
    mix = root / "MixMHCpred"
    checkout(mix, "https://github.com/GfellerLab/MixMHCpred.git",
             "0a7f9b9e20d1cf02236f4a0a90d16735be879b38")
    # Keep the launcher's parent at the checkout root: tests also use the
    # official sequence-alignment fixture shipped beside the executable.
    entry = mix / "mhctools-test-launcher"
    launcher(entry, [mix / "MixMHCpred"],
             'export PATH=%s:"$PATH"\n' % shlex.quote(str(runtime.parent)))
    config.update(MIXMHCPRED_PATH=str(entry), MIXMHCPRED_V3_PATH=str(entry))
    prime = root / "PRIME"
    checkout(prime, "https://github.com/GfellerLab/PRIME.git",
             "7b18d4e11042141e7102f7c69be2b0e03d138dab")
    if platform.system() == "Linux":
        # Upstream ships a macOS binary. Build in a separate runtime tree so
        # the pinned source checkout stays clean and setup remains repeatable.
        runtime_tree = root / "PRIME-linux"
        library = runtime_tree / "lib"
        library.mkdir(parents=True, exist_ok=True)
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


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("groups", nargs="+", choices=["half-life", "recognition", "gfeller", "legacy"])
    parser.add_argument("--python", default="python3.11", help="Python 3.11 interpreter for isolated runtimes")
    parser.add_argument("--root", type=Path, default=ROOT / "env/test-backends")
    parser.add_argument("--accept-license", action="store_true",
                        help="Accept upstream academic/non-commercial terms (recognition/gfeller)")
    args = parser.parse_args()
    if set(args.groups) & {"recognition", "gfeller"} and not args.accept_license:
        parser.error("recognition/gfeller require --accept-license; see docs/testing.md")
    root = args.root.expanduser().resolve()
    root.mkdir(parents=True, exist_ok=True)
    state = root / "config.json"
    config = json.loads(state.read_text()) if state.exists() else {}
    functions = {"half-life": half_life, "recognition": recognition,
                 "gfeller": gfeller, "legacy": legacy}
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
