"""Standalone JSON sidecar; also runs in the legacy Python 3.8 environment.

Only standard-library imports occur during asset inspection. This file must
not import mhctools: the prediction interpreter need not install the host app.
"""

import hashlib
import importlib.metadata
import importlib.util
import json
from pathlib import Path
import platform
import socket
import sys


def _sha256(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def runtime_identity(human_only, model_type):
    """Describe the actual interpreter's assets without importing inference."""
    spec = importlib.util.find_spec("pepsickle")
    if spec is None or spec.origin is None:
        raise ImportError("pepsickle is not installed")
    package_dir = Path(spec.origin).resolve().parent
    weights_path = package_dir / (
        "model.joblib" if model_type == "in-vitro" else "trained_model_dict.pickle")
    inference_path = package_dir / "model_functions.py"
    features_path = package_dir / "sequence_featurization_tools.py"
    if not weights_path.is_file():
        raise RuntimeError("Installed pepsickle is missing model weights at %s" % weights_path)
    packages = {}
    for name in ("pepsickle", "numpy", "scipy", "scikit-learn", "torch", "joblib", "biopython"):
        try:
            packages[name] = importlib.metadata.version(name)
        except importlib.metadata.PackageNotFoundError:
            packages[name] = None
    return {
        "package_version": packages["pepsickle"],
        "weights_path": str(weights_path), "weights_sha256": _sha256(weights_path),
        "inference_path": str(inference_path), "inference_sha256": _sha256(inference_path),
        "features_path": str(features_path), "features_sha256": _sha256(features_path),
        "python": platform.python_version(), "executable": sys.executable,
        "platform": platform.system(), "machine": platform.machine(), "packages": packages,
        "model_key": "gradient_boosting" if model_type == "in-vitro" else "+".join(
            "%s_%s_%s_mod" % (
                "human" if human_only else "all_mammal",
                "epitope" if model_type == "epitope" else "20S_digestion", component)
            for component in ("sequence", "motif")),
    }


def _network_disabled(*args, **kwargs):
    raise RuntimeError("Network access is disabled during Pepsickle inference")


def main():
    request = json.load(sys.stdin)
    identity = runtime_identity(request["human_only"], request["model_type"])
    if request.get("operation") == "identity":
        json.dump({"identity": identity}, sys.stdout)
        return
    # Inference uses installed assets only. The Docker launcher additionally
    # disables networking at the container boundary.
    socket.socket.connect = _network_disabled
    socket.create_connection = _network_disabled
    socket.getaddrinfo = _network_disabled
    from pepsickle.model_functions import (
        initialize_epitope_model, initialize_digestion_model,
        initialize_digestion_gb_model, predict_protein_cleavage_locations,
    )
    initialize = {"epitope": initialize_epitope_model,
                  "in-vitro": initialize_digestion_gb_model,
                  "in-vitro-2": initialize_digestion_model}[request["model_type"]]
    model = initialize(human_only=request["human_only"])
    results = {}
    for sequence in request["sequences"]:
        predictions = predict_protein_cleavage_locations(
            sequence, model, mod_type=request["model_type"],
            proteasome_type=request["proteasome_type"] or "C",
            threshold=request["threshold"])
        results[sequence] = [float(entry[2]) for entry in predictions]
    json.dump({"identity": identity, "results": results}, sys.stdout, allow_nan=False)


if __name__ == "__main__":
    main()
