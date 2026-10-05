# Licensed under the Apache License, Version 2.0 (the "License");
# https://www.apache.org/licenses/LICENSE-2.0

"""Prediction-only bridge to the pinned, separately installed CleaveNet source."""

import json
import os
from pathlib import Path
import sys


def main():
    root, input_path, output_path = map(Path, sys.argv[1:])
    # Upstream locates its vocabulary splits relative to the working directory.
    # Refuse missing splits so DataLoader never falls back to creating a dataset.
    for name in ("X_all.csv", "y_all.csv"):
        if not (root / "splits/kukreja" / name).is_file():
            raise FileNotFoundError("Missing CleaveNet vocabulary split: %s" % name)
    sys.path.insert(0, str(root))
    os.chdir(root)
    import numpy as np
    import tensorflow as tf
    if tf.__version__ != "2.18.0":
        raise RuntimeError("CleaveNet requires the reviewed TensorFlow 2.18.0 runtime")
    import cleavenet
    from cleavenet.utils import mmps

    sequences = json.loads(input_path.read_text())
    loader = cleavenet.data.DataLoader(
        str(root / "data/kukreja.csv"), seed=0, task="regression",
        model="transformer", test_split=0, dataset="kukreja")
    tokens = cleavenet.data.tokenize_sequences(sequences, loader)
    tokens = np.stack([
        np.append(np.array(loader.char2idx[loader.CLS]), item)
        for item in tokens])
    outputs = []
    for index in range(5):
        model = cleavenet.models.load_predictor_model(
            "transformer", str(root / ("weights/transformer_%d/model.h5" % index)),
            batch_size=len(sequences), mask_zero=True)
        outputs.append(np.asarray(model(tokens, training=False)))
    means, deviations = cleavenet.analysis.confidence_score(np.stack(outputs), mmps)
    Path(output_path).write_text(json.dumps({
        "sequences": sequences, "enzymes": mmps,
        "means": means.tolist(), "ensemble_sd": deviations.tolist(),
        "runtime": {"tensorflow": tf.__version__, "python": sys.version.split()[0]},
    }, allow_nan=False))


if __name__ == "__main__":
    main()
