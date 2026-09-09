# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0

"""Isolated runtime bridge to user-provided PlifePred2 and Pfeature checkouts.

Both upstreams are GPLv3, so nothing is vendored and neither is imported into
the mhctools interpreter. This script runs in whatever interpreter owns them
(``PLIFEPRED2_PYTHON``) and does two things:

1. Shells out to Pfeature's ``pfeature_comp.py`` for the quasi-sequence-order
   descriptor. That script reads ``Data/Schneider-Wrede.csv`` and
   ``Data/Grantham.csv`` by relative path, so it must run with its own
   directory as the working directory.
2. Loads PlifePred2's natural-peptide RandomForest and predicts.

The model is loaded with joblib, which unpickles; the caller is responsible for
pointing at checkouts it trusts.
"""

from argparse import ArgumentParser
import csv
from pathlib import Path
import subprocess
import sys
import tempfile


def _parse_args():
    parser = ArgumentParser()
    parser.add_argument("--pfeature-home", required=True)
    parser.add_argument("--model", required=True)
    parser.add_argument("--input", required=True)
    parser.add_argument("--output", required=True)
    return parser.parse_args()


def _run_pfeature(pfeature_home, fasta_path, features_path):
    """Run Pfeature's QSO job, returning nothing but writing *features_path*."""
    script = Path(pfeature_home) / "pfeature_comp.py"
    command = [
        sys.executable,
        # runpy rather than the script path so that the Pfeature directory is
        # not prepended to sys.path in this process's child.
        "-c",
        ("import runpy,sys; sys.argv=sys.argv[1:]; "
         "runpy.run_path(sys.argv[0],run_name='__main__')"),
        str(script),
        "-i", str(fasta_path),
        "-o", str(features_path),
        "-j", "QSO",
    ]
    process = subprocess.run(
        command,
        cwd=str(pfeature_home),
        capture_output=True,
        text=True,
    )
    if process.returncode != 0:
        raise RuntimeError(
            "Pfeature QSO extraction failed (exit %d):\n%s" % (
                process.returncode,
                (process.stderr or process.stdout).strip()))


def main():
    args = _parse_args()

    import joblib
    import pandas as pd

    with open(args.input, newline="") as handle:
        rows = list(csv.DictReader(handle))

    with tempfile.TemporaryDirectory(prefix="mhctools_pfeature_") as tmp:
        fasta_path = Path(tmp) / "input.fasta"
        features_path = Path(tmp) / "qso.csv"
        with open(fasta_path, "w") as handle:
            for row in rows:
                handle.write(">%s\n%s\n" % (row["__mhctools_id"], row["peptide"]))
        _run_pfeature(args.pfeature_home, fasta_path, features_path)
        features = pd.read_csv(features_path)

    if len(features) != len(rows):
        raise RuntimeError(
            "Pfeature returned %d feature rows for %d peptides"
            % (len(features), len(rows)))

    model = joblib.load(args.model)
    expected = list(model.feature_names_in_)
    missing = [name for name in expected if name not in features.columns]
    if missing:
        raise RuntimeError(
            "Pfeature output is missing feature(s) the model requires: %s. "
            "Is the Pfeature checkout the version PlifePred2 was built "
            "against?" % missing[:10])
    # Raw model output is log10(half-life in seconds); the wrapper converts.
    predictions = model.predict(features[expected])

    with open(args.output, "w", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=["__mhctools_id", "peptide", "log10_seconds"])
        writer.writeheader()
        for row, prediction in zip(rows, predictions):
            writer.writerow({
                "__mhctools_id": row["__mhctools_id"],
                "peptide": row["peptide"],
                "log10_seconds": float(prediction),
            })


if __name__ == "__main__":
    main()
