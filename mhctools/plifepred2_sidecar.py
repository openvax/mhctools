# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0

"""Isolated runtime bridge to user-provided PlifePred2 1.1.0 installations.

Upstream PlifePred2 and its bundled Pfeature script are GPLv3, so nothing is
vendored and neither is imported into the mhctools interpreter. This script
runs in whatever interpreter owns them (``PLIFEPRED2_PYTHON``) and does two
things:

1. Shells out to Pfeature's ``pfeature_comp.py`` for the 13 sequence
   descriptor groups used by PlifePred2 1.1.0. Because that script reads
   ``Data/`` relative to its own path and writes scratch files into its
   working directory, each invocation gets an isolated working directory
   containing only the verified script and reviewed ``Data/`` resources.
2. Selects the 100 features from ``sequence_features.tsv``, sets the five
   natural-peptide modification flags (``d``, ``ct``, ``nt``, ``cyc``, ``ptm``)
   to zero, loads ``Sequence_model.sav`` with joblib, and predicts.

The model is loaded with joblib, which unpickles; the caller is responsible for
pointing at installations it trusts.
"""

from argparse import ArgumentParser
import csv
from functools import lru_cache
import importlib.util
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile
import warnings


_PFEATURE_DATA_FILES = (
    "Schneider-Wrede.csv",
    "Grantham.csv",
    "PhysicoChemical.csv",
    "aaind.txt",
    "aaindex.csv",
    "AAIndexNames.csv",
    "data",
    "aa_attr_group.csv",
)

_SEQUENCE_FEATURE_TYPES = (
    "AAC,DPC,PCP,RRI,DDR,SEP,SER,SPC,CTC,CeTD,PAAC,APAAC,QSO")
_MOD_FLAG_COLUMNS = ("d", "ct", "nt", "cyc", "ptm")
_EXPECTED_SCIKIT_LEARN = "1.7.2"


@lru_cache(maxsize=1)
def _offline_bootstrap():
    """Load ``_OFFLINE_BOOTSTRAP`` from sibling ``optional_backend.py``.

    This sidecar runs in an isolated interpreter where ``mhctools`` is not
    installed, so loading the pure-stdlib sibling module by path reuses the
    single bootstrap definition without importing the ``mhctools`` package.
    """
    spec = importlib.util.spec_from_file_location(
        "_mhctools_optional_backend_bootstrap",
        Path(__file__).with_name("optional_backend.py"),
    )
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module._OFFLINE_BOOTSTRAP


def _parse_args():
    parser = ArgumentParser()
    parser.add_argument("--pfeature-home", required=True)
    parser.add_argument("--model", required=True)
    parser.add_argument("--features", required=True)
    parser.add_argument("--input", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--pfeature-timeout", type=float)
    return parser.parse_args()


def _require_sklearn_version(actual):
    if actual != _EXPECTED_SCIKIT_LEARN:
        raise RuntimeError(
            "PlifePred2 1.1.0 requires scikit-learn %s for verified runtime "
            "reproduction; got %s" % (_EXPECTED_SCIKIT_LEARN, actual))


def _run_pfeature(pfeature_home, fasta_path, features_path, timeout=None):
    """Run Pfeature's sequence descriptor job, writing *features_path*.

    Upstream's script resolves ``Data/*.csv`` relative to ``__file__`` and
    scatters scratch files into ``-wd`` and the current working directory.
    Running directly against the installation would risk leaking intermediates
    or concurrent collisions.

    Each extraction therefore gets a private workspace with a copy of
    ``pfeature_comp.py`` and only the reviewed ``Data/`` resources, and runs
    with both ``-wd`` and ``cwd`` pointed at that scratch directory. The
    installation is only ever read, cleanup can only reach this invocation's
    own scratch files, and concurrent calls against one shared installation
    cannot collide.
    """
    pfeature_home = Path(pfeature_home)
    script = (pfeature_home / "pfeature_comp.py").resolve()
    data_source = pfeature_home / "Data"
    if not data_source.is_dir():
        raise RuntimeError(
            "Pfeature installation %r has no Data directory" % str(pfeature_home))

    with tempfile.TemporaryDirectory(prefix="mhctools_pfeature_run_") as run_dir:
        # Copy only the script and resources in the reviewed inventory.
        # A symlinked tree would let upstream's cleanup and writes reach the
        # installation; a whole-tree copy would expose unverified files.
        run_root = Path(run_dir)
        run_script = run_root / "pfeature_comp.py"
        shutil.copy2(script, run_script)
        runtime_data = run_root / "Data"
        runtime_data.mkdir()
        for name in _PFEATURE_DATA_FILES:
            shutil.copy2(data_source / name, runtime_data / name)
        command = [
            sys.executable,
            "-c",
            _offline_bootstrap(),
            str(run_script),
            "-i", str(Path(fasta_path).resolve()),
            "-o", str(Path(features_path).resolve()),
            "-j", _SEQUENCE_FEATURE_TYPES,
            "-wd", str(run_root),
        ]
        try:
            process = subprocess.run(
                command,
                cwd=run_dir,
                capture_output=True,
                text=True,
                timeout=timeout,
            )
        except subprocess.TimeoutExpired as error:
            raise RuntimeError(
                "Pfeature descriptor extraction timed out after %s seconds"
                % timeout) from error
        if process.returncode != 0:
            raise RuntimeError(
                "Pfeature descriptor extraction failed (exit %d):\n%s" % (
                    process.returncode,
                    (process.stderr or process.stdout).strip()))


def main():
    args = _parse_args()

    import joblib
    import pandas as pd
    import sklearn
    from sklearn.exceptions import InconsistentVersionWarning

    _require_sklearn_version(sklearn.__version__)

    with open(args.input, newline="") as handle:
        rows = list(csv.DictReader(handle))

    selected_descriptors = [
        str(column).strip()
        for column in pd.read_csv(args.features, sep="\t").columns
        if str(column).strip()
    ]
    if not selected_descriptors:
        raise RuntimeError(
            "No feature names found in %r" % args.features)

    with tempfile.TemporaryDirectory(prefix="mhctools_pfeature_") as tmp:
        fasta_path = Path(tmp) / "input.fasta"
        features_path = Path(tmp) / "descriptors.csv"
        with open(fasta_path, "w") as handle:
            for row in rows:
                handle.write(">%s\n%s\n" % (row["__mhctools_id"], row["peptide"]))
        _run_pfeature(
            args.pfeature_home,
            fasta_path,
            features_path,
            timeout=args.pfeature_timeout)
        features = pd.read_csv(features_path)

    if len(features) != len(rows):
        raise RuntimeError(
            "Pfeature returned %d feature rows for %d peptides"
            % (len(features), len(rows)))

    missing_selected = [
        name for name in selected_descriptors if name not in features.columns]
    if missing_selected:
        raise RuntimeError(
            "Pfeature output is missing feature(s) listed in sequence_features.tsv: %s"
            % missing_selected[:10])

    matrix = features.loc[:, selected_descriptors].copy()
    for column in _MOD_FLAG_COLUMNS:
        matrix[column] = 0

    with warnings.catch_warnings():
        # Upstream serialized Sequence_model.sav with scikit-learn 1.6.1 while
        # plifepred2 1.1.0 pins scikit-learn==1.7.2.
        warnings.simplefilter("ignore", category=InconsistentVersionWarning)
        model = joblib.load(args.model)
    expected = list(model.feature_names_in_)
    missing = [name for name in expected if name not in matrix.columns]
    if missing:
        raise RuntimeError(
            "Feature matrix is missing feature(s) the model requires: %s. "
            "Is the Pfeature installation the version PlifePred2 1.1.0 was "
            "built against?" % missing[:10])
    # Raw model output is log10(t); the wrapper optionally converts assuming seconds.
    predictions = model.predict(matrix[expected])

    with open(args.output, "w", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=[
                "__mhctools_id", "peptide", "log10_seconds",
                "scikit_learn_version",
            ])
        writer.writeheader()
        for row, prediction in zip(rows, predictions):
            writer.writerow({
                "__mhctools_id": row["__mhctools_id"],
                "peptide": row["peptide"],
                "log10_seconds": float(prediction),
                "scikit_learn_version": sklearn.__version__,
            })


if __name__ == "__main__":
    main()
