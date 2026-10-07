# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0

"""Isolated runtime bridge to a user-provided PeptiVerse snapshot.

Runs in the interpreter that owns PeptiVerse's dependency stack (torch,
transformers, xgboost and ESM2) so none of it has to be importable from the
mhctools environment.

The caller selects one exact sequence endpoint: log-trained half-life or the
released CPU SVM CPP classifier. The two unused SMILES embedders are replaced
before construction. Half-life remains the default for existing callers.
"""

from argparse import ArgumentParser
import csv
import json
from pathlib import Path
import sys


def _parse_args():
    parser = ArgumentParser()
    parser.add_argument("--home", required=True)
    parser.add_argument("--esm-home", required=True)
    parser.add_argument("--manifest", required=True)
    parser.add_argument("--input", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--device", default="")
    parser.add_argument("--uncertainty", action="store_true")
    parser.add_argument("--endpoint", choices=("half-life", "cpp"), default="half-life")
    return parser.parse_args()


def _verify_tokenized_length(embedder, peptide):
    """Check the exact residue mask that the pinned WTEmbedder uses."""
    tokens = embedder._tokenize([peptide])
    valid = embedder._valid_mask(tokens["input_ids"], tokens["attention_mask"])
    encoded_length = int(valid.sum().item())
    if encoded_length != len(peptide):
        raise RuntimeError(
            "PeptiVerse tokenizer retained %d of %d residues; refusing to "
            "score truncated or altered input" % (encoded_length, len(peptide)))


def main():
    args = _parse_args()
    home = Path(args.home).resolve()
    sys.path.insert(0, str(home))

    if args.endpoint == "cpp":
        # Verify the serialized estimator's supported version BEFORE loading
        # the joblib payload. Cross-version unpickling is not a fallback.
        import sklearn
        if sklearn.__version__ != "1.7.2":
            raise RuntimeError(
                "PeptiVerse CPP requires scikit-learn==1.7.2 in the isolated "
                "interpreter; got %s" % sklearn.__version__)

    # Imported only inside this subprocess; these modules belong to the
    # separately licensed upstream snapshot, not to mhctools.
    import inference

    class _UnusedEmbedder:
        def __init__(self, *args, **kwargs):
            pass

    # Upstream eagerly constructs these two large remote models even when the
    # manifest contains no SMILES endpoint. They cannot affect sequence
    # inference, so prevent both the load and any attempted hub access.
    inference.SMILESEmbedder = _UnusedEmbedder
    inference.ChemBERTaEmbedder = _UnusedEmbedder

    predictor = inference.PeptiVersePredictor(
        manifest_path=args.manifest,
        classifier_weight_root=str(home),
        esm_name=str(Path(args.esm_home).resolve()),
        device=args.device or None,
    )

    cpp = args.endpoint == "cpp"
    prop_key = "permeability_penetrance" if cpp else "halflife"
    meta = predictor.meta.get((prop_key, "wt"))
    expected_artifact = (home / "training_classifiers" /
        ("permeability_penetrance/svm_gpu_wt/best_model.joblib" if cpp else
         "half_life/transformer_wt_log/best_model.pt")).resolve()
    if not meta or Path(meta.get("artifact", "")).resolve() != expected_artifact:
        raise RuntimeError(
            "PeptiVerse loaded %r instead of exact artifact %s"
            % (None if not meta else meta.get("artifact"), expected_artifact))
    expected_meta = ("svm_gpu", "wt", "joblib") if cpp else (
        "transformer_wt_log", "wt", "torch_ckpt")
    if (meta.get("model_name"), meta.get("emb_tag"), meta.get("kind")) != expected_meta:
        raise RuntimeError("Unexpected PeptiVerse model metadata: %r" % meta)
    if cpp:
        from sklearn.svm import SVC
        model = predictor.models[(prop_key, "wt")]
        if (not isinstance(model, SVC) or model.classes_.tolist() != [0, 1] or
                not model.probability or meta.get("threshold") != 0.5493 or
                meta.get("task_type", "").lower() != "classifier"):
            raise RuntimeError("Unexpected PeptiVerse CPP estimator/classes/threshold")

    with open(args.input, newline="") as handle:
        rows = list(csv.DictReader(handle))

    results = []
    for row in rows:
        _verify_tokenized_length(predictor.wt_embedder, row["peptide"])
        # `predict_property` takes (prop_key, col, input_str) positionally --
        # upstream's README examples pass `mode=`, which is not a parameter of
        # the shipped signature and raises TypeError.
        out = predictor.predict_property(
            prop_key,
            "wt",
            row["peptide"],
            uncertainty=args.uncertainty)
        record = {
            "__mhctools_id": row["__mhctools_id"],
            "peptide": row["peptide"],
            # Half-life is already inverse-transformed to hours upstream;
            # CPP is the native positive-class score. Transform neither here.
            "score" if cpp else "hours": out["score"],
            "emb_tag": out.get("emb_tag", ""),
            "uncertainty": (
                "" if out.get("uncertainty") is None
                else json.dumps(out["uncertainty"])),
            "uncertainty_type": out.get("uncertainty_type", ""),
        }
        if cpp:
            record.update(label=out["label"], threshold=out["threshold"])
        results.append(record)

    fieldnames = [
        "__mhctools_id", "peptide", "score" if cpp else "hours",
        "emb_tag", "uncertainty", "uncertainty_type",
    ]
    if cpp:
        fieldnames.extend(("label", "threshold"))
    with open(args.output, "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(results)


if __name__ == "__main__":
    main()
