# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0

"""Isolated runtime bridge to a user-provided PeptiVerse snapshot.

Runs in the interpreter that owns PeptiVerse's dependency stack (torch,
transformers, xgboost, the ESM2 / PeptideCLM / ChemBERTa embedding models) so
none of it has to be importable from the mhctools environment.

Only the half-life endpoint is exposed. The manifest written by the caller lists
``Halflife`` and nothing else, so ``PeptiVersePredictor`` loads one model
instead of the nine it would otherwise pull in for hemolysis, solubility,
permeability, toxicity, non-fouling and binding affinity.
"""

from argparse import ArgumentParser
import csv
import json
from pathlib import Path
import sys


def _parse_args():
    parser = ArgumentParser()
    parser.add_argument("--home", required=True)
    parser.add_argument("--manifest", required=True)
    parser.add_argument("--input", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--device", default="")
    parser.add_argument("--uncertainty", action="store_true")
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

    # Imported only inside this subprocess; these modules belong to the
    # separately licensed upstream snapshot, not to mhctools.
    from inference import PeptiVersePredictor

    predictor = PeptiVersePredictor(
        manifest_path=args.manifest,
        classifier_weight_root=str(home),
        device=args.device or None,
    )

    with open(args.input, newline="") as handle:
        rows = list(csv.DictReader(handle))

    results = []
    for row in rows:
        _verify_tokenized_length(predictor.wt_embedder, row["peptide"])
        # `predict_property` takes (prop_key, col, input_str) positionally --
        # upstream's README examples pass `mode=`, which is not a parameter of
        # the shipped signature and raises TypeError.
        out = predictor.predict_property(
            "halflife",
            "wt",
            row["peptide"],
            uncertainty=args.uncertainty)
        results.append({
            "__mhctools_id": row["__mhctools_id"],
            "peptide": row["peptide"],
            # `inference.py` applies np.expm1 to the log-trained wt model
            # itself, so this is already hours. Do not transform again.
            "hours": out["score"],
            "emb_tag": out.get("emb_tag", ""),
            "uncertainty": (
                "" if out.get("uncertainty") is None
                else json.dumps(out["uncertainty"])),
            "uncertainty_type": out.get("uncertainty_type", ""),
        })

    fieldnames = [
        "__mhctools_id", "peptide", "hours",
        "emb_tag", "uncertainty", "uncertainty_type",
    ]
    with open(args.output, "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(results)


if __name__ == "__main__":
    main()
