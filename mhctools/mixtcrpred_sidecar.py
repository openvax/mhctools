# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0

"""Isolated runtime bridge to a user-licensed MixTCRpred checkout."""

from argparse import ArgumentParser
import os
from pathlib import Path
import pickle
import sys


def _parse_args():
    parser = ArgumentParser()
    parser.add_argument("--home", required=True)
    parser.add_argument("--model", required=True)
    parser.add_argument("--peptide", required=True)
    parser.add_argument("--checkpoint", required=True)
    parser.add_argument("--input", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--batch-size", type=int, required=True)
    parser.add_argument("--host", required=True)
    return parser.parse_args()


def _clean_derived(value):
    return str(value).replace("X", "") if value == value else ""


def main():
    args = _parse_args()
    home = Path(args.home).resolve()
    sys.path.insert(0, str(home))

    # Imports deliberately happen only inside this subprocess. The modules
    # below belong to the separately licensed upstream checkout.
    import pandas as pd
    import torch
    from torch.utils.data import DataLoader
    from src.dataloaders import db_transformer
    from src.models import TransformerPredictor_AB_cdr123
    from src.utils import compute_perc_rank

    frame = pd.read_csv(args.input, keep_default_na=False)
    dataset = db_transformer(
        frame,
        padding=[len(args.peptide), 20, 20, 10],
        istrain=False,
        host=args.host,
    )
    if len(dataset) != len(frame):
        raise RuntimeError(
            "MixTCRpred's input encoder removed one or more validated rows")

    loader = DataLoader(
        dataset,
        batch_size=args.batch_size,
        shuffle=False,
        num_workers=0,
    )
    # The wrapper sets TORCH_FORCE_NO_WEIGHTS_ONLY_LOAD for compatibility with
    # both the original Lightning 1.6 runtime and current PL. Managed files are
    # pinned/checksummed; callers must trust explicit or user-managed paths.
    model = TransformerPredictor_AB_cdr123.load_from_checkpoint(
        checkpoint_path=args.checkpoint,
        map_location=torch.device("cpu"),
    )
    model.eval()
    scores = []
    encoded_names = []
    with torch.inference_mode():
        for names, inputs, _ in loader:
            encoded_names.extend(names)
            scores.extend(model(inputs, mask=True)[:, 1].cpu().tolist())

    anchor_path = home / "pretrained_models" / "anchors_perc_rank.pickle"
    with anchor_path.open("rb") as input_file:
        anchors = pickle.load(input_file)
    ranks = compute_perc_rank(args.model, anchors, scores)

    if len(scores) != len(frame) or len(encoded_names) != len(frame):
        raise RuntimeError("MixTCRpred returned an unexpected number of rows")

    corrected = []
    for encoded in encoded_names:
        parts = encoded.split("_")
        if len(parts) != 7:
            raise RuntimeError("Unexpected MixTCRpred encoded TCR identifier")
        corrected.append(parts[1:])

    output = pd.DataFrame({
        "__mhctools_id": frame["__mhctools_id"].astype(int),
        "score": scores,
        "percentile_rank": ranks,
        "trav_corrected": [_clean_derived(value) for value in dataset.TRAV],
        "traj_corrected": [_clean_derived(value) for value in dataset.TRAJ],
        "trbv_corrected": [_clean_derived(value) for value in dataset.TRBV],
        "trbj_corrected": [_clean_derived(value) for value in dataset.TRBJ],
        "cdr1a_derived": [_clean_derived(row[1]) for row in corrected],
        "cdr2a_derived": [_clean_derived(row[2]) for row in corrected],
        "cdr1b_derived": [_clean_derived(row[4]) for row in corrected],
        "cdr2b_derived": [_clean_derived(row[5]) for row in corrected],
    })
    alpha_error = (
        (output["cdr1a_derived"] == "") |
        (output["cdr2a_derived"] == ""))
    beta_error = (
        (output["cdr1b_derived"] == "") |
        (output["cdr2b_derived"] == ""))
    output["warning"] = "-"
    output.loc[alpha_error, "warning"] = "Error TRAV"
    output.loc[beta_error, "warning"] = "Error TRBV"
    output.loc[alpha_error & beta_error, "warning"] = "Error TRAV-TRBV"
    output.to_csv(args.output, index=False)


if __name__ == "__main__":
    # Avoid user-site packages leaking into an explicitly selected sidecar env.
    os.environ.setdefault("PYTHONNOUSERSITE", "1")
    main()
