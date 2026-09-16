#!/usr/bin/env python3
"""Reproduce the osteosarc.com vaccine-sequence cleavage analysis.

The script deliberately keeps native model outputs separate.  A score from
one model is never averaged with, ranked against, or relabelled as a score
from another model.  Motif matches are recorded as recognition evidence,
not as cleavage probabilities.
"""

from __future__ import annotations

import argparse
import csv
from datetime import datetime
import hashlib
import importlib.metadata
import json
import math
from pathlib import Path
import platform
import subprocess
import tempfile
from typing import Any, Iterable

import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import numpy as np
import pandas as pd

from mhctools import (
    NetCleave_I,
    NetCleave_II,
    Pepsickle,
    cleavage_models,
    get_cleavage_model,
)
from mhctools.eramer_cleavage import ERAMERCleavage
from mhctools.netchop import NetChop


CANONICAL_AA = frozenset("ACDEFGHIKLMNPQRSTVWY")
SLP_VACCINES = frozenset(("JLF V1", "JLF V2", "JLF V3", "CeGaT"))
THRESHOLD = 0.5
NETCHOP_IMAGE = "i386/debian@sha256:75efd55b326373cf69989912388c0d50c5390638af7378d2fedc3aeb9d100e46"
EXTRACELLULAR_MOTIF_MODELS = [
    "ace-dipeptidyl",
    "mme-hydrophobic",
    "cpb2-basic",
    "cpn-basic",
    "app2-xp",
    "fap-dipeptidyl",
    "fap-endo-gp",
    "enpep-acidic",
    "anpep-ala",
]
INTRACELLULAR_ER_MOTIF_MODELS = [
    "app1-xp",
    "tpp2-tripeptidyl",
    "npepps-n-terminal",
    "dpp8-xp-xa",
    "dpp9-xp-xa",
    "prep-pro",
    "erap2-basic",
]


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def git_value(repo: Path, format_string: str) -> str:
    return subprocess.check_output(
        ["git", "-C", str(repo), "log", "-1", f"--format={format_string}"],
        text=True,
    ).strip()


def write_csv(path: Path, rows: Iterable[dict[str, Any]], columns: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def extract_inventory(
    variants_path: Path,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    variants = json.loads(variants_path.read_text(encoding="utf-8"))
    inventory: list[dict[str, Any]] = []
    missing: list[dict[str, Any]] = []

    for variant in variants:
        peptides = variant.get("vaccine_peptides") or []
        disclosed_by_vaccine: dict[str, list[str]] = {}
        for index, peptide in enumerate(peptides, start=1):
            vaccines = peptide.get("in_vaccines") or []
            if not vaccines:
                continue
            sequence = peptide["sequence"].replace(" ", "").upper()
            if not sequence or not set(sequence) <= CANONICAL_AA:
                raise ValueError(
                    f"Non-canonical sequence in {variant['id']} peptide {index}: {sequence!r}"
                )
            for vaccine in vaccines:
                disclosed_by_vaccine.setdefault(vaccine, []).append(sequence)

            if vaccines == ["mRNA"]:
                sequence_type = (
                    "mrna_minimal_epitope"
                    if peptide.get("is_mrna_minimal_epitope")
                    else "mrna_encoded_context"
                )
            elif any(vaccine in SLP_VACCINES for vaccine in vaccines):
                sequence_type = "synthetic_long_peptide"
            else:
                sequence_type = "other_disclosed_vaccine_sequence"

            inventory.append(
                {
                    "sequence_record_id": f"{variant['id']}:vaccine-peptide-{index}",
                    "variant_id": variant["id"],
                    "gene": variant.get("gene"),
                    "protein_change": variant.get("protein_change"),
                    "sequence": sequence,
                    "length": len(sequence),
                    "sequence_sha256": hashlib.sha256(
                        sequence.encode("ascii")
                    ).hexdigest(),
                    "sequence_type": sequence_type,
                    "vaccines": ";".join(vaccines),
                    "is_mrna_minimal_epitope": bool(
                        peptide.get("is_mrna_minimal_epitope")
                    ),
                    "minimal_epitope_offset": peptide.get("minimal_epitope_offset"),
                    "source_url": f"https://osteosarc.com/variant/{variant['id']}/",
                }
            )

        for vaccine, selected in (variant.get("vaccines") or {}).items():
            if selected and not disclosed_by_vaccine.get(vaccine):
                missing.append(
                    {
                        "variant_id": variant["id"],
                        "gene": variant.get("gene"),
                        "protein_change": variant.get("protein_change"),
                        "vaccine": vaccine,
                        "reason": "target listed but no sequence is disclosed for this vaccine",
                        "source_url": f"https://osteosarc.com/variant/{variant['id']}/",
                    }
                )

    inventory_df = pd.DataFrame(inventory).sort_values(
        ["sequence_type", "gene", "sequence_record_id"]
    )
    missing_df = pd.DataFrame(missing).sort_values(["vaccine", "gene", "variant_id"])

    conflicts: list[dict[str, Any]] = []
    for sequence, group in inventory_df.groupby("sequence", sort=True):
        variant_ids = sorted(group["variant_id"].unique())
        if len(variant_ids) <= 1:
            continue
        conflicts.append(
            {
                "sequence": sequence,
                "length": len(sequence),
                "sequence_sha256": hashlib.sha256(sequence.encode("ascii")).hexdigest(),
                "variant_count": len(variant_ids),
                "variant_ids": ";".join(variant_ids),
                "genes": ";".join(sorted(group["gene"].unique())),
                "vaccines": ";".join(
                    sorted(
                        {
                            vaccine
                            for value in group["vaccines"]
                            for vaccine in value.split(";")
                        }
                    )
                ),
                "warning": "identical disclosed sequence is assigned to distinct variants",
            }
        )
    conflicts_df = pd.DataFrame(conflicts)
    return inventory_df, missing_df, conflicts_df


def slp_records(inventory_df: pd.DataFrame) -> pd.DataFrame:
    return inventory_df.loc[
        inventory_df["sequence_type"] == "synthetic_long_peptide"
    ].copy()


def common_prediction_row(record: pd.Series, bond: int) -> dict[str, Any]:
    sequence = record["sequence"]
    return {
        "sequence_record_id": record["sequence_record_id"],
        "variant_id": record["variant_id"],
        "gene": record["gene"],
        "protein_change": record["protein_change"],
        "vaccines": record["vaccines"],
        "sequence": sequence,
        "length": len(sequence),
        "bond": bond,
        "left_residue": sequence[bond - 1],
        "right_residue": sequence[bond],
        "bond_label": f"{sequence[bond - 1]}{bond}|{sequence[bond]}{bond + 1}",
    }


def pepsickle_rows(
    records: pd.DataFrame,
) -> tuple[list[dict[str, Any]], list[dict[str, Any]]]:
    rows: list[dict[str, Any]] = []
    models: list[dict[str, Any]] = []
    unique_sequences = list(dict.fromkeys(records["sequence"]))
    for human_only, model_name in (
        (False, "pepsickle-in-vivo-all-mammal"),
        (True, "pepsickle-in-vivo-human-only"),
    ):
        predictor = Pepsickle(human_only=human_only, isolate_subprocess=True)
        predictions = predictor.cleavage_probs_many(unique_sequences)
        models.append(
            {
                "model": model_name,
                "family": "Pepsickle",
                "biological_context": "constitutive proteasome; in-vivo epitope-trained model",
                "evidence_type": "quantitative_model",
                "input_definition": "one score after each residue; terminal score excluded from internal-bond analysis",
                "score_units": "native 0-1 cleavage score",
                "display_threshold": THRESHOLD,
                "threshold_basis": "common display threshold; not calibrated on these SLPs",
                "included": True,
                "exclusion_reason": "",
                "reference": "https://doi.org/10.1093/bioinformatics/btab628",
            }
        )
        for _, record in records.iterrows():
            scores = predictions[record["sequence"]]
            if len(scores) != len(record["sequence"]):
                raise RuntimeError(f"{model_name} returned the wrong score count")
            for bond in range(1, len(record["sequence"])):
                score = float(scores[bond - 1])
                row = common_prediction_row(record, bond)
                row.update(
                    {
                        "model": model_name,
                        "family": "Pepsickle",
                        "biological_context": "proteasome/cytosol",
                        "score": score,
                        "score_units": "native 0-1 cleavage score",
                        "display_threshold": THRESHOLD,
                        "above_display_threshold": score >= THRESHOLD,
                        "assessable": True,
                        "unsupported_reason": "",
                    }
                )
                rows.append(row)
    return rows, models


def run_netchop_docker(
    sequences: list[str], netchop_dir: Path, model_variant: int
) -> list[list[float]]:
    if model_variant not in (0, 1):
        raise ValueError("NetChop model variant must be 0 (Cterm) or 1 (20S)")
    with tempfile.TemporaryDirectory(prefix="osteosarc_netchop_") as tmp:
        tmp_path = Path(tmp)
        fasta_path = tmp_path / "sequences.fasta"
        with fasta_path.open("w", encoding="ascii") as handle:
            for index, sequence in enumerate(sequences):
                handle.write(f">slp_{index}\n{sequence}\n")
        command = [
            "docker",
            "run",
            "--rm",
            "--platform",
            "linux/386",
            "-e",
            "NETCHOP=/netchop",
            "-e",
            "TMPDIR=/tmp",
            "-v",
            f"{netchop_dir.resolve()}:/netchop:ro",
            "-v",
            f"{tmp_path.resolve()}:/work:ro",
            NETCHOP_IMAGE,
            "/netchop/bin/netChop",
            "-v",
            str(model_variant),
            "/work/sequences.fasta",
        ]
        completed = subprocess.run(command, capture_output=True, check=True)
    parsed = NetChop.parse_netchop(completed.stdout)
    if len(parsed) != len(sequences):
        raise RuntimeError(
            f"NetChop returned {len(parsed)} sequences for {len(sequences)} inputs"
        )
    return parsed


def netchop_rows(
    records: pd.DataFrame, netchop_dir: Path
) -> tuple[list[dict[str, Any]], list[dict[str, Any]]]:
    rows: list[dict[str, Any]] = []
    models: list[dict[str, Any]] = []
    unique_sequences = list(dict.fromkeys(records["sequence"]))
    for variant, suffix, context in (
        (0, "cterm-3.0", "C-terminal epitope-trained proteasome model"),
        (1, "20s-3.0", "in-vitro 20S proteasome model"),
    ):
        model_name = f"netchop-3.1-{suffix}"
        predictions = dict(
            zip(
                unique_sequences,
                run_netchop_docker(unique_sequences, netchop_dir, variant),
            )
        )
        models.append(
            {
                "model": model_name,
                "family": "NetChop",
                "biological_context": context,
                "evidence_type": "quantitative_model",
                "input_definition": "one score after each residue; terminal score excluded from internal-bond analysis",
                "score_units": "native 0-1 cleavage score",
                "display_threshold": THRESHOLD,
                "threshold_basis": "NetChop default threshold",
                "included": True,
                "exclusion_reason": "",
                "reference": "https://doi.org/10.1007/s00251-005-0781-7",
            }
        )
        for _, record in records.iterrows():
            scores = predictions[record["sequence"]]
            if len(scores) != len(record["sequence"]):
                raise RuntimeError(f"{model_name} returned the wrong score count")
            for bond in range(1, len(record["sequence"])):
                score = float(scores[bond - 1])
                row = common_prediction_row(record, bond)
                row.update(
                    {
                        "model": model_name,
                        "family": "NetChop",
                        "biological_context": "proteasome/cytosol",
                        "score": score,
                        "score_units": "native 0-1 cleavage score",
                        "display_threshold": THRESHOLD,
                        "above_display_threshold": score >= THRESHOLD,
                        "assessable": True,
                        "unsupported_reason": "",
                    }
                )
                rows.append(row)
    return rows, models


def netcleave_rows(
    records: pd.DataFrame, netcleave_dir: Path
) -> tuple[list[dict[str, Any]], list[dict[str, Any]]]:
    rows: list[dict[str, Any]] = []
    models: list[dict[str, Any]] = []
    for mhc_class, predictor, epitope_length, model_name, context in (
        (
            "I",
            NetCleave_I(netcleave_path=str(netcleave_dir)),
            8,
            "netcleave-i-hla",
            "MHC-I C-terminal processing; proteasome-associated",
        ),
        (
            "II",
            NetCleave_II(netcleave_path=str(netcleave_dir)),
            13,
            "netcleave-ii-hla",
            "MHC-II C-terminal processing; endolysosomal-associated",
        ),
    ):
        requests: list[tuple[pd.Series, int, str, str]] = []
        peptides: list[str] = []
        c_flanks: list[str] = []
        for _, record in records.iterrows():
            sequence = record["sequence"]
            for bond in range(epitope_length, len(sequence) - 2):
                peptide = sequence[bond - epitope_length : bond]
                c_flank = sequence[bond : bond + 3]
                requests.append((record, bond, peptide, c_flank))
                peptides.append(peptide)
                c_flanks.append(c_flank)
        predictions = predictor.predict(peptides, c_flanks=c_flanks)
        if len(predictions) != len(requests):
            raise RuntimeError(f"{model_name} returned the wrong result count")
        for (record, bond, _peptide, _c_flank), result in zip(requests, predictions):
            if not result.preds:
                raise RuntimeError(f"{model_name} omitted an in-domain request")
            score = float(result.preds[0].score)
            row = common_prediction_row(record, bond)
            row.update(
                {
                    "model": model_name,
                    "family": "NetCleave",
                    "biological_context": (
                        "proteasome/cytosol" if mhc_class == "I" else "endolysosome"
                    ),
                    "score": score,
                    "score_units": "native neural-network output (0-1)",
                    "display_threshold": THRESHOLD,
                    "above_display_threshold": score >= THRESHOLD,
                    "assessable": True,
                    "unsupported_reason": "",
                }
            )
            rows.append(row)
        models.append(
            {
                "model": model_name,
                "family": "NetCleave",
                "biological_context": context,
                "evidence_type": "quantitative_model",
                "input_definition": (
                    f"{epitope_length}-residue peptide ending at candidate bond "
                    "plus three downstream residues"
                ),
                "score_units": "native neural-network output (0-1)",
                "display_threshold": THRESHOLD,
                "threshold_basis": (
                    "display threshold only; class-II model has much weaker published discrimination"
                ),
                "included": True,
                "exclusion_reason": "",
                "reference": "https://doi.org/10.1038/s41598-021-92632-y",
            }
        )
    return rows, models


def peptidase_rows(
    records: pd.DataFrame, eramer_dir: Path
) -> tuple[list[dict[str, Any]], list[dict[str, Any]], list[dict[str, Any]]]:
    quantitative: list[dict[str, Any]] = []
    motifs: list[dict[str, Any]] = []
    model_rows: list[dict[str, Any]] = []

    built_in_models = list(cleavage_models(include_optional=False))
    selected = [
        model
        for model in built_in_models
        if model.evidence in ("motif_rule", "quantitative_model")
    ]
    for metadata in built_in_models:
        included = metadata in selected
        model_rows.append(
            {
                "model": metadata.name,
                "family": "mhctools peptidase panel",
                "biological_context": ";".join(metadata.compartments),
                "evidence_type": metadata.evidence,
                "input_definition": "intact SLP with free N and C termini",
                "score_units": metadata.score_units or "matched/not_matched",
                "display_threshold": "",
                "threshold_basis": "",
                "included": included,
                "exclusion_reason": (
                    ""
                    if included
                    else "exact-sequence source reference; not extrapolated"
                ),
                "reference": ";".join(metadata.references),
            }
        )

    for metadata in selected:
        predictor = get_cleavage_model(
            metadata.name,
            enzyme_state="active" if metadata.name == "cpb2-basic" else None,
        )
        for _, record in records.iterrows():
            result = predictor.predict(record["sequence"])
            if metadata.evidence == "quantitative_model":
                for site in result.sites:
                    row = common_prediction_row(record, site.bond)
                    row.update(
                        {
                            "model": metadata.name,
                            "family": "mhctools peptidase panel",
                            "biological_context": ";".join(metadata.compartments),
                            "score": float(site.score),
                            "score_units": metadata.score_units,
                            "display_threshold": math.nan,
                            "above_display_threshold": None,
                            "assessable": True,
                            "unsupported_reason": "",
                        }
                    )
                    quantitative.append(row)
                if result.unsupported_reason:
                    quantitative.append(
                        {
                            **common_prediction_row(record, 1),
                            "bond": math.nan,
                            "left_residue": "",
                            "right_residue": "",
                            "bond_label": "",
                            "model": metadata.name,
                            "family": "mhctools peptidase panel",
                            "biological_context": ";".join(metadata.compartments),
                            "score": math.nan,
                            "score_units": metadata.score_units,
                            "display_threshold": math.nan,
                            "above_display_threshold": None,
                            "assessable": False,
                            "unsupported_reason": result.unsupported_reason,
                        }
                    )
            else:
                if result.unsupported_reason:
                    motifs.append(
                        {
                            "sequence_record_id": record["sequence_record_id"],
                            "variant_id": record["variant_id"],
                            "gene": record["gene"],
                            "protein_change": record["protein_change"],
                            "vaccines": record["vaccines"],
                            "sequence": record["sequence"],
                            "length": len(record["sequence"]),
                            "model": metadata.name,
                            "enzyme": metadata.enzyme,
                            "compartments": ";".join(metadata.compartments),
                            "motif_strictness": metadata.motif_strictness,
                            "status": "unsupported",
                            "bond": math.nan,
                            "bond_label": "",
                            "reason": "",
                            "unsupported_reason": result.unsupported_reason,
                        }
                    )
                for site in result.sites:
                    bond = common_prediction_row(record, site.bond)
                    motifs.append(
                        {
                            **bond,
                            "model": metadata.name,
                            "enzyme": metadata.enzyme,
                            "compartments": ";".join(metadata.compartments),
                            "motif_strictness": metadata.motif_strictness,
                            "status": site.status,
                            "reason": site.reason,
                            "unsupported_reason": "",
                        }
                    )

    eramer = ERAMERCleavage(eramer_home=str(eramer_dir))
    metadata = eramer.model
    model_rows.append(
        {
            "model": metadata.name,
            "family": "ERAMER",
            "biological_context": "ERAP1 trimming in the endoplasmic reticulum",
            "evidence_type": metadata.evidence,
            "input_definition": "intact 9-16-residue SLP; first N-terminal trimming step",
            "score_units": metadata.score_units,
            "display_threshold": "",
            "threshold_basis": "native score has no validated binary threshold",
            "included": True,
            "exclusion_reason": "",
            "reference": ";".join(metadata.references),
        }
    )
    for _, record in records.iterrows():
        result = eramer.predict(record["sequence"])
        if result.sites:
            site = result.sites[0]
            row = common_prediction_row(record, site.bond)
            row.update(
                {
                    "model": metadata.name,
                    "family": "ERAMER",
                    "biological_context": "ERAP1/ER",
                    "score": float(site.score),
                    "score_units": metadata.score_units,
                    "display_threshold": math.nan,
                    "above_display_threshold": None,
                    "assessable": True,
                    "unsupported_reason": "",
                }
            )
            quantitative.append(row)
        else:
            quantitative.append(
                {
                    **common_prediction_row(record, 1),
                    "bond": math.nan,
                    "left_residue": "",
                    "right_residue": "",
                    "bond_label": "",
                    "model": metadata.name,
                    "family": "ERAMER",
                    "biological_context": "ERAP1/ER",
                    "score": math.nan,
                    "score_units": metadata.score_units,
                    "display_threshold": math.nan,
                    "above_display_threshold": None,
                    "assessable": False,
                    "unsupported_reason": result.unsupported_reason,
                }
            )
    return quantitative, motifs, model_rows


def summarize_models(
    records: pd.DataFrame,
    quantitative_df: pd.DataFrame,
    motifs_df: pd.DataFrame,
    model_catalog_df: pd.DataFrame,
) -> pd.DataFrame:
    rows: list[dict[str, Any]] = []
    included_models = model_catalog_df.loc[model_catalog_df["included"], "model"]
    for _, record in records.iterrows():
        for model in included_models:
            q = quantitative_df.loc[
                (quantitative_df["sequence_record_id"] == record["sequence_record_id"])
                & (quantitative_df["model"] == model)
                & (quantitative_df["assessable"] == True)  # noqa: E712
            ]
            m = motifs_df.loc[
                (motifs_df["sequence_record_id"] == record["sequence_record_id"])
                & (motifs_df["model"] == model)
            ]
            catalog = model_catalog_df.loc[model_catalog_df["model"] == model].iloc[0]
            if not q.empty:
                thresholded = q["display_threshold"].notna().any()
                hit_count = (
                    int(q["above_display_threshold"].eq(True).sum())
                    if thresholded
                    else math.nan
                )
                assessed = len(q)
                rows.append(
                    {
                        "sequence_record_id": record["sequence_record_id"],
                        "variant_id": record["variant_id"],
                        "gene": record["gene"],
                        "protein_change": record["protein_change"],
                        "vaccines": record["vaccines"],
                        "sequence": record["sequence"],
                        "length": record["length"],
                        "model": model,
                        "family": catalog["family"],
                        "evidence_type": catalog["evidence_type"],
                        "assessed_sites": assessed,
                        "candidate_sites": hit_count,
                        "candidate_fraction": hit_count / assessed
                        if thresholded
                        else math.nan,
                        "matched_sites": math.nan,
                        "mean_score": q["score"].mean(),
                        "max_score": q["score"].max(),
                        "unsupported_reason": "",
                    }
                )
            elif not m.empty:
                assessed = int((m["status"] != "unsupported").sum())
                matched = int((m["status"] == "matched").sum())
                unsupported = ";".join(
                    sorted(
                        set(m.loc[m["status"] == "unsupported", "unsupported_reason"])
                    )
                )
                rows.append(
                    {
                        "sequence_record_id": record["sequence_record_id"],
                        "variant_id": record["variant_id"],
                        "gene": record["gene"],
                        "protein_change": record["protein_change"],
                        "vaccines": record["vaccines"],
                        "sequence": record["sequence"],
                        "length": record["length"],
                        "model": model,
                        "family": catalog["family"],
                        "evidence_type": catalog["evidence_type"],
                        "assessed_sites": assessed,
                        "candidate_sites": math.nan,
                        "candidate_fraction": math.nan,
                        "matched_sites": matched if assessed else math.nan,
                        "mean_score": math.nan,
                        "max_score": math.nan,
                        "unsupported_reason": unsupported,
                    }
                )
            else:
                unsupported = quantitative_df.loc[
                    (
                        quantitative_df["sequence_record_id"]
                        == record["sequence_record_id"]
                    )
                    & (quantitative_df["model"] == model),
                    "unsupported_reason",
                ]
                rows.append(
                    {
                        "sequence_record_id": record["sequence_record_id"],
                        "variant_id": record["variant_id"],
                        "gene": record["gene"],
                        "protein_change": record["protein_change"],
                        "vaccines": record["vaccines"],
                        "sequence": record["sequence"],
                        "length": record["length"],
                        "model": model,
                        "family": catalog["family"],
                        "evidence_type": catalog["evidence_type"],
                        "assessed_sites": 0,
                        "candidate_sites": math.nan,
                        "candidate_fraction": math.nan,
                        "matched_sites": math.nan,
                        "mean_score": math.nan,
                        "max_score": math.nan,
                        "unsupported_reason": ";".join(
                            sorted(set(unsupported.dropna()))
                        ),
                    }
                )
    return pd.DataFrame(rows)


def predictor_matrix(records: pd.DataFrame, summary_df: pd.DataFrame) -> pd.DataFrame:
    """Return one row per SLP with a self-describing column per model."""
    identity_columns = [
        "sequence_record_id",
        "variant_id",
        "gene",
        "protein_change",
        "vaccines",
        "sequence",
        "length",
    ]
    matrix = records[identity_columns].copy().set_index("sequence_record_id")
    for model in summary_df["model"].drop_duplicates():
        model_rows = summary_df.loc[summary_df["model"] == model].set_index(
            "sequence_record_id"
        )
        evidence_type = model_rows["evidence_type"].iloc[0]
        if evidence_type == "motif_rule":
            metric = "matched_sites"
        elif model_rows["candidate_fraction"].notna().any():
            metric = "fraction_ge_0_5"
            model_rows = model_rows.rename(
                columns={"candidate_fraction": "fraction_ge_0_5"}
            )
        else:
            metric = "native_score"
            model_rows = model_rows.rename(columns={"max_score": "native_score"})
        matrix[f"{model}__{metric}"] = model_rows[metric]
    return matrix.reset_index()


def validate_analysis(
    inventory_df: pd.DataFrame,
    records: pd.DataFrame,
    quantitative_df: pd.DataFrame,
    motifs_df: pd.DataFrame,
    model_catalog_df: pd.DataFrame,
    summary_df: pd.DataFrame,
) -> None:
    """Fail before publishing internally inconsistent analysis artifacts."""
    if inventory_df.empty or records.empty:
        raise RuntimeError("The source inventory or SLP subset is empty")
    if not records["sequence_record_id"].is_unique:
        raise RuntimeError("SLP sequence record IDs are not unique")
    if not inventory_df["sequence"].map(lambda value: set(value) <= CANONICAL_AA).all():
        raise RuntimeError("The inventory contains a non-canonical sequence")

    included_models = model_catalog_df.loc[model_catalog_df["included"], "model"]
    expected_pairs = len(records) * len(included_models)
    if (
        len(summary_df) != expected_pairs
        or summary_df.duplicated(["sequence_record_id", "model"]).any()
    ):
        raise RuntimeError("The model summary is not a complete SLP-by-model grid")

    assessable_quantitative = quantitative_df.loc[quantitative_df["assessable"]]
    invalid_quantitative_bond = (assessable_quantitative["bond"] < 1) | (
        assessable_quantitative["bond"] >= assessable_quantitative["length"]
    )
    if invalid_quantitative_bond.any():
        raise RuntimeError("An assessable quantitative result has an invalid bond")
    if not np.isfinite(assessable_quantitative["score"]).all():
        raise RuntimeError("An assessable quantitative result has a non-finite score")
    thresholded = assessable_quantitative.loc[
        assessable_quantitative["display_threshold"].notna()
    ]
    if not thresholded["score"].between(0, 1).all():
        raise RuntimeError("A thresholded score falls outside its documented 0-1 scale")

    assessed_motifs = motifs_df.loc[motifs_df["status"] != "unsupported"]
    invalid_motif_bond = (assessed_motifs["bond"] < 1) | (
        assessed_motifs["bond"] >= assessed_motifs["length"]
    )
    if invalid_motif_bond.any():
        raise RuntimeError("An assessed motif result has an invalid bond")
    if not motifs_df["status"].isin(("matched", "not_matched", "unsupported")).all():
        raise RuntimeError("A motif result has an unknown status")

    candidate_fractions = summary_df["candidate_fraction"].dropna()
    if not candidate_fractions.between(0, 1).all():
        raise RuntimeError("A summary candidate fraction falls outside 0-1")


def row_labels(records: pd.DataFrame) -> dict[str, str]:
    labels: dict[str, str] = {}
    seen: dict[str, int] = {}
    for _, record in records.iterrows():
        vaccines = record["vaccines"].replace(";", ", ")
        base = (
            f"{record['gene']} {record['protein_change'] or ''} · "
            f"{vaccines} · {record['length']} aa"
        )
        seen[base] = seen.get(base, 0) + 1
        labels[record["sequence_record_id"]] = (
            base if seen[base] == 1 else f"{base} [{seen[base]}]"
        )
    return labels


def save_figure(
    fig: plt.Figure,
    stem: Path,
    pdf_pages: PdfPages,
    figure_number: int,
    figure_count: int,
) -> None:
    fig.text(
        0.995,
        0.005,
        f"Osteosarc vaccine cleavage · Figure {figure_number} of {figure_count}",
        ha="right",
        va="bottom",
        fontsize=7,
        color="#555555",
    )
    fig.savefig(stem.with_suffix(".png"), dpi=180, bbox_inches="tight")
    pdf_pages.savefig(fig, bbox_inches="tight")
    plt.close(fig)


def plot_quantitative_heatmap(
    records: pd.DataFrame,
    summary_df: pd.DataFrame,
    figures_dir: Path,
    pdf_pages: PdfPages,
) -> None:
    models = [
        "pepsickle-in-vivo-all-mammal",
        "pepsickle-in-vivo-human-only",
        "netchop-3.1-cterm-3.0",
        "netchop-3.1-20s-3.0",
        "netcleave-i-hla",
        "netcleave-ii-hla",
    ]
    labels = row_labels(records)
    matrix = summary_df.loc[summary_df["model"].isin(models)].pivot(
        index="sequence_record_id", columns="model", values="candidate_fraction"
    )
    matrix = matrix.reindex(index=records["sequence_record_id"], columns=models)
    display_labels = [
        "Pepsickle\nall-mammal",
        "Pepsickle\nhuman-only",
        "NetChop\nCterm 3.0",
        "NetChop\n20S 3.0",
        "NetCleave I\npan-HLA",
        "NetCleave II\npan-HLA",
    ]
    fig, ax = plt.subplots(figsize=(16, 10.5))
    masked = np.ma.masked_invalid(matrix.to_numpy(dtype=float))
    cmap = plt.get_cmap("magma").copy()
    cmap.set_bad("#d9d9d9")
    image = ax.imshow(masked, aspect="auto", vmin=0, vmax=1, cmap=cmap)
    ax.set_xticks(range(len(models)), labels=display_labels)
    ax.tick_params(axis="x", labelsize=9, pad=6)
    ax.set_yticks(
        range(len(matrix)), labels=[labels[index] for index in matrix.index], fontsize=8
    )
    ax.set_title(
        "Fraction of assessable internal bonds at or above 0.5",
        fontsize=14,
        pad=34,
    )
    ax.text(
        2,
        -2.2,
        "Proteasome / MHC-I processing models",
        ha="center",
        va="center",
        fontsize=9,
        fontweight="bold",
    )
    ax.text(
        5,
        -2.2,
        "MHC-II processing model",
        ha="center",
        va="center",
        fontsize=9,
        fontweight="bold",
    )
    ax.axvline(4.5, color="white", linewidth=4)
    ax.axvline(4.5, color="#444444", linewidth=0.8)
    ax.set_xlabel("Model (native scores; fractions are not cross-model probabilities)")
    ax.set_ylabel("Disclosed synthetic long-peptide record")
    for y in range(masked.shape[0]):
        for x in range(masked.shape[1]):
            if masked.mask[y, x]:
                ax.text(
                    x,
                    y,
                    "NA",
                    ha="center",
                    va="center",
                    fontsize=6,
                    color="#555555",
                )
            else:
                ax.text(
                    x,
                    y,
                    f"{masked[y, x]:.2f}",
                    ha="center",
                    va="center",
                    fontsize=6,
                    color="white" if masked[y, x] > 0.48 else "black",
                )
    colorbar = fig.colorbar(image, ax=ax, label="within-model fraction")
    colorbar.ax.text(
        0.5,
        -0.035,
        "Gray = not assessable",
        transform=colorbar.ax.transAxes,
        ha="center",
        va="top",
        fontsize=7,
    )
    save_figure(
        fig,
        figures_dir / "slp_quantitative_hit_fraction",
        pdf_pages,
        figure_number=2,
        figure_count=4,
    )


def plot_motif_heatmap(
    records: pd.DataFrame,
    summary_df: pd.DataFrame,
    figures_dir: Path,
    pdf_pages: PdfPages,
    models: list[str],
    title: str,
    stem: str,
    figure_number: int,
) -> None:
    labels = row_labels(records)
    matrix = summary_df.loc[summary_df["model"].isin(models)].pivot(
        index="sequence_record_id", columns="model", values="matched_sites"
    )
    matrix = matrix.reindex(index=records["sequence_record_id"], columns=models)
    fig, ax = plt.subplots(figsize=(16, 10.5))
    masked = np.ma.masked_invalid(matrix.to_numpy(dtype=float))
    max_count = max(1.0, float(np.nanmax(matrix.to_numpy(dtype=float))))
    cmap = plt.get_cmap("viridis").copy()
    cmap.set_bad("#d9d9d9")
    image = ax.imshow(masked, aspect="auto", vmin=0, vmax=max_count, cmap=cmap)
    ax.set_xticks(
        range(len(models)), labels=models, rotation=45, ha="right", fontsize=8
    )
    ax.set_yticks(
        range(len(matrix)), labels=[labels[index] for index in matrix.index], fontsize=8
    )
    ax.set_title(title, fontsize=14, pad=28)
    ax.text(
        0.5,
        1.01,
        "Compartment grouping is for readability; see model_catalog.csv for full context.",
        transform=ax.transAxes,
        ha="center",
        va="bottom",
        fontsize=7,
        color="#555555",
    )
    ax.set_xlabel("Curated recognition rule (a match is not a cleavage probability)")
    ax.set_ylabel("Disclosed synthetic long-peptide record")
    for y in range(masked.shape[0]):
        for x in range(masked.shape[1]):
            if masked.mask[y, x]:
                ax.text(
                    x,
                    y,
                    "NA",
                    ha="center",
                    va="center",
                    fontsize=6,
                    color="#555555",
                )
            elif masked[y, x] > 0:
                normalized = masked[y, x] / max_count
                ax.text(
                    x,
                    y,
                    str(int(masked[y, x])),
                    ha="center",
                    va="center",
                    fontsize=7,
                    color="white" if normalized < 0.5 else "black",
                )
    colorbar = fig.colorbar(image, ax=ax, label="matched sites")
    colorbar.ax.text(
        0.5,
        -0.035,
        "Gray = unsupported",
        transform=colorbar.ax.transAxes,
        ha="center",
        va="top",
        fontsize=7,
    )
    save_figure(
        fig,
        figures_dir / stem,
        pdf_pages,
        figure_number=figure_number,
        figure_count=4,
    )


def plot_score_distributions(
    quantitative_df: pd.DataFrame, figures_dir: Path, pdf_pages: PdfPages
) -> None:
    models = [
        "pepsickle-in-vivo-all-mammal",
        "netchop-3.1-cterm-3.0",
        "netcleave-i-hla",
        "pepsickle-in-vivo-human-only",
        "netchop-3.1-20s-3.0",
        "netcleave-ii-hla",
    ]
    fig, axes = plt.subplots(2, 3, figsize=(15, 9), sharex=True, sharey=True)
    for ax, model in zip(axes.flat, models):
        values = quantitative_df.loc[
            (quantitative_df["model"] == model) & quantitative_df["assessable"],
            "score",
        ].dropna()
        weights = np.full(len(values), 100.0 / len(values))
        ax.hist(
            values,
            bins=np.linspace(0, 1, 31),
            weights=weights,
            color="#3366a6",
            alpha=0.85,
        )
        ax.axvline(
            THRESHOLD,
            color="#b33a3a",
            linestyle="--",
            linewidth=1,
            label="0.5 display threshold",
        )
        ax.set_title(f"{model}\nn={len(values):,} assessable bonds", fontsize=10)
        ax.set_xlabel("native score")
        ax.set_ylabel("assessable bonds per bin (%)")
        ax.tick_params(axis="x", labelbottom=True)
    fig.suptitle(
        "Native output distributions (panels are not calibrated to one another)",
        fontsize=14,
    )
    handles, labels = axes.flat[0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="upper center", bbox_to_anchor=(0.5, 0.945))
    fig.tight_layout(rect=(0, 0.02, 1, 0.91))
    save_figure(
        fig,
        figures_dir / "quantitative_score_distributions",
        pdf_pages,
        figure_number=1,
        figure_count=4,
    )


def correlations(quantitative_df: pd.DataFrame) -> pd.DataFrame:
    thresholded = quantitative_df.loc[
        quantitative_df["display_threshold"].notna() & quantitative_df["assessable"]
    ].copy()
    wide = thresholded.pivot_table(
        index=["sequence_record_id", "bond"], columns="model", values="score"
    )
    return wide.corr(method="spearman", min_periods=10)


def model_hashes(
    netchop_dir: Path, netcleave_dir: Path, eramer_dir: Path
) -> dict[str, Any]:
    import pepsickle

    pepsickle_root = Path(pepsickle.__file__).parent
    pepsickle_files = [
        pepsickle_root / "model.joblib",
        pepsickle_root / "trained_model_dict.pickle",
        pepsickle_root / "in-vitro_mammal" / "model.joblib",
        pepsickle_root / "model_functions.py",
        pepsickle_root / "sequence_featurization_tools.py",
    ]
    netchop_files = [netchop_dir / "bin" / "netChop"] + sorted(
        (netchop_dir / "data").rglob("*")
    )
    netcleave_files = [
        netcleave_dir
        / "data/models/I_mass-spectrometry_HLA/I_mass-spectrometry_HLA_model.h5",
        netcleave_dir
        / "data/models/II_mass-spectrometry_HLA/II_mass-spectrometry_HLA_model.h5",
    ]
    return {
        "pepsickle": {
            "package_version": importlib.metadata.version("pepsickle"),
            "files": {
                str(path.relative_to(pepsickle_root)): sha256_file(path)
                for path in pepsickle_files
            },
        },
        "netchop": {
            "version": "3.1",
            "docker_runtime_image": NETCHOP_IMAGE,
            "files": {
                str(path.relative_to(netchop_dir)): sha256_file(path)
                for path in netchop_files
                if path.is_file()
            },
        },
        "netcleave": {
            "repository": "https://github.com/BSC-CNS-EAPM/NetCleave",
            "commit": git_value(netcleave_dir, "%H"),
            "files": {
                str(path.relative_to(netcleave_dir)): sha256_file(path)
                for path in netcleave_files
            },
        },
        "eramer": {
            "repository": "https://github.com/aalokaily/ERAMER",
            "commit": git_value(eramer_dir, "%H"),
            "files": {"PWM.xlsx": sha256_file(eramer_dir / "PWM.xlsx")},
        },
    }


def write_report(
    path: Path,
    inventory_df: pd.DataFrame,
    records: pd.DataFrame,
    missing_df: pd.DataFrame,
    conflicts_df: pd.DataFrame,
    summary_df: pd.DataFrame,
    source_commit: str,
    source_date: str,
) -> None:
    thresholded = summary_df.loc[summary_df["candidate_fraction"].notna()]
    by_model = (
        thresholded.groupby("model")
        .agg(
            records_assessed=("sequence_record_id", "count"),
            median_candidate_fraction=("candidate_fraction", "median"),
            min_candidate_fraction=("candidate_fraction", "min"),
            max_candidate_fraction=("candidate_fraction", "max"),
        )
        .reset_index()
    )
    motif = summary_df.loc[summary_df["evidence_type"] == "motif_rule"]
    motif_hits = (
        motif.groupby("model")["matched_sites"]
        .sum(min_count=1)
        .sort_values(ascending=False)
    )
    unthresholded = summary_df.loc[
        (summary_df["evidence_type"] == "quantitative_model")
        & summary_df["candidate_fraction"].isna()
        & summary_df["mean_score"].notna()
    ]
    unthresholded_by_model = (
        unthresholded.groupby("model")
        .agg(
            records_assessed=("sequence_record_id", "count"),
            min_score=("mean_score", "min"),
            max_score=("mean_score", "max"),
        )
        .reset_index()
    )
    lines = [
        "# Osteosarc vaccine cleavage analysis",
        "",
        f"Source snapshot: osteosarc.com repository `{source_commit}` ({source_date}).",
        "All inference was local; no peptide sequence was uploaded.",
        "",
        "## Coverage",
        "",
        f"- {len(inventory_df)} disclosed sequence records ({inventory_df['sequence'].nunique()} unique sequences).",
        f"- {len(records)} disclosed synthetic-long-peptide records ({records['sequence'].nunique()} unique sequences).",
        f"- {len(missing_df)} vaccine-target assignments have no vaccine-specific sequence disclosed on the site.",
        f"- {len(conflicts_df)} identical-sequence/across-variant provenance conflict was detected.",
        "",
        "The inventory distinguishes mRNA encoded contexts, displayed mRNA minimal epitopes, and SLPs. "
        "Only SLP records are included in the cleavage tables and figures.",
        "The compact [SLP-by-predictor matrix](tables/slp_predictor_matrix.csv) and exact "
        "[bond-level scores](tables/slp_quantitative_bond_scores.csv) are provided separately.",
        "",
        "## Quantitative model output distribution",
        "",
        "The 0.5 cutoff is used only as a within-model display threshold. It is not a calibrated "
        "probability of SLP degradation, and fractions must not be compared as if the models shared a scale.",
        "",
        "| Model | SLP records assessed | Median fraction ≥0.5 | Range |",
        "| --- | ---: | ---: | ---: |",
    ]
    for _, row in by_model.iterrows():
        lines.append(
            f"| {row['model']} | {int(row['records_assessed'])} | "
            f"{row['median_candidate_fraction']:.3f} | "
            f"{row['min_candidate_fraction']:.3f}–{row['max_candidate_fraction']:.3f} |"
        )
    lines.extend(
        [
            "",
            "![Within-model candidate-site fractions](figures/slp_quantitative_hit_fraction.png)",
            "",
            "![Native quantitative score distributions](figures/quantitative_score_distributions.png)",
            "",
            "NetCleave-I uses an 8-residue peptide ending at each candidate bond plus three "
            "downstream residues. NetCleave-II uses a 13-residue ending peptide plus the same "
            "three-residue downstream context; bonds lacking that context are not assessed.",
            "",
            "## Quantitative peptidase outputs without a binary threshold",
            "",
            "These native scores have no validated common cutoff and are therefore not included in "
            "the ≥0.5 figure. Their scales are model-specific.",
            "",
            "| Model | SLP records assessed | Native-score range |",
            "| --- | ---: | ---: |",
        ]
    )
    for _, row in unthresholded_by_model.iterrows():
        lines.append(
            f"| {row['model']} | {int(row['records_assessed'])} | "
            f"{row['min_score']:.4f}–{row['max_score']:.4f} |"
        )
    lines.extend(
        [
            "",
            "## Motif matches",
            "",
            "Counts below are matched recognition sites across all disclosed SLP records. Required, "
            "preferred, and permissive rules have different meanings; see `model_catalog.csv` and "
            "the project cleavage guide before interpreting a match.",
            "",
            "| Motif model | Total matched sites |",
            "| --- | ---: |",
        ]
    )
    for model, count in motif_hits.items():
        lines.append(f"| {model} | {int(count)} |")
    lines.extend(
        [
            "",
            "![Extracellular and plasma peptidase motif matches](figures/slp_peptidase_motifs_extracellular.png)",
            "",
            "![Intracellular and ER peptidase motif matches](figures/slp_peptidase_motifs_intracellular_er.png)",
            "",
            "## Interpretation boundary",
            "",
            "- Proteasome results are conditional on cytosolic access. An injected SLP ordinarily begins "
            "outside the cytosol; these scores do not predict uptake or cross-presentation.",
            "- NetCleave-II is an MHC-II C-terminal-processing model, not a named-cathepsin model. Its "
            "published class-II discrimination is much weaker than its class-I result.",
            "- Peptidase motif matches describe partial recognition rules on an intact peptide with free "
            "termini. They do not model abundance, activation, competition, ordered digestion, kinetics, "
            "formulation, structure, or MHC protection.",
            "- ERAMER is evaluated only where the intact SLP is inside its documented 9–16-residue domain. "
            "It is an ERAP1 trimming score, not a whole-SLP degradation score.",
            "- Missing sequences and the cross-variant duplicate are left visible. No sequence was guessed, "
            "corrected, or reassigned.",
            "",
        ]
    )
    path.write_text("\n".join(lines), encoding="utf-8")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--osteosarc-repo", type=Path, required=True)
    parser.add_argument("--netchop-dir", type=Path, required=True)
    parser.add_argument("--netcleave-dir", type=Path, required=True)
    parser.add_argument("--eramer-dir", type=Path, required=True)
    parser.add_argument(
        "--output-dir", type=Path, default=Path(__file__).resolve().parent / "results"
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    output_dir = args.output_dir.resolve()
    tables_dir = output_dir / "tables"
    figures_dir = output_dir / "figures"
    tables_dir.mkdir(parents=True, exist_ok=True)
    figures_dir.mkdir(parents=True, exist_ok=True)

    variants_path = args.osteosarc_repo / "src/data/variants.json"
    source_commit = git_value(args.osteosarc_repo, "%H")
    source_date = git_value(args.osteosarc_repo, "%cI")
    inventory_df, missing_df, conflicts_df = extract_inventory(variants_path)
    records = slp_records(inventory_df)

    inventory_df.to_csv(tables_dir / "vaccine_sequence_inventory.csv", index=False)
    missing_df.to_csv(tables_dir / "undisclosed_vaccine_sequences.csv", index=False)
    conflicts_df.to_csv(tables_dir / "sequence_provenance_conflicts.csv", index=False)

    quantitative_rows: list[dict[str, Any]] = []
    model_rows: list[dict[str, Any]] = []
    rows, models = pepsickle_rows(records)
    quantitative_rows.extend(rows)
    model_rows.extend(models)
    rows, models = netchop_rows(records, args.netchop_dir)
    quantitative_rows.extend(rows)
    model_rows.extend(models)
    rows, models = netcleave_rows(records, args.netcleave_dir)
    quantitative_rows.extend(rows)
    model_rows.extend(models)
    rows, motif_rows, models = peptidase_rows(records, args.eramer_dir)
    quantitative_rows.extend(rows)
    model_rows.extend(models)

    quantitative_df = pd.DataFrame(quantitative_rows)
    motifs_df = pd.DataFrame(motif_rows)
    model_catalog_df = pd.DataFrame(model_rows).drop_duplicates("model", keep="last")
    summary_df = summarize_models(records, quantitative_df, motifs_df, model_catalog_df)
    validate_analysis(
        inventory_df,
        records,
        quantitative_df,
        motifs_df,
        model_catalog_df,
        summary_df,
    )

    quantitative_df.to_csv(tables_dir / "slp_quantitative_bond_scores.csv", index=False)
    motifs_df.to_csv(tables_dir / "slp_motif_assessments.csv", index=False)
    motifs_df.loc[motifs_df["status"] == "matched"].to_csv(
        tables_dir / "slp_motif_matches.csv", index=False
    )
    model_catalog_df.to_csv(tables_dir / "model_catalog.csv", index=False)
    summary_df.to_csv(tables_dir / "slp_model_summary.csv", index=False)
    predictor_matrix(records, summary_df).to_csv(
        tables_dir / "slp_predictor_matrix.csv", index=False
    )
    correlations(quantitative_df).to_csv(
        tables_dir / "quantitative_model_spearman_correlations.csv"
    )

    included_motifs = set(
        model_catalog_df.loc[
            model_catalog_df["included"]
            & (model_catalog_df["evidence_type"] == "motif_rule"),
            "model",
        ]
    )
    grouped_motifs = set(EXTRACELLULAR_MOTIF_MODELS) | set(
        INTRACELLULAR_ER_MOTIF_MODELS
    )
    if included_motifs != grouped_motifs:
        raise RuntimeError(
            "Motif figure groups do not cover the included models exactly: "
            f"missing={sorted(included_motifs - grouped_motifs)}, "
            f"extra={sorted(grouped_motifs - included_motifs)}"
        )

    old_motif_figure = figures_dir / "slp_peptidase_motif_matches.png"
    old_motif_figure.unlink(missing_ok=True)
    pdf_metadata = {
        "Title": "Osteosarc vaccine SLP cleavage prediction figures",
        "Author": "mhctools",
        "Subject": "Model-specific cleavage scores and peptidase motif matches",
        "Keywords": "osteosarcoma vaccine SLP cleavage proteasome peptidase",
        "CreationDate": datetime.fromisoformat(source_date),
        "ModDate": datetime.fromisoformat(source_date),
    }
    with PdfPages(output_dir / "all-figures.pdf", metadata=pdf_metadata) as pdf_pages:
        plot_score_distributions(quantitative_df, figures_dir, pdf_pages)
        plot_quantitative_heatmap(records, summary_df, figures_dir, pdf_pages)
        plot_motif_heatmap(
            records,
            summary_df,
            figures_dir,
            pdf_pages,
            EXTRACELLULAR_MOTIF_MODELS,
            "Extracellular and plasma peptidase motif matches",
            "slp_peptidase_motifs_extracellular",
            figure_number=3,
        )
        plot_motif_heatmap(
            records,
            summary_df,
            figures_dir,
            pdf_pages,
            INTRACELLULAR_ER_MOTIF_MODELS,
            "Intracellular and ER peptidase motif matches",
            "slp_peptidase_motifs_intracellular_er",
            figure_number=4,
        )
    write_report(
        output_dir / "REPORT.md",
        inventory_df,
        records,
        missing_df,
        conflicts_df,
        summary_df,
        source_commit,
        source_date,
    )

    provenance = {
        "schema_version": 1,
        "source": {
            "site": "https://osteosarc.com/",
            "repository": "https://gitlab.com/slowkow/osteosarc.com",
            "commit": source_commit,
            "commit_date": source_date,
            "variants_json_sha256": sha256_file(variants_path),
        },
        "analysis": {
            "mhctools_version": __import__("mhctools").__version__,
            "mhctools_commit": git_value(Path(__file__).resolve().parents[2], "%H"),
            "analysis_script_sha256": sha256_file(Path(__file__)),
            "display_threshold": THRESHOLD,
            "termini_assumption": "free N and C termini",
            "sequence_uploads": "none; all inference ran locally",
        },
        "runtime": {
            "python": platform.python_version(),
            "packages": {
                name: importlib.metadata.version(name)
                for name in (
                    "keras",
                    "matplotlib",
                    "numpy",
                    "pandas",
                    "scikit-learn",
                    "tensorflow",
                )
            },
        },
        "models": model_hashes(args.netchop_dir, args.netcleave_dir, args.eramer_dir),
    }
    (output_dir / "provenance.json").write_text(
        json.dumps(provenance, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )

    manifest: dict[str, str] = {}
    for path in sorted(output_dir.rglob("*")):
        if path.is_file() and path.name != "SHA256SUMS.json":
            manifest[str(path.relative_to(output_dir))] = sha256_file(path)
    (output_dir / "SHA256SUMS.json").write_text(
        json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )


if __name__ == "__main__":
    main()
