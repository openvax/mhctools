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
from scipy.cluster.hierarchy import leaves_list, linkage
from scipy.spatial.distance import squareform

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
QUANTITATIVE_SITE_MODELS = [
    "netchop-3.1-20s-3.0",
    "pepsickle-in-vivo-human-only",
    "pepsickle-in-vivo-all-mammal",
    "netchop-3.1-cterm-3.0",
    "netcleave-i-hla",
    "netcleave-ii-hla",
]
MODEL_DISPLAY_NAMES = {
    "netchop-3.1-20s-3.0": "NetChop 20S",
    "pepsickle-in-vivo-human-only": "Pepsickle human",
    "pepsickle-in-vivo-all-mammal": "Pepsickle all-mammal",
    "netchop-3.1-cterm-3.0": "NetChop Cterm",
    "netcleave-i-hla": "NetCleave I",
    "netcleave-ii-hla": "NetCleave II",
}
MOTIF_DISPLAY_NAMES = {
    "ace-dipeptidyl": "ACE - dipeptidyl",
    "mme-hydrophobic": "MME - hydrophobic",
    "cpb2-basic": "CPB2 - basic C-term",
    "cpn-basic": "CPN1 - basic C-term",
    "app2-xp": "XPNPEP2 - X|Pro",
    "fap-dipeptidyl": "FAP - dipeptidyl",
    "fap-endo-gp": "FAP - Gly|Pro",
    "enpep-acidic": "ENPEP - acidic N-term",
    "anpep-ala": "ANPEP - Ala N-term",
    "app1-xp": "XPNPEP1 - X|Pro",
    "tpp2-tripeptidyl": "TPP2 - tripeptidyl",
    "npepps-n-terminal": "NPEPPS - N-term",
    "dpp8-xp-xa": "DPP8 - X-Pro|X",
    "dpp9-xp-xa": "DPP9 - X-Pro|X",
    "prep-pro": "PREP - Pro-associated",
    "erap2-basic": "ERAP2 - basic N-term",
}
MAX_BONDS_PER_ATLAS_PAGE = 27


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


def save_pdf_page(
    fig: plt.Figure,
    pdf_pages: PdfPages,
    page_number: int,
    page_count: int,
    png_path: Path | None = None,
) -> None:
    fig.text(
        0.995,
        0.005,
        f"Osteosarc vaccine cleavage atlas - page {page_number} of {page_count}",
        ha="right",
        va="bottom",
        fontsize=7,
        color="#555555",
    )
    if png_path is not None:
        fig.savefig(png_path, dpi=180, bbox_inches="tight")
    pdf_pages.savefig(fig, bbox_inches="tight")
    plt.close(fig)


def clustered_orders(
    records: pd.DataFrame,
    quantitative_df: pd.DataFrame,
    summary_df: pd.DataFrame,
) -> tuple[list[str], list[str], pd.DataFrame, pd.DataFrame]:
    """Cluster predictors by bond scores and SLPs by threshold-hit profiles."""
    site_scores = quantitative_df.loc[
        quantitative_df["model"].isin(QUANTITATIVE_SITE_MODELS)
        & quantitative_df["assessable"]
    ].pivot_table(
        index=["sequence_record_id", "bond"], columns="model", values="score"
    )
    correlations_df = site_scores.corr(method="spearman", min_periods=10).reindex(
        index=QUANTITATIVE_SITE_MODELS, columns=QUANTITATIVE_SITE_MODELS
    )
    distances = np.clip(1.0 - correlations_df.to_numpy(dtype=float), 0.0, 2.0)
    distances = (distances + distances.T) / 2.0
    np.fill_diagonal(distances, 0.0)
    model_tree = linkage(
        squareform(distances, checks=False), method="average", optimal_ordering=True
    )
    model_order = [
        QUANTITATIVE_SITE_MODELS[index] for index in leaves_list(model_tree)
    ]

    profiles = summary_df.loc[
        summary_df["model"].isin(QUANTITATIVE_SITE_MODELS)
    ].pivot(
        index="sequence_record_id", columns="model", values="candidate_fraction"
    )
    profiles = profiles.reindex(
        index=records["sequence_record_id"], columns=QUANTITATIVE_SITE_MODELS
    )
    clustering_values = profiles.copy()
    for column in clustering_values:
        clustering_values[column] = clustering_values[column].fillna(
            clustering_values[column].median()
        )
        standard_deviation = clustering_values[column].std(ddof=0)
        if standard_deviation > 0:
            clustering_values[column] = (
                clustering_values[column] - clustering_values[column].mean()
            ) / standard_deviation
    record_tree = linkage(
        clustering_values.to_numpy(dtype=float),
        method="average",
        metric="euclidean",
        optimal_ordering=True,
    )
    record_order = [profiles.index[index] for index in leaves_list(record_tree)]
    return record_order, model_order, profiles, correlations_df


def plot_agreement_overview(
    records: pd.DataFrame,
    quantitative_df: pd.DataFrame,
    summary_df: pd.DataFrame,
    figures_dir: Path,
    pdf_pages: PdfPages,
    page_count: int,
) -> tuple[list[str], list[str]]:
    record_order, model_order, profiles, correlations_df = clustered_orders(
        records, quantitative_df, summary_df
    )
    labels = row_labels(records)
    ordered_profiles = profiles.reindex(index=record_order, columns=model_order)
    ordered_correlations = correlations_df.reindex(
        index=model_order, columns=model_order
    )
    fig = plt.figure(figsize=(16, 10.5))
    grid = fig.add_gridspec(
        1,
        2,
        width_ratios=(1.55, 1),
        left=0.18,
        right=0.95,
        top=0.87,
        bottom=0.16,
        wspace=0.48,
    )
    ax_profiles = fig.add_subplot(grid[0, 0])
    profile_cmap = plt.get_cmap("magma").copy()
    profile_cmap.set_bad("#d9d9d9")
    profile_image = ax_profiles.imshow(
        np.ma.masked_invalid(ordered_profiles.to_numpy(dtype=float)),
        aspect="auto",
        vmin=0,
        vmax=1,
        cmap=profile_cmap,
    )
    ax_profiles.set_xticks(
        range(len(model_order)),
        labels=[MODEL_DISPLAY_NAMES[model] for model in model_order],
        rotation=35,
        ha="right",
        fontsize=8,
    )
    ax_profiles.set_yticks(
        range(len(record_order)),
        labels=[labels[record_id] for record_id in record_order],
        fontsize=6.5,
    )
    ax_profiles.set_title("SLPs clustered by within-model candidate-site fraction")
    fig.colorbar(
        profile_image,
        ax=ax_profiles,
        fraction=0.025,
        pad=0.02,
        label="fraction of assessable bonds >= 0.5",
    )

    ax_correlations = fig.add_subplot(grid[0, 1])
    correlation_image = ax_correlations.imshow(
        ordered_correlations.to_numpy(dtype=float),
        cmap="coolwarm",
        vmin=-1,
        vmax=1,
    )
    correlation_labels = [MODEL_DISPLAY_NAMES[model] for model in model_order]
    ax_correlations.set_xticks(
        range(len(model_order)),
        labels=correlation_labels,
        rotation=35,
        ha="right",
        fontsize=8,
    )
    ax_correlations.set_yticks(
        range(len(model_order)), labels=correlation_labels, fontsize=8
    )
    ax_correlations.set_title("Predictor agreement at shared peptide bonds")
    for row in range(len(model_order)):
        for column in range(len(model_order)):
            value = ordered_correlations.iloc[row, column]
            ax_correlations.text(
                column,
                row,
                f"{value:.2f}",
                ha="center",
                va="center",
                fontsize=8,
                color="white" if abs(value) > 0.55 else "black",
            )
    fig.colorbar(
        correlation_image,
        ax=ax_correlations,
        fraction=0.046,
        pad=0.04,
        label="Spearman rho",
    )
    fig.suptitle(
        "Agreement map for the vaccine SLP cleavage atlas",
        fontsize=16,
        fontweight="bold",
        y=0.955,
    )
    fig.text(
        0.5,
        0.06,
        "Clustering organizes the atlas; it does not create an ensemble probability. "
        "Sequence pages follow the left-panel SLP order. "
        "The left panel uses a common display threshold on native 0-1 outputs. The right "
        "panel correlates native scores only where both models assess the same bond.",
        ha="center",
        va="center",
        fontsize=9,
        color="#333333",
        wrap=True,
    )
    save_pdf_page(
        fig,
        pdf_pages,
        page_number=1,
        page_count=page_count,
        png_path=figures_dir / "predictor_agreement_and_slp_clusters.png",
    )
    return record_order, model_order


def sequence_segments(sequence: str) -> list[tuple[int, int]]:
    """Return inclusive 1-based bond ranges that cover a sequence."""
    last_bond = len(sequence) - 1
    return [
        (start, min(start + MAX_BONDS_PER_ATLAS_PAGE - 1, last_bond))
        for start in range(1, last_bond + 1, MAX_BONDS_PER_ATLAS_PAGE)
    ]


def _score_matrix(
    record_id: str,
    quantitative_df: pd.DataFrame,
    models: list[str],
    bonds: list[int],
) -> np.ndarray:
    matrix = np.full((len(models), len(bonds)), np.nan)
    bond_indices = {bond: index for index, bond in enumerate(bonds)}
    model_indices = {model: index for index, model in enumerate(models)}
    subset = quantitative_df.loc[
        (quantitative_df["sequence_record_id"] == record_id)
        & quantitative_df["model"].isin(models)
        & quantitative_df["assessable"]
    ]
    for row in subset.itertuples():
        bond = int(row.bond)
        if bond in bond_indices:
            matrix[model_indices[row.model], bond_indices[bond]] = float(row.score)
    return matrix


def _motif_matrix(
    record_id: str,
    motifs_df: pd.DataFrame,
    models: list[str],
    bonds: list[int],
) -> np.ndarray:
    matrix = np.full((len(models), len(bonds)), np.nan)
    bond_indices = {bond: index for index, bond in enumerate(bonds)}
    model_indices = {model: index for index, model in enumerate(models)}
    subset = motifs_df.loc[
        (motifs_df["sequence_record_id"] == record_id)
        & motifs_df["model"].isin(models)
        & (motifs_df["status"] != "unsupported")
    ]
    for row in subset.itertuples():
        bond = int(row.bond)
        if bond in bond_indices:
            matrix[model_indices[row.model], bond_indices[bond]] = (
                1.0 if row.status == "matched" else 0.0
            )
    return matrix


def _draw_grid(ax: plt.Axes, rows: int, columns: int) -> None:
    ax.set_xticks(np.arange(-0.5, columns, 1), minor=True)
    ax.set_yticks(np.arange(-0.5, rows, 1), minor=True)
    ax.grid(which="minor", color="#dddddd", linewidth=0.35)
    ax.tick_params(which="minor", bottom=False, left=False)


def plot_sequence_atlas_page(
    record: pd.Series,
    quantitative_df: pd.DataFrame,
    motifs_df: pd.DataFrame,
    pdf_pages: PdfPages,
    page_number: int,
    page_count: int,
    segment: tuple[int, int],
    segment_number: int,
    segment_count: int,
) -> None:
    sequence = record["sequence"]
    start_bond, end_bond = segment
    bonds = list(range(start_bond, end_bond + 1))
    record_id = record["sequence_record_id"]
    score_matrix = _score_matrix(
        record_id, quantitative_df, QUANTITATIVE_SITE_MODELS, bonds
    )
    motif_models = EXTRACELLULAR_MOTIF_MODELS + INTRACELLULAR_ER_MOTIF_MODELS
    motif_matrix = _motif_matrix(record_id, motifs_df, motif_models, bonds)

    fig = plt.figure(figsize=(16, 10.5))
    grid = fig.add_gridspec(
        5,
        1,
        height_ratios=(0.72, 0.48, 2.05, 0.64, 4.45),
        left=0.18,
        right=0.94,
        top=0.87,
        bottom=0.105,
        hspace=0.22,
    )
    x_positions = np.arange(len(bonds))
    bond_labels = [
        f"{bond}\n{sequence[bond - 1]}|{sequence[bond]}" for bond in bonds
    ]

    ax_sequence = fig.add_subplot(grid[0, 0])
    ax_sequence.set_xlim(-0.5, len(bonds) - 0.5)
    ax_sequence.set_ylim(0, 1)
    ax_sequence.axis("off")
    for x, bond in zip(x_positions, bonds):
        ax_sequence.text(
            x,
            0.55,
            f"{sequence[bond - 1]}|{sequence[bond]}",
            ha="center",
            va="center",
            family="monospace",
            fontsize=8,
            fontweight="bold",
        )
        ax_sequence.text(
            x,
            0.08,
            str(bond),
            ha="center",
            va="bottom",
            fontsize=6.5,
            color="#555555",
        )
    ax_sequence.text(
        -0.012,
        0.55,
        "bond",
        transform=ax_sequence.transAxes,
        ha="right",
        va="center",
        fontsize=8,
        color="#555555",
    )

    ax_agreement = fig.add_subplot(grid[1, 0])
    assessed = np.sum(~np.isnan(score_matrix), axis=0)
    hits = np.sum(score_matrix >= THRESHOLD, axis=0)
    agreement = np.divide(
        hits,
        assessed,
        out=np.full_like(hits, np.nan, dtype=float),
        where=assessed > 0,
    )[None, :]
    agreement_cmap = plt.get_cmap("Blues").copy()
    agreement_cmap.set_bad("#d9d9d9")
    ax_agreement.imshow(agreement, aspect="auto", vmin=0, vmax=1, cmap=agreement_cmap)
    for column, (hit_count, assessed_count) in enumerate(zip(hits, assessed)):
        label = "NA" if assessed_count == 0 else f"{hit_count}/{assessed_count}"
        ax_agreement.text(
            column,
            0,
            label,
            ha="center",
            va="center",
            fontsize=6.5,
            color="white" if assessed_count and hit_count / assessed_count > 0.55 else "black",
        )
    ax_agreement.set_yticks([0], labels=["models >= 0.5"])
    ax_agreement.set_xticks([])
    ax_agreement.tick_params(axis="y", labelsize=8)
    _draw_grid(ax_agreement, 1, len(bonds))

    ax_scores = fig.add_subplot(grid[2, 0])
    score_cmap = plt.get_cmap("magma").copy()
    score_cmap.set_bad("#d9d9d9")
    score_image = ax_scores.imshow(
        np.ma.masked_invalid(score_matrix),
        aspect="auto",
        vmin=0,
        vmax=1,
        cmap=score_cmap,
    )
    ax_scores.set_yticks(
        range(len(QUANTITATIVE_SITE_MODELS)),
        labels=[MODEL_DISPLAY_NAMES[model] for model in QUANTITATIVE_SITE_MODELS],
        fontsize=8,
    )
    ax_scores.set_xticks([])
    ax_scores.set_ylabel("native 0-1 score", fontsize=8, labelpad=34)
    for row in range(score_matrix.shape[0]):
        for column in range(score_matrix.shape[1]):
            value = score_matrix[row, column]
            if np.isnan(value):
                continue
            ax_scores.text(
                column,
                row,
                f"{value:.2f}",
                ha="center",
                va="center",
                fontsize=5.8,
                color="white" if value < 0.62 else "black",
            )
    _draw_grid(ax_scores, len(QUANTITATIVE_SITE_MODELS), len(bonds))
    colorbar = fig.colorbar(score_image, ax=ax_scores, fraction=0.018, pad=0.012)
    colorbar.ax.tick_params(labelsize=7)
    colorbar.set_label("native model output", fontsize=8)

    ax_terminal = fig.add_subplot(grid[3, 0])
    terminal_models = ["dpp4-qpisa", "eramer-step"]
    terminal_labels = ["DPP4 native score", "ERAP1 ERAMER score"]
    ax_terminal.set_xlim(-0.5, len(bonds) - 0.5)
    ax_terminal.set_ylim(1.5, -0.5)
    ax_terminal.set_yticks(range(2), labels=terminal_labels, fontsize=8)
    ax_terminal.set_xticks([])
    ax_terminal.set_facecolor("#f7f7f7")
    terminal_subset = quantitative_df.loc[
        (quantitative_df["sequence_record_id"] == record_id)
        & quantitative_df["model"].isin(terminal_models)
        & quantitative_df["assessable"]
    ]
    for row in terminal_subset.itertuples():
        bond = int(row.bond)
        if bond not in bonds:
            continue
        y = terminal_models.index(row.model)
        x = bonds.index(bond)
        ax_terminal.scatter(
            [x], [y], marker="D", s=80, color="#f0a202", edgecolor="#3b2a00", zorder=2
        )
        ax_terminal.text(
            x + 0.35,
            y,
            f"{float(row.score):.3f}",
            ha="left",
            va="center",
            fontsize=7,
            color="#3b2a00",
        )
    ax_terminal.text(
        1.003,
        0.5,
        "incompatible native scales; no cutoff",
        transform=ax_terminal.transAxes,
        ha="left",
        va="center",
        fontsize=7,
        color="#555555",
    )
    _draw_grid(ax_terminal, 2, len(bonds))

    ax_motifs = fig.add_subplot(grid[4, 0])
    motif_colors = np.empty((*motif_matrix.shape, 4))
    motif_colors[:] = (0.85, 0.85, 0.85, 1.0)
    assessed_cells = ~np.isnan(motif_matrix)
    motif_colors[assessed_cells] = (1.0, 1.0, 1.0, 1.0)
    for row in range(len(motif_models)):
        matched = motif_matrix[row] == 1
        motif_colors[row, matched] = (
            (0.0, 0.48, 0.52, 1.0)
            if row < len(EXTRACELLULAR_MOTIF_MODELS)
            else (0.43, 0.24, 0.62, 1.0)
        )
    ax_motifs.imshow(motif_colors, aspect="auto")
    ax_motifs.set_yticks(
        range(len(motif_models)),
        labels=[MOTIF_DISPLAY_NAMES[model] for model in motif_models],
        fontsize=7.2,
    )
    for label_index, label in enumerate(ax_motifs.get_yticklabels()):
        label.set_color(
            "#006d73"
            if label_index < len(EXTRACELLULAR_MOTIF_MODELS)
            else "#60408a"
        )
    ax_motifs.set_xticks(x_positions, labels=bond_labels, fontsize=6.5)
    ax_motifs.tick_params(axis="x", pad=3)
    ax_motifs.set_xlabel("cleavage after numbered left residue", fontsize=8)
    ax_motifs.axhline(
        len(EXTRACELLULAR_MOTIF_MODELS) - 0.5,
        color="#333333",
        linewidth=1.4,
    )
    for row in range(motif_matrix.shape[0]):
        for column in range(motif_matrix.shape[1]):
            if motif_matrix[row, column] == 1:
                ax_motifs.text(
                    column,
                    row,
                    "●",
                    ha="center",
                    va="center",
                    fontsize=6,
                    color="white",
                )
    _draw_grid(ax_motifs, len(motif_models), len(bonds))
    ax_motifs.text(
        1.003,
        0.74,
        "extracellular / plasma",
        transform=ax_motifs.transAxes,
        ha="left",
        va="center",
        fontsize=7,
        color="#006d73",
        rotation=90,
    )
    ax_motifs.text(
        1.003,
        0.22,
        "cytosol / ER",
        transform=ax_motifs.transAxes,
        ha="left",
        va="center",
        fontsize=7,
        color="#60408a",
        rotation=90,
    )

    vaccines = record["vaccines"].replace(";", ", ")
    continuation = (
        f" - segment {segment_number}/{segment_count}, bonds {start_bond}-{end_bond}"
        if segment_count > 1
        else ""
    )
    fig.suptitle(
        f"{record['gene']} {record['protein_change'] or ''} | {vaccines} | "
        f"{record['length']} aa{continuation}",
        fontsize=15,
        fontweight="bold",
        y=0.965,
    )
    displayed_sequence = sequence[start_bond - 1 : end_bond + 1]
    sequence_prefix = (
        f"residues {start_bond}-{end_bond + 1}: " if segment_count > 1 else ""
    )
    fig.text(
        0.5,
        0.913,
        sequence_prefix + " ".join(displayed_sequence),
        ha="center",
        va="center",
        family="monospace",
        fontsize=9,
        color="#222222",
    )
    fig.text(
        0.5,
        0.047,
        "Figure legend. Scores are model-native and not calibrated across rows; 0.5 is a display "
        "threshold only. Concurrence is hits/assessable models, not a cleavage probability. Filled "
        "motif cells are partial recognition-rule matches; white is assessed/no match and gray is "
        "not assessed. Proteasome and ER evidence is conditional on intracellular access.",
        ha="center",
        va="center",
        fontsize=8,
        color="#333333",
        wrap=True,
    )
    save_pdf_page(fig, pdf_pages, page_number, page_count)


def render_figures(
    output_dir: Path,
    records: pd.DataFrame,
    quantitative_df: pd.DataFrame,
    motifs_df: pd.DataFrame,
    summary_df: pd.DataFrame,
    source_date: str,
) -> pd.DataFrame:
    """Render the clustered overview and bond-aligned sequence atlas."""
    figures_dir = output_dir / "figures"
    figures_dir.mkdir(parents=True, exist_ok=True)
    page_count = 1 + sum(
        len(sequence_segments(record["sequence"])) for _, record in records.iterrows()
    )
    pdf_metadata = {
        "Title": "Osteosarc vaccine SLP bond-level cleavage atlas",
        "Author": "mhctools",
        "Subject": "Sequence-aligned cleavage scores and peptidase motif evidence",
        "Keywords": "osteosarcoma vaccine SLP cleavage proteasome peptidase",
        "CreationDate": datetime.fromisoformat(source_date),
        "ModDate": datetime.fromisoformat(source_date),
    }
    atlas_rows: list[dict[str, Any]] = []
    with PdfPages(output_dir / "all-figures.pdf", metadata=pdf_metadata) as pdf_pages:
        record_order, model_order = plot_agreement_overview(
            records,
            quantitative_df,
            summary_df,
            figures_dir,
            pdf_pages,
            page_count,
        )
        records_by_id = records.set_index("sequence_record_id", drop=False)
        page_number = 2
        for atlas_order, record_id in enumerate(record_order, start=1):
            record = records_by_id.loc[record_id]
            segments = sequence_segments(record["sequence"])
            first_page = page_number
            for segment_number, segment in enumerate(segments, start=1):
                plot_sequence_atlas_page(
                    record,
                    quantitative_df,
                    motifs_df,
                    pdf_pages,
                    page_number,
                    page_count,
                    segment,
                    segment_number,
                    len(segments),
                )
                page_number += 1
            atlas_rows.append(
                {
                    "atlas_order": atlas_order,
                    "sequence_record_id": record_id,
                    "gene": record["gene"],
                    "protein_change": record["protein_change"],
                    "vaccines": record["vaccines"],
                    "sequence": record["sequence"],
                    "length": record["length"],
                    "pdf_page_start": first_page,
                    "pdf_page_end": page_number - 1,
                }
            )
        if page_number - 1 != page_count:
            raise RuntimeError(
                f"Rendered {page_number - 1} pages but expected {page_count}"
            )

    for obsolete_name in (
        "quantitative_score_distributions.png",
        "slp_quantitative_hit_fraction.png",
        "slp_peptidase_motifs_extracellular.png",
        "slp_peptidase_motifs_intracellular_er.png",
    ):
        (figures_dir / obsolete_name).unlink(missing_ok=True)
    atlas_order_df = pd.DataFrame(atlas_rows)
    atlas_order_df.attrs["model_cluster_order"] = model_order
    return atlas_order_df


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
        "## Sequence cleavage atlas",
        "",
        "The [complete PDF atlas](all-figures.pdf) places every quantitative score and motif-rule "
        "match at its exact peptide bond. It contains one clustered overview followed by one page "
        "per SLP (with the 80-aa outlier split across three continuation pages). "
        "[Atlas order and PDF page numbers](tables/atlas_sequence_order.csv) are provided for navigation.",
        "",
        "![Predictor agreement and clustered SLP order](figures/predictor_agreement_and_slp_clusters.png)",
        "",
        "The overview clusters predictors using Spearman correlation of native scores at shared, "
        "assessable bonds. It clusters SLPs using standardized within-model fractions above the "
        "0.5 display threshold. Clustering is organizational only and is not an ensemble model.",
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

    atlas_order_df = render_figures(
        output_dir,
        records,
        quantitative_df,
        motifs_df,
        summary_df,
        source_date,
    )
    atlas_order_df.to_csv(tables_dir / "atlas_sequence_order.csv", index=False)
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
