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
import os
from pathlib import Path
import platform
import re
import subprocess
import tempfile
from typing import Any, Iterable

import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.patches import FancyBboxPatch, Rectangle
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
from mhctools.netchop import NETCHOP_CONTAINER_IMAGE, NetChop


CANONICAL_AA = frozenset("ACDEFGHIKLMNPQRSTVWY")
SLP_VACCINES = frozenset(("JLF V1", "JLF V2", "JLF V3", "CeGaT"))
THRESHOLD = 0.5
MHC_I_DISPLAY_RANK = 2.0
MHC_II_DISPLAY_RANK = 5.0
PDF_FILENAME = "mhctools-all-figures.pdf"
MANUSCRIPT_PDF_FILENAME = "mhctools-manuscript-figures.pdf"
MHC_I_ALLELES = (
    "HLA-A*01:01",
    "HLA-B*08:01",
    "HLA-B*27:05",
    "HLA-C*01:02",
    "HLA-C*07:01",
)
# These are the class-II combinations already named in the osteosarc source's
# candidate-prediction fields. The HLA table itself is not phased, so do not
# manufacture additional alpha/beta pairings from it.
MHC_II_ALLELES = (
    "HLA-DPA1*01:03-DPB1*04:01",
    "HLA-DQA1*04:01-DQB1*04:02",
    "HLA-DQA1*05:01-DQB1*02:01",
    "HLA-DRB1*03:01",
    "HLA-DRB1*08:01",
)
NETCHOP_IMAGE = NETCHOP_CONTAINER_IMAGE
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
MHC_I_CLEAVAGE_MODELS = QUANTITATIVE_SITE_MODELS[:5]
FIGURE_CYTOSOL_MODELS = [
    "netchop-3.1-20s-3.0",
    "netcleave-i-hla",
    "pepsickle-in-vivo-human-only",
    "netchop-3.1-cterm-3.0",
]
FIGURE_QUANTITATIVE_MODELS = FIGURE_CYTOSOL_MODELS + ["netcleave-ii-hla"]
FIGURE_RED_SUPPORT_REQUIRED = 3
MODEL_DISPLAY_NAMES = {
    "netchop-3.1-20s-3.0": "NetChop 20S",
    "pepsickle-in-vivo-human-only": "Pepsickle human",
    "pepsickle-in-vivo-all-mammal": "Pepsickle all-mammal",
    "netchop-3.1-cterm-3.0": "NetChop Cterm",
    "netcleave-i-hla": "NetCleave I",
    "netcleave-ii-hla": "NetCleave II",
}
ATLAS_MODEL_LABELS = {
    "netchop-3.1-20s-3.0": "NetChop 20S (in vitro)",
    "pepsickle-in-vivo-human-only": "Pepsickle human (in vivo)",
    "netchop-3.1-cterm-3.0": "NetChop Cterm (ligand)",
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
                    "minimal_epitope": variant.get("minimal_epitope"),
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


def peptide_windows(
    sequence: str, lengths: Iterable[int]
) -> list[tuple[int, int, str]]:
    """Return 1-based inclusive peptide windows for the requested lengths."""
    return [
        (start + 1, start + length, sequence[start : start + length])
        for length in lengths
        if length <= len(sequence)
        for start in range(len(sequence) - length + 1)
    ]


def _minimal_epitope_bounds(record: pd.Series) -> tuple[int, int] | None:
    epitope = record.get("minimal_epitope")
    offset = record.get("minimal_epitope_offset")
    if not epitope or pd.isna(offset):
        return None
    start = int(offset) + 1
    end = start + len(epitope) - 1
    if record["sequence"][start - 1 : end] != epitope:
        raise ValueError(
            f"Minimal epitope offset does not match {record['sequence_record_id']}: "
            f"expected {epitope!r} at {start}-{end}"
        )
    return start, end


def _window_overlaps_minimal(record: pd.Series, start: int, end: int) -> bool:
    bounds = _minimal_epitope_bounds(record)
    return bool(bounds and start <= bounds[1] and end >= bounds[0])


def mhcflurry_ligand_rows(
    records: pd.DataFrame,
) -> tuple[list[dict[str, Any]], dict[str, Any]]:
    """Predict class-I ligand windows locally with MHCflurry presentation."""
    from mhcflurry import Class1PresentationPredictor
    from mhcflurry.downloads import (
        get_current_release,
        get_default_class1_models_dir,
    )

    predictor = Class1PresentationPredictor.load()
    metadata: list[tuple[pd.Series, int, int]] = []
    peptides: list[str] = []
    n_flanks: list[str] = []
    c_flanks: list[str] = []
    for _, record in records.iterrows():
        sequence = record["sequence"]
        for start, end, peptide in peptide_windows(sequence, range(8, 12)):
            metadata.append((record, start, end))
            peptides.append(peptide)
            n_flanks.append(sequence[max(0, start - 16) : start - 1])
            c_flanks.append(sequence[end : min(len(sequence), end + 15)])

    sample_to_allele = {
        f"class_i_{i}": [allele] for i, allele in enumerate(MHC_I_ALLELES)
    }
    frame = predictor.predict(
        peptides=peptides,
        alleles=sample_to_allele,
        n_flanks=n_flanks,
        c_flanks=c_flanks,
        include_affinity_percentile=True,
        verbose=0,
    )
    rows: list[dict[str, Any]] = []
    for result in frame.itertuples(index=False):
        record, start, end = metadata[int(result.peptide_num)]
        rank = float(result.presentation_percentile)
        rows.append(
            {
                "sequence_record_id": record["sequence_record_id"],
                "gene": record["gene"],
                "mhc_class": "I",
                "predictor": "MHCflurry",
                "model_version": f"models_class1_pan/{get_current_release()}",
                "allele": result.best_allele,
                "peptide": result.peptide,
                "start": start,
                "end": end,
                "length": end - start + 1,
                "score": float(result.presentation_score),
                "percentile_rank": rank,
                "rank_threshold": MHC_I_DISPLAY_RANK,
                "display_candidate": rank <= MHC_I_DISPLAY_RANK,
                "affinity_nM": float(result.affinity),
                "affinity_percentile": float(result.affinity_percentile),
                "processing_score": float(result.processing_score),
                "binding_core": "",
                "n_flank": result.n_flank,
                "c_flank": result.c_flank,
                "overlaps_disclosed_minimal_epitope": _window_overlaps_minimal(
                    record, start, end
                ),
            }
        )
    return rows, {
        "package_version": importlib.metadata.version("mhcflurry"),
        "model_release": get_current_release(),
        "model_provenance": predictor.provenance_string,
        "models_path": str(get_default_class1_models_dir()),
    }


def _parse_netmhciipan_rows(stdout: str) -> list[dict[str, Any]]:
    """Parse the stable tabular fields needed from NetMHCIIpan 4.3 output."""
    rows: list[dict[str, Any]] = []
    for line in stdout.splitlines():
        fields = line.split()
        if len(fields) < 11 or not fields[0].isdigit():
            continue
        try:
            score = float(fields[8])
            rank = float(fields[9])
        except ValueError:
            continue
        rows.append(
            {
                "allele": fields[1],
                "peptide": fields[2],
                "binding_core": fields[4],
                "score": score,
                "percentile_rank": rank,
            }
        )
    if not rows:
        tail = "\n".join(stdout.splitlines()[-20:])
        raise ValueError(f"No NetMHCIIpan 4.3 predictions parsed. Output tail:\n{tail}")
    return rows


def netmhciipan_cli_allele(allele: str) -> str:
    """Convert normalized human class-II notation to NetMHCIIpan CLI notation."""
    value = allele.removeprefix("HLA-")
    if value.startswith("DRB"):
        return value.replace("*", "_").replace(":", "")
    return "HLA-" + value.replace("*", "").replace(":", "")


def netmhciipan_ligand_rows(
    records: pd.DataFrame, netmhciipan_path: Path
) -> tuple[list[dict[str, Any]], dict[str, Any]]:
    """Predict class-II ligand windows locally with NetMHCIIpan 4.3 EL mode."""
    occurrences: dict[str, list[tuple[pd.Series, int, int]]] = {}
    for _, record in records.iterrows():
        for start, end, peptide in peptide_windows(record["sequence"], range(13, 22)):
            occurrences.setdefault(peptide, []).append((record, start, end))

    resolved = netmhciipan_path.resolve()
    installation_root = resolved.parent
    bundle_root = installation_root.parent
    environment = os.environ.copy()
    environment["NETMHC_BUNDLE_HOME"] = str(bundle_root)
    environment["NETMHC_BUNDLE_TMPDIR"] = tempfile.gettempdir()
    cli_to_normalized = {
        netmhciipan_cli_allele(allele): allele for allele in MHC_II_ALLELES
    }
    with tempfile.NamedTemporaryFile("w", suffix=".txt", encoding="ascii") as handle:
        handle.write("\n".join(occurrences))
        handle.flush()
        completed = subprocess.run(
            [
                str(netmhciipan_path),
                "-f",
                handle.name,
                "-inptype",
                "1",
                "-a",
                ",".join(cli_to_normalized),
                "-BA",
            ],
            check=True,
            capture_output=True,
            text=True,
            env=environment,
        )
    # The tcsh launcher writes its platform line to stdout while some portable
    # builds emit the predictor table on stderr. Preserve and parse both.
    predictor_output = completed.stdout + "\n" + completed.stderr
    parsed = _parse_netmhciipan_rows(predictor_output)
    version_line = next(
        (
            line.lstrip("# ")
            for line in predictor_output.splitlines()
            if "version 4.3" in line
        ),
        "NetMHCIIpan version 4.3",
    )
    rows: list[dict[str, Any]] = []
    for result in parsed:
        rank = result["percentile_rank"]
        for record, start, end in occurrences[result["peptide"]]:
            rows.append(
                {
                    "sequence_record_id": record["sequence_record_id"],
                    "gene": record["gene"],
                    "mhc_class": "II",
                    "predictor": "NetMHCIIpan",
                    "model_version": version_line,
                    "allele": cli_to_normalized[result["allele"]],
                    "peptide": result["peptide"],
                    "start": start,
                    "end": end,
                    "length": end - start + 1,
                    "score": result["score"],
                    "percentile_rank": rank,
                    "rank_threshold": MHC_II_DISPLAY_RANK,
                    "display_candidate": rank <= MHC_II_DISPLAY_RANK,
                    "affinity_nM": np.nan,
                    "affinity_percentile": np.nan,
                    "processing_score": np.nan,
                    "binding_core": result["binding_core"],
                    "n_flank": "",
                    "c_flank": "",
                    "overlaps_disclosed_minimal_epitope": _window_overlaps_minimal(
                        record, start, end
                    ),
                }
            )
    return rows, {
        "model_version": version_line,
        "program_path": str(resolved),
        "installation_root": str(installation_root),
    }


def annotate_ligand_cleavage_exposure(
    ligand_df: pd.DataFrame, quantitative_df: pd.DataFrame
) -> pd.DataFrame:
    """Annotate ligand spans with bond-level evidence; never aggregate to probability."""
    output = ligand_df.copy()
    lookup: dict[tuple[str, str, int], tuple[int, int]] = {}
    for (record_id, model, bond), group in quantitative_df.loc[
        quantitative_df["assessable"]
    ].groupby(["sequence_record_id", "model", "bond"]):
        lookup[(record_id, model, int(bond))] = (
            int((group["score"] >= THRESHOLD).sum()),
            len(group),
        )

    internal_values: list[str] = []
    n_values: list[str] = []
    c_values: list[str] = []
    for row in output.itertuples(index=False):
        models = FIGURE_CYTOSOL_MODELS if row.mhc_class == "I" else ["netcleave-ii-hla"]

        def support(bond: int) -> tuple[int, int]:
            values = [
                lookup.get((row.sequence_record_id, model, bond)) for model in models
            ]
            values = [value for value in values if value is not None]
            return sum(value[0] for value in values), sum(value[1] for value in values)

        internal: list[str] = []
        for bond in range(int(row.start), int(row.end)):
            hits, assessed = support(bond)
            required = FIGURE_RED_SUPPORT_REQUIRED if row.mhc_class == "I" else 1
            all_display_models_assessed = assessed == len(models)
            if all_display_models_assessed and hits >= required:
                internal.append(f"{bond}({hits}/{assessed})")
        internal_values.append(";".join(internal))
        n_hits, n_assessed = support(int(row.start) - 1)
        c_hits, c_assessed = support(int(row.end))
        n_values.append("" if not n_assessed else f"{n_hits}/{n_assessed}")
        c_values.append("" if not c_assessed else f"{c_hits}/{c_assessed}")
    output["internal_candidate_cleavage_bonds"] = internal_values
    output["n_boundary_support"] = n_values
    output["c_boundary_support"] = c_values
    return output


def conservative_cut_support(
    quantitative_df: pd.DataFrame,
    record_id: str,
    bonds: Iterable[int],
    models: list[str],
    required_hits: int,
) -> dict[int, tuple[int, int]]:
    """Return fully assessed bonds meeting a conservative support-count rule."""
    requested_bonds = {int(bond) for bond in bonds}
    subset = quantitative_df.loc[
        (quantitative_df["sequence_record_id"] == record_id)
        & quantitative_df["model"].isin(models)
        & quantitative_df["bond"].isin(requested_bonds)
        & quantitative_df["assessable"],
        ["model", "bond", "score"],
    ]
    # A model gets one vote per bond even if an upstream adapter emitted a
    # duplicate row. The atlas must never turn duplicate output into support.
    unique_scores = subset.groupby(["bond", "model"], as_index=False)["score"].max()
    supported: dict[int, tuple[int, int]] = {}
    for bond, group in unique_scores.groupby("bond"):
        assessed = int(group["model"].nunique())
        hits = int((group["score"] >= THRESHOLD).sum())
        if assessed == len(models) and hits >= required_hits:
            supported[int(bond)] = (hits, assessed)
    return supported


def vulnerable_bond_table(
    records: pd.DataFrame, quantitative_df: pd.DataFrame, motifs_df: pd.DataFrame
) -> pd.DataFrame:
    """Return context-separated multi-model/rule support without probability claims."""
    record_lookup = records.set_index("sequence_record_id")
    rows: list[dict[str, Any]] = []

    def append_row(
        record_id: str, bond: int, context: str, models: list[str], rule: str
    ) -> None:
        record = record_lookup.loc[record_id]
        bounds = _minimal_epitope_bounds(record)
        rows.append(
            {
                "sequence_record_id": record_id,
                "gene": record["gene"],
                "bond": bond,
                "bond_label": f"{record['sequence'][bond - 1]}|{record['sequence'][bond]}",
                "biological_context": context,
                "support_count": len(models),
                "supporting_models_or_rules": ";".join(sorted(models)),
                "in_disclosed_minimal_epitope": bool(
                    bounds and bounds[0] <= bond < bounds[1]
                ),
                "selection_rule": rule,
                "interpretation": "support count only; not a cleavage probability",
            }
        )

    proteasome = quantitative_df.loc[
        quantitative_df["assessable"]
        & quantitative_df["model"].isin(MHC_I_CLEAVAGE_MODELS)
        & (quantitative_df["score"] >= THRESHOLD)
    ]
    for (record_id, bond), group in proteasome.groupby(["sequence_record_id", "bond"]):
        models = sorted(group["model"].unique())
        if len(models) >= 3:
            append_row(
                record_id,
                int(bond),
                "cytosolic/proteasome and MHC-I processing",
                models,
                ">=3 of 5 native 0-1 models at the 0.5 display threshold",
            )

    for context, models in (
        ("extracellular/plasma recognition motifs", EXTRACELLULAR_MOTIF_MODELS),
        ("cytosol/ER recognition motifs", INTRACELLULAR_ER_MOTIF_MODELS),
    ):
        matched = motifs_df.loc[
            (motifs_df["status"] == "matched") & motifs_df["model"].isin(models)
        ]
        for (record_id, bond), group in matched.groupby(["sequence_record_id", "bond"]):
            supporters = sorted(group["model"].unique())
            if len(supporters) >= 2:
                append_row(
                    record_id,
                    int(bond),
                    context,
                    supporters,
                    ">=2 distinct curated recognition rules match the same bond",
                )
    return pd.DataFrame(rows).sort_values(
        ["sequence_record_id", "bond", "biological_context"]
    )


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
    predictor = NetChop(
        execution="container",
        netchop_dir=netchop_dir,
        model_variant=model_variant,
    )
    predictions = predictor.cleavage_probs_many(sequences)
    return [predictions[sequence] for sequence in sequences]


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
    standalone_pdf_path: Path | None = None,
    standalone_png_path: Path | None = None,
    standalone_title: str | None = None,
    footer_label: str = "Osteosarc vaccine cleavage atlas",
) -> None:
    if standalone_pdf_path is not None or standalone_png_path is not None:
        metadata = {
            "Title": standalone_title or "Osteosarc vaccine SLP cleavage map",
            "Author": "mhctools",
            "Subject": "Sequence-aligned cleavage and MHC ligand predictions",
        }
        if standalone_pdf_path is not None:
            fig.savefig(
                standalone_pdf_path,
                bbox_inches="tight",
                metadata=metadata,
            )
        if standalone_png_path is not None:
            fig.savefig(
                standalone_png_path,
                dpi=300,
                bbox_inches="tight",
                metadata=metadata,
            )
    fig.text(
        0.995,
        0.005,
        f"{footer_label} - page {page_number} of {page_count}",
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
    """Cluster the deliberately reduced figure model set and its SLP profiles."""
    site_scores = quantitative_df.loc[
        quantitative_df["model"].isin(FIGURE_QUANTITATIVE_MODELS)
        & quantitative_df["assessable"]
    ].pivot_table(index=["sequence_record_id", "bond"], columns="model", values="score")
    correlations_df = site_scores.corr(method="spearman", min_periods=10).reindex(
        index=FIGURE_QUANTITATIVE_MODELS, columns=FIGURE_QUANTITATIVE_MODELS
    )
    distances = np.clip(1.0 - correlations_df.to_numpy(dtype=float), 0.0, 2.0)
    distances = (distances + distances.T) / 2.0
    np.fill_diagonal(distances, 0.0)
    model_tree = linkage(
        squareform(distances, checks=False), method="average", optimal_ordering=True
    )
    model_order = [
        FIGURE_QUANTITATIVE_MODELS[index] for index in leaves_list(model_tree)
    ]

    profiles = summary_df.loc[
        summary_df["model"].isin(FIGURE_QUANTITATIVE_MODELS)
    ].pivot(index="sequence_record_id", columns="model", values="candidate_fraction")
    profiles = profiles.reindex(
        index=records["sequence_record_id"], columns=FIGURE_QUANTITATIVE_MODELS
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
        "The figure shows human-only Pepsickle, both biologically distinct NetChop modes, NetCleave-I, and NetCleave-II; "
        "the near-redundant all-mammal Pepsickle track remains in the tables. The left panel uses a common display threshold on native 0-1 outputs. The right "
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


def continuous_score_runs(
    bonds: list[int], scores: np.ndarray
) -> list[tuple[np.ndarray, np.ndarray]]:
    """Return exact, contiguous bond-score runs without spanning missing values."""
    if len(bonds) != len(scores):
        raise ValueError("bonds and scores must have the same length")
    runs: list[tuple[np.ndarray, np.ndarray]] = []
    run_bonds: list[int] = []
    run_scores: list[float] = []
    for bond, score in zip(bonds, scores):
        if np.isfinite(score):
            if run_bonds and bond != run_bonds[-1] + 1:
                runs.append(
                    (np.asarray(run_bonds, dtype=float) + 0.5, np.asarray(run_scores))
                )
                run_bonds = []
                run_scores = []
            run_bonds.append(bond)
            run_scores.append(float(score))
        elif run_bonds:
            runs.append(
                (np.asarray(run_bonds, dtype=float) + 0.5, np.asarray(run_scores))
            )
            run_bonds = []
            run_scores = []
    if run_bonds:
        runs.append((np.asarray(run_bonds, dtype=float) + 0.5, np.asarray(run_scores)))
    return runs


def standalone_map_stem(
    atlas_order: int,
    record_id: str,
    segment_number: int,
    segment_count: int,
) -> str:
    """Return a stable, filesystem-safe stem for an individual map page."""
    slug = re.sub(r"[^a-z0-9]+", "-", record_id.lower()).strip("-")
    segment = (
        f"-segment-{segment_number}-of-{segment_count}" if segment_count > 1 else ""
    )
    return f"{atlas_order:02d}-{slug}{segment}"


def merge_residue_spans(spans: list[tuple[int, int]]) -> list[tuple[int, int]]:
    """Merge overlapping or adjacent inclusive residue spans."""
    merged: list[list[int]] = []
    for start, end in sorted(spans):
        if start > end:
            raise ValueError(f"invalid residue span: {start}-{end}")
        if merged and start <= merged[-1][1] + 1:
            merged[-1][1] = max(merged[-1][1], end)
        else:
            merged.append([start, end])
    return [(start, end) for start, end in merged]


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


def segment_ligand_candidates(
    ligand_df: pd.DataFrame,
    record_id: str,
    mhc_class: str,
    residue_start: int,
    residue_end: int,
) -> pd.DataFrame:
    """Select display candidates owned by one half-open sequence segment."""
    candidates = ligand_df.loc[
        (ligand_df["sequence_record_id"] == record_id)
        & (ligand_df["mhc_class"] == mhc_class)
        & ligand_df["display_candidate"]
    ].copy()
    candidates["midpoint"] = (candidates["start"] + candidates["end"]) / 2
    candidates = candidates.loc[
        (candidates["midpoint"] >= residue_start)
        & (candidates["midpoint"] < residue_end)
    ]
    return candidates.sort_values(
        ["percentile_rank", "allele", "start", "end"], kind="stable"
    ).drop_duplicates(
        ["allele", "start", "end"], keep="first"
    )


def select_ligand_candidates_for_display(
    candidates: pd.DataFrame,
    lane_count: int = 3,
    max_candidates: int = 6,
) -> pd.DataFrame:
    """Select informative non-overlapping ligand windows for fixed display lanes.

    Native percentile rank remains the only within-predictor strength ordering.
    The passes prevent a single allele from consuming every available lane:
    first retain the strongest window overlapping a disclosed intended epitope,
    then represent distinct alleles, and finally fill unused capacity by rank.
    """
    if candidates.empty:
        return candidates.assign(
            display_lane=pd.Series(dtype="int64"),
            selection_reason=pd.Series(dtype="object"),
        )
    required = {
        "allele",
        "start",
        "end",
        "percentile_rank",
        "overlaps_disclosed_minimal_epitope",
    }
    missing = required - set(candidates.columns)
    if missing:
        raise ValueError(f"ligand candidates missing columns: {sorted(missing)}")
    if lane_count < 1 or max_candidates < 1:
        raise ValueError("lane_count and max_candidates must be positive")

    ordered = candidates.sort_values(
        ["percentile_rank", "allele", "start", "end"], kind="stable"
    ).reset_index(drop=True)
    selected: list[tuple[int, str]] = []
    selected_indices: set[int] = set()
    used_spans: set[tuple[int, int]] = set()
    represented_alleles: set[str] = set()

    def lane_assignment(indices: list[int]) -> dict[int, int] | None:
        lane_ends = [-math.inf] * lane_count
        assignments: dict[int, int] = {}
        spatial_order = ordered.loc[indices].sort_values(
            ["start", "end", "percentile_rank", "allele"], kind="stable"
        )
        for index, row in spatial_order.iterrows():
            lane = next(
                (
                    lane_index
                    for lane_index, lane_end in enumerate(lane_ends)
                    if float(row["start"]) > lane_end + 0.5
                ),
                None,
            )
            if lane is None:
                return None
            lane_ends[lane] = float(row["end"])
            assignments[int(index)] = lane
        return assignments

    def try_add(index: int, reason: str) -> bool:
        if index in selected_indices or len(selected) >= max_candidates:
            return False
        row = ordered.loc[index]
        span = (int(row["start"]), int(row["end"]))
        if span in used_spans:
            return False
        if lane_assignment([item[0] for item in selected] + [index]) is None:
            return False
        selected_indices.add(index)
        used_spans.add(span)
        represented_alleles.add(str(row["allele"]))
        selected.append((index, reason))
        return True

    intended = ordered.index[
        ordered["overlaps_disclosed_minimal_epitope"].fillna(False).astype(bool)
    ]
    if len(intended):
        try_add(int(intended[0]), "intended_epitope")

    best_by_allele = (
        ordered.groupby("allele", sort=False)["percentile_rank"].min().sort_values()
    )
    for allele in best_by_allele.index:
        if str(allele) in represented_alleles:
            continue
        allele_rows = ordered.loc[ordered["allele"] == allele].copy()
        allele_rows["_intended"] = (
            allele_rows["overlaps_disclosed_minimal_epitope"]
            .fillna(False)
            .astype(bool)
        )
        allele_rows = allele_rows.sort_values(
            ["_intended", "percentile_rank", "start", "end"],
            ascending=[False, True, True, True],
            kind="stable",
        )
        for index in allele_rows.index:
            if try_add(int(index), "allele_representative"):
                break

    fill = ordered.copy()
    fill["_intended"] = (
        fill["overlaps_disclosed_minimal_epitope"].fillna(False).astype(bool)
    )
    fill = fill.sort_values(
        ["_intended", "percentile_rank", "allele", "start", "end"],
        ascending=[False, True, True, True, True],
        kind="stable",
    )
    for index in fill.index:
        try_add(int(index), "rank_fill")
        if len(selected) >= max_candidates:
            break

    assignments = lane_assignment([item[0] for item in selected])
    if assignments is None:
        raise RuntimeError("selected ligand windows cannot be assigned to display lanes")
    result: list[dict[str, Any]] = []
    for index, reason in selected:
        values = ordered.loc[index].to_dict()
        values["display_lane"] = assignments[index]
        values["selection_reason"] = reason
        result.append(values)
    return pd.DataFrame(result).reset_index(drop=True)


def mhc_display_selections(
    records: pd.DataFrame, ligand_df: pd.DataFrame
) -> pd.DataFrame:
    """Record every ligand window selected for a rendered sequence segment."""
    rows: list[dict[str, Any]] = []
    for record in records.itertuples(index=False):
        segments = sequence_segments(record.sequence)
        for segment_number, (start_bond, end_bond) in enumerate(segments, start=1):
            residue_start, residue_end = start_bond, end_bond + 1
            for mhc_class in ("I", "II"):
                eligible = segment_ligand_candidates(
                    ligand_df,
                    record.sequence_record_id,
                    mhc_class,
                    residue_start,
                    residue_end,
                )
                selected = select_ligand_candidates_for_display(eligible)
                for candidate in selected.to_dict("records"):
                    rows.append(
                        {
                            "sequence_record_id": record.sequence_record_id,
                            "gene": record.gene,
                            "mhc_class": mhc_class,
                            "segment_number": segment_number,
                            "segment_count": len(segments),
                            "bond_start": start_bond,
                            "bond_end": end_bond,
                            "eligible_candidate_count": len(eligible),
                            **candidate,
                        }
                    )
    return pd.DataFrame(rows)


def plot_sequence_atlas_page(
    record: pd.Series,
    quantitative_df: pd.DataFrame,
    motifs_df: pd.DataFrame,
    ligand_df: pd.DataFrame,
    pdf_pages: PdfPages,
    page_number: int,
    page_count: int,
    segment: tuple[int, int],
    segment_number: int,
    segment_count: int,
    standalone_pdf_path: Path | None,
    standalone_png_path: Path | None,
    panel_label: str | None = None,
    footer_label: str = "Osteosarc vaccine cleavage atlas",
) -> None:
    sequence = record["sequence"]
    start_bond, end_bond = segment
    bonds = list(range(start_bond, end_bond + 1))
    record_id = record["sequence_record_id"]
    residue_start, residue_end = start_bond, end_bond + 1
    fig = plt.figure(figsize=(16, 10.5))
    ax = fig.add_axes([0.12, 0.15, 0.835, 0.71])
    ax.set_xlim(residue_start - 0.8, residue_end + 0.8)
    ax.set_ylim(-4.55, 5.95)
    ax.axis("off")

    # The sequence is the coordinate system, not merely another annotation.
    ax.add_patch(
        Rectangle(
            (residue_start - 0.47, -0.48),
            residue_end - residue_start + 0.94,
            0.96,
            facecolor="#f4f6f8",
            edgecolor="none",
            zorder=-1,
        )
    )
    for bond in bonds:
        ax.plot(
            [bond + 0.5, bond + 0.5],
            [-0.16, 0.16],
            color="#9da7b0",
            linewidth=0.55,
            alpha=0.65,
            zorder=2,
        )

    # Disclosed intended minimal epitope: gold behind the actual residue letters.
    minimal_bounds = _minimal_epitope_bounds(record)
    if minimal_bounds:
        left = max(residue_start, minimal_bounds[0])
        right = min(residue_end, minimal_bounds[1])
        if left <= right:
            ax.add_patch(
                Rectangle(
                    (left - 0.47, -0.58),
                    right - left + 0.94,
                    1.16,
                    facecolor="#fff0b3",
                    edgecolor="#c78c00",
                    linewidth=2.0,
                    zorder=0,
                )
            )

    # Large sequence strip. Bonds live at half-integer x coordinates.
    for position in range(residue_start, residue_end + 1):
        ax.text(
            position,
            0,
            sequence[position - 1],
            ha="center",
            va="center",
            family="monospace",
            fontsize=21,
            fontweight="bold",
            color="#18222d",
            zorder=5,
        )
        if position == residue_start or position == residue_end or position % 5 == 0:
            ax.text(
                position,
                -0.76,
                str(position),
                ha="center",
                va="top",
                fontsize=8,
                color="#56616c",
            )

    score_subset = quantitative_df.loc[
        (quantitative_df["sequence_record_id"] == record_id)
        & quantitative_df["assessable"]
        & quantitative_df["bond"].isin(bonds)
    ]
    model_colors = {
        "netchop-3.1-20s-3.0": "#0072b2",
        "pepsickle-in-vivo-human-only": "#009e73",
        "netchop-3.1-cterm-3.0": "#d55e00",
        "netcleave-i-hla": "#cc79a7",
        "netcleave-ii-hla": "#6f4aa8",
    }
    model_y = {
        model: 2.42 + index * 0.78 for index, model in enumerate(FIGURE_CYTOSOL_MODELS)
    }
    model_y["netcleave-ii-hla"] = -2.08
    model_scores = _score_matrix(
        record_id, quantitative_df, FIGURE_QUANTITATIVE_MODELS, bonds
    )
    for model_index, model in enumerate(FIGURE_QUANTITATIVE_MODELS):
        y = model_y[model]
        direction = 1 if model in FIGURE_CYTOSOL_MODELS else -1
        amplitude = 0.64 if direction > 0 else 0.56
        ax.plot(
            [residue_start - 0.45, residue_end + 0.45],
            [y, y],
            color="#c9d0d6",
            linewidth=0.7,
            zorder=0,
        )
        ax.plot(
            [residue_start - 0.45, residue_end + 0.45],
            [y + direction * amplitude * THRESHOLD] * 2,
            color=model_colors[model],
            linewidth=0.55,
            alpha=0.35,
            linestyle="--",
            zorder=0,
        )
        ax.text(
            residue_start - 0.68,
            y,
            f"{ATLAS_MODEL_LABELS[model]}  native 0-1",
            ha="right",
            va="center",
            fontsize=8.2,
            color=model_colors[model],
        )
        ax.text(
            residue_end + 0.55,
            y + direction * amplitude * THRESHOLD,
            "0.5 display",
            ha="left",
            va="center",
            fontsize=6.8,
            color=model_colors[model],
        )
        scores = model_scores[model_index]
        for x, values in continuous_score_runs(bonds, scores):
            endpoints = y + direction * amplitude * values
            ax.fill_between(
                x,
                y,
                endpoints,
                color=model_colors[model],
                alpha=0.12,
                linewidth=0,
                zorder=1,
            )
            ax.plot(
                x,
                endpoints,
                color=model_colors[model],
                linewidth=1.8,
                solid_joinstyle="round",
                solid_capstyle="round",
                zorder=2,
            )
            ax.scatter(
                x,
                endpoints,
                s=7,
                color=model_colors[model],
                alpha=0.48,
                linewidth=0,
                zorder=3,
            )
            strong = values >= THRESHOLD
            ax.scatter(
                x[strong],
                endpoints[strong],
                s=25,
                color=model_colors[model],
                edgecolor="white",
                linewidth=0.45,
                zorder=4,
            )

    ax.text(
        residue_start - 0.68,
        5.65,
        "INTRACELLULAR / CLASS-I PROCESSING - four complementary score tracks",
        ha="left",
        va="center",
        fontsize=8.6,
        fontweight="bold",
        color="#46515b",
    )
    ax.text(
        residue_end + 0.45,
        5.65,
        "higher native score upward",
        ha="right",
        va="center",
        fontsize=7.2,
        color="#66717b",
    )

    # A bond mark and faint guide connect conservative support in the four
    # distant traces to the exact inter-residue coordinate. Four beads encode
    # the vote without placing a text label on top of an MHC candidate bar.
    cut_support = conservative_cut_support(
        quantitative_df,
        record_id,
        bonds,
        FIGURE_CYTOSOL_MODELS,
        FIGURE_RED_SUPPORT_REQUIRED,
    )
    for bond, (hits, assessed) in cut_support.items():
        x = bond + 0.5
        ax.plot(
            [x, x],
            [0.48, 5.42],
            color="#b3212d",
            linewidth=0.65,
            alpha=0.18,
            linestyle=(0, (2, 3)),
            zorder=0,
        )
        ax.plot([x, x], [-0.47, 0.47], color="#b3212d", linewidth=2.5, zorder=6)
        bead_y = np.linspace(-0.27, 0.27, assessed)
        ax.scatter(
            [x] * assessed,
            bead_y,
            s=17,
            facecolors=[
                "#b3212d" if index < hits else "white" for index in range(assessed)
            ],
            edgecolors="#8f1520",
            linewidths=0.65,
            zorder=7,
        )

    # Native peptidase outputs without a validated common threshold.
    terminal_y = {"eramer-step": 1.98, "dpp4-qpisa": -3.18}
    terminal_label = {
        "eramer-step": "ERAP1 / ERAMER",
        "dpp4-qpisa": "DPP4 native (intact SLP)",
    }
    terminal_color = {"eramer-step": "#6a4492", "dpp4-qpisa": "#8a5a00"}
    for model, y in terminal_y.items():
        ax.text(
            residue_start - 0.68,
            y,
            terminal_label[model],
            ha="right",
            va="center",
            fontsize=8.2,
            color=terminal_color[model],
        )
        for result in score_subset.loc[score_subset["model"] == model].itertuples():
            x = int(result.bond) + 0.5
            ax.scatter(
                [x],
                [y],
                marker="D",
                s=38,
                color=terminal_color[model],
                edgecolor="#49323f" if model == "eramer-step" else "#6a4600",
            )
            ax.text(x + 0.08, y + 0.12, f"{float(result.score):.2g}", fontsize=6.5)

    def short_allele(value: str) -> str:
        value = value.replace("HLA-", "").replace("DRA1*01:01-", "")
        if "-" in value and value.startswith(("DPA1", "DQA1")):
            alpha, beta = value.split("-", 1)
            locus = "DP" if alpha.startswith("DPA1") else "DQ"
            return f"{locus}{alpha.split('*')[1]}/{beta.split('*')[1]}"
        if value.startswith("DRB1*"):
            return "DR" + value.split("*", 1)[1]
        return value.replace("*", "")

    def draw_ligands(mhc_class: str, lane_y: list[float], color: str) -> None:
        eligible_candidates = segment_ligand_candidates(
            ligand_df,
            record_id,
            mhc_class,
            residue_start,
            residue_end,
        )
        selected_candidates = select_ligand_candidates_for_display(eligible_candidates)
        if selected_candidates.empty:
            segment_suffix = " in this segment" if segment_count > 1 else ""
            ax.text(
                residue_start - 0.68,
                lane_y[0],
                f"MHC-{mhc_class}: none <= "
                f"{MHC_I_DISPLAY_RANK if mhc_class == 'I' else MHC_II_DISPLAY_RANK:g}% "
                f"rank{segment_suffix}",
                ha="right",
                va="center",
                fontsize=8,
                color="#68717a",
            )
            return
        ax.text(
            residue_start - 0.68,
            sum(lane_y) / len(lane_y),
            f"MHC-{mhc_class} candidate ligands",
            ha="right",
            va="center",
            fontsize=8.2,
            color=color,
        )
        drawn_spans: list[tuple[int, int]] = []
        for candidate in selected_candidates.itertuples(index=False):
            lane = int(candidate.display_lane)
            drawn_spans.append(
                (
                    max(int(candidate.start), residue_start),
                    min(int(candidate.end), residue_end),
                )
            )
            y = lane_y[lane]
            left = max(candidate.start - 0.44, residue_start - 0.55)
            right = min(candidate.end + 0.44, residue_end + 0.55)
            patch = FancyBboxPatch(
                (left, y - 0.14),
                right - left,
                0.28,
                boxstyle="round,pad=0.02,rounding_size=0.10",
                facecolor="white",
                edgecolor=color,
                linewidth=1.1,
                alpha=0.94,
                zorder=3,
            )
            ax.add_patch(patch)
            label = f"{short_allele(candidate.allele)} {candidate.percentile_rank:.2g}%"
            ax.text(
                (left + right) / 2,
                y,
                label,
                ha="center",
                va="center",
                fontsize=6.3,
                color=color,
                fontweight="bold",
                clip_on=True,
                zorder=4,
            )
            for item in str(candidate.internal_candidate_cleavage_bonds).split(";"):
                if not item or item == "nan":
                    continue
                bond = int(item.split("(", 1)[0])
                if residue_start <= bond <= residue_end:
                    ax.plot(
                        [bond + 0.5, bond + 0.5],
                        [y - 0.17, y + 0.17],
                        color="#b3212d",
                        linewidth=1.8,
                        zorder=4,
                    )
        band_y = 0.0 if mhc_class == "I" else -0.48
        for start, end in merge_residue_spans(drawn_spans):
            ax.add_patch(
                Rectangle(
                    (start - 0.47, band_y),
                    end - start + 0.94,
                    0.48,
                    facecolor=color,
                    edgecolor="none",
                    alpha=0.10,
                    zorder=1,
                )
            )
        omitted = len(eligible_candidates) - len(selected_candidates)
        if omitted > 0:
            omitted_y = max(lane_y) + 0.34 if mhc_class == "I" else min(lane_y) - 0.12
            ax.text(
                residue_end + 0.45,
                omitted_y,
                f"+{omitted} more in CSV",
                ha="right",
                va="center",
                fontsize=6.6,
                color=color,
            )

    draw_ligands("I", [0.68, 1.03, 1.38], "#2878b5")
    draw_ligands("II", [-0.72, -1.07, -1.42], "#5b3d91")

    ax.text(
        residue_start - 0.68,
        -1.78,
        "ENDOLYSOSOMAL / CLASS-II PROCESSING - C-terminal score track",
        ha="left",
        va="center",
        fontsize=8.3,
        fontweight="bold",
        color="#46515b",
    )
    ax.text(
        residue_end + 0.45,
        -1.78,
        "higher native score downward",
        ha="right",
        va="center",
        fontsize=7.2,
        color="#66717b",
    )

    motif_abbreviation = {
        model: MOTIF_DISPLAY_NAMES[model].split(" - ", 1)[0]
        for model in MOTIF_DISPLAY_NAMES
    }
    matched = motifs_df.loc[
        (motifs_df["sequence_record_id"] == record_id)
        & (motifs_df["status"] == "matched")
        & motifs_df["bond"].isin(bonds)
    ]
    for models, y, label, color in (
        (INTRACELLULAR_ER_MOTIF_MODELS, 1.66, "cytosol / ER motifs", "#6a4492"),
        (EXTRACELLULAR_MOTIF_MODELS, -3.82, "serum / extracellular motifs", "#007b83"),
    ):
        ax.text(
            residue_start - 0.68,
            y,
            label,
            ha="right",
            va="center",
            fontsize=8.2,
            color=color,
        )
        for bond, group in matched.loc[matched["model"].isin(models)].groupby("bond"):
            names = "+".join(
                motif_abbreviation[model] for model in sorted(group["model"].unique())
            )
            x = int(bond) + 0.5
            marker = "^" if y > 0 else "v"
            label_y = y + 0.16 if y > 0 else y - 0.16
            va = "bottom" if y > 0 else "top"
            ax.scatter([x], [y], marker=marker, s=38, color=color, zorder=3)
            ax.text(
                x,
                label_y,
                names,
                rotation=45,
                ha="right",
                va=va,
                fontsize=6.1,
                color=color,
            )

    ax.text(
        residue_start - 0.68,
        -2.78,
        "SERUM / EXTRACELLULAR - native scores and recognition motifs",
        ha="left",
        va="center",
        fontsize=8.6,
        fontweight="bold",
        color="#46515b",
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
        fontsize=16,
        fontweight="bold",
        y=0.965,
    )
    if panel_label:
        fig.text(
            0.018,
            0.965,
            panel_label,
            ha="left",
            va="top",
            fontsize=22,
            fontweight="bold",
            color="#18222d",
        )
    epitope_label = (
        f"disclosed minimal epitope {record['minimal_epitope']} at "
        f"{minimal_bounds[0]}-{minimal_bounds[1]}"
        if minimal_bounds
        else "minimal epitope/window not disclosed for this SLP"
    )
    fig.text(
        0.5,
        0.913,
        epitope_label,
        ha="center",
        va="center",
        fontsize=9.5,
        color="#8a6300" if minimal_bounds else "#626b73",
    )
    fig.text(
        0.5,
        0.086,
        "SCORES  Native 0-1 bond scores; lines join adjacent assessed bonds only (no smoothing). Dashed line = 0.5 display threshold; gaps = unassessed.",
        ha="center",
        va="center",
        fontsize=8.8,
        color="#333333",
    )
    fig.text(
        0.5,
        0.057,
        "RED BOND MARKS  All four intracellular tracks assess the bond; >=3 reach 0.5. Four beads show support (filled = hit, open = miss), not probability.",
        ha="center",
        va="center",
        fontsize=8.8,
        color="#333333",
    )
    fig.text(
        0.5,
        0.028,
        "MHC WINDOWS  <=2%/<=5% rank; display favors intended overlap, then allele diversity, then rank. Internal red ticks are pre-binding cut evidence; all predictions remain in CSV.",
        ha="center",
        va="center",
        fontsize=8.8,
        color="#333333",
    )
    standalone_title = (
        f"{record['gene']} {record['protein_change'] or ''} vaccine SLP cleavage map"
    )
    save_pdf_page(
        fig,
        pdf_pages,
        page_number,
        page_count,
        standalone_pdf_path=standalone_pdf_path,
        standalone_png_path=standalone_png_path,
        standalone_title=standalone_title,
        footer_label=footer_label,
    )


def render_figures(
    output_dir: Path,
    records: pd.DataFrame,
    quantitative_df: pd.DataFrame,
    motifs_df: pd.DataFrame,
    ligand_df: pd.DataFrame,
    summary_df: pd.DataFrame,
    generated_at: datetime,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Render the clustered overview and bond-aligned sequence atlas."""
    figures_dir = output_dir / "figures"
    figures_dir.mkdir(parents=True, exist_ok=True)
    map_dir = figures_dir / "slp-maps"
    map_dir.mkdir(parents=True, exist_ok=True)
    page_count = 1 + sum(
        len(sequence_segments(record["sequence"])) for _, record in records.iterrows()
    )
    # Matplotlib 3.11 mis-encodes negative UTC offsets in PDF date strings
    # (for example, -04:00 becomes -20:00). Preserve the correct local clock
    # time as a timezone-naive PDF date; provenance.json retains the exact
    # timezone-aware generation timestamp.
    pdf_generated_at = generated_at.replace(tzinfo=None)
    pdf_metadata = {
        "Title": "Osteosarc vaccine SLP bond-level cleavage atlas",
        "Author": "mhctools",
        "Subject": "Sequence-aligned cleavage scores and peptidase motif evidence",
        "Keywords": "osteosarcoma vaccine SLP cleavage proteasome peptidase",
        "CreationDate": pdf_generated_at,
        "ModDate": pdf_generated_at,
    }
    atlas_rows: list[dict[str, Any]] = []
    map_rows: list[dict[str, Any]] = []
    with PdfPages(output_dir / PDF_FILENAME, metadata=pdf_metadata) as pdf_pages:
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
                map_stem = standalone_map_stem(
                    atlas_order,
                    record_id,
                    segment_number,
                    len(segments),
                )
                standalone_pdf_path = map_dir / f"{map_stem}.pdf"
                standalone_png_path = map_dir / f"{map_stem}.png"
                plot_sequence_atlas_page(
                    record,
                    quantitative_df,
                    motifs_df,
                    ligand_df,
                    pdf_pages,
                    page_number,
                    page_count,
                    segment,
                    segment_number,
                    len(segments),
                    standalone_pdf_path,
                    standalone_png_path,
                )
                map_rows.append(
                    {
                        "atlas_order": atlas_order,
                        "sequence_record_id": record_id,
                        "gene": record["gene"],
                        "protein_change": record["protein_change"],
                        "segment_number": segment_number,
                        "segment_count": len(segments),
                        "bond_start": segment[0],
                        "bond_end": segment[1],
                        "atlas_pdf_page": page_number,
                        "pdf_file": str(standalone_pdf_path.relative_to(output_dir)),
                        "png_file": str(standalone_png_path.relative_to(output_dir)),
                    }
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
    return atlas_order_df, pd.DataFrame(map_rows)


def select_manuscript_records(
    records: pd.DataFrame,
    ligand_df: pd.DataFrame,
    vulnerable_df: pd.DataFrame,
) -> pd.DataFrame:
    """Choose four complementary, reproducible examples for a compact figure set."""
    metrics: list[dict[str, Any]] = []
    for record in records.itertuples(index=False):
        if len(sequence_segments(record.sequence)) != 1:
            continue
        ligands = ligand_df.loc[
            (ligand_df["sequence_record_id"] == record.sequence_record_id)
            & ligand_df["display_candidate"]
        ]
        vulnerable = vulnerable_df.loc[
            (vulnerable_df["sequence_record_id"] == record.sequence_record_id)
            & vulnerable_df["biological_context"].str.startswith(
                "cytosolic/proteasome", na=False
            )
        ]
        internal = vulnerable.loc[
            vulnerable["in_disclosed_minimal_epitope"].fillna(False).astype(bool)
        ]
        class_i = ligands.loc[ligands["mhc_class"] == "I"]
        class_ii = ligands.loc[ligands["mhc_class"] == "II"]
        intended = ligands.loc[
            ligands["overlaps_disclosed_minimal_epitope"].fillna(False).astype(bool)
        ]
        metrics.append(
            {
                "sequence_record_id": record.sequence_record_id,
                "gene": record.gene,
                "protein_change": record.protein_change,
                "length": record.length,
                "has_disclosed_minimal_epitope": pd.notna(record.minimal_epitope)
                and bool(str(record.minimal_epitope)),
                "intended_epitope_internal_conservative_cuts": len(internal),
                "intended_epitope_overlapping_mhc_candidates": len(intended),
                "mhc_i_candidate_count": len(class_i),
                "mhc_i_allele_count": class_i["allele"].nunique(),
                "mhc_ii_candidate_count": len(class_ii),
                "mhc_ii_allele_count": class_ii["allele"].nunique(),
            }
        )
    candidates = pd.DataFrame(metrics)
    if len(candidates) < 4:
        raise RuntimeError("At least four single-page SLP records are required")

    selected: list[dict[str, Any]] = []
    used: set[str] = set()
    used_genes: set[str] = set()

    def add_best(pool: pd.DataFrame, columns: list[str], reason: str) -> None:
        pool = pool.loc[
            ~pool["sequence_record_id"].isin(used) & ~pool["gene"].isin(used_genes)
        ]
        if pool.empty:
            return
        ordered = pool.sort_values(
            columns + ["sequence_record_id"],
            ascending=[False] * len(columns) + [True],
            kind="stable",
        )
        row = ordered.iloc[0].to_dict()
        row["selection_reason"] = reason
        used.add(str(row["sequence_record_id"]))
        used_genes.add(str(row["gene"]))
        selected.append(row)

    with_minimal = candidates.loc[candidates["has_disclosed_minimal_epitope"]]
    add_best(
        with_minimal,
        [
            "intended_epitope_internal_conservative_cuts",
            "intended_epitope_overlapping_mhc_candidates",
        ],
        "most conservative cut evidence inside a disclosed intended epitope",
    )
    add_best(
        with_minimal.loc[
            with_minimal["intended_epitope_internal_conservative_cuts"] == 0
        ],
        ["intended_epitope_overlapping_mhc_candidates"],
        "no conservative internal cut and most MHC windows overlapping the intended epitope",
    )
    add_best(
        with_minimal,
        ["mhc_i_candidate_count", "mhc_i_allele_count"],
        "most eligible class-I windows, then greatest class-I allele diversity",
    )
    add_best(
        with_minimal,
        ["mhc_ii_candidate_count", "mhc_ii_allele_count"],
        "most eligible class-II windows, then greatest class-II allele diversity",
    )
    while len(selected) < 4:
        before = len(selected)
        add_best(
            candidates,
            ["intended_epitope_overlapping_mhc_candidates"],
            "highest remaining intended-epitope MHC-window coverage",
        )
        if len(selected) == before:
            break
    if len(selected) != 4:
        raise RuntimeError(f"Selected {len(selected)} manuscript records, expected 4")
    result = pd.DataFrame(selected)
    result.insert(0, "panel", list("ABCD"))
    return result


def render_manuscript_figures(
    output_dir: Path,
    records: pd.DataFrame,
    quantitative_df: pd.DataFrame,
    motifs_df: pd.DataFrame,
    ligand_df: pd.DataFrame,
    selection_df: pd.DataFrame,
    generated_at: datetime,
) -> None:
    """Render the objective four-example manuscript subset as full-size panels."""
    metadata = {
        "Title": "Osteosarc vaccine SLP cleavage manuscript figure set",
        "Author": "mhctools",
        "Subject": "Four complementary sequence-aligned cleavage maps",
        "Keywords": "osteosarcoma vaccine SLP cleavage MHC manuscript",
        "CreationDate": generated_at.replace(tzinfo=None),
        "ModDate": generated_at.replace(tzinfo=None),
    }
    records_by_id = records.set_index("sequence_record_id", drop=False)
    with PdfPages(output_dir / MANUSCRIPT_PDF_FILENAME, metadata=metadata) as pdf_pages:
        for page_number, selected in enumerate(
            selection_df.itertuples(index=False), start=1
        ):
            record = records_by_id.loc[selected.sequence_record_id]
            segments = sequence_segments(record["sequence"])
            if len(segments) != 1:
                raise RuntimeError("Manuscript selection must fit on one full-size page")
            plot_sequence_atlas_page(
                record,
                quantitative_df,
                motifs_df,
                ligand_df,
                pdf_pages,
                page_number,
                len(selection_df),
                segments[0],
                1,
                1,
                None,
                None,
                panel_label=selected.panel,
                footer_label="Osteosarc manuscript figure set",
            )


def write_manuscript_caption(path: Path, selection_df: pd.DataFrame) -> None:
    """Write a self-contained caption and the objective panel-selection audit."""
    lines = [
        "# Manuscript figure caption",
        "",
        "**Figure. Sequence-aligned cleavage and MHC-ligand evidence for four disclosed "
        "osteosarcoma vaccine SLPs.** Each full-size panel centers the vaccine sequence; "
        "gold marks the disclosed intended minimal epitope. Intracellular/class-I "
        "processing evidence is above the sequence, and endolysosomal/class-II plus "
        "serum/extracellular evidence is below it. Lines connect adjacent exact native "
        "bond scores without smoothing. Red sequence marks require all four displayed "
        "intracellular models to assess a bond and at least three to reach the 0.5 display "
        "threshold; filled/open beads show the support count and are not probabilities. "
        "Blue and purple bars are selected MHC-I and MHC-II ligand windows at <=2% and "
        "<=5% native percentile rank, respectively. Their selection favors intended-"
        "epitope overlap, then distinct alleles, then native rank; all predictions remain "
        "in the accompanying CSV.",
        "",
        "Panels were chosen by declared complementary criteria, with different genes and "
        "a disclosed intended epitope on every page:",
        "",
    ]
    for row in selection_df.itertuples(index=False):
        lines.append(
            f"- **{row.panel}, {row.gene} {row.protein_change}.** {row.selection_reason.capitalize()}. "
            f"Audit metrics: {row.intended_epitope_internal_conservative_cuts} conservative "
            "internal intended-epitope cut site(s), "
            f"{row.intended_epitope_overlapping_mhc_candidates} eligible MHC window(s) "
            "overlapping the intended epitope, "
            f"{row.mhc_i_candidate_count} class-I and {row.mhc_ii_candidate_count} class-II "
            "eligible windows overall."
        )
    lines.extend(
        [
            "",
            "The 0.5 score cutoff is a common visualization threshold, not a calibrated "
            "cleavage probability. Predictors have different training data and biological "
            "scope; agreement does not establish in-vivo degradation, uptake, presentation, "
            "or MHC protection.",
            "",
        ]
    )
    path.write_text("\n".join(lines), encoding="utf-8")


def correlations(quantitative_df: pd.DataFrame) -> pd.DataFrame:
    thresholded = quantitative_df.loc[
        quantitative_df["display_threshold"].notna() & quantitative_df["assessable"]
    ].copy()
    wide = thresholded.pivot_table(
        index=["sequence_record_id", "bond"], columns="model", values="score"
    )
    return wide.corr(method="spearman", min_periods=10)


def validate_hla_inputs(hla_path: Path, variants_path: Path) -> None:
    """Require every scored allele/combination to be disclosed by the source."""
    hla = pd.read_csv(hla_path, sep="\t")
    normal = hla.loc[hla["sample"] == "normal"]
    disclosed = {f"HLA-{allele}" for allele in normal["allele"]}
    for allele in MHC_I_ALLELES:
        if allele not in disclosed:
            raise ValueError(
                f"Class-I prediction allele not found in source HLA table: {allele}"
            )

    variants = json.loads(variants_path.read_text(encoding="utf-8"))
    source_pairs = {
        value
        for variant in variants
        for value in [(variant.get("peptides") or {}).get("pvac25_best_allele")]
        if value and ("DPA1" in value or "DQA1" in value or "DRB1" in value)
    }
    for allele in MHC_II_ALLELES:
        source_name = allele.removeprefix("HLA-")
        if source_name not in source_pairs:
            raise ValueError(
                "Class-II prediction combination is not named in source candidate fields: "
                f"{allele}"
            )


def model_file_inventory(
    mhcflurry_metadata: dict[str, Any], netmhciipan_metadata: dict[str, Any]
) -> pd.DataFrame:
    """Inventory every file under the two local MHC model installations."""
    roots = {
        "mhcflurry": Path(mhcflurry_metadata["models_path"]),
        "netmhciipan": Path(netmhciipan_metadata["installation_root"]),
    }
    rows: list[dict[str, Any]] = []
    for model, root in roots.items():
        for path in sorted(root.rglob("*")):
            if path.is_file():
                rows.append(
                    {
                        "model": model,
                        "relative_path": str(path.relative_to(root)),
                        "size_bytes": path.stat().st_size,
                        "sha256": sha256_file(path),
                    }
                )
    return pd.DataFrame(rows)


def inventory_digest(inventory_df: pd.DataFrame, model: str) -> str:
    digest = hashlib.sha256()
    subset = inventory_df.loc[inventory_df["model"] == model].sort_values(
        "relative_path"
    )
    for row in subset.itertuples(index=False):
        digest.update(f"{row.relative_path}\0{row.size_bytes}\0{row.sha256}\n".encode())
    return digest.hexdigest()


def model_hashes(
    netchop_dir: Path,
    netcleave_dir: Path,
    eramer_dir: Path,
    mhcflurry_metadata: dict[str, Any],
    netmhciipan_metadata: dict[str, Any],
    mhc_inventory_df: pd.DataFrame,
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
        "mhcflurry": {
            **mhcflurry_metadata,
            "inventory_file": "tables/mhc_model_file_inventory.csv",
            "file_count": int((mhc_inventory_df["model"] == "mhcflurry").sum()),
            "inventory_sha256": inventory_digest(mhc_inventory_df, "mhcflurry"),
        },
        "netmhciipan": {
            **netmhciipan_metadata,
            "inventory_file": "tables/mhc_model_file_inventory.csv",
            "file_count": int((mhc_inventory_df["model"] == "netmhciipan").sum()),
            "inventory_sha256": inventory_digest(mhc_inventory_df, "netmhciipan"),
        },
    }


def write_report(
    path: Path,
    inventory_df: pd.DataFrame,
    records: pd.DataFrame,
    missing_df: pd.DataFrame,
    conflicts_df: pd.DataFrame,
    summary_df: pd.DataFrame,
    ligand_df: pd.DataFrame,
    vulnerable_df: pd.DataFrame,
    source_commit: str,
    source_date: str,
    generated_at: datetime,
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
        f"Generated locally at {generated_at.isoformat()}.",
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
        "The exact [MHC ligand-window predictions](tables/slp_mhc_ligand_predictions.csv), "
        "[context-separated vulnerable bonds](tables/slp_vulnerable_bonds.csv), and full "
        "[MHC model file inventory](tables/mhc_model_file_inventory.csv) are also retained.",
        "",
        "## Sequence cleavage atlas",
        "",
        f"The [complete PDF atlas]({PDF_FILENAME}) makes each amino-acid sequence the central axis, "
        "with continuous exact-score profiles, motif flags, disclosed minimal-epitope spans, and candidate "
        "MHC ligand windows hugging the sequence at exact residues and bonds. Adjacent assessed scores are joined for readability "
        "without smoothing; unassessed gaps remain open. It contains one clustered overview followed by one page "
        "per SLP (with the 80-aa outlier split across three continuation pages). "
        "[Atlas order and PDF page numbers](tables/atlas_sequence_order.csv) are provided for navigation. "
        "Every map page is also available as a vector PDF and 300 dpi PNG, indexed in "
        "[the individual-map export table](tables/slp_map_exports.csv).",
        f"A compact [four-panel manuscript figure set]({MANUSCRIPT_PDF_FILENAME}) retains "
        "full-size readable maps rather than shrinking them into a crowded contact sheet. "
        "Its examples are selected by four declared criteria—internal intended-epitope cut evidence, "
        "apparent preservation with overlapping MHC windows, class-I density, and class-II density—"
        "recorded in [the manuscript selection table](tables/manuscript_figure_selection.csv). "
        "A [self-contained caption](MANUSCRIPT_CAPTION.md) defines every visual encoding and "
        "lists the panel-specific selection audit.",
        "The atlas shows four intracellular quantitative tracks: human-only Pepsickle, NetChop Cterm, NetChop 20S, "
        "and NetCleave-I. Cterm is ligand-trained and emphasizes candidate MHC-I boundaries, whereas 20S is retained "
        "as a distinct in-vitro proteasome view rather than a substitute for Cterm. The human-only Pepsickle model is "
        "species-matched but experimental and trained on less data than its all-mammal counterpart. The near-redundant "
        "all-mammal Pepsickle output remains available in the exact-score and summary tables.",
        "",
        "![Predictor agreement and clustered SLP order](figures/predictor_agreement_and_slp_clusters.png)",
        "",
        "The overview clusters predictors using Spearman correlation of native scores at shared, "
        "assessable bonds. It clusters SLPs using standardized within-model fractions above the "
        "0.5 display threshold. Clustering is organizational only and is not an ensemble model.",
        "",
        "## MHC ligand and integrity overlays",
        "",
        f"Local inference produced {len(ligand_df):,} peptide/allele predictions. The atlas displays "
        f"MHC-I windows at <= {MHC_I_DISPLAY_RANK:g}% MHCflurry presentation rank and MHC-II windows "
        f"at <= {MHC_II_DISPLAY_RANK:g}% NetMHCIIpan EL rank. The table retains every prediction, "
        "including those outside the display cutoffs.",
        "MHC-I predictions use all five disclosed classical class-I alleles. MHC-II predictions use "
        "only alpha/beta combinations already named by the osteosarc source; the unphased HLA table "
        "is not used to invent additional combinations. All inference ran locally.",
        "Within each map segment and class, display selection first retains the strongest window "
        "overlapping a disclosed intended epitope, then tries to represent distinct alleles, and "
        "finally fills remaining non-overlapping lane capacity by native percentile rank. Exact "
        "selected rows and reasons are recorded in `slp_mhc_display_selection.csv`; every raw "
        "prediction remains in `slp_mhc_ligand_predictions.csv`.",
        "For class I, a red sequence mark or ligand-window tick requires all four displayed intracellular "
        "tracks to assess the bond and at least three to reach the 0.5 display threshold. Four beads on the sequence "
        "mark encode that support count (filled = hit, open = miss), not a probability. Red ticks inside a ligand bar are pre-binding internal cleavage evidence in the relevant "
        "processing view. They do not establish that a bound pMHC complex will be cleaved or protected; "
        "binding occupancy and timing are not modeled.",
        "",
        "## Conservative vulnerable bonds",
        "",
        f"The vulnerable-bond table contains {len(vulnerable_df)} context-separated rows. A proteasome/MHC-I "
        "row requires at least three of five native 0-1 models at the 0.5 display threshold. A motif row "
        "requires at least two distinct recognition rules in the same broad biological context. These are "
        "support counts, not calibrated probabilities, and unlike contexts are never combined.",
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
    parser.add_argument("--netmhciipan-path", type=Path, required=True)
    parser.add_argument(
        "--hla-table",
        type=Path,
        help="Defaults to scripts/dragen/tables/hla.tsv inside --osteosarc-repo",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path(__file__).resolve().parent / "results",
        help="Base directory; each run creates a date/time-stamped child directory",
    )
    return parser.parse_args()


def timestamped_output_dir(base_dir: Path, generated_at: datetime) -> Path:
    """Return a collision-resistant, filesystem-safe directory for one run."""
    if generated_at.tzinfo is None or generated_at.utcoffset() is None:
        raise ValueError("generated_at must include a timezone")
    stamp = generated_at.strftime("%Y-%m-%dT%H%M%S-%f%z")
    return base_dir / stamp


def main() -> None:
    args = parse_args()
    generated_at = datetime.now().astimezone()
    output_dir = timestamped_output_dir(args.output_dir.resolve(), generated_at)
    tables_dir = output_dir / "tables"
    figures_dir = output_dir / "figures"
    tables_dir.mkdir(parents=True, exist_ok=False)
    figures_dir.mkdir(parents=True, exist_ok=False)

    variants_path = args.osteosarc_repo / "src/data/variants.json"
    hla_path = args.hla_table or args.osteosarc_repo / "scripts/dragen/tables/hla.tsv"
    source_commit = git_value(args.osteosarc_repo, "%H")
    source_date = git_value(args.osteosarc_repo, "%cI")
    validate_hla_inputs(hla_path, variants_path)
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

    class_i_rows, mhcflurry_metadata = mhcflurry_ligand_rows(records)
    class_ii_rows, netmhciipan_metadata = netmhciipan_ligand_rows(
        records, args.netmhciipan_path
    )
    ligand_df = annotate_ligand_cleavage_exposure(
        pd.DataFrame(class_i_rows + class_ii_rows), quantitative_df
    )
    ligand_df.to_csv(tables_dir / "slp_mhc_ligand_predictions.csv", index=False)
    vulnerable_df = vulnerable_bond_table(records, quantitative_df, motifs_df)
    vulnerable_df.to_csv(tables_dir / "slp_vulnerable_bonds.csv", index=False)
    display_selection_df = mhc_display_selections(records, ligand_df)
    display_selection_df.to_csv(
        tables_dir / "slp_mhc_display_selection.csv", index=False
    )
    manuscript_selection_df = select_manuscript_records(
        records, ligand_df, vulnerable_df
    )
    manuscript_selection_df.to_csv(
        tables_dir / "manuscript_figure_selection.csv", index=False
    )
    mhc_inventory_df = model_file_inventory(mhcflurry_metadata, netmhciipan_metadata)
    mhc_inventory_df.to_csv(tables_dir / "mhc_model_file_inventory.csv", index=False)

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

    atlas_order_df, map_exports_df = render_figures(
        output_dir,
        records,
        quantitative_df,
        motifs_df,
        ligand_df,
        summary_df,
        generated_at,
    )
    atlas_order_df.to_csv(tables_dir / "atlas_sequence_order.csv", index=False)
    map_exports_df.to_csv(tables_dir / "slp_map_exports.csv", index=False)
    render_manuscript_figures(
        output_dir,
        records,
        quantitative_df,
        motifs_df,
        ligand_df,
        manuscript_selection_df,
        generated_at,
    )
    write_manuscript_caption(
        output_dir / "MANUSCRIPT_CAPTION.md", manuscript_selection_df
    )
    write_report(
        output_dir / "REPORT.md",
        inventory_df,
        records,
        missing_df,
        conflicts_df,
        summary_df,
        ligand_df,
        vulnerable_df,
        source_commit,
        source_date,
        generated_at,
    )

    provenance = {
        "schema_version": 1,
        "source": {
            "site": "https://osteosarc.com/",
            "repository": "https://gitlab.com/slowkow/osteosarc.com",
            "commit": source_commit,
            "commit_date": source_date,
            "variants_json_sha256": sha256_file(variants_path),
            "hla_table_sha256": sha256_file(hla_path),
            "mhc_i_alleles": MHC_I_ALLELES,
            "mhc_ii_alleles": MHC_II_ALLELES,
            "mhc_ii_pairing_policy": (
                "only combinations already named in osteosarc source candidate fields"
            ),
        },
        "analysis": {
            "mhctools_version": __import__("mhctools").__version__,
            "mhctools_commit": git_value(Path(__file__).resolve().parents[2], "%H"),
            "mhctools_worktree_dirty": subprocess.run(
                ["git", "diff", "--quiet", "HEAD", "--"],
                cwd=Path(__file__).resolve().parents[2],
                check=False,
            ).returncode
            != 0,
            "analysis_script_sha256": sha256_file(Path(__file__)),
            "generated_at": generated_at.isoformat(),
            "output_directory": output_dir.name,
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
        "models": model_hashes(
            args.netchop_dir,
            args.netcleave_dir,
            args.eramer_dir,
            mhcflurry_metadata,
            netmhciipan_metadata,
            mhc_inventory_df,
        ),
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
    print(output_dir)


if __name__ == "__main__":
    main()
