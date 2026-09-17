#!/usr/bin/env python3
"""Re-render a timestamped atlas from a frozen, checksummed prediction run."""

import argparse
from copy import deepcopy
from datetime import datetime
import importlib.metadata
import json
from pathlib import Path
import platform
import shutil
import subprocess

import pandas as pd

from analyze import (
    MANUSCRIPT_PDF_FILENAME,
    PDF_FILENAME,
    git_value,
    mhc_display_selections,
    render_figures,
    render_manuscript_figures,
    select_manuscript_records,
    sha256_file,
    slp_records,
    timestamped_output_dir,
    write_manuscript_caption,
    write_report,
)


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--source-run", type=Path, required=True,
        help="Prior timestamped run whose exact prediction tables are reused",
    )
    parser.add_argument(
        "--output-dir", type=Path,
        default=Path(__file__).resolve().parent / "results",
        help="Base directory for the new date/time-stamped run",
    )
    return parser.parse_args()


def _verified_source_manifest(source_run):
    manifest_path = source_run / "SHA256SUMS.json"
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    failures = []
    for relative, expected in manifest.items():
        path = source_run / relative
        if not path.is_file() or sha256_file(path) != expected:
            failures.append(relative)
    if failures:
        raise RuntimeError(
            "Source run checksum verification failed: %s" % ", ".join(failures)
        )
    return manifest


def main():
    args = parse_args()
    source_run = args.source_run.resolve()
    source_manifest = _verified_source_manifest(source_run)
    source_provenance = json.loads(
        (source_run / "provenance.json").read_text(encoding="utf-8")
    )
    generated_at = datetime.now().astimezone()
    output_dir = timestamped_output_dir(args.output_dir.resolve(), generated_at)
    output_dir.mkdir(parents=True, exist_ok=False)
    shutil.copytree(source_run / "tables", output_dir / "tables")
    tables = output_dir / "tables"

    inventory_df = pd.read_csv(tables / "vaccine_sequence_inventory.csv")
    missing_df = pd.read_csv(tables / "undisclosed_vaccine_sequences.csv")
    conflicts_df = pd.read_csv(tables / "sequence_provenance_conflicts.csv")
    quantitative_df = pd.read_csv(tables / "slp_quantitative_bond_scores.csv")
    motifs_df = pd.read_csv(tables / "slp_motif_assessments.csv")
    summary_df = pd.read_csv(tables / "slp_model_summary.csv")
    ligand_df = pd.read_csv(tables / "slp_mhc_ligand_predictions.csv")
    vulnerable_df = pd.read_csv(tables / "slp_vulnerable_bonds.csv")
    records = slp_records(inventory_df)

    display_selection_df = mhc_display_selections(records, ligand_df)
    display_selection_df.to_csv(
        tables / "slp_mhc_display_selection.csv", index=False
    )
    manuscript_selection_df = select_manuscript_records(
        records, ligand_df, vulnerable_df
    )
    manuscript_selection_df.to_csv(
        tables / "manuscript_figure_selection.csv", index=False
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
    atlas_order_df.to_csv(tables / "atlas_sequence_order.csv", index=False)
    map_exports_df.to_csv(tables / "slp_map_exports.csv", index=False)
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
        source_provenance["source"]["commit"],
        source_provenance["source"]["commit_date"],
        generated_at,
    )

    provenance = deepcopy(source_provenance)
    source_analysis = deepcopy(source_provenance["analysis"])
    repo = Path(__file__).resolve().parents[2]
    provenance["analysis"].update({
        "mhctools_version": __import__("mhctools").__version__,
        "mhctools_commit": git_value(repo, "%H"),
        "mhctools_worktree_dirty": subprocess.run(
            ["git", "diff", "--quiet", "HEAD", "--"], cwd=repo, check=False
        ).returncode != 0,
        "analysis_script_sha256": sha256_file(Path(__file__).with_name("analyze.py")),
        "rerender_script_sha256": sha256_file(Path(__file__)),
        "generated_at": generated_at.isoformat(),
        "output_directory": output_dir.name,
        "prediction_source_run": source_run.name,
        "prediction_source_manifest_sha256": sha256_file(
            source_run / "SHA256SUMS.json"
        ),
        "prediction_source_files_verified": len(source_manifest),
        "prediction_execution": "not repeated; verified frozen tables reused",
        "prediction_source_analysis": source_analysis,
    })
    provenance["runtime"] = {
        "python": platform.python_version(),
        "packages": {
            name: importlib.metadata.version(name)
            for name in ("matplotlib", "numpy", "pandas", "scipy")
        },
    }
    (output_dir / "provenance.json").write_text(
        json.dumps(provenance, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )

    manifest = {}
    for path in sorted(output_dir.rglob("*")):
        if path.is_file() and path.name != "SHA256SUMS.json":
            manifest[str(path.relative_to(output_dir))] = sha256_file(path)
    (output_dir / "SHA256SUMS.json").write_text(
        json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    if not (output_dir / PDF_FILENAME).is_file():
        raise RuntimeError("Missing %s" % PDF_FILENAME)
    if not (output_dir / MANUSCRIPT_PDF_FILENAME).is_file():
        raise RuntimeError("Missing %s" % MANUSCRIPT_PDF_FILENAME)
    print(output_dir)


if __name__ == "__main__":
    main()
