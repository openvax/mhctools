#!/usr/bin/env python3
"""Run native enzyme-specific evidence without inferring serum survival."""

import argparse
from datetime import datetime
import hashlib
import json
from pathlib import Path
import shutil
import subprocess
import sys

from mechanism_scorecards import (
    core_bounds, enriched_assessments, mechanism_catalog, read_csv, sha256,
    verified_source, write_csv,
)
from mhctools import CleaveNet, CleavageInput, PeptideInput, Pepsickle
from mhctools.cleavenet import ENZYMES


SERUM_MODELS = (
    "dpp4-qpisa", "ace-dipeptidyl", "cpn-basic", "cpb2-basic", "app2-xp",
    "fap-endo-gp", "fap-dipeptidyl", "mme-hydrophobic", "anpep-ala", "enpep-acidic",
)


def vaccine_records(source):
    """Keep disclosed intact SLP occurrences, including identical sequences."""
    records = [r for r in read_csv(source / "tables/vaccine_sequence_inventory.csv")
               if r["sequence_type"] == "synthetic_long_peptide"]
    ids = [r["sequence_record_id"] for r in records]
    if not records or len(set(ids)) != len(ids):
        raise ValueError("Expected distinct vaccine records")
    for record in records:
        sequence = record["sequence"]
        if (int(record["length"]) != len(sequence) or
                hashlib.sha256(sequence.encode()).hexdigest() != record["sequence_sha256"]):
            raise ValueError("Source sequence identity mismatch")
        core_bounds(record)
    return records


def window_inputs(records):
    """Every overlapping ten-mer; coordinates never designate a cut."""
    return [PeptideInput(
        record["sequence"][start:start + 10],
        occurrence_id=record["sequence_record_id"] + ":window:" + str(start),
        source_sequence_name=record["sequence_record_id"], source_start=start)
        for record in records for start in range(len(record["sequence"]) - 9)]


def window_rows(records, inputs, results):
    """Preserve exact native means/spreads and source-window/core overlap."""
    parents = {r["sequence_record_id"]: r for r in records}
    if len(inputs) != len(results):
        raise ValueError("Missing CleaveNet window results")
    rows = []
    for item, result in zip(inputs, results):
        if result.peptide_input != item or tuple(s.enzyme for s in result.scores) != ENZYMES:
            raise ValueError("CleaveNet occurrence or enzyme identity mismatch")
        record = parents[item.source_sequence_name]
        bounds = core_bounds(record)
        end = item.source_start + len(item.sequence)
        overlap = (max(0, min(end, bounds[1]) - max(item.source_start, bounds[0]))
                   if bounds else "")
        for score in result.scores:
            rows.append(dict(
                sequence_record_id=item.source_sequence_name, gene=record["gene"],
                sequence_sha256=record["sequence_sha256"],
                window_id=item.occurrence_id, window_start=item.source_start,
                window_end=end, window_sequence=item.sequence, enzyme=score.enzyme,
                z_score=score.z_score, ensemble_sd=score.ensemble_sd,
                source_core_overlap_residues=overlap,
                endpoint="whole-substrate relative cleavage Z-score",
                interpretation="conditional ten-mer substrate; exact cut unknown",
                cache_key=result.cache_key))
    return rows


def bond_rows(records, results):
    """Validate canonical sites and annotate strict source-core interior."""
    parents = {r["sequence_record_id"]: r for r in records}
    rows = []
    for result in results:
        record = parents[result.peptide.source_id]
        if result.peptide.sequence != record["sequence"] or result.peptide.source_start != 0:
            raise ValueError("Cleavage input identity mismatch")
        bounds = core_bounds(record)
        for site in result.sites:
            rows.append(dict(
                sequence_record_id=record["sequence_record_id"], gene=record["gene"],
                sequence_sha256=record["sequence_sha256"], model=result.model.name,
                enzyme=result.model.enzyme, bond=site.bond,
                left_aa=record["sequence"][site.bond - 1], right_aa=record["sequence"][site.bond],
                score=site.score, score_units=result.model.score_units,
                inside_source_minimal_epitope=(bounds[0] < site.bond < bounds[1] if bounds else ""),
                reason=site.reason))
    return rows


def run(source, output_root):
    manifest = verified_source(source)
    records = vaccine_records(source)
    output = output_root / (datetime.now().astimezone().strftime("%Y-%m-%dT%H%M%S-%f%z") + "-stability")
    output.mkdir(parents=True, exist_ok=False)
    shutil.copyfile(Path(__file__).with_name("stability_models.json"), output / "model_availability.json")
    write_csv(output / "vaccine_records.csv", records)
    inputs = window_inputs(records)
    print("Running CleaveNet: %d windows x 18 MMPs" % len(inputs), flush=True)
    mmp_results = CleaveNet(subprocess_timeout=600).predict(inputs)
    write_csv(output / "mmp_window_scores.csv", window_rows(records, inputs, mmp_results))
    # Inventory/runtime is identical for the batch; retain it once, not per window.
    mmp_native = [r.to_dict() for r in mmp_results]
    inventory = mmp_native[0]["inventory"]
    runtime = mmp_native[0]["runtime"]
    for result in mmp_native:
        if result.pop("inventory") != inventory or result.pop("runtime") != runtime:
            raise ValueError("CleaveNet batch changed runtime identity")
    peptides = [CleavageInput(r["sequence"], n_term="unknown", c_term="unknown",
                             source_id=r["sequence_record_id"]) for r in records]
    peps_results = []
    for proteasome in ("C", "I"):
        print("Running human Pepsickle digestion " + proteasome, flush=True)
        peps_results.extend(Pepsickle(
            human_only=True, model_type="in-vitro-2", proteasome_type=proteasome,
            isolate_subprocess=True).predict_cleavage_many(peptides))
    write_csv(output / "proteasome_bond_scores.csv", bond_rows(records, peps_results))
    from mhctools.itcell_cleavage import ITCELL_MODELS, ITCellCleavage
    cat_results = []
    # H assesses only the initial, exposed N-terminal bond; no cascade simulation.
    cat_inputs = [CleavageInput(p.sequence, source_id=p.source_id) for p in peptides]
    for settings in ITCELL_MODELS.values():
        predictor = ITCellCleavage(**settings)
        cat_results.extend(predictor.predict(p) for p in cat_inputs)
    write_csv(output / "cathepsin_bond_scores.csv", bond_rows(records, cat_results))
    catalog = mechanism_catalog(read_csv(source / "tables/model_catalog.csv"), "human")
    assessments = enriched_assessments(
        records, read_csv(source / "tables/slp_quantitative_bond_scores.csv"),
        read_csv(source / "tables/slp_motif_assessments.csv"), catalog)
    # Keep every frozen extracellular model; selections are in the renderer only.
    write_csv(output / "existing_enzyme_assessments.csv", assessments)
    write_csv(output / "existing_model_catalog.csv", catalog)
    native = dict(cleavenet=dict(inventory=inventory, runtime=runtime, results=mmp_native),
                  proteasome=[r.to_dict() for r in peps_results],
                  cathepsin=[r.to_dict() for r in cat_results])
    (output / "native_predictions.json").write_text(json.dumps(native, indent=2, allow_nan=False) + "\n")
    shutil.copyfile(source / "provenance.json", output / "source_provenance.json")
    provenance = dict(
        generated_at=datetime.now().astimezone().isoformat(), source=str(source),
        source_manifest_sha256=sha256(source / "SHA256SUMS.json"), source_manifest=manifest,
        source_git_commit=subprocess.check_output(["git", "rev-parse", "HEAD"], text=True).strip(),
        runtime_python=sys.version, scripts={Path(__file__).name: sha256(__file__),
            "itcell_cleavage.py": sha256(Path(__file__).parents[2] / "mhctools/itcell_cleavage.py"),
            "itcell_profiles.json": sha256(Path(__file__).parents[2] / "mhctools/data/itcell_profiles.json"),
            "stability_models.json": sha256(Path(__file__).with_name("stability_models.json"))},
        record_count=len(records), unique_sequence_count=len({r["sequence"] for r in records}),
        mmp_window_count=len(inputs), mmp_head_count=len(ENZYMES),
        terminal_chemistry="unestablished; free canonical termini assumed for terminal cathepsin H and frozen enzyme panel",
        organism="human; human-only Pepsickle is experimental, preference does not establish superiority",
        core_scope="source minimal/candidate epitopes; not all vaccine targets are established",
        window_scope="overlapping ten-mer susceptibility; windows are not claimed released fragments",
        limitations="No dose, concentration, uptake, fragment cascade, half-life or combined serum loss inferred")
    (output / "provenance.json").write_text(json.dumps(provenance, indent=2) + "\n")
    print(output, flush=True)
    return output


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("source", type=Path)
    parser.add_argument("--output-root", type=Path, default=Path(__file__).parent / "local_results")
    args = parser.parse_args()
    run(args.source.resolve(), args.output_root.resolve())


if __name__ == "__main__":
    main()
