#!/usr/bin/env python3
"""Replay a source-backed SLP case study; never infer negative cleavage labels."""

import argparse
import hashlib
import json
from pathlib import Path

from mhctools import predict_cleavage_batch, write_cleavage_batch
from mhctools import __version__


ROOT = Path(__file__).resolve().parents[1]
FIXTURE = ROOT / "tests/data/wada2018/figure2ac.json"
TRAINING_REVISION = "c448c4db81925afad78477e74a7d25e0209d3bce"
TRAINING_MAPS_SHA256 = "0c63172f5b877d6a06bcafa5027ef51d04c5947b6374966371d6e90564090543"


def product_bonds(value):
    """Map detected product boundaries, excluding retained parent endpoints."""
    sequence = value["sequence"]
    bonds, ids = {}, set()
    for product in value["observed_products"]:
        start, end = product["start"], product["end"]
        if (product["id"] in ids or not 0 <= start < end <= len(sequence) or
                sequence[start:end] != product["sequence"]):
            raise ValueError("Invalid source product %s" % product["id"])
        ids.add(product["id"])
        for boundary, bond in (("N", start), ("C", end)):
            if 0 < bond < len(sequence):
                bonds.setdefault(bond, []).append(dict(
                    product_id=product["id"], boundary=boundary,
                    first_detected_hours=product["first_detected_hours"]))
    return bonds


def context_window(sequence, bond):
    # Seven residues centered on P1 (the residue before the bond), matching
    # the GB model. Include terminal padding in exact-context overlap tests.
    padded = "***" + sequence + "***"
    return padded[bond - 1:bond + 6]


def audit_training(inputs, snapshot):
    """Audit all raw source windows, a conservative superset of fitted windows."""
    directory = Path(snapshot) / "data/raw/digestion_map_files"
    files = sorted(path for path in directory.iterdir() if path.is_file())
    inventory = [[path.name, hashlib.sha256(path.read_bytes()).hexdigest()] for path in files]
    digest = hashlib.sha256(json.dumps(inventory, separators=(",", ":")).encode()).hexdigest()
    if digest != TRAINING_MAPS_SHA256:
        raise ValueError("Training-map inventory differs from the pinned author snapshot")
    sequences, studies, unresolved_study_files = set(), set(), []
    for path in files:
        parts = path.read_text().split(">")
        references = [line.split("=", 1)[1].strip() for line in parts[0].splitlines()
                      if line.startswith("# DOI =")]
        studies.update(references)
        if not references or "?" in references:
            unresolved_study_files.append(path.name)
        for block in parts[1::2]:
            sequence = "".join(block.split())
            if not sequence or set(sequence) - set("ACDEFGHIKLMNPQRSTVWY"):
                raise ValueError("Unrecognized training source sequence")
            sequences.add(sequence)
    windows = {context_window(sequence, bond) for sequence in sequences
               for bond in range(1, len(sequence) + 1)}
    return dict(
        source="https://github.com/pdxgx/pepsickle-paper/tree/" + TRAINING_REVISION,
        inventory_sha256=digest, files=len(files), source_sequences=len(sequences),
        studies=sorted(studies), wada_study_present="10.1371/journal.pone.0199249" in studies,
        unresolved_study_files=unresolved_study_files,
        inputs=[dict(id=value["id"], full_sequence_overlap=value["sequence"] in sequences,
                     observed_bond_context_overlaps=[bond for bond in sorted(product_bonds(value))
                         if context_window(value["sequence"], bond) in windows]) for value in inputs],
        limitations="Raw-map study/exact-sequence/7-residue-context audit only; no homology-family "
        "audit or reconstruction of the original fitted/validation partitions. The publication "
        "labels Wada held out; this does not certify complete training independence.")


def evaluate(request):
    report = predict_cleavage_batch(request["inputs"], request["scenarios"], raise_on_error=True)
    report["experimental_source"] = request["source"]
    report["experimental_assay"] = request["assay"]
    report["evaluation_provenance"] = dict(
        mhctools_version=__version__,
        script_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        request_sha256=hashlib.sha256(json.dumps(request, sort_keys=True).encode()).hexdigest())
    inputs = {value["id"]: value for value in request["inputs"]}
    observed = []
    for row in report["assessments"]:
        sites = {site["bond"]: site for site in row["result"]["sites"]}
        observed.append(dict(
            input_id=row["input_id"], model=row["model"],
            observed_bonds=[dict(bond=bond, products=products,
                                 prediction=sites.get(bond))
                            for bond, products in sorted(product_bonds(inputs[row["input_id"]]).items())]))
    report["experimental_product_overlays"] = observed
    report["validation_scope"] = (
        "Two complete figure panels, 47 detected products from two reordered constructs. "
        "Native site scores are compared with observed internal product boundaries. "
        "No negative labels, cleavage kinetics, product yields, AUC, survival, or presentation "
        "accuracy are inferred. First detection of a product does not date each cleavage event.")
    return report


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", required=True, help="Output batch evidence JSON")
    parser.add_argument("--html", help="Optional batch HTML report")
    parser.add_argument("--training-snapshot", type=Path,
                        help="Extracted pepsickle-paper snapshot at " + TRAINING_REVISION)
    args = parser.parse_args()
    request = json.loads(FIXTURE.read_text())
    report = evaluate(request)
    if args.training_snapshot:
        report["training_audit"] = audit_training(request["inputs"], args.training_snapshot)
    else:
        report["training_audit"] = {"status": "not_run"}
    write_cleavage_batch(report, args.out, html_path=args.html)


if __name__ == "__main__":
    main()
