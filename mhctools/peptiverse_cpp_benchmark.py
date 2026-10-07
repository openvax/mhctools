"""Audit and reproduce released CPP source labels, without inferring assays.

Consumes a local, hash-verified metadata file. No datasets, weights, family
assignments, experimental chemistry or biological calibration are generated.
"""

from collections import Counter, defaultdict
import csv
import hashlib
from io import StringIO
import math
from pathlib import Path
import re

from .benchmark import AssayMeasurement, BenchmarkPrediction, ModelLineage, evaluate_benchmark
from .peptide_input import PeptideInput
from .peptiverse import UPSTREAM_REVISION
from .peptiverse_cpp import CPP_THRESHOLD
from .pred import Kind


CPP_METADATA_PATH = "training_data_cleaned/permeability_penetrance/permeability_meta_with_split.csv"
CPP_METADATA_SHA256 = "c924f1b92fac14f2007afb1b4b3641047896219a17a5780beb865ee0f4b35ec8"
CPP_METADATA_URL = ("https://huggingface.co/ChatterjeeLab/PeptiVerse/resolve/" +
                    UPSTREAM_REVISION + "/" + CPP_METADATA_PATH)
_DATASET = "PeptiVerse canonical CPP source metadata@" + UPSTREAM_REVISION
_MODEL = "peptiverse-cpp"
_UNKNOWN_FIELDS = ("study", "assay", "chemical_form", "species", "matrix", "cell_type",
                   "family", "concentration", "exposure_time", "internalization_readout")
_DOMAINS = (
    {"name": "beyond released training lengths", "min_length": 62},
    {"name": "primary human DC cytosolic delivery", "species": "Homo sapiens",
     "cell_type": "primary dendritic cell", "endpoint": "cytosolic_delivery"},
    {"name": "productive antigen presentation", "endpoint": "antigen_presentation"},
)
_NOTICE = (
    "Released source-label reproduction, not independent external validation. "
    "The publication reports validation-based model selection. Unknown original "
    "assays and chemical forms remain unknown; a CPP label or P(CPP) is not a "
    "physical uptake fraction, primary-DC delivery or antigen presentation. "
    "Exact-sequence nonoverlap does not certify study/family independence."
)


def _parse_metadata(text):
    reader = csv.DictReader(StringIO(text))
    if reader.fieldnames != ["sequence", "label", "id", "split"]:
        raise ValueError("CPP metadata requires exact sequence,label,id,split columns")
    rows, ids = [], set()
    for line, raw in enumerate(reader, 2):
        if (set(raw) != set(reader.fieldnames) or
                not all(isinstance(v, str) and v for v in raw.values()) or
                not re.fullmatch(r"seq_\d+", raw["id"]) or
                raw["id"] in ids or raw["label"] not in ("0", "1") or
                raw["split"] not in ("train", "val") or
                set(raw["sequence"]) - set("ACDEFGHIKLMNPQRSTVWY")):
            raise ValueError("Invalid CPP metadata or duplicate source ID at CSV line %d" % line)
        ids.add(raw["id"])
        rows.append(dict(raw, label=int(raw["label"])))
    if {row["split"] for row in rows} != {"train", "val"}:
        raise ValueError("CPP metadata requires nonempty train and val splits")
    return tuple(rows)


def load_cpp_metadata(path):
    """Read only the exact released source file; fail before loading models."""
    content = Path(path).read_bytes()
    if hashlib.sha256(content).hexdigest() != CPP_METADATA_SHA256:
        raise ValueError("CPP metadata SHA-256 does not match the pinned source snapshot")
    return _parse_metadata(content.decode("utf-8"))


def _audit(rows):
    labels = defaultdict(set)
    for row in rows:
        labels[row["sequence"]].add(row["label"])
    splits = {}
    for split in ("train", "val"):
        selected = [r for r in rows if r["split"] == split]
        lengths = Counter(len(r["sequence"]) for r in selected)
        splits[split] = dict(
            records=len(selected), unique_sequences=len({r["sequence"] for r in selected}),
            labels=dict(sorted(Counter(str(r["label"]) for r in selected).items())),
            min_length=min(lengths), max_length=max(lengths),
            length_counts=dict(sorted(lengths.items())),
            records_with_unknown_metadata={field: len(selected) for field in _UNKNOWN_FIELDS})
    train = {r["sequence"] for r in rows if r["split"] == "train"}
    val = {r["sequence"] for r in rows if r["split"] == "val"}
    return dict(source=CPP_METADATA_URL, sha256=CPP_METADATA_SHA256,
                source_revision=UPSTREAM_REVISION, records=len(rows), splits=splits,
                duplicate_source_ids=len(rows) - len({r["id"] for r in rows}),
                duplicate_sequences=len(rows) - len(labels),
                conflicting_label_sequences=sum(len(values) > 1 for values in labels.values()),
                exact_sequence_overlap_count=len(train & val),
                study_overlap="unknown", family_overlap="unknown", chemical_form_overlap="unknown",
                original_assay_metadata="unknown", notice=_NOTICE)


def _score_rows(selected, predictor):
    if not selected:
        return [], None
    try:
        if predictor is None:
            from .peptiverse_cpp import PeptiVerseCPP
            predictor = PeptiVerseCPP(device="cpu")
        inputs = [PeptideInput(r["sequence"], occurrence_id=r["id"]) for r in selected]
        results = predictor.predict(inputs, on_unsupported="record")
        if len(results) != len(inputs):
            raise ValueError("CPP benchmark predictor lost source records")
        predictions = []
        for row, item, result in zip(selected, inputs, results):
            if len(result.preds) != 1:
                raise ValueError("CPP benchmark requires one prediction per source record")
            pred = result.preds[0]
            context = pred.measurement_context
            if (pred.kind != Kind.cpp_classification or pred.peptide_input != item or
                    pred.peptide != row["sequence"] or pred.value is not None or context is None):
                raise ValueError("CPP benchmark prediction changed source identity or endpoint")
            kwargs = dict(measurement_id=row["id"], model=_MODEL,
                          endpoint="cpp_classification", units="binary", scale="probability")
            if context.status == "unsupported":
                predictions.append(BenchmarkPrediction(**kwargs, status="unsupported",
                    reason=context.detail or "Native CPP input unsupported"))
            elif (context.status != "available" or pred.score is None or
                  not math.isfinite(pred.score) or not 0 <= pred.score <= 1 or
                  context.score_semantics != "native SVC positive-class P(CPP)" or
                  context.class_label != ("CPP" if pred.score >= CPP_THRESHOLD else "non-CPP")):
                raise ValueError("CPP benchmark received invalid native class scores or labels")
            else:
                predictions.append(BenchmarkPrediction(**kwargs, value=pred.score))
        return predictions, predictor.artifact_inventory.to_dict()
    except Exception as error:
        # A sidecar or identity failure invalidates its batch. Keep every
        # selected source row rather than publishing a truncated success.
        return [BenchmarkPrediction(r["id"], _MODEL, "cpp_classification", "binary",
            status="failed", scale="probability", reason=str(error) or type(error).__name__)
            for r in selected], None


def evaluate_cpp_metadata(path, predict=False, source_ids=None, predictor=None):
    """Audit all splits and optionally reproduce CPP scores on validation IDs.

    The complete validation denominator remains in the report even for a
    partial cohort or runtime failure. ``predictor`` is an optional existing
    CPP adapter; prediction is always opt-in. Experimental chemistry is
    unknown, distinct from the canonical/free-terminus model-input assumption.
    """
    rows = load_cpp_metadata(path)
    val = [r for r in rows if r["split"] == "val"]
    if source_ids is not None:
        source_ids = tuple(source_ids)
        if (not source_ids or len(set(source_ids)) != len(source_ids) or
                set(source_ids) - {r["id"] for r in val}):
            raise ValueError("Select unique, nonempty validation source IDs only")
    if (source_ids is not None or predictor is not None) and not predict:
        raise ValueError("Cohort selection and predictor require explicit prediction")
    selected = [r for r in val if source_ids is None or r["id"] in source_ids] if predict else []
    predictions, inventory = _score_rows(selected, predictor)
    selected_ids = {r["id"] for r in selected}
    predictions.extend(BenchmarkPrediction(r["id"], _MODEL, "cpp_classification", "binary",
        status="not_assessed", scale="probability",
        reason="Outside selected source cohort" if predict else "Audit only; prediction not requested")
        for r in val if r["id"] not in selected_ids)
    measurements = [AssayMeasurement(
        measurement_id=r["id"], source_measurement_id=r["id"], source=CPP_METADATA_URL,
        dataset=_DATASET, study="unknown", assay="unknown", sequence=r["sequence"],
        chemistry="unknown", endpoint="cpp_classification", units="binary",
        species="unknown", matrix="unknown", value=r["label"], split="reference") for r in val]
    lineage = ModelLineage(_MODEL, CPP_METADATA_URL, provenance="partial", datasets=(_DATASET,),
        sequences=tuple(r["sequence"] for r in rows if r["split"] == "train"),
        notes="Released train membership; study, family, cell and chemical-form inventories unavailable.")
    benchmark = evaluate_benchmark(measurements, predictions, [lineage],
                                   evaluation="reproduction", requested_domains=_DOMAINS)
    by_id = {r["id"]: r for r in val}
    scored = [p for p in predictions if p.status == "scored"]
    confusion = {key: 0 for key in ("tp", "tn", "fp", "fn")}
    for pred in scored:
        observed, decision = by_id[pred.measurement_id]["label"], int(pred.value >= CPP_THRESHOLD)
        confusion[{(1, 1): "tp", (0, 0): "tn", (0, 1): "fp", (1, 0): "fn"}[observed, decision]] += 1
    from . import __version__
    return dict(schema_version=1, mhctools_version=__version__,
        evaluator_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(), source_audit=_audit(rows),
        cohort=dict(validation_records=len(val), selected_records=len(selected),
            selected_source_ids=[r["id"] for r in selected],
            full_validation_cohort=predict and len(selected) == len(val),
            status_counts=dict(Counter(p.status for p in predictions)),
            scored_records=len(scored), native_threshold=CPP_THRESHOLD,
            native_threshold_confusion=confusion if scored else None),
        model_input_assumption="Canonical L sequence with free termini for source-model reproduction; "
            "original measured chemical form is unknown",
        artifact_inventory=inventory, benchmark=benchmark, notice=_NOTICE)
