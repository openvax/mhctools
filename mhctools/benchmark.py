"""Assay-aware measurement evaluation and explicit provenance audits.

No unit conversion, homolog clustering, unreported-negative construction or
model training is performed. User-supplied family assignments are retained.
"""

from collections import Counter, defaultdict
from dataclasses import asdict, dataclass
from importlib.resources import files
import json
import math
from typing import Optional, Tuple


ENDPOINTS = frozenset((
    "site_cleavage", "substrate_depletion", "serum_half_life", "plasma_half_life",
    "whole_blood_half_life", "systemic_half_life", "systemic_clearance",
    "distribution_volume", "fluorescence_uptake", "cpp_classification",
    "cytosolic_delivery", "antigen_presentation"))


def _finite(value):
    return isinstance(value, (int, float)) and not isinstance(value, bool) and math.isfinite(value)


@dataclass(frozen=True)
class AssayMeasurement:
    """One source measurement; repetitions keep distinct measurement IDs.

    Chemistry is an exact descriptive identifier, independent of sequence.
    Unknown conditions must be explicit, rather than assumed physiological.
    """

    measurement_id: str
    source_measurement_id: str
    source: str
    dataset: str
    study: str
    assay: str
    sequence: str
    chemistry: str
    endpoint: str
    units: str
    species: str
    matrix: str
    value: Optional[float]
    split: str = "test"
    family: Optional[str] = None
    cell_type: Optional[str] = None
    enzyme: Optional[str] = None
    bond: Optional[int] = None
    censoring: str = "none"
    conditions: Tuple[Tuple[str, str], ...] = ()

    def __post_init__(self):
        for key in ("measurement_id", "source_measurement_id", "source", "dataset", "study",
                    "assay", "sequence", "chemistry", "units", "species", "matrix"):
            if not isinstance(getattr(self, key), str) or not getattr(self, key).strip():
                raise ValueError("Measurement requires explicit %s" % key)
        if not self.source.startswith(("https://", "http://")):
            raise ValueError("Measurement source must be a source URL")
        if self.endpoint not in ENDPOINTS:
            raise ValueError("Unknown measurement endpoint %r" % self.endpoint)
        if self.split not in ("train", "test", "reference"):
            raise ValueError("split must be train, test or reference")
        if self.censoring not in ("none", "less_than", "greater_than", "unknown"):
            raise ValueError("Unknown censoring status")
        if self.value is not None and not _finite(self.value):
            raise ValueError("Measurement value must be finite or None")
        if self.units == "binary" and self.value not in (0, 1, None):
            raise ValueError("Binary observations must be 0, 1 or None")
        if self.endpoint == "site_cleavage":
            if self.value not in (0, 1, None) or self.units != "binary":
                raise ValueError("Site cleavage requires binary observations or unknown")
            if (not isinstance(self.bond, int) or isinstance(self.bond, bool) or
                    not 1 <= self.bond < len(self.sequence) or not self.enzyme):
                raise ValueError("Site cleavage requires enzyme and an internal peptide bond")
        conditions = self.conditions.items() if isinstance(self.conditions, dict) else self.conditions
        conditions = tuple(sorted(tuple(pair) for pair in conditions))
        if any(len(p) != 2 or not all(isinstance(v, str) for v in p) for p in conditions):
            raise ValueError("Conditions must be string key/value pairs")
        if len(dict(conditions)) != len(conditions):
            raise ValueError("Duplicate condition keys")
        object.__setattr__(self, "conditions", conditions)


@dataclass(frozen=True)
class BenchmarkPrediction:
    """A prediction keyed to one exact measurement, with explicit output scale."""

    measurement_id: str
    model: str
    endpoint: str
    units: str
    value: Optional[float] = None
    status: str = "scored"
    scale: str = "native"
    reason: Optional[str] = None
    interval: Optional[Tuple[float, float]] = None
    interval_level: Optional[float] = None
    interval_target: Optional[str] = None

    def __post_init__(self):
        if not self.measurement_id or not self.model or self.endpoint not in ENDPOINTS:
            raise ValueError("Prediction requires measurement, model and known endpoint")
        if self.status not in ("scored", "unsupported", "failed", "not_assessed"):
            raise ValueError("Unknown prediction status")
        if self.scale not in ("native", "probability", "decision"):
            raise ValueError("Unknown prediction scale")
        if self.status == "scored":
            if not _finite(self.value):
                raise ValueError("Scored predictions require finite values")
            if self.scale == "probability" and not 0 <= self.value <= 1:
                raise ValueError("Probability must be within 0..1")
            if self.scale == "decision" and self.value not in (0, 1):
                raise ValueError("Decision must be 0 or 1")
        elif self.value is not None or not self.reason:
            raise ValueError("Unscored predictions require a reason and no value")
        if self.interval is not None:
            interval = tuple(self.interval)
            if (self.status != "scored" or len(interval) != 2 or
                    not all(_finite(v) for v in interval) or interval[0] > interval[1] or
                    not _finite(self.interval_level) or not 0 < self.interval_level < 1 or
                    self.interval_target != "individual_measurement"):
                raise ValueError("Intervals require ordered bounds, level and individual_measurement target")
            object.__setattr__(self, "interval", interval)
        elif self.interval_level is not None or self.interval_target is not None:
            raise ValueError("Interval metadata requires bounds")


@dataclass(frozen=True)
class ModelLineage:
    """Declared model training provenance; complete must be justified by source.

    Empty collections with unknown/partial provenance never establish absence
    of overlap. Families are externally curated, not guessed by this module.
    """

    model: str
    source: str
    provenance: str = "unknown"
    datasets: Tuple[str, ...] = ()
    studies: Tuple[str, ...] = ()
    assays: Tuple[str, ...] = ()
    sequences: Tuple[str, ...] = ()
    chemical_forms: Tuple[Tuple[str, str], ...] = ()
    families: Tuple[str, ...] = ()
    cell_types: Tuple[str, ...] = ()
    notes: str = ""

    def __post_init__(self):
        if self.provenance not in ("unknown", "partial", "complete"):
            raise ValueError("Unknown training provenance state")
        if not self.model or not self.source:
            raise ValueError("Lineage requires model and source")
        if self.provenance == "complete" and not (self.sequences and self.studies and self.families):
            raise ValueError("Complete lineage requires sequence, study and family inventories")
        for key in ("datasets", "studies", "assays", "sequences", "families", "cell_types"):
            object.__setattr__(self, key, tuple(getattr(self, key)))
        object.__setattr__(self, "chemical_forms", tuple(tuple(x) for x in self.chemical_forms))


def _overlaps(row, lineage):
    flags = []
    for attr, values in (("dataset", lineage.datasets), ("study", lineage.studies),
                         ("assay", lineage.assays), ("sequence", lineage.sequences),
                         ("family", lineage.families), ("cell_type", lineage.cell_types)):
        if getattr(row, attr) is not None and getattr(row, attr) in values:
            flags.append(attr)
    if (row.sequence, row.chemistry) in lineage.chemical_forms:
        flags.append("chemical_form")
    return flags


def _metrics(pairs, scale):
    if not pairs:
        return None
    truth = [m.value for m, p in pairs]
    predicted = [p.value for m, p in pairs]
    if pairs[0][0].units == "binary":
        if any(v not in (0, 1) for v in truth):
            return None
        if scale == "decision":
            return {"tp": sum(t == 1 and p == 1 for t, p in zip(truth, predicted)),
                    "tn": sum(t == 0 and p == 0 for t, p in zip(truth, predicted)),
                    "fp": sum(t == 0 and p == 1 for t, p in zip(truth, predicted)),
                    "fn": sum(t == 1 and p == 0 for t, p in zip(truth, predicted))}
        if scale == "probability":
            return {"brier": sum((t - p) ** 2 for t, p in zip(truth, predicted)) / len(pairs)}
        return None
    errors = [p - t for t, p in zip(truth, predicted)]
    return {"mae": sum(abs(e) for e in errors) / len(errors),
            "rmse": math.sqrt(sum(e * e for e in errors) / len(errors)),
            "bias": sum(errors) / len(errors)}


def evaluate_benchmark(measurements, predictions, lineages=(), evaluation="external_validation", requested_domains=()):
    """Evaluate comparable records, report overlap and retain every failure.

    Metrics are descriptive within exact assay/domain/output-scale strata.
    They never certify independence where training or family data are absent.
    """
    if evaluation not in ("external_validation", "reproduction"):
        raise ValueError("Unknown evaluation purpose")
    measurements = tuple(measurements)
    by_id = {m.measurement_id: m for m in measurements}
    if len(by_id) != len(measurements):
        raise ValueError("Duplicate measurement IDs")
    lineages = tuple(lineages)
    lineage_map = {x.model: x for x in lineages}
    if len(lineage_map) != len(lineages):
        raise ValueError("Duplicate model lineage")
    predictions = tuple(predictions)
    prediction_map = {(p.measurement_id, p.model): p for p in predictions}
    if len(prediction_map) != len(predictions):
        raise ValueError("Duplicate measurement/model predictions")
    if any(p.measurement_id not in by_id for p in predictions):
        raise ValueError("Prediction references an unknown measurement")
    models = sorted(set(p.model for p in predictions) | set(lineage_map))
    if not models:
        raise ValueError("At least one prediction model or lineage is required")
    training = [m for m in measurements if m.split == "train"]
    partition = ModelLineage("partition", "supplied benchmark train partition", "partial",
        studies=tuple(m.study for m in training), assays=tuple(m.assay for m in training),
        sequences=tuple(m.sequence for m in training), families=tuple(m.family for m in training if m.family),
        chemical_forms=tuple((m.sequence, m.chemistry) for m in training),
        cell_types=tuple(m.cell_type for m in training if m.cell_type))
    rows = []
    groups = defaultdict(list)
    for model in models:
        lineage = lineage_map.get(model, ModelLineage(model, "unknown"))
        for m in measurements:
            if m.split == "train":
                continue
            p = prediction_map.get((m.measurement_id, model))
            flags = _overlaps(m, lineage)
            partition_flags = _overlaps(m, partition) if m.split == "test" else []
            independent = (lineage.provenance == "complete" and m.family is not None and
                           not flags and not partition_flags and m.split == "test")
            reason = None
            if p is None:
                reason = "missing_prediction"
            elif p.status != "scored":
                reason = p.status
            elif m.value is None or m.censoring != "none":
                reason = "unknown_or_censored_observation"
            elif evaluation == "external_validation" and m.split != "test":
                reason = "reference_record"
            elif p.endpoint != m.endpoint or p.units != m.units:
                reason = "incompatible_endpoint_or_units"
            elif m.units == "binary" and p.scale not in ("decision", "probability"):
                reason = "native_score_has_no_binary_calibration"
            elif m.units != "binary" and p.scale != "native":
                reason = "incompatible_output_scale"
            row = {"measurement": asdict(m), "model": model,
                   "prediction": asdict(p) if p else None, "exclusion": reason,
                   "training_overlap": flags, "partition_overlap": partition_flags,
                   "training_provenance": lineage.provenance,
                   "family_audit": "available" if m.family is not None else "unknown",
                   "independence_established": independent}
            rows.append(row)
            key = (model, m.endpoint, m.units, m.species, m.matrix, m.cell_type,
                   m.assay, m.conditions, m.enzyme, p.scale if p else "unknown", m.study)
            groups[key].append((m, p, row))
    summaries = []
    for key, values in groups.items():
        comparable = [(m, p) for m, p, r in values if r["exclusion"] is None]
        independent_count = sum(r["independence_established"] for m, p, r in values)
        intervals = defaultdict(list)
        for m, p in comparable:
            if p.interval is not None:
                intervals[p.interval_level].append(p.interval[0] <= m.value <= p.interval[1])
        summaries.append({
            "model": key[0], "endpoint": key[1], "units": key[2], "species": key[3],
            "matrix": key[4], "cell_type": key[5], "assay": key[6],
            "conditions": dict(key[7]), "enzyme": key[8], "scale": key[9], "study": key[10],
            "measurement_count": len(values), "comparable_count": len(comparable),
            "unique_sequences": len({m.sequence for m, p, r in values}),
            "unique_source_measurements": len({(m.source, m.source_measurement_id) for m, p, r in values}),
            "unique_chemical_forms": len({(m.sequence, m.chemistry) for m, p, r in values}),
            "length_counts": dict(sorted(Counter(len(m.sequence) for m, p, r in values).items())),
            "exclusions": dict(Counter(r["exclusion"] for m, p, r in values if r["exclusion"])),
            "independent_measurement_count": independent_count,
            "claim": ("reproduction" if evaluation == "reproduction" else
                      "audited_external" if independent_count == len(values) else "unverified_external"),
            "descriptive_metrics": _metrics(comparable, key[9]),
            "prediction_interval_coverage": [{"nominal_level": level, "n": len(hits),
                "coverage": sum(hits) / len(hits)} for level, hits in sorted(intervals.items())]})
    domains = []
    for domain in requested_domains:
        allowed = {"name", "endpoint", "species", "matrix", "cell_type", "min_length", "max_length"}
        if not domain.get("name") or set(domain) - allowed:
            raise ValueError("Requested domains require a name and recognized constraints")
        for bound in ("min_length", "max_length"):
            if bound in domain and (not isinstance(domain[bound], int) or isinstance(domain[bound], bool) or domain[bound] < 1):
                raise ValueError("Domain length bounds must be positive integers")
        if domain.get("min_length", 1) > domain.get("max_length", float("inf")):
            raise ValueError("Inverted domain length bounds")
        def matches(m):
            return (all(getattr(m, key) == value for key, value in domain.items()
                        if key in ("endpoint", "species", "matrix", "cell_type")) and
                    domain.get("min_length", 1) <= len(m.sequence) <= domain.get("max_length", float("inf")))
        domain_rows = [m for m in measurements if m.split != "train" and matches(m)]
        domain_ids = {m.measurement_id for m in domain_rows}
        domains.append({"requested": dict(domain), "measurement_count": len(domain_rows),
                        "unique_sequences": len({m.sequence for m in domain_rows}),
                        "comparable_prediction_count": sum(r["exclusion"] is None and r["measurement"]["measurement_id"] in domain_ids for r in rows),
                        "evidence": "records_available" if domain_rows else "no_evidence_in_supplied_data"})
    return {"schema_version": 1, "evaluation": evaluation, "requested_domains": domains, "groups": summaries, "records": rows,
            "lineages": [asdict(x) for x in lineages],
            "limitations": "Repeated measurements are not independent peptides. Missing family/training data cannot establish external validation. No serum or delivery calibration is inferred from motif rules."}


def model_lineage_inventory():
    """Read the curated model/data relationships and unresolved provenance."""
    return json.loads(files("mhctools").joinpath("data/model_lineage.json").read_text())


def predict_cleavage_measurements(measurements, models):
    """Run site models against explicit measurements; unassessed bonds abstain.

    Only the exact canonical chemistry identifiers below can be represented
    by CleavageInput; all other chemistry remains an unsupported record.
    """
    from .cleavage import CleavageInput
    from .peptidases import get_cleavage_model
    chemistry = {"linear_L_free": ("free", "free"),
                 "linear_L_N_acetylated": ("acetylated", "free"),
                 "linear_L_C_amidated": ("free", "amidated")}
    predictions = []
    for name in models:
        try:
            predictor = get_cleavage_model(name)
        except (OSError, ImportError) as error:
            predictor = None
            failure = str(error)
        for m in measurements:
            kwargs = dict(measurement_id=m.measurement_id, model=name, endpoint=m.endpoint, units=m.units)
            if predictor is None:
                predictions.append(BenchmarkPrediction(**kwargs, status="failed", reason=failure))
                continue
            if m.endpoint != "site_cleavage" or m.enzyme != predictor.model.enzyme:
                predictions.append(BenchmarkPrediction(**kwargs, status="not_assessed", reason="Different endpoint or enzyme"))
                continue
            if m.chemistry not in chemistry:
                predictions.append(BenchmarkPrediction(**kwargs, status="unsupported", reason="Unrepresented chemistry"))
                continue
            try:
                n_term, c_term = chemistry[m.chemistry]
                selected = (get_cleavage_model(name, enzyme_state=dict(m.conditions).get("enzyme_state"))
                            if name == "cpb2-basic" else predictor)
                result = selected.predict(CleavageInput(m.sequence, n_term, c_term))
                site = next((s for s in result.sites if s.bond == m.bond), None)
                if site is None:
                    predictions.append(BenchmarkPrediction(**kwargs, status="not_assessed",
                        reason=result.unsupported_reason or "Bond outside model topology"))
                elif site.status == "scored":
                    predictions.append(BenchmarkPrediction(m.measurement_id, name,
                        "substrate_depletion" if name == "dpp4-qpisa" else "site_cleavage",
                        result.model.score_units, site.score))
                else:
                    predictions.append(BenchmarkPrediction(**kwargs, value=int(site.status == "matched"), scale="decision"))
            except (ValueError, TypeError, RuntimeError, OSError, ImportError) as error:
                predictions.append(BenchmarkPrediction(**kwargs, status="failed", reason=str(error)))
    return tuple(predictions)
