# Copyright (c) 2026 Mount Sinai School of Medicine
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0

"""Exact peptide chemical form and descriptive prediction context."""

from collections.abc import Mapping
from dataclasses import asdict, dataclass, fields
import hashlib
import json
import math


CANONICAL_AMINO_ACIDS = frozenset("ACDEFGHIKLMNPQRSTVWY")
N_TERM_VALUES = frozenset(("free", "acetylated", "unknown"))
C_TERM_VALUES = frozenset(("free", "amidated", "unknown"))


def _optional_text(name, value):
    if value is not None and (not isinstance(value, str) or not value.strip()):
        raise ValueError(f"{name} must be a nonempty string or None")


def _identity_sha256(payload):
    encoded = json.dumps(
        payload, sort_keys=True, separators=(",", ":")).encode("utf-8")
    return hashlib.sha256(encoded).hexdigest()


@dataclass(frozen=True)
class PeptideContext:
    """Versioned administration, study, assay, and cellular context.

    Every field is descriptive. ``None`` means unknown; this class does not
    infer patient-specific properties, physiological defaults, or a preferred
    formulation, route, or modification.
    """

    administration_route: str | None = None
    formulation: str | None = None
    cargo: str | None = None
    administered_material: str | None = None
    released_material: str | None = None
    measured_analyte: str | None = None
    study_id: str | None = None
    assay_species: str | None = None
    matrix: str | None = None
    cell_type: str | None = None
    cell_subtype: str | None = None
    maturation_state: str | None = None
    readout: str | None = None
    timepoint: float | None = None
    time_unit: str | None = None
    conditions: tuple[tuple[str, str], ...] = ()
    schema_version: int = 1

    def __post_init__(self):
        if self.schema_version != 1:
            raise ValueError(
                f"Unsupported PeptideContext schema_version "
                f"{self.schema_version!r}")
        for name in (
                "administration_route", "formulation", "cargo",
                "administered_material", "released_material",
                "measured_analyte", "study_id", "assay_species", "matrix",
                "cell_type", "cell_subtype", "maturation_state", "readout",
                "time_unit"):
            _optional_text(name, getattr(self, name))
        if ((self.timepoint is None) != (self.time_unit is None)):
            raise ValueError("timepoint and time_unit must be provided together")
        if (self.timepoint is not None and
                (isinstance(self.timepoint, bool) or
                 not isinstance(self.timepoint, (int, float)) or
                 not math.isfinite(self.timepoint))):
            raise ValueError("timepoint must be a finite number or None")
        conditions = (
            self.conditions.items()
            if isinstance(self.conditions, Mapping) else self.conditions)
        conditions = tuple(sorted(tuple(pair) for pair in conditions))
        if any(
                len(pair) != 2 or
                not all(isinstance(value, str) and value for value in pair)
                for pair in conditions):
            raise ValueError("conditions must be nonempty string key/value pairs")
        if len(dict(conditions)) != len(conditions):
            raise ValueError("conditions keys must be unique")
        object.__setattr__(self, "conditions", conditions)

    def to_dict(self):
        """Return a JSON-friendly representation."""
        return asdict(self)

    @classmethod
    def from_dict(cls, value):
        """Restore a context while ignoring unknown same-version fields."""
        if isinstance(value, cls):
            return value
        valid = {field.name for field in fields(cls)}
        return cls(**{key: item for key, item in value.items() if key in valid})


@dataclass(frozen=True)
class PeptideInput:
    """Canonical L-peptide chemical form plus provenance and context.

    ``attachments`` contains ``(site, identity)`` pairs. Sites and identities
    are deliberately descriptive strings: adapters may support a documented
    subset, and must reject every other form rather than dropping it. Source
    occurrence provenance is kept separate from the chemical/context identity
    used for inference caching.
    """

    sequence: str
    n_term: str = "free"
    c_term: str = "free"
    attachments: tuple[tuple[str, str], ...] = ()
    occurrence_id: str | None = None
    source_sequence_name: str | None = None
    source_start: int = 0
    source_gene: str | None = None
    source_species: str | None = None
    context: PeptideContext | None = None
    schema_version: int = 1

    def __post_init__(self):
        if self.schema_version != 1:
            raise ValueError(
                f"Unsupported PeptideInput schema_version "
                f"{self.schema_version!r}")
        if not isinstance(self.sequence, str) or not self.sequence:
            raise ValueError(
                "Empty peptide: expected a canonical uppercase L-peptide")
        if (self.sequence != self.sequence.upper() and
                set(self.sequence.upper()) <= CANONICAL_AMINO_ACIDS):
            raise ValueError("PeptideInput sequence must be uppercase")
        invalid = set(self.sequence) - CANONICAL_AMINO_ACIDS
        if invalid:
            raise ValueError(
                "non-standard residues cannot describe an unmodified "
                f"sequence: {''.join(sorted(invalid))}")
        if self.n_term not in N_TERM_VALUES:
            raise ValueError(
                f"n_term must be one of {sorted(N_TERM_VALUES)}, got "
                f"{self.n_term!r}")
        if self.c_term not in C_TERM_VALUES:
            raise ValueError(
                f"c_term must be one of {sorted(C_TERM_VALUES)}, got "
                f"{self.c_term!r}")
        attachments = (
            self.attachments.items()
            if isinstance(self.attachments, Mapping) else self.attachments)
        attachments = tuple(tuple(pair) for pair in attachments)
        if any(
                len(pair) != 2 or
                not all(isinstance(value, str) and value for value in pair)
                for pair in attachments):
            raise ValueError(
                "attachments must be nonempty (site, identity) string pairs")
        if len(set(attachments)) != len(attachments):
            raise ValueError("attachments must be unique")
        object.__setattr__(self, "attachments", tuple(sorted(attachments)))
        if (not isinstance(self.source_start, int) or
                isinstance(self.source_start, bool) or self.source_start < 0):
            raise ValueError("source_start must be a nonnegative integer")
        for name in (
                "occurrence_id", "source_sequence_name", "source_gene",
                "source_species"):
            _optional_text(name, getattr(self, name))
        if isinstance(self.context, Mapping):
            object.__setattr__(
                self, "context", PeptideContext.from_dict(self.context))
        elif self.context is not None and not isinstance(
                self.context, PeptideContext):
            raise TypeError("context must be PeptideContext, a mapping, or None")

    @property
    def chemical_identity_sha256(self):
        """Identity of sequence, termini, and attachments only."""
        return _identity_sha256({
            "sequence": self.sequence,
            "n_term": self.n_term,
            "c_term": self.c_term,
            "attachments": self.attachments,
        })

    @property
    def inference_identity_sha256(self):
        """Chemical and model-input context identity, excluding occurrence."""
        return _identity_sha256({
            "schema": "mhctools.PeptideInput.inference.v1",
            "chemical_identity_sha256": self.chemical_identity_sha256,
            "context": (
                self.context.to_dict()
                if self.context else PeptideContext().to_dict()),
        })

    @property
    def record_identity_sha256(self):
        """Full record identity including occurrence and source provenance."""
        return _identity_sha256(self.to_dict())

    def to_dict(self):
        """Return a JSON-friendly representation."""
        result = asdict(self)
        result["context"] = self.context.to_dict() if self.context else None
        return result

    @classmethod
    def from_dict(cls, value):
        """Restore an input while ignoring unknown same-version fields."""
        if isinstance(value, cls):
            return value
        valid = {field.name for field in fields(cls)}
        values = {
            key: item for key, item in value.items() if key in valid}
        if values.get("context") is not None:
            values["context"] = PeptideContext.from_dict(values["context"])
        return cls(**values)


def coerce_peptide_inputs(values):
    """Return exact inputs; strings remain the legacy natural-peptide shorthand."""
    if isinstance(values, (str, PeptideInput)):
        values = [values]
    result = []
    for value in values:
        if isinstance(value, PeptideInput):
            result.append(value)
        elif isinstance(value, str):
            result.append(PeptideInput(value.strip().upper()))
        else:
            raise TypeError("Expected PeptideInput or peptide string")
    return result


def sequence_only_chemistry_error(peptide_input):
    """Reason a sequence-only natural-peptide model cannot score an input."""
    if peptide_input.n_term != "free":
        return f"unsupported N-terminal chemistry: {peptide_input.n_term}"
    if peptide_input.c_term != "free":
        return f"unsupported C-terminal chemistry: {peptide_input.c_term}"
    if peptide_input.attachments:
        return "attachments are unsupported by this sequence-only model"
    return None
