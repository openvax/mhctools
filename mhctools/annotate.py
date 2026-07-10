# Copyright (c) 2016. Mount Sinai School of Medicine
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

"""Annotate a table of peptides (and optional alleles) with predictor scores.

Downstream evaluation workflows often start from an annotated benchmark
table (columns like ``sample_id``, ``hit``, ``peptide``, and genotype/allele
information) rather than a bare peptide list. :func:`annotate_table` runs one
or more predictors over such a table and appends one score column per
predictor, choosing the best allele per row, while preserving every input
column unchanged.

The core function is I/O-free and works on any :class:`pandas.DataFrame`; the
``mhctools predict-table`` subcommand (see :mod:`mhctools.cli.annotate_table`)
wraps it with CSV reading/writing.

Design notes
------------
* Each predictor runs **once** over the union of all ``(peptide, allele)``
  pairs in the table (leveraging the predictors' allele batching), then a
  per-row lookup picks the best allele. Predictors are not re-invoked per row.
* "Best" comes from :func:`mhctools.pred.best_direction`: ``score`` is
  higher-better, ``affinity`` and ``percentile_rank`` are lower-better.
* Row alleles are normalized with
  :func:`mhctools.allele_normalization.normalize_allele_name_or_raw` so a
  cell written ``A0201`` matches a prediction emitted as ``HLA-A*02:01``, and
  exotic un-normalizable alleles still round-trip (see #220). Peptides are
  upper-cased and stripped on both sides of the join for the same reason.
"""

from __future__ import annotations

import re
from dataclasses import dataclass, field as _dc_field
from typing import Optional

import pandas as pd

from .allele_normalization import normalize_allele_name_or_raw
from .pred import Kind, reduce_op


# Whitespace, comma, or semicolon separate multiple alleles packed into one
# cell (e.g. a genotype column "HLA-A*02:01 HLA-B*07:02" or "A0201,B0702").
_DEFAULT_ALLELE_SEP = re.compile(r"[\s,;]+")


# Output-field token -> (kind, Prediction attribute). ``kind`` is ``None`` when
# the field's "best" direction is kind-independent (``score`` is always
# higher-better, ``percentile_rank`` always lower-better); ``affinity`` maps to
# the ``value`` of ``pMHC_affinity`` predictions (IC50 nM, lower-better).
_OUTPUT_FIELDS = {
    "affinity": (Kind.pMHC_affinity, "value"),
    "score": (None, "score"),
    "percentile_rank": (None, "percentile_rank"),
    # convenience aliases
    "rank": (None, "percentile_rank"),
    "value": (Kind.pMHC_affinity, "value"),
    "presentation": (Kind.pMHC_presentation, "score"),
    "processing": (Kind.antigen_processing, "score"),
    # NetMHCstabpan reports half-life (Thalf) in `score`, not `value`
    # (`parse_netmhcstabpan` sets no ic50), so read it from there.
    "stability": (Kind.pMHC_stability, "score"),
    "immunogenicity": (Kind.immunogenicity, "score"),
    "tap_transport": (Kind.tap_transport, "score"),
    "erap_trimming": (Kind.erap_trimming, "score"),
}


def output_field_tokens():
    """Sorted list of accepted output-field tokens (for CLI help / errors)."""
    return sorted(_OUTPUT_FIELDS)


@dataclass
class AnnotationSpec:
    """One predictor mapped to one output column.

    Parameters
    ----------
    predictor : predictor instance, or callable ``alleles -> predictor``
        A built predictor (used as-is), or a factory that builds one given
        the union of alleles in the table. Commandline predictors validate
        their alleles at construction, so passing a factory lets the predictor
        be built with exactly the alleles the table needs.
    output_column : str
        Name of the appended score column.
    field : str
        Which prediction field to write, one of :func:`output_field_tokens`
        (``"affinity"``, ``"score"``, ``"percentile_rank"``, ...). Determines
        both the value written and the best-allele direction.
    add_best_allele : bool
        If True (default) and the table has an allele column, also append a
        provenance column naming the allele that won each row.
    best_allele_column : str, optional
        Name of the provenance column. Defaults to
        ``f"{output_column}_best_allele"``.
    """

    predictor: object
    output_column: str
    field: str = "affinity"
    add_best_allele: bool = True
    best_allele_column: Optional[str] = None
    kind: Optional[str] = _dc_field(default=None, init=False)
    prediction_field: str = _dc_field(default="", init=False)

    def __post_init__(self):
        if self.field not in _OUTPUT_FIELDS:
            raise ValueError(
                "Unknown output field %r. Available: %s"
                % (self.field, ", ".join(output_field_tokens())))
        self.kind, self.prediction_field = _OUTPUT_FIELDS[self.field]

    def resolved_best_allele_column(self):
        return self.best_allele_column or ("%s_best_allele" % self.output_column)

    def build_predictor(self, alleles):
        """Return a predictor instance for the given union of alleles.

        A callable that is not itself a predictor (no ``predict`` method) is
        treated as a factory and called with ``alleles``; anything else is
        assumed to already be a built predictor and returned unchanged.
        """
        predictor = self.predictor
        if callable(predictor) and not hasattr(predictor, "predict"):
            return predictor(alleles)
        return predictor

    def direction_op(self):
        """``max`` or ``min`` callable for reducing candidate predictions."""
        return reduce_op(self.kind, self.prediction_field)


def parse_annotation_spec(token):
    """Parse a CLI predictor spec ``NAME:OUTPUT_COLUMN:FIELD``.

    ``OUTPUT_COLUMN`` and ``FIELD`` are optional; they default to
    ``f"{NAME}_{FIELD}"`` and ``"affinity"`` respectively. ``NAME`` must be a
    key in :data:`mhctools.cli.args.mhc_predictors`.

    Examples
    --------
    ``"netmhcpan42-ba:netmhcpan4.2.ba:affinity"``,
    ``"netmhcpan42-el:netmhcpan4.2.el:score"``, ``"mixmhcpred"``.
    """
    # Import the registry lazily to keep the core module free of the CLI's
    # (heavier) import graph and avoid an import cycle.
    from .cli.args import mhc_predictors, _cls_accepts

    parts = [p.strip() for p in str(token).split(":")]
    name = parts[0].lower()
    field = parts[2].lower() if len(parts) > 2 and parts[2] else "affinity"
    output_column = parts[1] if len(parts) > 1 and parts[1] else "%s_%s" % (name, field)

    if name not in mhc_predictors:
        raise ValueError(
            "Unknown predictor %r. Available: %s"
            % (name, ", ".join(sorted(mhc_predictors.keys()))))
    if field not in _OUTPUT_FIELDS:
        raise ValueError(
            "Unknown output field %r in spec %r. Available: %s"
            % (field, token, ", ".join(output_field_tokens())))

    cls = mhc_predictors[name]

    def factory(alleles):
        if alleles is not None and _cls_accepts(cls, "alleles"):
            return cls(alleles=alleles)
        return cls()

    return AnnotationSpec(
        predictor=factory, output_column=output_column, field=field)


def _normalize_peptide(peptide):
    """Canonicalize a peptide for matching: strip whitespace, uppercase.

    Predictions are joined back to rows by exact peptide string, and some
    predictors uppercase/strip the peptides they echo. Normalizing both the
    values sent to the predictor and the per-row lookup key the same way keeps
    a table written ``siinfekl`` (or with stray whitespace) from silently
    missing every prediction. Amino-acid sequences are canonically uppercase,
    so this never loses information.
    """
    return str(peptide).strip().upper()


def _split_alleles(cell, allele_sep):
    """Split one table cell into a list of normalized allele names."""
    if cell is None or (isinstance(cell, float) and pd.isna(cell)):
        return []
    tokens = [t for t in allele_sep.split(str(cell).strip()) if t]
    # Normalize + de-dup while preserving order.
    seen = {}
    for token in tokens:
        normalized = normalize_allele_name_or_raw(token)
        seen.setdefault(normalized, None)
    return list(seen)


def _build_lookup(results, spec):
    """Index one predictor's results by (peptide, allele) and by peptide.

    Returns ``(by_pair, by_peptide)`` where ``by_pair`` maps
    ``(peptide, allele) -> Prediction`` (allele-bearing predictions) and
    ``by_peptide`` maps ``peptide -> Prediction`` (allele-free predictions,
    e.g. processing predictors). Predictions whose target field is ``None``
    are skipped.

    A kind-agnostic token (``score``/``percentile_rank``/``rank``, i.e.
    ``spec.kind is None``) matches every kind, so if the predictor emits more
    than one kind for the same key the result would be an arbitrary
    last-one-wins pick. That is raised as an error rather than silently
    resolved — the caller should use a kind-specific token instead.
    """
    by_pair = {}
    by_peptide = {}
    for peptide_result in results:
        for pred in peptide_result.filter(kind=spec.kind):
            if getattr(pred, spec.prediction_field) is None:
                continue
            target, key = (
                (by_pair, (pred.peptide, pred.allele)) if pred.allele
                else (by_peptide, pred.peptide))
            existing = target.get(key)
            if (spec.kind is None and existing is not None
                    and existing.kind != pred.kind):
                raise ValueError(
                    "Ambiguous %r field: predictor emits multiple kinds "
                    "(%s, %s) for peptide %r allele %r. Use a kind-specific "
                    "token (e.g. affinity, presentation, immunogenicity) "
                    "instead of %r."
                    % (spec.field, existing.kind, pred.kind,
                       pred.peptide, pred.allele, spec.field))
            target[key] = pred
    return by_pair, by_peptide


def annotate_table(
        table,
        specs,
        peptide_column="peptide",
        allele_column=None,
        allele_sep=_DEFAULT_ALLELE_SEP,
        overwrite=False):
    """Append predictor score columns to a table, best-allele per row.

    Parameters
    ----------
    table : pandas.DataFrame
        Input table; returned unmodified with new columns appended to a copy.
    specs : sequence of AnnotationSpec
        One entry per output column.
    peptide_column : str
        Column holding the peptide sequence for each row.
    allele_column : str, optional
        Column holding the row's allele(s). May contain several alleles per
        cell (see ``allele_sep``). If omitted, predictors run allele-free and
        no best-allele provenance column is written.
    allele_sep : compiled regex
        Splits multiple alleles packed into one cell. Defaults to splitting on
        whitespace, commas, and semicolons.
    overwrite : bool
        If False (default), raise when an output column already exists.

    Returns
    -------
    pandas.DataFrame
        A copy of ``table`` with one numeric column per spec appended (plus a
        ``<output_column>_best_allele`` column when an allele column is given).
        Rows with no usable prediction get ``NaN`` (and ``None`` best allele).

    Notes
    -----
    Every allele in the table must be supported by each predictor;
    commandline predictors raise ``UnsupportedAllele`` at construction
    otherwise, matching mhctools' behavior elsewhere.
    """
    if not isinstance(table, pd.DataFrame):
        raise TypeError(
            "annotate_table expects a pandas DataFrame, got %s; use the "
            "predict-table CLI to read a file." % type(table).__name__)

    df = table.copy()

    if peptide_column not in df.columns:
        raise KeyError(
            "peptide column %r not found; columns are: %s"
            % (peptide_column, list(df.columns)))
    if allele_column is not None and allele_column not in df.columns:
        raise KeyError(
            "allele column %r not found; columns are: %s"
            % (allele_column, list(df.columns)))

    specs = list(specs)
    # Fail fast on any output-column collision before running predictors:
    # against existing input columns (unless overwrite), and against columns
    # any other spec plans to write (always — two specs writing the same
    # column would silently clobber each other, which overwrite can't fix).
    planned_columns = set()
    for spec in specs:
        new_columns = [spec.output_column]
        if spec.add_best_allele and allele_column is not None:
            new_columns.append(spec.resolved_best_allele_column())
        for column in new_columns:
            if column in df.columns and not overwrite:
                raise ValueError(
                    "output column %r already exists; pass overwrite=True to "
                    "replace it" % column)
            if column in planned_columns:
                raise ValueError(
                    "output column %r is produced by more than one spec" % column)
            planned_columns.add(column)

    peptides = [_normalize_peptide(p) for p in df[peptide_column]]
    if allele_column is not None:
        row_alleles = [_split_alleles(cell, allele_sep) for cell in df[allele_column]]
    else:
        row_alleles = [[] for _ in peptides]

    union_peptides = list(dict.fromkeys(peptides))
    union_alleles = sorted({a for alleles in row_alleles for a in alleles})

    for spec in specs:
        predictor = spec.build_predictor(union_alleles or None)
        results = predictor.predict(union_peptides)
        by_pair, by_peptide = _build_lookup(results, spec)
        reducer = spec.direction_op()
        field = spec.prediction_field

        values = []
        best_alleles = []
        for peptide, alleles in zip(peptides, row_alleles):
            candidates = [by_pair[(peptide, a)] for a in alleles
                          if (peptide, a) in by_pair]
            if not candidates:
                # Allele-free predictors (e.g. processing) index by peptide
                # only; fall back to that lookup.
                allele_free = by_peptide.get(peptide)
                if allele_free is not None:
                    candidates = [allele_free]
            if candidates:
                best = reducer(candidates, key=lambda p: getattr(p, field))
                values.append(getattr(best, field))
                best_alleles.append(best.allele or None)
            else:
                values.append(float("nan"))
                best_alleles.append(None)

        df[spec.output_column] = values
        if spec.add_best_allele and allele_column is not None:
            df[spec.resolved_best_allele_column()] = best_alleles

    return df
