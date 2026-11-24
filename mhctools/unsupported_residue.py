# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#       http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

"""Typed exception for peptide or sequence inputs with unsupported residues."""

from __future__ import annotations

from collections.abc import Iterable, Sequence


def format_unsupported_residue_message(
    *,
    peptides: tuple[str, ...] = (),
    residues: tuple[str, ...] = (),
    field: str | None = None,
    predictor_name: str | None = None,
) -> str:
    """Build a human-readable message when a caller omits an explicit one."""
    prefix = f"{predictor_name}: " if predictor_name else ""
    label = field or "peptide"
    parts = [f"{prefix}Unsupported residue(s)"]
    if residues:
        parts.append(", ".join(repr(residue) for residue in residues))
    if peptides:
        if len(peptides) == 1:
            parts.append(f"in {label} {peptides[0]!r}")
        else:
            sample = ", ".join(repr(peptide) for peptide in peptides[:5])
            suffix = ", ..." if len(peptides) > 5 else ""
            parts.append(f"in {len(peptides)} {label} input(s): [{sample}{suffix}]")
    return " ".join(parts)


class UnsupportedResidueError(ValueError):
    """Raised when an input sequence contains residues unsupported by a model.

    Subclasses :class:`ValueError` so existing callers that catch
    ``ValueError`` continue to work, while callers performing batch recovery
    can distinguish unsupported-residue rejections from configuration,
    length, allele, or malformed-output errors by exception type.

    Attributes
    ----------
    peptides : tuple of str
        Distinct offending input sequences in order of first appearance.
    residues : tuple of str
        Sorted distinct unsupported residue/character strings, when known.
    input_indices : tuple of int
        Zero-based indices of offending inputs in the caller's batch, when
        known.
    field : str or None
        Name of the sequence field that contained the offending residue(s)
        (for example ``"peptide"``, ``"n_flank"``, ``"cdr3a"``), if known.
    predictor_name : str or None
        Predictor or backend name that rejected the residue(s), if known.
    """

    def __init__(
        self,
        message: str | None = None,
        *,
        peptides: Iterable[str] | str = (),
        residues: Iterable[str] | str = (),
        input_indices: Iterable[int] = (),
        field: str | None = None,
        predictor_name: str | None = None,
    ) -> None:
        if isinstance(peptides, str):
            normalized_peptides = (peptides,)
        else:
            normalized_peptides = tuple(dict.fromkeys(peptides))

        if isinstance(residues, str):
            normalized_residues = tuple(sorted(set(residues)))
        else:
            normalized_residues = tuple(sorted(set(residues)))

        normalized_indices = tuple(dict.fromkeys(int(i) for i in input_indices))

        self.peptides: tuple[str, ...] = normalized_peptides
        self.residues: tuple[str, ...] = normalized_residues
        self.input_indices: tuple[int, ...] = normalized_indices
        self.field: str | None = field
        self.predictor_name: str | None = predictor_name

        resolved_message = (
            message
            if message is not None
            else format_unsupported_residue_message(
                peptides=self.peptides,
                residues=self.residues,
                field=self.field,
                predictor_name=self.predictor_name,
            )
        )
        super().__init__(resolved_message)

    @property
    def peptide(self) -> str | None:
        """First offending input sequence, or ``None`` if unavailable."""
        return self.peptides[0] if self.peptides else None

    def __reduce__(self):
        return (
            self.__class__,
            (str(self),),
            {
                "peptides": self.peptides,
                "residues": self.residues,
                "input_indices": self.input_indices,
                "field": self.field,
                "predictor_name": self.predictor_name,
            },
        )


def remap_unsupported_residue_indices(
    error: UnsupportedResidueError,
    inputs: Sequence[str],
    *,
    field: str = "peptide",
) -> UnsupportedResidueError:
    """Remap offending indices onto the caller's original sequence batch.

    When a predictor wrapper deduplicates inputs before scoring, any
    :class:`UnsupportedResidueError` raised on the deduplicated list carries
    indices into the deduplicated list rather than the caller's original
    batch. This function produces a copy whose ``input_indices`` point to
    every matching occurrence in ``inputs``.
    """
    if not error.peptides or error.field not in (None, field):
        return error
    offending = set(error.peptides)
    remapped = tuple(
        index for index, value in enumerate(inputs) if value in offending
    )
    if not remapped or remapped == error.input_indices:
        return error
    return UnsupportedResidueError(
        str(error),
        peptides=error.peptides,
        residues=error.residues,
        input_indices=remapped,
        field=error.field,
        predictor_name=error.predictor_name,
    )
