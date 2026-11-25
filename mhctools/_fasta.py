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

"""Lightweight standard-library FASTA readers and writers via ``datacache.fasta``."""

from datacache.fasta import (
    iter_fasta_records as _iter_fasta_records,
    parse_fasta_header_id as _parse_fasta_header_id,
    read_fasta as _read_fasta,
    write_fasta_records as _write_fasta_records,
)


def read_fasta_records(source, *, reject_pre_header=False):
    """Yield ``(header, sequence)`` pairs from a FASTA path or text stream.

    Blank lines and semicolon-prefixed comment lines are ignored. Full header
    text after ``>`` is preserved (with trailing whitespace stripped) and
    wrapped sequence lines are concatenated. Supports plain text and
    ``.gz``/``.gzip`` files when given a path.
    """
    try:
        for record in _iter_fasta_records(
            source,
            ignore_comments=True,
            allow_empty_header=True,
            allow_pre_header=not reject_pre_header,
        ):
            yield record.header.rstrip(), record.sequence
    except ValueError as error:
        if reject_pre_header and "before the first header" in str(error):
            raise ValueError(
                "Sequence data appears before the first FASTA header"
            ) from error
        raise


def write_fasta_records(target, records):
    """Write ``(header, sequence)`` pairs to a FASTA path or text stream."""
    _write_fasta_records(target, records)


def _parse_header_id(header):
    """Extract the leading identifier from a FASTA header (after ``>``)."""
    return _parse_fasta_header_id(header, split_pipe=True)


def parse_fasta_dictionary(fasta_path):
    """Parse a plain or gzipped FASTA file into ``{identifier: sequence}``.

    Identifiers are taken from the first space-delimited token of each header
    (and the first ``|``-delimited field for GENCODE-style pipe headers).
    Records without an identifier or sequence are skipped with a warning.
    """
    return _read_fasta(
        fasta_path,
        duplicate_policy="last",
        skip_empty_sequence=True,
        ignore_comments=True,
        allow_empty_header=True,
        allow_pre_header=True,
        split_pipe=True,
    )
