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

"""Lightweight standard-library FASTA readers and writers."""

from contextlib import nullcontext
import gzip
import logging
import os

logger = logging.getLogger(__name__)


def _open_text_read(source):
    if isinstance(source, (str, os.PathLike)):
        path = os.fspath(source)
        if path.lower().endswith((".gz", ".gzip")):
            return gzip.open(path, "rt", encoding="utf-8")
        return open(path, "r", encoding="utf-8")
    return nullcontext(source)


def _open_text_write(target):
    if isinstance(target, (str, os.PathLike)):
        return open(os.fspath(target), "w", encoding="utf-8")
    return nullcontext(target)


def read_fasta_records(source, *, reject_pre_header=False):
    """Yield ``(header, sequence)`` pairs from a FASTA path or text stream.

    Blank lines and semicolon-prefixed comment lines are ignored. Full header
    text after ``>`` is preserved and wrapped sequence lines are concatenated.
    Supports plain text and ``.gz``/``.gzip`` files when given a path.
    """
    header = None
    sequence_lines = []
    with _open_text_read(source) as handle:
        for raw_line in handle:
            line = raw_line.strip()
            if not line or line.startswith(";"):
                continue
            if line.startswith(">"):
                if header is not None:
                    yield header, "".join(sequence_lines)
                header = line[1:]
                sequence_lines = []
            elif header is None:
                if reject_pre_header:
                    raise ValueError(
                        "Sequence data appears before the first FASTA header")
            else:
                sequence_lines.append(line)
    if header is not None:
        yield header, "".join(sequence_lines)


def write_fasta_records(target, records):
    """Write ``(header, sequence)`` pairs to a FASTA path or text stream."""
    with _open_text_write(target) as handle:
        for header, sequence in records:
            if not isinstance(header, str):
                raise TypeError("FASTA headers must be strings")
            if not isinstance(sequence, str):
                raise TypeError("FASTA sequences must be strings")
            clean_header = header.replace("\n", "").strip()
            handle.write(">%s\n%s\n" % (clean_header, sequence))


def _parse_header_id(header):
    """Extract the leading identifier from a FASTA header (after ``>``)."""
    if not header:
        raise ValueError("No identifier on FASTA line")
    identifier = header.split(" ", 1)[0]
    if "|" in identifier:
        identifier = identifier.split("|", 1)[0]
    return identifier


def parse_fasta_dictionary(fasta_path):
    """Parse a plain or gzipped FASTA file into ``{identifier: sequence}``.

    Identifiers are taken from the first space-delimited token of each header
    (and the first ``|``-delimited field for GENCODE-style pipe headers).
    Records without an identifier or sequence are skipped with a warning.
    """
    records = {}
    for header, sequence in read_fasta_records(fasta_path):
        identifier = _parse_header_id(header)
        if not identifier:
            logger.warning("Unable to parse ID from header line: >%s", header)
            continue
        if not sequence:
            logger.warning("No sequence data for '%s'", identifier)
            continue
        records[identifier] = sequence
    return records
