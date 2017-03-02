# Copyright (c) 2014-2017. Mount Sinai School of Medicine
#
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

from __future__ import print_function, division, absolute_import
import tempfile
from math import ceil

def _create_input_files(
        n_sequences,
        max_sequences_per_file=None,
        suffix=".txt"):
    """
    Create input files that we're going to fill with n_sequences entries.
    If max_sequences_per_file is specified and is less than n_sequences then
    multiple files will be returned.
    """
    if not max_sequences_per_file:
        n_files = 1
    else:
        n_files = int(ceil(n_sequences / max_sequences_per_file))

    return [
        tempfile.NamedTemporaryFile(
            "w",
            prefix="input_file_%d" % (i + 1),
            suffix=suffix,
            delete=False)
        for i in range(n_files)
    ]

def _close_and_return_file_names(files):
    """
    Close each file handle and return their names.
    """
    file_names = []
    for f in files:
        file_names.append(f.name)
        f.close()
    return file_names

def create_input_fasta_files(fasta_dictionary, max_sequences_per_file=None):
    """
    Return the name of closed files which have to be manually deleted,
    and a dictionary from FASTA IDs to sequence records.

    If max_sequences_per_file is not provided, one file is created.
    If max_sequences_per_file is provided, each file will only hold (at most)
    max_sequences_per_file records.
    """
    n_fasta_records = len(fasta_dictionary)
    input_files = _create_input_files(
        n_sequences=n_fasta_records,
        max_sequences_per_file=max_sequences_per_file,
        suffix=".fasta")

    sequence_key_mapping = {}
    file_counter = 0
    for i, (original_key, seq) in enumerate(fasta_dictionary.items()):
        unique_id = str(i)
        # make a nicer string representation for Variants and Effects
        # if those are the keys we're given
        if hasattr(original_key, 'short_description'):
            unsanitized_string = original_key.short_description
        else:
            unsanitized_string = str(original_key)
        sanitized = "".join([
            c if c.isalnum() else "_"
            for c in unsanitized_string
        ])
        # this looks crazy but NetMHCpan seems to require key lengths
        # of 15 or shorter, but since I still want the generated FASTA
        # file to be vaguely readable I'm combining a truncation of the
        # original key with a unique indentifier
        key = sanitized[: 14 - len(unique_id)] + "_" + unique_id
        sequence_key_mapping[key] = original_key
        input_files[file_counter].write(">%s\n%s" % (key, seq))

        # Don't add a newline after the last record
        if i + 1 < n_fasta_records:
            # If we are splitting files, don't add a newline as the
            # last line
            if max_sequences_per_file:
                if (i + 1) % max_sequences_per_file != 0:
                    input_files[file_counter].write("\n")
                else:
                    file_counter += 1
            else:
                input_files[file_counter].write("\n")
    input_file_names = _close_and_return_file_names(input_files)
    return input_file_names, sequence_key_mapping

def create_input_peptides_files(peptides, max_peptides_per_file=None):
    """
    Creates one or more files containing one peptide per line,
    returns names of files.
    """
    n_peptides = len(peptides)
    if not max_peptides_per_file:
        max_peptides_per_file = n_peptides
    input_files = _create_input_files(
        n_sequences=n_peptides,
        max_sequences_per_file=max_peptides_per_file,
        suffix=".txt")
    for i, p in enumerate(peptides):
        f = input_files[i // max_peptides_per_file]
        f.write("%s\n" % p)
    return _close_and_return_file_names(input_files)
