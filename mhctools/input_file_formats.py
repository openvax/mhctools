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

def make_writable_tempfile(prefix_number, prefix_name, suffix):
    prefix = "input_file_%d_%s" % (prefix_number, prefix_name)
    return tempfile.NamedTemporaryFile(
        "w",
        prefix=prefix,
        suffix=suffix,
        delete=False)

def create_input_peptides_files(
        peptides,
        max_peptides_per_file=None,
        group_by_length=False):
    """
    Creates one or more files containing one peptide per line,
    returns names of files.
    """
    if group_by_length:
        peptide_lengths = {len(p) for p in peptides}
        peptide_groups = {l: [] for l in peptide_lengths}
        for p in peptides:
            peptide_groups[len(p)].append(p)
    else:
        peptide_groups = {"": peptides}

    file_names = []
    for key, group in peptide_groups.items():
        n_peptides = len(group)
        if not max_peptides_per_file:
            max_peptides_per_file = n_peptides
        input_file = None
        for i, p in enumerate(group):
            if i % max_peptides_per_file == 0:
                if input_file is not None:
                    file_names.append(input_file.name)
                    input_file.close()
                input_file = make_writable_tempfile(
                    prefix_number=i // max_peptides_per_file,
                    prefix_name=key,
                    suffix=".txt")
            input_file.write("%s\n" % p)
        if input_file is not None:
            file_names.append(input_file.name)
            input_file.close()
    return file_names
