# Copyright (c) 2014. Mount Sinai School of Medicine
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

from varcode import Variant, MutationEffect

from .epitope_collection_builder import EpitopeCollectionBuilder

def create_input_fasta_file(fasta_dictionary):
    """
    Turn peptide entries from a dataframe into a FASTA file.
    If mutation_window_size is an integer >0 then only use subsequence
    around mutated residues.

    Return the name of closed file which has to be manually deleted,
    and a dictionary from FASTA IDs to peptide records.
    """
    input_file = tempfile.NamedTemporaryFile(
        "w", prefix="peptide", delete=False)
    n_fasta_records = len(fasta_dictionary)
    sequence_key_mapping = {}
    for i, (original_key, seq) in enumerate(fasta_dictionary.items()):
        unique_id = str(i)
        # make a nicer string representation for Variants and Effects
        # if those are the keys we're given
        if isinstance(original_key, (Variant, MutationEffect)):
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
        input_file.write(">%s\n%s" % (key, seq))
        # newline unless at end of file
        if i + 1 < n_fasta_records:
            input_file.write("\n")
    input_file.close()
    return input_file.name, sequence_key_mapping

def parse_netmhc_stdout(
        netmhc_output,
        fasta_dictionary,
        prediction_method_name="netmhc",
        sequence_key_mapping=None):
    """
    Parse the output format for NetMHC predictors, which looks like:

    # Affinity Threshold for Strong binding peptides  50.000',
    # Affinity Threshold for Weak binding peptides 500.000',
    # Rank Threshold for Strong binding peptides   0.500',
    # Rank Threshold for Weak binding peptides   2.000',
    ----------------------------------------------------------------------------
    pos  HLA  peptide  Identity 1-log50k(aff) Affinity(nM)    %Rank  BindLevel
    ----------------------------------------------------------------------------
      0  HLA-A*02:03    QQQQQYFPE   id0         0.024     38534.25   50.00
      1  HLA-A*02:03    QQQQYFPEI   id0         0.278      2461.53   15.00
      2  HLA-A*02:03    QQQYFPEIT   id0         0.078     21511.53   50.00
      3  HLA-A*02:03    QQYFPEITH   id0         0.041     32176.84   50.00
      4  HLA-A*02:03    QYFPEITHI   id0         0.085     19847.09   32.00
      5  HLA-A*02:03    YFPEITHII   id0         0.231      4123.85   15.00
      6  HLA-A*02:03    FPEITHIII   id0         0.060     26134.28   50.00
      7  HLA-A*02:03    PEITHIIIA   id0         0.034     34524.63   50.00
      8  HLA-A*02:03    EITHIIIAS   id0         0.076     21974.48   50.00
      9  HLA-A*02:03    ITHIIIASS   id0         0.170      7934.26   32.00
     10  HLA-A*02:03    THIIIASSS   id0         0.040     32361.18   50.00
     11  HLA-A*02:03    HIIIASSSL   id0         0.515       189.74    4.00 <= WB
    """

    builder = EpitopeCollectionBuilder(
        fasta_dictionary=fasta_dictionary,
        prediction_method_name=prediction_method_name)

    lines = netmhc_output.split("\n")
    lines = [l.strip() for l in lines]
    # remove empty lines
    lines = [l for l in lines if len(l) > 0]
    # remove comments
    lines = [l for l in lines if not l.startswith("#")]
    for line in lines:
        fields = line.split()
        n_required_fields = 7
        if len(fields) >= n_required_fields:
            pos, allele, peptide, key, log_affinity, ic50, rank = \
                fields[:n_required_fields]
            try:
                pos = int(pos)
                allele = str(allele)
                peptide = str(peptide)
                key = str(key)
                log_affinity = float(log_affinity)
                ic50 = float(ic50)
                rank = float(rank)
            except:
                # if position or affinity values can't be parsed,
                # then skip this line
                continue
            if sequence_key_mapping:
                original_key = sequence_key_mapping[key]
            else:
                # if sequence_key_mapping isn't provided then let's assume it's the
                # identity function
                original_key = key
            builder.add_binding_prediction(
                source_sequence_key=original_key,
                offset=pos,
                peptide=peptide,
                allele=allele,
                ic50=ic50,
                rank=rank,
                log_ic50=log_affinity)
    return builder.get_collection()

def parse_xls_file(
        xls_contents,
        fasta_dictionary,
        prediction_method_name,
        sequence_key_mapping=None):
    """
    XLS is a wacky output format used by NetMHCpan and NetMHCcons
    for peptide binding predictions.

    First line of XLS file format has HLA alleles
    and second line has fields like:
        ['Pos', 'Peptide', 'ID',
         '1-log50k', 'nM', 'Rank',
         '1-log50k', 'nM', 'Rank',
         '1-log50k', 'nM', 'Rank',
         ...'Ave', 'NB']
    """
    lines = [
        line.split("\t")
        for line in xls_contents.split("\n")
        if len(line) > 0
    ]

    if len(lines) == 0:
        raise ValueError("Empty XLS file for %s result" % (
            prediction_method_name,))

    # top line of XLS file has alleles
    alleles = [x for x in lines[0] if len(x) > 0]
    # skip alleles and column headers
    lines = lines[2:]

    builder = EpitopeCollectionBuilder(
        fasta_dictionary=fasta_dictionary,
        prediction_method_name=prediction_method_name)

    for fields in lines:
        pos = int(fields[0])
        epitope = fields[1]
        key = fields[2]
        for i, allele in enumerate(alleles):
            # we start at an offset of 3 to skip the allele-invariant
            # pos, epitope, identifier columns
            # each allele has three columns: log IC50, IC50, rank
            offset = 3 + 3 * i
            log_ic50 = float(fields[offset])
            ic50 = float(fields[offset + 1])
            rank = float(fields[offset + 2])
            if sequence_key_mapping:
                original_key = sequence_key_mapping[key]
            else:
                original_key = key
            builder.add_binding_prediction(
                source_sequence_key=original_key,
                offset=pos,
                peptide=epitope,
                allele=allele,
                ic50=ic50,
                rank=rank,
                log_ic50=log_ic50)
    return builder.get_collection()
