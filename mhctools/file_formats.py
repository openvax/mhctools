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
from math import ceil

from varcode import Variant, MutationEffect

from .epitope_collection_builder import EpitopeCollectionBuilder

def create_input_fasta_files(fasta_dictionary, max_file_records=None):
    """
    Turn peptide entries from a dataframe into one or more FASTA files.
    If mutation_window_size is an integer >0 then only use subsequence
    around mutated residues.

    Return the name of closed files which have to be manually deleted,
    and a dictionary from FASTA IDs to peptide records.

    If max_file_records is not provided, one file is created.
    If max_file_records is provided, each file will only hold (at most)
    max_file_records records.
    """
    n_fasta_records = len(fasta_dictionary)
    input_files = []
    # A file for every max_file_records records
    i_range = [0] if not max_file_records else range(
        int(ceil(n_fasta_records / max_file_records)))
    input_files = [tempfile.NamedTemporaryFile(
        "w", prefix="peptide", delete=False) for i in i_range]

    sequence_key_mapping = {}
    file_counter = 0
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
        input_files[file_counter].write(">%s\n%s" % (key, seq))
        # Don't add a newline after the last record
        if i + 1 < n_fasta_records:
            # If we are splitting files, don't add a newline as the
            # last line
            if max_file_records:
                if (i + 1) % max_file_records != 0:
                    input_files[file_counter].write("\n")
                else:
                    file_counter += 1
            else:
                input_files[file_counter].write("\n")
    input_file_names = []
    for input_file in input_files:
        input_file_names.append(input_file.name)
        input_file.close()
    return input_file_names, sequence_key_mapping


def split_stdout_lines(stdout):
    """
    Given the standard output from NetMHC/NetMHCpan/NetMHCcons tools,
    drop all {comments, lines of hyphens, empty lines} and split the
    remaining lines by whitespace.
    """
    for l in stdout.split("\n"):
        l = l.strip()
        if len(l) > 0 and not l.startswith("#"):
            yield l.split()

def parse_netmhc3_stdout(
        stdout,
        fasta_dictionary,
        prediction_method_name="netmhc",
        sequence_key_mapping=None):
    """
    Parse the output format for NetMHC 3.x, which looks like:

    ----------------------------------------------------------------------------------------------------
    pos    peptide      logscore affinity(nM) Bind Level    Protein Name     Allele
    ----------------------------------------------------------------------------------------------------
    0  SIINKFELL         0.437          441         WB              A1 HLA-A02:01
    --------------------------------------------------------------------------------------------------
    0  SIINKFFFQ         0.206         5411                         A2 HLA-A02:01
    1  IINKFFFQQ         0.128        12544                         A2 HLA-A02:01
    2  INKFFFQQQ         0.046        30406                         A2 HLA-A02:01
    3  NKFFFQQQQ         0.050        29197                         A2 HLA-A02:01
    --------------------------------------------------------------------------------------------------
    """

    builder = EpitopeCollectionBuilder(
        fasta_dictionary=fasta_dictionary,
        prediction_method_name=prediction_method_name)

    for fields in split_stdout_lines(stdout):
        if len(fields) >= 6:
            pos, peptide, log_affinity, ic50 = fields[:4]
            # annoyingly, the space between "affinity" and "Protein Name" may
            # have "WB" for weak binders and "SB" for strong binders. Couldn't
            # they at least have left those at the end?
            #
            # WARNING: if the sequence key is called "SB" or "WB" then those
            # lines will be ignored.
            #
            # TODO: use NetMHC's XLS output instead? Strangely, it seems
            # different from the format parsed below.
            if fields[4] == "WB" or fields[4] == "SB":
                if len(fields) < 7:
                    continue
                key, allele = fields[5:7]
            else:
                key, allele = fields[4:6]
            try:
                pos = int(pos)
                allele = str(allele)
                peptide = str(peptide)
                key = str(key)
                log_affinity = float(log_affinity)
                ic50 = float(ic50)
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
                log_ic50=log_affinity,
                rank=0.0)
    return builder.get_collection()

def parse_netmhc4_stdout(
        stdout,
        fasta_dictionary,
        prediction_method_name="netmhc",
        sequence_key_mapping=None):
    """
    # Peptide length 9
    # Rank Threshold for Strong binding peptides   0.500
    # Rank Threshold for Weak binding peptides   2.000
    -----------------------------------------------------------------------------------
      pos          HLA      peptide         Core Offset  I_pos  I_len  D_pos  D_len        iCore        Identity 1-log50k(aff) Affinity(nM)    %Rank  BindLevel
    -----------------------------------------------------------------------------------
        0    HLA-A0201    TMDKSELVQ    TMDKSELVQ      0      0      0      0      0    TMDKSELVQ 143B_BOVIN_P293         0.051     28676.59    43.00
        1    HLA-A0201    MDKSELVQK    MDKSELVQK      0      0      0      0      0    MDKSELVQK 143B_BOVIN_P293         0.030     36155.15    70.00
        2    HLA-A0201    DKSELVQKA    DKSELVQKA      0      0      0      0      0    DKSELVQKA 143B_BOVIN_P293         0.030     36188.42    70.00
        3    HLA-A0201    KSELVQKAK    KSELVQKAK      0      0      0      0      0    KSELVQKAK 143B_BOVIN_P293         0.032     35203.22    65.00
        4    HLA-A0201    SELVQKAKL    SELVQKAKL      0      0      0      0      0    SELVQKAKL 143B_BOVIN_P293         0.031     35670.99    65.00
        5    HLA-A0201    ELVQKAKLA    ELVQKAKLA      0      0      0      0      0    ELVQKAKLA 143B_BOVIN_P293         0.080     21113.07    29.00
        6    HLA-A0201    LVQKAKLAE    LVQKAKLAE      0      0      0      0      0    LVQKAKLAE 143B_BOVIN_P293         0.027     37257.56    75.00
        7    HLA-A0201    VQKAKLAEQ    VQKAKLAEQ      0      0      0      0      0    VQKAKLAEQ 143B_BOVIN_P293         0.040     32404.62    55.00
      219    HLA-A0201    QLLRDNLTL    QLLRDNLTL      0      0      0      0      0    QLLRDNLTL 143B_BOVIN_P293         0.527       167.10     1.50 <= WB
    -----------------------------------------------------------------------------------
    """
    builder = EpitopeCollectionBuilder(
        fasta_dictionary=fasta_dictionary,
        prediction_method_name=prediction_method_name)

    n_fields = 14
    for fields in split_stdout_lines(stdout):
        if len(fields) < n_fields:
            continue

        pos, allele, peptide, core, offset, i_pos, i_len, d_pos, d_len, i_core, key, log_affinity, ic50, rank = fields[:n_fields]
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

def parse_netmhcpan_stdout(
        stdout,
        fasta_dictionary,
        prediction_method_name="netmhcpan",
        sequence_key_mapping=None,
        contains_class2_columns=False):
    """
    Parse the output format for NetMHCpan, NetMHCIIpan* and NetMHCcons, which looks like:

     * netMHCIIpan has two extra fields

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

    # netMHCIIpan has some extra fields
    n_required_fields = 9 if contains_class2_columns else 7
    for fields in split_stdout_lines(stdout):
        if len(fields) >= n_required_fields:
            if contains_class2_columns:
                pos, allele, peptide, key, Pos, Core, log_affinity, ic50, rank = (
                    fields[:n_required_fields])
            else:
                pos, allele, peptide, key, log_affinity, ic50, rank = fields[:n_required_fields]
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

def parse_netmhciipan_stdout(
        stdout,
        fasta_dictionary,
        prediction_method_name,
        sequence_key_mapping=None):
    return parse_netmhcpan_stdout(
        stdout=stdout,
        fasta_dictionary=fasta_dictionary,
        prediction_method_name=prediction_method_name,
        sequence_key_mapping=sequence_key_mapping,
        contains_class2_columns=True)

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
