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
import logging

from mhcnames import parse_classi_or_classii_allele_name
from typechecks import require_integer

from .common import check_sequence_dictionary
from .base_commandline_predictor import BaseCommandlinePredictor
from .parsing import parse_netmhciipan_stdout
from .input_file_formats import create_input_fasta_files

logger = logging.getLogger(__name__)

class NetMHCIIpan(BaseCommandlinePredictor):
    def __init__(
            self,
            alleles,
            default_peptide_lengths=[15, 16, 17, 18, 19, 20],
            program_name="netMHCIIpan",
            max_sequences_per_fasta_file=1000,
            process_limit=-1):
        BaseCommandlinePredictor.__init__(
            self,
            program_name=program_name,
            alleles=alleles,
            default_peptide_lengths=default_peptide_lengths,
            parse_output_fn=parse_netmhciipan_stdout,
            supported_alleles_flag="-list",
            input_file_flag="-f",
            allele_flag="-a",
            peptide_mode_flags=["-inptype", "1"],
            length_flag="-length",
            tempdir_flag="-tdir",
            process_limit=process_limit)
        require_integer(
            max_sequences_per_fasta_file,
            "Maximum number of protein sequences per FASTA file")
        self.max_sequences_per_fasta_file = max_sequences_per_fasta_file

    def prepare_allele_name(self, allele_name):
        """
        netMHCIIpan has some unique requirements for allele formats,
        expecting the following forms:
         - DRB1_0101 (for non-alpha/beta pairs)
         - HLA-DQA10501-DQB10636 (for alpha and beta pairs)

        Other than human class II alleles, the only other alleles that
        netMHCIIpan accepts are the following mouse alleles:
         - H-2-IAb
         - H-2-IAd
        """
        parsed_alleles = parse_classi_or_classii_allele_name(allele_name)
        if len(parsed_alleles) == 1:
            allele = parsed_alleles[0]
            if allele.species == "H-2":
                return "%s-%s%s" % (
                    allele.species,
                    allele.gene,
                    allele.allele_code)
            else:
                # assume that we're dealing with a human DRB allele
                # which NetMHCIIpan treats differently because there is
                # supposedly no population diversity in the DR-alpha gene
                return "%s_%s%s" % (
                    allele.gene,
                    allele.allele_family,
                    allele.allele_code)
        else:
            alpha, beta = parsed_alleles
            return "HLA-%s%s%s-%s%s%s" % (
                alpha.gene,
                alpha.allele_family,
                alpha.allele_code,
                beta.gene,
                beta.allele_family,
                beta.allele_code)

    def predict_subsequences_using_fasta_mode(
            self,
            sequence_dict,
            peptide_lengths=None):
        """
        Predict MHC ligands from collection of protein sequences.

        Runs multiple predictors simultaneously, split across:
            1) multiple input files, each containing
                self.max_sequences_per_fasta_file
            2) every allele + peptide length combination

        Parameters
        ----------
        sequence_dict : dict
            Dictionary whose keys are expected to be unique identifiers for
            each protein sequence and whose values are amino acid sequences.

        peptide_lengths : list of int
            List of peptide lengths for which to make predictions. If omitted
            then use the default_peptide_lengths property of the predictor.

        Returns list of BindingPrediction objects
        """
        peptide_lengths = self._check_peptide_lengths(peptide_lengths)
        sequence_dict = check_sequence_dictionary(sequence_dict)

        input_filenames, sequence_key_mapping = create_input_fasta_files(
            sequence_dict,
            max_sequences_per_file=self.max_sequences_per_fasta_file)
        commands = {}
        dirs = []
        for i, input_filename in enumerate(input_filenames):
            for j, allele in enumerate(self.alleles):
                for length in peptide_lengths:
                    if self.tempdir_flag:
                        temp_dirname = tempfile.mkdtemp(
                            prefix="tmp_%s_length_%d" % (
                                self.program_name, length))
                        logger.debug(
                            "Created temporary directory %s for allele %s, length %d",
                            temp_dirname,
                            allele,
                            length)
                        dirs.append(temp_dirname)
                    else:
                        temp_dirname = None
                    output_file = tempfile.NamedTemporaryFile(
                        "w",
                        prefix="%s_output_%d_%d" % (self.program_name, i, j),
                        delete=False)
                    commands[output_file] = self._build_command(
                        input_filename=input_filename,
                        allele=allele,
                        length=length,
                        temp_dirname=temp_dirname)
        return self._run_commands_and_collect_predictions(
            commands=commands,
            input_filenames=input_filenames,
            temp_dir_list=dirs,
            sequence_key_mapping=sequence_key_mapping)
