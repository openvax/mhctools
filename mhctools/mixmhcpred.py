# Copyright (c) 2019. Mount Sinai School of Medicine
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

import pandas as pd
from tempfile import mkdtemp, NamedTemporaryFile
from os.path import join, exists
from os import remove

from mhcnames import normalize_allele_name

from .base_predictor import BasePredictor
from .binding_prediction import BindingPrediction
from .binding_prediction_collection import BindingPredictionCollection
from .process_helpers import run_command
from .cleanup_context import CleanupFiles

class MixMHCpred(BasePredictor):
    def __init__(
            self,
            alleles,
            default_peptide_lengths=[9],
            program_name="MixMHCpred",
            exclude_peptides_with_cysteine=False):
        """
        Wrapper for MixMHCpred

        Parameters
        ----------
        alleles : list of str
        
        default_peptide_lengths : list of int
        
        program_name : str
        
        exclude_peptides_with_cysteine : bool
            If True then drop peptides which contain 'C' from predictions, 
            default is True.    
        """
        BasePredictor.__init__(
            self,
            alleles=alleles,
            default_peptide_lengths=default_peptide_lengths,
            min_peptide_length=8,
            max_peptide_length=14,
            allow_X_in_peptides=False,
            allow_lowercase_in_peptides=False)
        self.program_name = program_name
        self.exclude_peptides_with_cysteine = exclude_peptides_with_cysteine
        if self.exclude_peptides_with_cysteine:
            self.extra_commandline_args = ["-c"]
        else:
            self.extra_commandline_args = []

    def predict_peptides(self, peptides):
        """

        Parameters
        ----------
        peptides : list of str


        Returns
        -------
        list of BindingPrediction
        """
        self._check_peptide_inputs(peptides)
        results = []
        for allele in self.alleles:

            temp_dir = mkdtemp(prefix="mhctools", suffix="mixmhcpred")
            input_file_path = join(temp_dir, "mixmhcpred_inputs.txt")
            output_file_path = join(temp_dir, "mixmhcpred_outputs.txt")

            with open(input_file_path, "w") as f:
                for i, p in enumerate(peptides):
                    f.write(p)
                    if i < len(peptides) - 1:
                        f.write("\n")
            with CleanupFiles(
                    filenames=[input_file_path, output_file_path],
                    directories=[temp_dir]):
                with NamedTemporaryFile(prefix="MixMHCpred_stdout", mode="w", delete=False) as stdout_file:
                    stdout_file_name = stdout_file.name
                    run_command([
                        self.program_name,
                        "-i", input_file_path,
                        "-o", output_file_path,
                        "-a", normalize_allele_name(allele)] + self.extra_commandline_args,
                        suppress_stderr=False,
                        redirect_stdout_file=stdout_file)
                if exists(output_file_path):
                    results.extend(parse_mixmhcpred_results(output_file_path))
                else:
                    with open(stdout_file_name, "r") as f:
                        stdout = f.read().strip()
                    raise ValueError(
                        "MixMHCpred failed on allele '%s' with stdout '%s'" % (allele, stdout))
                remove(stdout_file_name)
        return BindingPredictionCollection(results)

def parse_mixmhcpred_results(filename):
    """
    Parses output files of MixMHCpred that are expected to look like:

        Peptide  Score_bestAllele BestAllele  %Rank_bestAllele  Score_A0201  %Rank_A0201
        MLDDFSAGA          0.182093      A0201               0.3     0.182093          0.3
        SPEGEETII         -0.655341      A0201              51.0    -0.655341         51.0
        ILDRIITNA          0.203906      A0201               0.3     0.203906          0.3

    Parameters
    ----------
    filename : str
    
    Returns list of BindingPrediction
    """
    df = pd.read_csv(filename, comment="#", sep="\t")
    binding_predictions = []
    for peptide, allele, score, pr in zip(
            df["Peptide"],
            df["BestAllele"],
            df["Score_bestAllele"],
            df["%Rank_bestAllele"]):
        binding_predictions.append(BindingPrediction(
            peptide=peptide,
            allele=normalize_allele_name(allele),
            score=score,
            percentile_rank=pr,
            prediction_method_name="mixmhcpred"))
    return binding_predictions
