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

from .base_predictor import BasePredictor

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
        # create input file
        # create output file
        # run MixMHCpred
        # delete files
        results = parse_mixmhcpred_results(output_filename)
        return results

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
    for _, row in df.iterrows():
        binding_predictions.append(BindingPrediction(
            source_sequence_name=original_key,
            offset=offset,
            peptide=peptide,
            allele=normalize_allele_name(allele),
            affinity=ic50,
            elution_score=elution_score,
            percentile_rank=rank,
            log_affinity=log_ic50,
            prediction_method_name=prediction_method_name))
    return binding_predictions