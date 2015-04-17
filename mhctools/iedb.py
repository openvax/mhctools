# Copyright (c) 2014. Mount Sinai School of Medicine
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

import io
import urllib2
import urllib
import logging

import pandas as pd

from .base_predictor import BasePredictor
from .binding_measure import ic50_nM
from .binding_prediction import BindingPrediction
from .epitope_collection import EpitopeCollection
from .common import seq_to_str
"""
A note about prediction methods, copied from the IEDB website:

The prediction method list box allows choosing from a number of MHC class I
binding prediction methods:
- Artificial neural network (ANN) (is that just NetMHC?)
- Stabilized matrix method (SMM)
- SMM with a Peptide:MHC Binding Energy Covariance matrix (SMMPMBEC),
- NetMHCpan
- NetMHCcons
Excluded because IEDB results lack unique IC50 & Percentile Rank columns:
- Scoring Matrices from Combinatorial Peptide Libraries (Comblib_Sidney2008),
- Consensus,
"""

VALID_CLASS_I_METHODS = [
    "netmhccons",
    "netmhcpan",
    "ann",
    "smmpmbec",
    "smm",
    "pickpocket"
]

VALID_CLASS_II_METHODS = [
    "NetMHCIIpan",
    "nn_align",
    "smm_align",
]

def _parse_iedb_response(response):
    """Take the binding predictions returned by IEDB's web API
    and parse them into a DataFrame

    Expect response to look like:
    allele  seq_num start   end length  peptide ic50    rank
    HLA-A*01:01 1   2   10  9   LYNTVATLY   2145.70 3.7
    HLA-A*01:01 1   5   13  9   TVATLYCVH   2216.49 3.9
    HLA-A*01:01 1   7   15  9   ATLYCVHQR   2635.42 5.1
    HLA-A*01:01 1   4   12  9   NTVATLYCV   6829.04 20
    HLA-A*01:01 1   1   9   9   SLYNTVATL   8032.38 24
    HLA-A*01:01 1   8   16  9   TLYCVHQRI   8853.90 26
    HLA-A*01:01 1   3   11  9   YNTVATLYC   9865.62 29
    HLA-A*01:01 1   6   14  9   VATLYCVHQ   27575.71    58
    HLA-A*01:01 1   10  18  9   YCVHQRIDV   48929.64    74
    HLA-A*01:01 1   9   17  9   LYCVHQRID   50000.00    75
    """
    if len(response) == 0:
        raise ValueError("Empty response from IEDB!")
    df = pd.read_csv(io.BytesIO(response), delim_whitespace=True, header=0)
    if len(df) == 0:
        raise ValueError(
            "No binding predictions in response from IEDB: %s" % (response,))
    required_columns = [
        "allele",
        "peptide",
        "ic50",
        "rank",
        "start",
        "end",
    ]
    for column in required_columns:
        if column not in df.columns:
            raise ValueError(
                "Response from IEDB is missing '%s' column: %s" % (
                    column,
                    df.ix[0],))
    return df

def _query_iedb(request_values, url):
    """
    Call into IEDB's web API for MHC binding prediction using request dictionary
    with fields:
        - "method"
        - "length"
        - "sequence_text"
        - "allele"

    Parse the response into a DataFrame.
    """
    data = urllib.urlencode(request_values)
    req = urllib2.Request(url, data)
    response = urllib2.urlopen(req).read()
    return _parse_iedb_response(response)


class IedbBasePredictor(BasePredictor):

    def __init__(
            self,
            alleles,
            epitope_lengths,
            prediction_method,
            url):
        BasePredictor.__init__(
            self,
            alleles=alleles,
            epitope_lengths=epitope_lengths)
        self.prediction_method = prediction_method

        if not isinstance(url, str):
            raise TypeError("Expected URL to be string, not %s : %s" % (
                url, type(url)))
        self.url = url

    def _get_iedb_request_params(self, sequence, allele):

        params = {
            "method": seq_to_str(self.prediction_method),
            "length": seq_to_str(self.epitope_lengths),
            "sequence_text": sequence,
            # have to repeat allele for each length
            "allele": ",".join([allele] * len(self.epitope_lengths)),
        }
        return params

    def predict(self, fasta_dictionary):
        """Given a dictionary mapping unique keys to amino acid sequences,
        run MHC binding predictions on all candidate epitopes extracted from
        sequences and return a EpitopeCollection.

        Parameters
        ----------
        fasta_dictionary : dict
            Mapping of protein identifiers to protein amino acid sequences
        """
        # take each mutated sequence in the dataframe
        # and general MHC binding scores for all k-mer substrings
        binding_predictions = []
        for key, amino_acid_sequence in fasta_dictionary.items():
            for allele in self.alleles:
                request = self._get_iedb_request_params(
                    amino_acid_sequence, allele)
                logging.info(
                    "Calling IEDB (%s) with request %s",
                    self.url,
                    request)
                response_df = _query_iedb(request, self.url)
                binding_predictions.extend([
                    BindingPrediction(
                        allele=row['allele'],
                        peptide=row['peptide'],
                        length=row['length'],
                        source_sequence=amino_acid_sequence,
                        source_sequence_key=key,
                        # IEDB's start is 1-based, need to subtract 1
                        base0_start=row['start'] - 1,
                        base_end=row['end'] - 1,
                        value=row['ic50'],
                        percentile_rank=row['rank'],
                        measure=ic50_nM,
                        prediction_method_name=self.prediction_method)
                    for (_, row)
                    in response_df.iterrows()
                ])
        return EpitopeCollection(binding_predictions)

class IedbMhc1(IedbBasePredictor):
    def __init__(
            self,
            alleles,
            epitope_lengths=[8, 9, 10, 11],
            prediction_method="netmhccons",
            url="http://tools-api.iedb.org/tools_api/mhci/"):
        if prediction_method not in VALID_CLASS_I_METHODS:
            raise ValueError(
                "Invalid IEDB MHC class I binding prediction method: %s" % (
                    prediction_method,))
        IedbBasePredictor.__init__(
            self,
            alleles=alleles,
            epitope_lengths=epitope_lengths,
            prediction_method=prediction_method,
            url=url)

class IedbMhc2(IedbBasePredictor):
    def __init__(self,
            alleles,
            epitope_lengths=[15, 16, 17, 18, 19, 20],
            prediction_method="NetMHCIIpan",
            url="http://tools-api.iedb.org/tools_api/mhcii/"):
        if prediction_method not in VALID_CLASS_II_METHODS:
            raise ValueError(
                "Invalid IEDB MHC class II binding prediction method: %s" % (
                    prediction_method,))
        IedbBasePredictor.__init__(
            self,
            alleles=alleles,
            # only epitope lengths of 15 currently supported by IEDB's web API
            epitope_lengths=epitope_lengths,
            prediction_method=prediction_method,
            url=url)