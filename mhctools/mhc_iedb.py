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
from .common import seq_to_str, convert_str
from .peptide_binding_measure import (
        IC50_FIELD_NAME, PERCENTILE_RANK_FIELD_NAME
)

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

    def predict(self, peptides):
        """
        Given a dataframe with long amino acid sequences in the
        'SourceSequence' field, return an augmented dataframe
        with shorter k-mers in the 'Epitope' column and several
        columns of MHC binding predictions with names such as 'percentile_rank'
        """
        # take each mutated sequence in the dataframe
        # and general MHC binding scores for all k-mer substrings
        responses = {}
        for i, peptide in enumerate(peptides):
            for allele in self.alleles:
                key = (peptide, allele)
                if key not in responses:
                    request = self._get_iedb_request_params(peptide, allele)
                    logging.info(
                        "Calling IEDB (%s) with request %s",
                        self.url,
                        request)
                    response_df = _query_iedb(request, self.url)
                    print("after query", response_df)
                    response_df.rename(
                        columns={
                            'peptide': 'Epitope',
                            'length': 'EpitopeLength',
                            'start': 'EpitopeStart',
                            'end': 'EpitopeEnd',
                            'allele': 'Allele',
                        },
                        inplace=True)
                    response_df['EpitopeStart'] -= 1
                    response_df['EpitopeEnd'] -= 1
                    print("after rename", response_df)
                    responses[key] = response_df
                else:
                    logging.info(
                        "Already made predictions for peptide %s with allele %s",
                        peptide,
                        allele)

        # concatenating the responses makes a MultiIndex with two columns
        # - SourceSequence
        # - index of epitope from that sequence's IEDB call
        #
        # ...when we reset the index, we turn these into two columns
        # named 'level_0', and 'level_1'. We want to rename the former
        # and delete the latter.
        responses = pd.concat(responses).reset_index()
        print("concat", responses)
        responses['SourceSequence'] = responses['level_0']
        del responses['level_0']
        del responses['level_1']
        print("dropped levels", responses)

        # IEDB has inclusive end positions, change to exclusive
        responses['EpitopeEnd'] += 1

        assert 'ann_rank' in responses, \
            "Missing 'ann_rank' in %s" % (responses.head(),)
        responses[PERCENTILE_RANK_FIELD_NAME] = responses['ann_rank']

        assert 'ann_ic50' in responses, \
            "Missing 'ann_ic50' in %s" % (responses.head(),)
        responses[IC50_FIELD_NAME] = responses['ann_ic50']

        # instead of just building up a new dataframe I'm expliciting
        # dropping fields here to document what other information is available
        drop_fields = (
            'seq_num',
            'method',
            'ann_ic50',
            'ann_rank',
            'consensus_percentile_rank',
            'smm_ic50',
            'smm_rank',
            'comblib_sidney2008_score',
            'comblib_sidney2008_rank'
        )
        for field in drop_fields:
            if field in responses:
                responses = responses.drop(field, axis=1)
        print("Dropped fields", responses)
        # some of the MHC scores come back as all NaN so drop them
        responses = responses.dropna(axis=1, how='all')
        return responses

class IedbMhc1(IedbBasePredictor):
    def __init__(
            self,
            alleles,
            epitope_lengths=[9],
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
            epitope_lengths=[15],
            prediction_method=prediction_method,
            url=url)
