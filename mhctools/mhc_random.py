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

import random

import pandas as pd

from peptide_binding_measure import IC50_FIELD_NAME, PERCENTILE_RANK_FIELD_NAME


def generate_scored_epitopes(mutated_regions, alleles = ['HLA-A*02:01']):
    records = []
    # if wer'e not running the MHC prediction then we have to manually
    # extract 9mer substrings
    for _, row in mutated_regions.iterrows():
        seq = row.SourceSequence
        epitope_length = 9
        for i in xrange(len(seq) - epitope_length + 1):
            for allele in alleles:
                record = {}
                record['Epitope'] = seq[i:i+epitope_length]
                record['EpitopeStart'] = i
                record['EpitopeEnd'] = i + epitope_length
                record['SourceSequence'] = seq
                record[PERCENTILE_RANK_FIELD_NAME] = random.randint(0,99)
                record[IC50_FIELD_NAME] = random.random() * 10000.0
                record['Allele'] = allele
                for k, v in row.iteritems():
                    if k not in record:
                        record[k] = v
                records.append(record)
    scored_epitopes = pd.DataFrame.from_records(records)
    return scored_epitopes
