import pandas as pd

from mhctools import NetMHCpan
from mhctools.common import normalize_hla_allele_name


DEFAULT_ALLELE = 'HLA-A*02:01'


# TODO(tavi) Test requires chr/pos/ref/alt to run.
# TODO(tavi) Reconcile MutationEnd vs. MutationStop.
def test_netmhc_pan():
    alleles = [normalize_hla_allele_name(DEFAULT_ALLELE)]
    pan_predictor = NetMHCpan(alleles)
    input_df = pd.DataFrame({
        'chr': ['1', '1', 'X', 'X'],
        'pos': [10, 10, 2000, 2000],
        'ref': ['A', 'A', '', ''],
        'alt': ['T', 'T', 'C', 'C'],
        'SourceSequence':  [
            'ASIINFKELA', 'ASIINFKELA',
            'ASILLLVFYW', 'ASILLLVFYW',
        ],
        'MutationStart': [3, 3, 5, 5],
        'MutationEnd': [4, 4, 6, 6],
        'GeneInfo': [None, None, None, None],
        'Gene': ['SMAD4', 'SMAD4', 'TP53', 'TP53'],
        'GeneMutationInfo': ['A>T', 'A>T', 'insC', 'insC'],
        'PeptideMutationInfo': ['L>I', 'L>I', 'fs', 'fs'],
        'TranscriptId': [
            'ENST00000528762', 'ENST00000528762',
            'ENST00000544455', 'ENST00000544455'
        ],
    })
    output = pan_predictor.predict(input_df)
