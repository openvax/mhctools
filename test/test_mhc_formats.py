from mhctools.file_formats import parse_netmhc_stdout
from mhctools.peptide_binding_measure import (
    IC50_FIELD_NAME,
    PERCENTILE_RANK_FIELD_NAME,
)

def test_mhc_stdout():
    s = """
    # Affinity Threshold for Strong binding peptides  50.000',
    # Affinity Threshold for Weak binding peptides 500.000',
    # Rank Threshold for Strong binding peptides   0.500',
    # Rank Threshold for Weak binding peptides   2.000',
    ---------------------------------------------------x
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
    class MutationEntry(object):
    	pass

    mutation_entry = MutationEntry()
    mutation_entry.SourceSequence = "QQQQQYFPEITHIIASSSL"
    mutation_entry.MutationStart = 2
    mutation_entry.MutationEnd = 3
    mutation_entry.GeneInfo = "TP53 missense"
    mutation_entry.Gene = "TP53"
    mutation_entry.GeneMutationInfo = "g.2 some mutation info"
    mutation_entry.PeptideMutationInfo = "p.2 T>Q"
    mutation_entry.TranscriptId = "TID0"
    mutation_entry.chr = 'X'
    mutation_entry.pos = 39393
    mutation_entry.ref = 'A'
    mutation_entry.alt = 'T'

    peptide_entries = {"id0": mutation_entry}

    rows = parse_netmhc_stdout(s, peptide_entries)

    assert len(rows) == 12

    for i in xrange(len(rows)):
        assert rows[i]['EpitopeStart'] == i
        assert rows[i]['Allele'] == 'HLA-A*02:03'

    assert rows[0][IC50_FIELD_NAME] == 38534.25
    assert rows[0][PERCENTILE_RANK_FIELD_NAME] == 50.00
