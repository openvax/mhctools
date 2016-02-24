from mhctools.file_formats import parse_netmhcpan_stdout, parse_netmhc3_stdout, parse_netmhc4_stdout

def test_netmhc3_stdout():
    """
    Test parsing of NetMHC output of predictions of HLA-A*02:01
    and HLA-A*02:03 for three epitopes:
      - CFTWNQMNL
      - SLYNTVATL
      - SLYNTVATF
    """

    netmhc_output = """
NetMHC version 3.4. 9mer predictions using Artificial Neural Networks - Direct. Allele HLA-A02:01.
Strong binder threshold  50 nM. Weak binder threshold score 500 nM



----------------------------------------------------------------------------------------------------
 pos    peptide      logscore affinity(nM) Bind Level    Protein Name     Allele
----------------------------------------------------------------------------------------------------
   0  CFTWNQMNL         0.085        19899                       seq4 HLA-A02:01
--------------------------------------------------------------------------------------------------
   0  SLYNTVATL         0.579           94         WB            seq5 HLA-A02:01
--------------------------------------------------------------------------------------------------
   0  SLYNTVATF         0.289         2186                       seq6 HLA-A02:01
--------------------------------------------------------------------------------------------------



NetMHC version 3.4. 9mer predictions using Artificial Neural Networks - Direct. Allele HLA-A02:03.
Strong binder threshold  50 nM. Weak binder threshold score 500 nM



----------------------------------------------------------------------------------------------------
 pos    peptide      logscore affinity(nM) Bind Level    Protein Name     Allele
----------------------------------------------------------------------------------------------------
   0  SB                0.113        14800                       seq4 HLA-A02:03
--------------------------------------------------------------------------------------------------
   0  SLYNTVATL         0.730           18         SB            seq5 HLA-A02:03
--------------------------------------------------------------------------------------------------
   0  SLYNTVATF         0.493          239         WB            seq6 HLA-A02:03
--------------------------------------------------------------------------------------------------
    """
    fasta_dictionary = {
      "seq4": "CFTWNQMNL",
      "seq5": "SLYNTVATL",
      "seq6": "SLYNTVATF",
    }
    epitope_collection = parse_netmhc3_stdout(
      netmhc_output,
      fasta_dictionary=fasta_dictionary,
      prediction_method_name="netmhc3")
    assert len(epitope_collection) == 2 * len(fasta_dictionary), \
      "Wrong number of binding predictions: %d (expected %d)" % (
        len(epitope_collection), 2 * len(fasta_dictionary))
    for i, entry in enumerate(epitope_collection):
        # make sure both allele's tables get parsed
        assert entry.allele in ('HLA-A*02:01', 'HLA-A*02:03')
        # expect the HIV epitope SLYNTVATL to be a predicted binder for both
        # alleles
        if entry.peptide == "SLYNTVATL":
            assert entry.value < 100

def test_netmhc4_stdout():
    netmhc_output = """
# NetMHC version 4.0

# Read 132 elements on pairlist /Users/tavi/drive/work/repos/cancer/n-b/netMHC-4.0/Darwin_x86_64/data/allelelist
# Input is in PEPTIDE format
# Rank Threshold for Strong binding peptides   0.500
# Rank Threshold for Weak binding peptides   2.000
-----------------------------------------------------------------------------------
  pos          HLA      peptide         Core Offset  I_pos  I_len  D_pos  D_len        iCore        Identity 1-log50k(aff) Affinity(nM)    %Rank  BindLevel
-----------------------------------------------------------------------------------
    0    HLA-A0201  AAAAAWYLWEV    AAAWYLWEV      0      0      0      1      2  AAAAAWYLWEV         SEQ_A           0.349      1147.39     4.50
    0    HLA-A0201    AEFGPWQTV    AEFGPWQTV      0      0      0      0      0    AEFGPWQTV         SEQ_B           0.129     12361.73    18.00
-----------------------------------------------------------------------------------

Protein PEPLIST. Allele HLA-A0201. Number of high binders 0. Number of weak binders 0. Number of peptides 10

-----------------------------------------------------------------------------------
# Rank Threshold for Strong binding peptides   0.500
# Rank Threshold for Weak binding peptides   2.000
-----------------------------------------------------------------------------------
  pos          HLA      peptide         Core Offset  I_pos  I_len  D_pos  D_len        iCore        Identity 1-log50k(aff) Affinity(nM)    %Rank  BindLevel
-----------------------------------------------------------------------------------
    0    HLA-A0202    AEFGPWQTV    AEFGPWQTV      0      0      0      0      0    AEFGPWQTV         SEQ_C           0.136     11437.51    23.00
  219    HLA-A0202    QLLRDNLTL    QLLRDNLTL      0      0      0      0      0    QLLRDNLTL         SEQ_D           0.527       167.10     1.50 <= WB
-----------------------------------------------------------------------------------

Protein PEPLIST. Allele HLA-A0202. Number of high binders 0. Number of weak binders 0. Number of peptides 10

-----------------------------------------------------------------------------------
"""
    fasta_dictionary = {
        "SEQ_A": "AAAAAWYLWEV",
        "SEQ_B": "AEFGPWQTV",
        "SEQ_C": "AEFGPWQTV",
        "SEQ_D": "QLLRDNLTL"
    }
    epitope_collection = parse_netmhc4_stdout(
      netmhc_output,
      fasta_dictionary=fasta_dictionary,
      prediction_method_name="netmhc4")
    assert len(epitope_collection) == len(fasta_dictionary), \
      "Wrong number of binding predictions: %d (expected %d)" % (
        len(epitope_collection), len(fasta_dictionary))
    for i, entry in enumerate(epitope_collection):
        # make sure both allele's tables get parsed
        assert entry.allele in ('HLA-A*02:01', 'HLA-A*02:02')
        assert 0 < entry.value < 50000
        # expect the epitope AEFGPWQTV to have high affinity for both
        # alleles
        if entry.peptide == "AEFGPWQTV":
            assert entry.value > 10000

def test_mhcpan_stdout():
    netmhcpan_output = """
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

    fasta_dictionary = {
      "id0": "QQQQQYFPEITHIIASSSL"
    }

    epitope_collection = parse_netmhcpan_stdout(
      netmhcpan_output,
      fasta_dictionary=fasta_dictionary,
      prediction_method_name="netmhcpan")

    assert len(epitope_collection) == 12

    for i, entry in enumerate(epitope_collection):
        assert entry.allele == 'HLA-A*02:03', \
            "Expected entry %s to have allele 'HLA-A*02:03'" % (entry,)
        if i == 0:
            # expect the epitopes to be sorted in increasing IC50
            assert entry.value == 189.74
            assert entry.percentile_rank == 4.00
