from mhctools.parsing import (
  parse_netmhcpan28_stdout,
  parse_netmhcpan3_stdout,
  parse_netmhc3_stdout,
  parse_netmhc4_stdout,
)

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
   0  CFTWNQMNL         0.113        14800                       seq4 HLA-A02:03
--------------------------------------------------------------------------------------------------
   0  SLYNTVATL         0.730           18         SB            seq5 HLA-A02:03
--------------------------------------------------------------------------------------------------
   0  SLYNTVATF         0.493          239         WB            seq6 HLA-A02:03
--------------------------------------------------------------------------------------------------
    """
    n_sequences = 3
    n_alleles = 2
    n_expected = n_alleles * n_sequences

    binding_predictions = parse_netmhc3_stdout(netmhc_output)
    assert len(binding_predictions) == n_expected, \
      "Wrong number of binding predictions: %d (expected %d)" % (
        len(binding_predictions), n_expected)
    for entry in binding_predictions:
        # make sure both allele's tables get parsed
        assert entry.allele in ('HLA-A*02:01', 'HLA-A*02:03'), entry
        # expect the HIV epitope SLYNTVATL to be a predicted binder for both
        # alleles
        if entry.peptide == "SLYNTVATL":
            assert entry.value < 100, entry

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
    n_sequences = 2
    n_alleles = 2
    n_expected = n_sequences * n_alleles
    binding_predictions = parse_netmhc4_stdout(netmhc_output)
    assert len(binding_predictions) == n_expected, \
      "Wrong number of binding predictions: %d (expected %d)" % (
        len(binding_predictions), n_expected)
    for entry in binding_predictions:
        # make sure both allele's tables get parsed
        assert entry.allele in ('HLA-A*02:01', 'HLA-A*02:02'), entry
        assert 0 < entry.value < 50000, entry
        # expect the epitope AEFGPWQTV to have high affinity for both
        # alleles
        if entry.peptide == "AEFGPWQTV":
            assert entry.value > 10000, entry

def test_mhcpan28_stdout():
    netmhcpan28_output = """
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
    binding_predictions = parse_netmhcpan28_stdout(netmhcpan28_output)
    assert len(binding_predictions) == 12, \
      "Expected 12 binding predictions but got %d" % (len(binding_predictions),)

    for entry in binding_predictions:
        assert entry.allele == 'HLA-A*02:03', \
            "Expected entry %s to have allele 'HLA-A*02:03'" % (entry,)
        if entry.peptide == "HIIIASSSL":
            # expect the epitopes to be sorted in increasing IC50
            assert entry.value == 189.74, entry
            assert entry.percentile_rank == 4.00, entry

def test_mhcpan3_stdout():
    netmhcpan3_output = """
    # Rank Threshold for Strong binding peptides   0.500
    # Rank Threshold for Weak binding peptides   2.000
    -----------------------------------------------------------------------------------
      Pos          HLA         Peptide       Core Of Gp Gl Ip Il        Icore        Identity   Score Aff(nM)   %Rank  BindLevel
    -----------------------------------------------------------------------------------
        1  HLA-B*18:01        QQQQQYFP  QQQQQYFP-  0  0  0  8  1     QQQQQYFP             id0 0.06456 24866.4   17.00
        2  HLA-B*18:01        QQQQYFPE  QQQQYFPE-  0  0  0  8  1     QQQQYFPE             id0 0.06446 24892.8   17.00
        3  HLA-B*18:01        QQQYFPEI  QQ-QYFPEI  0  0  0  2  1     QQQYFPEI             id0 0.06108 25819.2   18.00
        4  HLA-B*18:01        QQYFPEIT  QQYFPEIT-  0  0  0  8  1     QQYFPEIT             id0 0.04229 31642.1   29.00
        5  HLA-B*18:01        QYFPEITH  -QYFPEITH  0  0  0  0  1     QYFPEITH             id0 0.05316 28130.5   22.00
        6  HLA-B*18:01        YFPEITHI  Y-FPEITHI  0  0  0  1  1     YFPEITHI             id0 0.02576 37836.9   50.00
        7  HLA-B*18:01        FPEITHII  FP-EITHII  0  0  0  2  1     FPEITHII             id0 0.06199 25566.2   18.00
        8  HLA-B*18:01        PEITHIIA  PEITHIIA-  0  0  0  8  1     PEITHIIA             id0 0.06692 24239.3   16.00
        9  HLA-B*18:01        EITHIIAS  -EITHIIAS  0  0  0  0  1     EITHIIAS             id0 0.09323 18234.7   10.00
       10  HLA-B*18:01        ITHIIASS  ITHIIASS-  0  0  0  8  1     ITHIIASS             id0 0.01784 41223.5   70.00
       11  HLA-B*18:01        THIIASSS  THIIASSS-  0  0  0  8  1     THIIASSS             id0 0.03335 34856.1   38.00
       12  HLA-B*18:01        HIIASSSL  -HIIASSSL  0  0  0  0  1     HIIASSSL             id0 0.03049 35949.6   42.00
    """
    binding_predictions = parse_netmhcpan3_stdout(netmhcpan3_output)

    assert len(binding_predictions) == 12
    for entry in binding_predictions:
        assert entry.allele == 'HLA-B*18:01', \
            "Expected entry %s to have allele 'HLA-B*18:01'" % (entry,)
        if entry.peptide == "EITHIIAS":
            # expect the epitopes to be sorted in increasing IC50
            assert entry.value == 18234.7, entry
            assert entry.percentile_rank == 10.00, entry
