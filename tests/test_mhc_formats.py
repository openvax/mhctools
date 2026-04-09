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

from mhctools.parsing import (
  parse_netmhcpan28_stdout,
  parse_netmhcpan3_stdout,
  parse_netmhcpan_stdout,
  parse_netmhcpan_to_preds,
  parse_netmhc3_stdout,
  parse_netmhc4_stdout,
  parse_netmhciipan43_stdout,
  parse_netmhciipan4_stdout,
)
from mhctools.pred import Kind

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


def test_auto_detect_netmhcpan28():
    """Auto-detecting parser should produce identical results to parse_netmhcpan28_stdout."""
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
     11  HLA-A*02:03    HIIIASSSL   id0         0.515       189.74    4.00 <= WB
    """
    results = parse_netmhcpan_stdout(netmhcpan28_output)
    assert len(results) == 3
    for entry in results:
        assert entry.allele == 'HLA-A*02:03'
    last = [e for e in results if e.peptide == "HIIIASSSL"][0]
    assert last.value == 189.74
    assert last.percentile_rank == 4.00
    assert last.offset == 11


def test_auto_detect_netmhcpan3():
    """Auto-detecting parser should produce identical results to parse_netmhcpan3_stdout."""
    netmhcpan3_output = """
    # Rank Threshold for Strong binding peptides   0.500
    # Rank Threshold for Weak binding peptides   2.000
    -----------------------------------------------------------------------------------
      Pos          HLA         Peptide       Core Of Gp Gl Ip Il        Icore        Identity   Score Aff(nM)   %Rank  BindLevel
    -----------------------------------------------------------------------------------
        1  HLA-B*18:01        QQQQQYFP  QQQQQYFP-  0  0  0  8  1     QQQQQYFP             id0 0.06456 24866.4   17.00
        9  HLA-B*18:01        EITHIIAS  -EITHIIAS  0  0  0  0  1     EITHIIAS             id0 0.09323 18234.7   10.00
    """
    results = parse_netmhcpan_stdout(netmhcpan3_output)
    assert len(results) == 2
    for entry in results:
        assert entry.allele == 'HLA-B*18:01'
    first = results[0]
    assert first.peptide == "QQQQQYFP"
    assert first.offset == 0  # 1-based converted to 0-based
    assert first.percentile_rank == 17.00

    second = results[1]
    assert second.value == 18234.7
    assert second.percentile_rank == 10.00
    assert second.offset == 8  # 9 -> 8 (1-based to 0-based)


def test_auto_detect_netmhcpan4_ba():
    """Auto-detecting parser for NetMHCpan 4.0 binding affinity mode."""
    netmhcpan4_ba_output = """
# NetMHCpan version 4.0

# Input is in PEPTIDE format

# Make binding affinity predictions

HLA-A02:01 : Distance to training data  0.000 (using nearest neighbor HLA-A02:01)

# Rank Threshold for Strong binding peptides   0.500
# Rank Threshold for Weak binding peptides   2.000
-----------------------------------------------------------------------------------
  Pos          HLA         Peptide       Core Of Gp Gl Ip Il        Icore        Identity     Score Aff(nM)   %Rank  BindLevel
-----------------------------------------------------------------------------------
    1  HLA-A*02:01        SIINFEKL  SIINF-EKL  0  0  0  5  1     SIINFEKL         PEPLIST 0.1141340 14543.1 18.9860
-----------------------------------------------------------------------------------
"""
    results = parse_netmhcpan_stdout(netmhcpan4_ba_output)
    assert len(results) == 1
    entry = results[0]
    assert entry.allele == 'HLA-A*02:01'
    assert entry.peptide == "SIINFEKL"
    assert entry.offset == 0  # 1-based to 0-based
    assert abs(entry.score - 0.1141340) < 1e-6
    assert abs(entry.value - 14543.1) < 0.1
    assert abs(entry.percentile_rank - 18.9860) < 0.001


def test_auto_detect_netmhcpan4_el():
    """Auto-detecting parser for NetMHCpan 4.0 elution score mode (no Aff column)."""
    netmhcpan4_el_output = """
# NetMHCpan version 4.0

# Input is in PEPTIDE format

HLA-A02:01 : Distance to training data  0.000 (using nearest neighbor HLA-A02:01)

# Rank Threshold for Strong binding peptides   0.500
# Rank Threshold for Weak binding peptides   2.000
-----------------------------------------------------------------------------------
  Pos          HLA         Peptide       Core Of Gp Gl Ip Il        Icore        Identity     Score   %Rank  BindLevel
-----------------------------------------------------------------------------------
    1  HLA-A*02:01        SIINFEKL  SIINF-EKL  0  0  0  5  1     SIINFEKL         PEPLIST 0.3456780  5.1230
-----------------------------------------------------------------------------------
"""
    results = parse_netmhcpan_stdout(netmhcpan4_el_output)
    assert len(results) == 1
    entry = results[0]
    assert entry.allele == 'HLA-A*02:01'
    assert entry.peptide == "SIINFEKL"
    assert entry.offset == 0
    assert abs(entry.score - 0.3456780) < 1e-6
    assert entry.value is None  # no Aff(nM) column
    assert abs(entry.percentile_rank - 5.1230) < 0.001


def test_auto_detect_netmhcpan41_ba():
    """Auto-detecting parser for NetMHCpan 4.1 in binding affinity mode."""
    netmhcpan41_output = """
# NetMHCpan version 4.1b

# Input is in PEPTIDE format

# Make both EL and BA predictions

HLA-A02:01 : Distance to training data  0.000 (using nearest neighbor HLA-A02:01)

# Rank Threshold for Strong binding peptides   0.500
# Rank Threshold for Weak binding peptides   2.000
---------------------------------------------------------------------------------------------------------------------------
 Pos         MHC        Peptide      Core Of Gp Gl Ip Il        Icore        Identity  Score_EL %Rank_EL Score_BA %Rank_BA  Aff(nM) BindLevel
---------------------------------------------------------------------------------------------------------------------------
   1 HLA-A*02:01       SIINFEKL SII-NFEKL  0  0  0  3  1     SIINFEKL         PEPLIST 0.0100620    6.723 0.110414   20.171 15140.42
---------------------------------------------------------------------------------------------------------------------------
"""
    results = parse_netmhcpan_stdout(netmhcpan41_output, mode="binding_affinity")
    assert len(results) == 1
    entry = results[0]
    assert entry.allele == 'HLA-A*02:01'
    assert entry.peptide == "SIINFEKL"
    assert entry.offset == 0
    assert abs(entry.score - 0.110414) < 1e-5  # Score_BA
    assert abs(entry.percentile_rank - 20.171) < 0.01  # %Rank_BA
    assert abs(entry.value - 15140.42) < 0.1  # Aff(nM)


def test_auto_detect_netmhcpan41_el():
    """Auto-detecting parser for NetMHCpan 4.1 in elution score mode."""
    netmhcpan41_output = """
# NetMHCpan version 4.1b

# Make both EL and BA predictions

---------------------------------------------------------------------------------------------------------------------------
 Pos         MHC        Peptide      Core Of Gp Gl Ip Il        Icore        Identity  Score_EL %Rank_EL Score_BA %Rank_BA  Aff(nM) BindLevel
---------------------------------------------------------------------------------------------------------------------------
   1 HLA-A*02:01       SIINFEKL SII-NFEKL  0  0  0  3  1     SIINFEKL         PEPLIST 0.0100620    6.723 0.110414   20.171 15140.42
---------------------------------------------------------------------------------------------------------------------------
"""
    results = parse_netmhcpan_stdout(netmhcpan41_output, mode="elution_score")
    assert len(results) == 1
    entry = results[0]
    assert abs(entry.score - 0.0100620) < 1e-6  # Score_EL
    assert abs(entry.percentile_rank - 6.723) < 0.01  # %Rank_EL
    assert entry.value is None  # ic50 not used in EL mode


# --- Tests for parse_netmhcpan_to_preds (new-style Pred output) ---

def test_to_preds_netmhcpan28():
    """New-style parser returns Pred objects for NetMHCpan 2.8."""
    output = """
    ---------------------------------------------------x
    pos  HLA  peptide  Identity 1-log50k(aff) Affinity(nM)    %Rank  BindLevel
    ----------------------------------------------------------------------------
      0  HLA-A*02:03    QQQQQYFPE   id0         0.024     38534.25   50.00
     11  HLA-A*02:03    HIIIASSSL   id0         0.515       189.74    4.00 <= WB
    """
    preds = parse_netmhcpan_to_preds(output)
    assert len(preds) == 2
    for p in preds:
        assert p.kind == Kind.pMHC_affinity
        assert p.allele == "HLA-A*02:03"
    last = [p for p in preds if p.peptide == "HIIIASSSL"][0]
    assert abs(last.value - 189.74) < 0.01
    assert last.percentile_rank == 4.00
    assert last.score > 0.4  # 1 - log(189.74)/log(50000)


def test_to_preds_netmhcpan41_emits_both_kinds():
    """NetMHCpan 4.1 emits both affinity and presentation Preds per row."""
    output = """
# NetMHCpan version 4.1b

# Make both EL and BA predictions

---------------------------------------------------------------------------------------------------------------------------
 Pos         MHC        Peptide      Core Of Gp Gl Ip Il        Icore        Identity  Score_EL %Rank_EL Score_BA %Rank_BA  Aff(nM) BindLevel
---------------------------------------------------------------------------------------------------------------------------
   1 HLA-A*02:01       SIINFEKL SII-NFEKL  0  0  0  3  1     SIINFEKL         PEPLIST 0.0100620    6.723 0.110414   20.171 15140.42
---------------------------------------------------------------------------------------------------------------------------
"""
    preds = parse_netmhcpan_to_preds(output)
    assert len(preds) == 2  # one affinity + one presentation

    affinity = [p for p in preds if p.kind == Kind.pMHC_affinity]
    presentation = [p for p in preds if p.kind == Kind.pMHC_presentation]
    assert len(affinity) == 1
    assert len(presentation) == 1

    aff = affinity[0]
    assert aff.allele == "HLA-A*02:01"
    assert aff.peptide == "SIINFEKL"
    assert abs(aff.value - 15140.42) < 0.1  # IC50 in nM
    assert abs(aff.percentile_rank - 20.171) < 0.01
    assert aff.score > 0  # higher-is-better transform of IC50

    pres = presentation[0]
    assert abs(pres.score - 0.0100620) < 1e-6
    assert abs(pres.percentile_rank - 6.723) < 0.01
    assert pres.value is None  # presentation has no native-unit value


def test_to_preds_netmhcpan4_el():
    """NetMHCpan 4.0 EL mode (no Aff column) returns presentation Preds."""
    output = """
# NetMHCpan version 4.0
-----------------------------------------------------------------------------------
  Pos          HLA         Peptide       Core Of Gp Gl Ip Il        Icore        Identity     Score   %Rank  BindLevel
-----------------------------------------------------------------------------------
    1  HLA-A*02:01        SIINFEKL  SIINF-EKL  0  0  0  5  1     SIINFEKL         PEPLIST 0.3456780  5.1230
-----------------------------------------------------------------------------------
"""
    preds = parse_netmhcpan_to_preds(output)
    assert len(preds) == 1
    p = preds[0]
    assert p.kind == Kind.pMHC_presentation  # no Aff(nM) → EL mode
    assert abs(p.score - 0.3456780) < 1e-6


def test_to_preds_auto_detects_version():
    """Predictor version auto-detected from stdout."""
    output = """
# NetMHCpan version 4.1b
---------------------------------------------------------------------------------------------------------------------------
 Pos         MHC        Peptide      Core Of Gp Gl Ip Il        Icore        Identity  Score_EL %Rank_EL Score_BA %Rank_BA  Aff(nM) BindLevel
---------------------------------------------------------------------------------------------------------------------------
   1 HLA-A*02:01       SIINFEKL SII-NFEKL  0  0  0  3  1     SIINFEKL         PEPLIST 0.0100620    6.723 0.110414   20.171 15140.42
---------------------------------------------------------------------------------------------------------------------------
"""
    preds = parse_netmhcpan_to_preds(output)
    assert preds[0].predictor_version == "4.1b"


def test_to_preds_self_contained():
    """Each Pred carries full context — peptide, allele, source, offset."""
    output = """
    ---------------------------------------------------x
    pos  HLA  peptide  Identity 1-log50k(aff) Affinity(nM)    %Rank  BindLevel
    ----------------------------------------------------------------------------
      5  HLA-A*02:03    YFPEITHII   id0         0.231      4123.85   15.00
    """
    preds = parse_netmhcpan_to_preds(output, predictor_name="netMHCpan", predictor_version="2.8")
    p = preds[0]
    assert p.peptide == "YFPEITHII"
    assert p.allele == "HLA-A*02:03"
    assert p.source_sequence_name == "id0"
    assert p.offset == 5
    assert p.predictor_name == "netMHCpan"
    assert p.predictor_version == "2.8"


# ---------- NetMHCIIpan 4.3 ----------

def test_netmhciipan43_el_with_bind_level():
    """Rows with <=SB BindLevel should parse without error (GH-169)."""
    output = """
    --------------------------------------------------------------------------------------------------------------------------------------------
     Pos               MHC              Peptide   Of        Core  Core_Rel Inverted        Identity      Score_EL %Rank_EL  Exp_Bind      Score_BA %Rank_BA  Affinity(nM)  BindLevel
    --------------------------------------------------------------------------------------------------------------------------------------------
       1         DRB1_0101      GAATVAAGAATTAAG    4   VAAGAATTA     0.920        0        Sequence      0.164981     6.81     0.000      0.514261    17.98        191.63
       2         DRB1_0101      KSVPLEMLLINLTTI    4   LEMLLINLT     0.980        0        Sequence      0.807346     0.50     0.560      0.662093     4.95         38.71 <=SB
    """
    results = parse_netmhciipan43_stdout(output, mode="elution_score")
    assert len(results) == 2
    assert results[0].peptide == "GAATVAAGAATTAAG"
    assert abs(results[0].percentile_rank - 6.81) < 0.01
    assert results[1].peptide == "KSVPLEMLLINLTTI"
    assert abs(results[1].percentile_rank - 0.50) < 0.01


def test_netmhciipan43_ba_with_bind_level():
    """BA mode should parse rows with <=WB BindLevel (GH-169)."""
    output = """
    --------------------------------------------------------------------------------------------------------------------------------------------
     Pos               MHC              Peptide   Of        Core  Core_Rel Inverted        Identity      Score_EL %Rank_EL  Exp_Bind      Score_BA %Rank_BA  Affinity(nM)  BindLevel
    --------------------------------------------------------------------------------------------------------------------------------------------
       1         DRB1_0101      GAATVAAGAATTAAG    4   VAAGAATTA     0.920        0        Sequence      0.164981     6.81     0.000      0.514261    17.98        191.63 <=WB
    """
    results = parse_netmhciipan43_stdout(output, mode="binding_affinity")
    assert len(results) == 1
    assert abs(results[0].affinity - 191.63) < 0.01
    assert abs(results[0].percentile_rank - 17.98) < 0.01


def test_netmhciipan43_exp_bind_na():
    """Exp_Bind = NA should not crash float conversion (GH-169)."""
    output = """
    --------------------------------------------------------------------------------------------------------------------------------------------
     Pos               MHC              Peptide   Of        Core  Core_Rel Inverted        Identity      Score_EL %Rank_EL  Exp_Bind      Score_BA %Rank_BA  Affinity(nM)  BindLevel
    --------------------------------------------------------------------------------------------------------------------------------------------
       1         DRB1_0101      PAPAPSWPLSSSVPS    4   PSWPLSSSV     0.327        0        Sequence      0.000857    79.79       NA      0.327674    54.35       1442.91
    """
    results = parse_netmhciipan43_stdout(output, mode="elution_score")
    assert len(results) == 1
    assert abs(results[0].score - 0.000857) < 0.0001


def test_netmhciipan43_two_token_bind_level():
    """'<= SB' (with space) is two tokens — both should be stripped (GH-169)."""
    output = """
    --------------------------------------------------------------------------------------------------------------------------------------------
     Pos               MHC              Peptide   Of        Core  Core_Rel Inverted        Identity      Score_EL %Rank_EL  Exp_Bind      Score_BA %Rank_BA  Affinity(nM)  BindLevel
    --------------------------------------------------------------------------------------------------------------------------------------------
       1         DRB1_0101      KSVPLEMLLINLTTI    4   LEMLLINLT     0.980        0        Sequence      0.807346     0.50     0.560      0.662093     4.95         38.71 <= SB
    """
    results = parse_netmhciipan43_stdout(output, mode="binding_affinity")
    assert len(results) == 1
    assert abs(results[0].affinity - 38.71) < 0.01


def test_netmhciipan4_with_bind_level():
    """v4.0 layout with BindLevel should also parse (GH-169)."""
    output = """
    --------------------------------------------------------------------------------------------------------------------------------------------
     Pos           MHC              Peptide   Of        Core  Core_Rel        Identity      Score_EL %Rank_EL Exp_Bind      Score_BA  Affinity(nM) %Rank_BA  BindLevel
    --------------------------------------------------------------------------------------------------------------------------------------------
       1     DRB1_0101      PAPAPSWPLSSSVPS    4   PSWPLSSSV     0.327            test      0.000857    79.79       NA      0.327674       1442.91    54.35
       2     DRB1_0101      GAATVAAGAATTAAG    4   VAAGAATTA     0.920            test      0.800000     0.50     0.000      0.514261        191.63    17.98 <=SB
    """
    results = parse_netmhciipan4_stdout(output, mode="elution_score")
    assert len(results) == 2
    assert abs(results[1].score - 0.800000) < 0.0001
