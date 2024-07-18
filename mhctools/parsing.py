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

import numpy as np

from mhcnames import normalize_allele_name

from .binding_prediction import BindingPrediction


NETMHC_TOKENS = {
    "pos",
    "Pos",
    "Seq",
    "Number",
    "Protein",
    "Allele",
    "NetMHC",
    "Strong",
}

def check_stdout_error(stdout, program_name):
    if "ERROR" in stdout.upper():
        # if NetMHC* failed with an error then let's pull out the error
        # message line and raise an exception with it
        error_index = stdout.upper().index("ERROR")
        stdout_after_error = stdout[error_index:]
        error_line = stdout_after_error.split("\n")[0]
        raise ValueError("%s failed - %s" % (program_name, error_line))

def split_stdout_lines(stdout):
    """
    Given the standard output from NetMHC/NetMHCpan/NetMHCcons tools,
    drop all {comments, lines of hyphens, empty lines} and split the
    remaining lines by whitespace.
    """
    # all the NetMHC formats use lines full of dashes before any actual
    # binding results
    seen_dash = False
    for l in stdout.split("\n"):
        l = l.strip()
        # wait for a line like '----------' before trying to parse entries
        # have to include multiple dashes here since NetMHC 4.0 sometimes
        # gives negative positions in its "peptide" input mode
        if l.startswith("---"):
            seen_dash = True
            continue
        if not seen_dash:
            continue
        # ignore empty lines and comments
        if not l or l.startswith("#"):
            continue
        # beginning of headers in NetMHC
        if any(l.startswith(word) for word in NETMHC_TOKENS):
            continue
        yield l.split()


def clean_fields(fields, ignored_value_indices, transforms):
    """
    Sometimes, NetMHC* has fields that are only populated sometimes, which results
    in different count/indexing of the fields when that happens.

    We handle this by looking for particular strings at particular indices, and
    deleting them.

    Warning: this may result in unexpected behavior sometimes. For example, we
    ignore "SB" and "WB" for NetMHC 3.x output; which also means that any line
    with a key called SB or WB will be ignored.

    Also, sometimes NetMHC* will have fields that we want to modify in some
    consistent way, e.g. NetMHCpan3 has 1-based offsets and all other predictors
    have 0-based offsets (and we rely on 0-based offsets). We handle this using
    a map from field index to transform function.
    """
    cleaned_fields = []
    for i, field in enumerate(fields):
        if field in ignored_value_indices:
            ignored_index = ignored_value_indices[field]

            # Is the value we want to ignore at the index where we'd ignore it?
            if ignored_index == i:
                continue

        # transform this field if the index is in transforms, otherwise leave alone
        cleaned_field = transforms[i](field) if i in transforms else field
        cleaned_fields.append(cleaned_field)
    return cleaned_fields

def valid_affinity(x):
    """
    Check that an IC50 affinity value is valid.

    Parameters
    ----------
    x : float

    Returns
    -------
    bool
    """
    if x is None:
        return False
    if np.isnan(x) or np.isinf(x):
        return False
    return x >= 0

def parse_stdout(
        stdout,
        prediction_method_name,
        sequence_key_mapping,
        key_index,
        offset_index,
        peptide_index,
        allele_index,
        score_index,
        rank_index=None,
        ic50_index=None,

        ignored_value_indices={},
        transforms={}):
    """
    Generic function for parsing any NetMHC* output, given expected indices
    of values of interest.

    Parameters
    ----------
    stdout : str

    prediction_method_name : str

    key_index : int

    offset_index : int

    peptide_index : int

    allele_index : int

    score_index : int

    rank_index : int

    ic50_index : int

    sequence_key_mapping : dict
        Dictionary mapping sequence names (which might be hashes or truncated)
        to the sequence names which should be used in the parsed
        BindingPrediction objects

    ignored_value_indices : dict
        Map from values to the positions we'll ignore them at.
        See clean_fields.

    transforms  : dict
        Map from field index to a transform function to be applied to values in
        that field. See clean_fields.

    Returns BindingPredictionCollection
    """

    binding_predictions = []
    for fields in split_stdout_lines(stdout):
        fields = clean_fields(fields, ignored_value_indices, transforms)

        offset = int(fields[offset_index])
        peptide = str(fields[peptide_index])
        allele = str(fields[allele_index])

        if score_index is None:
            score = None
        else:
            score = float(fields[score_index])

        if rank_index is None:
            rank = None
        else:
            rank = float(fields[rank_index])

        if ic50_index is None:
            ic50 = None
        else:
            ic50 = float(fields[ic50_index])

        key = str(fields[key_index])
        if sequence_key_mapping:
            original_key = sequence_key_mapping[key]
        else:
            # if sequence_key_mapping isn't provided then let's assume it's the
            # identity function
            original_key = key

        # if we have a bad IC50 score we might still get a salvageable
        # log of the score. Strangely, this is necessary sometimes!
        if ic50_index is not None and (not valid_affinity(ic50)) and np.isfinite(score):
            # pylint: disable=invalid-unary-operand-type
            ic50 = 50000 ** (1 - score)

        binding_predictions.append(BindingPrediction(
            source_sequence_name=original_key,
            offset=offset,
            peptide=peptide,
            allele=normalize_allele_name(allele),
            score=score,
            affinity=ic50,
            percentile_rank=rank,
            prediction_method_name=prediction_method_name))
    return binding_predictions

def parse_netmhc3_stdout(
        stdout,
        prediction_method_name="netmhc3",
        sequence_key_mapping=None):
    """
    Parse the output format for NetMHC 3.x, which looks like:

    ----------------------------------------------------------------------------------------------------
    pos    peptide      logscore affinity(nM) Bind Level    Protein Name     Allele
    ----------------------------------------------------------------------------------------------------
    0  SIINKFELL         0.437          441         WB              A1 HLA-A02:01
    --------------------------------------------------------------------------------------------------
    0  SIINKFFFQ         0.206         5411                         A2 HLA-A02:01
    1  IINKFFFQQ         0.128        12544                         A2 HLA-A02:01
    2  INKFFFQQQ         0.046        30406                         A2 HLA-A02:01
    3  NKFFFQQQQ         0.050        29197                         A2 HLA-A02:01
    --------------------------------------------------------------------------------------------------
    """
    return parse_stdout(
        stdout=stdout,
        prediction_method_name=prediction_method_name,
        sequence_key_mapping=sequence_key_mapping,
        key_index=4,
        offset_index=0,
        peptide_index=1,
        allele_index=5,
        score_index=2,
        ic50_index=3,
        rank_index=None,
        ignored_value_indices={"WB": 4, "SB": 4})

def parse_netmhc4_stdout(
        stdout,
        prediction_method_name="netmhc4",
        sequence_key_mapping=None):
    """
    # Peptide length 9
    # Rank Threshold for Strong binding peptides   0.500
    # Rank Threshold for Weak binding peptides   2.000
    -----------------------------------------------------------------------------------
      pos          HLA      peptide         Core Offset  I_pos  I_len  D_pos  D_len        iCore        Identity 1-log50k(aff) Affinity(nM)    %Rank  BindLevel
    -----------------------------------------------------------------------------------
        0    HLA-A0201    TMDKSELVQ    TMDKSELVQ      0      0      0      0      0    TMDKSELVQ 143B_BOVIN_P293         0.051     28676.59    43.00
        1    HLA-A0201    MDKSELVQK    MDKSELVQK      0      0      0      0      0    MDKSELVQK 143B_BOVIN_P293         0.030     36155.15    70.00
        2    HLA-A0201    DKSELVQKA    DKSELVQKA      0      0      0      0      0    DKSELVQKA 143B_BOVIN_P293         0.030     36188.42    70.00
        3    HLA-A0201    KSELVQKAK    KSELVQKAK      0      0      0      0      0    KSELVQKAK 143B_BOVIN_P293         0.032     35203.22    65.00
        4    HLA-A0201    SELVQKAKL    SELVQKAKL      0      0      0      0      0    SELVQKAKL 143B_BOVIN_P293         0.031     35670.99    65.00
        5    HLA-A0201    ELVQKAKLA    ELVQKAKLA      0      0      0      0      0    ELVQKAKLA 143B_BOVIN_P293         0.080     21113.07    29.00
        6    HLA-A0201    LVQKAKLAE    LVQKAKLAE      0      0      0      0      0    LVQKAKLAE 143B_BOVIN_P293         0.027     37257.56    75.00
        7    HLA-A0201    VQKAKLAEQ    VQKAKLAEQ      0      0      0      0      0    VQKAKLAEQ 143B_BOVIN_P293         0.040     32404.62    55.00
      219    HLA-A0201    QLLRDNLTL    QLLRDNLTL      0      0      0      0      0    QLLRDNLTL 143B_BOVIN_P293         0.527       167.10     1.50 <= WB
    -----------------------------------------------------------------------------------
    """
    return parse_stdout(
        stdout=stdout,
        prediction_method_name=prediction_method_name,
        sequence_key_mapping=sequence_key_mapping,
        key_index=10,
        offset_index=0,
        peptide_index=2,
        allele_index=1,
        ic50_index=12,
        rank_index=13,
        score_index=11)

def parse_netmhc41_stdout(
        stdout,
        mode,
        prediction_method_name="netmhc41",
        sequence_key_mapping=None):
    """
    # NetMHCpan version 4.1b

    # Tmpdir made /tmp/netMHCpanGjkOCb
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

    Protein PEPLIST. Allele HLA-A*02:01. Number of high binders 0. Number of weak binders 0. Number of peptides 1
    """
    assert mode in ("binding_affinity", "elution_score")
    return parse_stdout(
        stdout=stdout,
        prediction_method_name=prediction_method_name,
        sequence_key_mapping=sequence_key_mapping,
        key_index=10,
        offset_index=0,
        peptide_index=2,
        allele_index=1,
        ic50_index=15 if mode == "binding_affinity" else None,
        rank_index=14 if mode == "binding_affinity" else 12,
        score_index=13 if mode == "binding_affinity" else 11)

def parse_netmhcpan28_stdout(
        stdout,
        prediction_method_name="netmhcpan",
        sequence_key_mapping=None):
    """
    # Affinity Threshold for Strong binding peptides  50.000',
    # Affinity Threshold for Weak binding peptides 500.000',
    # Rank Threshold for Strong binding peptides   0.500',
    # Rank Threshold for Weak binding peptides   2.000',
    ----------------------------------------------------------------------------
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
    check_stdout_error(stdout, "NetMHCpan-2.8")
    return parse_stdout(
        stdout=stdout,
        prediction_method_name=prediction_method_name,
        sequence_key_mapping=sequence_key_mapping,
        key_index=3,
        offset_index=0,
        peptide_index=2,
        allele_index=1,
        ic50_index=5,
        rank_index=6,
        score_index=4)

def parse_netmhcpan3_stdout(
        stdout,
        prediction_method_name="netmhcpan",
        sequence_key_mapping=None):
    """
    # Rank Threshold for Strong binding peptides   0.500
    # Rank Threshold for Weak binding peptides   2.000
    -----------------------------------------------------------------------------------
    Pos          HLA         Peptide       Core Of Gp Gl Ip Il        Icore        Identity   Score Aff(nM)   %Rank  BindLevel
    -----------------------------------------------------------------------------------
    1  HLA-B*18:01        MFCQLAKT  MFCQLAKT-  0  0  0  8  1     MFCQLAKT     sequence0_0 0.02864 36676.0   45.00
    2  HLA-B*18:01        FCQLAKTY  F-CQLAKTY  0  0  0  1  1     FCQLAKTY     sequence0_0 0.07993 21056.5   13.00
    """

    # the offset specified in "pos" (at index 0) is 1-based instead of 0-based. we adjust it to be
    # 0-based, as in all the other netmhc predictors supported by this library.
    transforms = {
        0: lambda x: int(x) - 1,
    }
    return parse_stdout(
        stdout=stdout,
        prediction_method_name=prediction_method_name,
        sequence_key_mapping=sequence_key_mapping,
        key_index=10,
        offset_index=0,
        peptide_index=2,
        allele_index=1,
        ic50_index=12,
        rank_index=13,
        score_index=11,
        transforms=transforms)


def parse_netmhcpan4_stdout(
        stdout,
        prediction_method_name="netmhcpan",
        sequence_key_mapping=None,
        mode="binding_affinity"):
    """
    # NetMHCpan version 4.0

    # Tmpdir made /var/folders/jc/fyrvcrcs3sb8g4mkdg6nl_t80000gp/T//netMHCpanuH3SvY
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

    Protein PEPLIST. Allele HLA-A*02:01. Number of high binders 0. Number of weak binders 0. Number of peptides 1
    """

    # the offset specified in "pos" (at index 0) is 1-based instead of 0-based. we adjust it to be
    # 0-based, as in all the other netmhc predictors supported by this library.
    transforms = {
        0: lambda x: int(x) - 1,
    }
    return parse_stdout(
        stdout=stdout,
        prediction_method_name=prediction_method_name,
        sequence_key_mapping=sequence_key_mapping,
        key_index=10,
        offset_index=0,
        peptide_index=2,
        allele_index=1,
        score_index=11,
        ic50_index=None if mode == "elution_score" else 12,
        rank_index=12 if mode == "elution_score" else 13,
        transforms=transforms)

def parse_netmhcpan41_stdout(
        stdout,
        prediction_method_name="netmhcpan",
        sequence_key_mapping=None,
        mode="binding_affinity"):
    """
	NetMHCpan version 4.1b
	# Rank Threshold for Strong binding peptides   0.500
	# Rank Threshold for Weak binding peptides   2.000
	---------------------------------------------------------------------------------------------------------------------------
	 Pos         MHC        Peptide      Core Of Gp Gl Ip Il        Icore        Identity  Score_EL %Rank_EL Score_BA %Rank_BA  Aff(nM) BindLevel
	---------------------------------------------------------------------------------------------------------------------------
	   1 HLA-A*03:01    GKSGGGRCGGG GKSGGGRGG  0  7  2  0  0  GKSGGGRCGGG            seq1 0.0000000  100.000 0.009240   95.346 45243.03
	---------------------------------------------------------------------------------------------------------------------------

	Protein seq1. Allele HLA-A*03:01. Number of high binders 0. Number of weak binders 0. Number of peptides 1

	-----------------------------------------------------------------------------------
    """

    # the offset specified in "pos" (at index 0) is 1-based instead of 0-based. we adjust it to be
    # 0-based, as in all the other netmhc predictors supported by this library.
    transforms = {
        0: lambda x: int(x) - 1,
    }
    return parse_stdout(
        stdout=stdout,
        prediction_method_name=prediction_method_name,
        sequence_key_mapping=sequence_key_mapping,
        key_index=10,
        offset_index=0,
        peptide_index=2,
        allele_index=1,
        score_index=11 if mode == "elution_score" else 13,
        ic50_index=None if mode == "elution_score" else 15,
        rank_index=12 if mode == "elution_score" else 14,
        transforms=transforms)


def parse_netmhccons_stdout(
        stdout,
        prediction_method_name="netmhccons",
        sequence_key_mapping=None):
    """
    # Affinity Threshold for Strong binding peptides  50.000',
    # Affinity Threshold for Weak binding peptides 500.000',
    # Rank Threshold for Strong binding peptides   0.500',
    # Rank Threshold for Weak binding peptides   2.000',
    ----------------------------------------------------------------------------
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
    return parse_stdout(
        stdout=stdout,
        prediction_method_name=prediction_method_name,
        sequence_key_mapping=sequence_key_mapping,
        key_index=3,
        offset_index=0,
        peptide_index=2,
        allele_index=1,
        ic50_index=5,
        rank_index=6,
        score_index=4)

def parse_netmhciipan_stdout(
        stdout,
        prediction_method_name="netmhciipan",
        sequence_key_mapping=None):
    """
    # Threshold for Strong binding peptides (IC50)	50.000 nM
    # Threshold for Weak binding peptides (IC50)	500.000 nM

    # Threshold for Strong binding peptides (%Rank)	0.5%
    # Threshold for Weak binding peptides (%Rank)	2%

    # Allele: DRB1_0301
    --------------------------------------------------------------------------------------------------------------------------------------------
       Seq          Allele              Peptide    Identity  Pos      Core  Core_Rel 1-log50k(aff)  Affinity(nM)  %Rank Exp_Bind  BindingLevel
    --------------------------------------------------------------------------------------------------------------------------------------------
         0         DRB1_0301      AGFKGEQGPKGEPG    Sequence    2    FKGEQGPKG 0.810         0.080      21036.68  50.00   9.999
         1         DRB1_0301     GELIGTLNAAKVPAD    Sequence    2    LIGTLNAAK 0.650         0.340       1268.50  32.00   9.999
         2         DRB1_0301    PEVIPMFSALSEGATP    Sequence    5    MFSALSEGA 0.385         0.180       7161.16  50.00   9.999
         3         DRB1_0301       PKYVKQNTLKLAT    Sequence    2    YVKQNTLKL 0.575         0.442        418.70   6.00   9.999   <=WB
         4         DRB1_0301     VGSDWRFLRGYHQYA    Sequence    0    VGSDWRFLR 0.575         0.466        322.07  10.00   9.999   <=WB
         5         DRB1_0301         XFVKQNAAALX    Sequence    2    VKQNAAALX 0.500         0.262       2939.20  15.00   9.999
         6         DRB1_0301     AAYSDQATPLLLSPR    Sequence    1    AYSDQATPL 0.395         0.291       2152.21  50.00   9.999
         7         DRB1_0301     PVSKMRMATPLLMQA    Sequence    4    MRMATPLLM 0.890         0.770         12.00   0.01   9.999   <=SB
         8         DRB1_0301        AYMRADAAAGGA    Sequence    2    MRADAAAGG 0.835         0.303       1887.87  15.00   9.999
         9         DRB1_0301       PKYVKQNTLKLAT    Sequence    2    YVKQNTLKL 0.575         0.442        418.70   6.00   9.999   <=WB
        10         DRB1_0301     ENPVVHFFKNIVTPR    Sequence    6    FFKNIVTPR 0.425         0.357       1049.04  32.00   9.999
    """
    check_stdout_error(stdout, "NetMHCIIpan")
    return parse_stdout(
        stdout=stdout,
        prediction_method_name=prediction_method_name,
        sequence_key_mapping=sequence_key_mapping,
        key_index=3,
        offset_index=0,
        peptide_index=2,
        allele_index=1,
        ic50_index=8,
        rank_index=9,
        score_index=7)

def parse_netmhciipan4_stdout(
        stdout,
        prediction_method_name="netmhciipan",
        sequence_key_mapping=None,
        mode="elution_score"):
    """
    # Threshold for Strong binding peptides (%Rank) 2%
    # Threshold for Weak binding peptides (%Rank)   10%

    # Allele: DRB1_0101
    --------------------------------------------------------------------------------------------------------------------------------------------
     Pos           MHC              Peptide   Of        Core  Core_Rel        Identity      Score_EL %Rank_EL Exp_Bind      Score_BA  Affinity(nM) %Rank_BA  BindLevel
    --------------------------------------------------------------------------------------------------------------------------------------------
       1     DRB1_0101      PAPAPSWPLSSSVPS    4   PSWPLSSSV     0.327            test      0.000857    79.79       NA      0.327674       1442.91    54.35
       2     DRB1_0101      APAPSWPLSSSVPSQ    3   PSWPLSSSV     0.333            test      0.001268    71.87       NA      0.346949       1171.30    50.15
       3     DRB1_0101      PAPSWPLSSSVPSQK    4   WPLSSSVPS     0.713            test      0.002836    54.45       NA      0.412004        579.40    36.66
       4     DRB1_0101      APSWPLSSSVPSQKT    3   WPLSSSVPS     0.773            test      0.003677    49.14       NA      0.448939        388.53    29.75
       5     DRB1_0101      PSWPLSSSVPSQKTY    2   WPLSSSVPS     0.407            test      0.001602    66.79       NA      0.470979        306.09    25.98
       6     DRB1_0101      SWPLSSSVPSQKTYQ    3   LSSSVPSQK     0.633            test      0.001671    65.82       NA      0.476222        289.21    25.07
       7     DRB1_0101      WPLSSSVPSQKTYQG    3   SSSVPSQKT     0.553            test      0.001697    65.45       NA      0.447217        395.83    30.05


    """
    check_stdout_error(stdout, "NetMHCIIpan")

    if mode not in ["elution_score", "binding_affinity"]:
        raise ValueError("Mode is %s but must be one of: elution_score, binding affinity" % mode)

    # the offset specified in "pos" (at index 0) is 1-based instead of 0-based. we adjust it to be
    # 0-based, as in all the other netmhc predictors supported by this library.
    transforms = {
        0: lambda x: int(x) - 1,
    }

    # we're running NetMHCIIpan 4 with -BA every time so both EL and BA are available, but only
    # return one of them depending on the input mode
    return parse_stdout(
        stdout=stdout,
        prediction_method_name=prediction_method_name,
        sequence_key_mapping=sequence_key_mapping,
        key_index=6,
        offset_index=0,
        peptide_index=2,
        allele_index=1,
        ic50_index=11 if mode == "binding_affinity" else None,
        rank_index=8 if mode == "elution_score" else 12,
        score_index=7 if mode == "elution_score" else 10,
        transforms=transforms)

def parse_netmhcstabpan(
        stdout,
        prediction_method_name="netmhcstabpan",
        sequence_key_mapping=None,
        mode=None):
    """
    # NetMHCstabpan version 1.0

    # Input is in PEPTIDE format

    HLA-A02:01 : Distance to traning data  0.000 (using nearest neighbor HLA-A02:01)

    # Rank Threshold for Strong binding peptides   0.500
    # Rank Threshold for Weak binding peptides   2.000
    -----------------------------------------------------------------------------------------------------
     pos      HLA         peptide         Identity       Pred     Thalf(h) %Rank_Stab BindLevel
    -----------------------------------------------------------------------------------------------------
        0  HLA-A*02:01   AAAAAAAAAA         PEPLIST      0.075       0.27      19.00
    -----------------------------------------------------------------------------------------------------

    Protein PEPLIST. Allele HLA-A*02:01. Number of high binders 0. Number of weak binders 0. Number of peptides 1

    -----------------------------------------------------------------------------------------------------
    """

    # the offset specified in "pos" (at index 0) is 1-based instead of 0-based. we adjust it to be
    # 0-based, as in all the other netmhc predictors supported by this library.
    transforms = {
        0: lambda x: int(x) - 1,
    }
    return parse_stdout(
        stdout=stdout,
        prediction_method_name=prediction_method_name,
        sequence_key_mapping=sequence_key_mapping,
        key_index=3,
        offset_index=0,
        peptide_index=2,
        allele_index=1,
        score_index=5,
        rank_index=6,
        transforms=transforms)

def parse_netmhciipan43_stdout(
        stdout,
        prediction_method_name="netmhciipan",
        sequence_key_mapping=None,
        mode="elution_score"):
    """
    # Threshold for Strong binding peptides (%Rank)  1.00%
    # Threshold for Weak binding peptides (%Rank)    5.00%

    # DRB1_0101 : Distance to training data  0.000 (using nearest neighbor DRB1_0101)

    # Allele: DRB1_0101
    --------------------------------------------------------------------------------------------------------------------------------------------
     Pos               MHC              Peptide   Of        Core  Core_Rel Inverted        Identity      Score_EL %Rank_EL  Exp_Bind      Score_BA %Rank_BA  Affinity(nM)  BindLevel
    --------------------------------------------------------------------------------------------------------------------------------------------
       1         DRB1_0101        AAAGAEAGKATTE    1   AAGAEAGKA     0.740        0        Sequence      0.000143    72.70     0.000      0.086339    95.13      19645.62
       2         DRB1_0101      AALAAAAGVPPADKY    2   LAAAAGVPP     0.950        0        Sequence      0.030661    15.86     0.300      0.616438     7.87         63.44
       3         DRB1_0101         EKPGNRNPYENL    3   GNRNPYENL     0.680        0        Sequence      0.000000   100.00     0.670      0.030746    98.27      35850.53
       4         DRB1_0101      STWLLKPGAGIMIFD    2   WLLKPGAGI     0.420        0        Sequence      0.027506    16.58     0.000      0.708645     2.85         23.39
       5         DRB1_0101      KSVPLEMLLINLTTI    4   LEMLLINLT     0.980        0        Sequence      0.007346    26.97     0.560      0.662093     4.95         38.71
       6         DRB1_0101      KTKEDLFGKKNLIPS    5   LFGKKNLIP     0.560        0        Sequence      0.000144    72.62     0.670      0.298919    55.60       1969.51
       7         DRB1_0101      KIYHKCDNACIGSIR    2   YHKCDNACI     0.830        0        Sequence      0.003393    34.40     0.940      0.407849    34.13        606.04
       8         DRB1_0101      PCLFMRTVSHVILHG    3   FMRTVSHVI     0.960        0        Sequence      0.047565    13.05     0.345      0.664932     4.80         37.54
       9         DRB1_0101      GAATVAAGAATTAAG    4   VAAGAATTA     0.920        0        Sequence      0.164981     6.81     0.000      0.514261    17.98        191.63
  """
    check_stdout_error(stdout, "NetMHCIIpan")

    if mode not in ["elution_score", "binding_affinity"]:
        raise ValueError("Mode is %s but must be one of: elution_score, binding affinity" % mode)

    # the offset specified in "pos" (at index 0) is 1-based instead of 0-based. we adjust it to be
    # 0-based, as in all the other netmhc predictors supported by this library.
    transforms = {
        0: lambda x: int(x) - 1,
    }

    # we're running NetMHCIIpan 4.3 with -BA every time so both EL and BA are available, but only
    # return one of them depending on the input mode
    return parse_stdout(
        stdout=stdout,
        prediction_method_name=prediction_method_name,
        sequence_key_mapping=sequence_key_mapping,
        key_index=7,
        offset_index=0,
        peptide_index=2,
        allele_index=1,
        ic50_index=13 if mode == "binding_affinity" else None,
        rank_index=9 if mode == "elution_score" else 12,
        score_index=8 if mode == "elution_score" else 11,
        transforms=transforms)
