from mhctools.file_formats import parse_netmhc_stdout


def test_mhc_stdout():
    netmhc_output = """
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

    epitope_collection = parse_netmhc_stdout(
      netmhc_output,
      fasta_dictionary=fasta_dictionary,
      prediction_method_name="netmhc")

    assert len(epitope_collection) == 12

    for i, entry in enumerate(epitope_collection):
        assert entry.allele == 'HLA-A*02:03', \
            "Expected entry %s to have allele 'HLA-A*02:03'" (entry,)
        if i == 0:
            # expect the epitopes to be sorted in increasing IC50
            assert entry.value == 189.74
            assert entry.percentile_rank == 4.00
