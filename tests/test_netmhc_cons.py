from mhctools import NetMHCcons
from mhcnames import normalize_allele_name


DEFAULT_ALLELE = 'HLA-A*02:01'

def test_netmhc_cons():
    alleles = [normalize_allele_name(DEFAULT_ALLELE)]
    cons_predictor = NetMHCcons(
        alleles=alleles,
        default_peptide_lengths=[9])
    sequence_dict = {
        "SMAD4-001": "ASIINFKELA",
        "TP53-001": "ASILLLVFYW"
    }
    binding_predictions = cons_predictor.predict_subsequences(
        sequence_dict=sequence_dict)

    assert len(binding_predictions) == 4, \
        "Expected 4 epitopes from %s" % (binding_predictions,)

def test_netmhc_cons_multiple_lengths():
    cons_predictor = NetMHCcons(alleles=["A6801"])
    binding_predictions = cons_predictor.predict_peptides(
        ["A" * 8, "A" * 9, "A" * 10, "A" * 11])
    assert len(binding_predictions) == 4, \
        "Expected 4 epitopes from %s" % (binding_predictions,)

def test_netmhc_cons_multiple_alleles():
    alleles = 'A*02:01,B*35:02'
    cons_predictor = NetMHCcons(
        alleles=alleles,
        default_peptide_lengths=[9])
    sequence_dict = {
        "SMAD4-001": "ASIINFKELA",
        "TP53-001": "ASILLLVFYW"
    }
    binding_predictions = cons_predictor.predict_subsequences(
        sequence_dict=sequence_dict)
    assert len(binding_predictions) == 8, \
        "Expected 4 binding predictions from %s" % (binding_predictions,)

def test_netmhc_cons_process_limits():
    alleles = [normalize_allele_name(DEFAULT_ALLELE)]
    sequence_dict = {
        "SMAD4-001": "ASIINFKELA",
        "TP53-001": "ASILLLVFYW",
        "SMAD4-002": "ASIINFKELS",
        "TP53-002": "ASILLLVFYS",
        "TP53-003": "ASILLLVFYT",
        "TP53-004": "ASILLLVFYG",
        "TP53-005": "ASILLLVFYG"
    }
    for process_limit in [1, 2, 10]:
        cons_predictor = NetMHCcons(
            alleles=alleles,
            default_peptide_lengths=[9],
            process_limit=process_limit
        )
        binding_predictions = cons_predictor.predict_subsequences(
            sequence_dict=sequence_dict)
        assert len(binding_predictions) == 14, \
            "Expected 14 binding predictions from but got %d: %s" % (
                len(binding_predictions),
                binding_predictions)

        source_names = [bp.source_sequence_name for bp in binding_predictions]
        for fasta_key in sequence_dict.keys():
            fasta_count = source_names.count(fasta_key)
            assert fasta_count == 2, \
                ("Expected each fasta key to appear twice, once for "
                 "each length, but saw %s %d time(s)" % (
                     fasta_key, fasta_count))

if __name__ == "__main__":
    test_netmhc_cons()
    test_netmhc_cons_process_limits()
