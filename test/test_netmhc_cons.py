from mhctools import NetMHCcons
from mhctools.alleles import normalize_allele_name


DEFAULT_ALLELE = 'HLA-A*02:01'

def test_netmhc_cons():
    alleles = [normalize_allele_name(DEFAULT_ALLELE)]
    cons_predictor = NetMHCcons(
        alleles=alleles,
        epitope_lengths=[9])
    fasta_dictionary = {
        "SMAD4-001": "ASIINFKELA",
        "TP53-001": "ASILLLVFYW"
    }
    epitope_collection = cons_predictor.predict(
        fasta_dictionary=fasta_dictionary)

    assert len(epitope_collection) == 4, \
        "Expected 4 epitopes from %s" % (epitope_collection,)

def test_netmhc_cons_multiple_alleles():
    alleles = 'A*02:01,B*35:02'
    cons_predictor = NetMHCcons(
        alleles=alleles,
        epitope_lengths=[9])
    fasta_dictionary = {
        "SMAD4-001": "ASIINFKELA",
        "TP53-001": "ASILLLVFYW"
    }
    epitope_collection = cons_predictor.predict(
        fasta_dictionary=fasta_dictionary)
    assert len(epitope_collection) == 8, \
        "Expected 4 epitopes from %s" % (epitope_collection,)

def test_netmhc_cons_chunking():
    alleles = [normalize_allele_name(DEFAULT_ALLELE)]
    fasta_dictionary = {
        "SMAD4-001": "ASIINFKELA",
        "TP53-001": "ASILLLVFYW",
        "SMAD4-002": "ASIINFKELS",
        "TP53-002": "ASILLLVFYS",
        "TP53-003": "ASILLLVFYT",
        "TP53-004": "ASILLLVFYG",
        "TP53-005": "ASILLLVFYG"
    }
    for max_file_records in [1, 5, 20]:
        for process_limit in [1, 2, 10]:
            cons_predictor = NetMHCcons(
                alleles=alleles,
                epitope_lengths=[9],
                max_file_records=max_file_records,
                process_limit=process_limit
            )
            epitope_collection = cons_predictor.predict(
                fasta_dictionary=fasta_dictionary)
            assert len(epitope_collection) == 14, \
                "Expected 14 epitopes from %s" % (epitope_collection,)
            source_keys = []
            for epitope in epitope_collection:
                source_keys.append(epitope.source_sequence_key)
            for fasta_key in fasta_dictionary.keys():
                fasta_count = source_keys.count(fasta_key)
                assert fasta_count == 2, \
                    ("Expected each fasta key to appear twice, once for "
                     "each length, but saw %s %d time(s)" % (
                         fasta_key, fasta_count))

if __name__ == "__main__":
    test_netmhc_cons()
    test_netmhc_cons_chunking()