from mhctools.cli.script import parse_args, run_predictor
from nose.tools import eq_
import tempfile
from os import remove

def test_peptides_without_subsequences():
    peptide = "SIINFEKLQY"
    args = parse_args([
        "--mhc-predictor", "netmhc",
        "--mhc-peptide-lengths", "9",
        "--sequence", peptide,
        "--mhc-alleles", "H-2-Kb"])
    binding_predictions = run_predictor(args)
    eq_(len(binding_predictions), 1, binding_predictions)
    eq_(binding_predictions[0].peptide, peptide)

def test_peptides_with_subsequences():
    peptide = "SIINFEKLQY"
    args = parse_args([
        "--mhc-predictor", "netmhc",
        "--mhc-peptide-lengths", "9",
        "--sequence", peptide,
        "--extract-subsequences",
        "--mhc-alleles", "H-2-Kb"])
    binding_predictions = sorted(run_predictor(args), key=lambda bp: bp.offset)
    eq_(len(binding_predictions), 2, binding_predictions)
    eq_(binding_predictions[0].peptide, peptide[:9])
    eq_(binding_predictions[1].peptide, peptide[1:10])

def test_peptides_file_without_subsequences():
    peptide = "SIINFEKLQY"
    with tempfile.NamedTemporaryFile(mode="w", delete=False) as f:
        f.write("%s\n" % peptide)

    args = parse_args([
        "--mhc-predictor", "netmhc",
        "--mhc-peptide-lengths", "9",
        "--input-peptides-file", f.name,
        "--mhc-alleles", "H-2-Kb"])
    binding_predictions = run_predictor(args)
    eq_(len(binding_predictions), 1, binding_predictions)
    eq_(binding_predictions[0].peptide, peptide)
    remove(f.name)

def test_peptides_file_with_subsequences():
    peptide = "SIINFEKLQY"
    with tempfile.NamedTemporaryFile(mode="w", delete=False) as f:
        f.write("%s\n" % peptide)

    args = parse_args([
        "--mhc-predictor", "netmhc",
        "--mhc-peptide-lengths", "9",
        "--input-peptides-file", f.name,
        "--extract-subsequences",
        "--mhc-alleles", "H-2-Kb"])
    binding_predictions = sorted(run_predictor(args), key=lambda bp: bp.offset)
    eq_(len(binding_predictions), 2, binding_predictions)
    eq_(binding_predictions[0].peptide, peptide[:9])
    eq_(binding_predictions[1].peptide, peptide[1:10])
    remove(f.name)

