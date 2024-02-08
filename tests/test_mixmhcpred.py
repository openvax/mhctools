from tempfile import NamedTemporaryFile
from os import remove
from mhctools.mixmhcpred import parse_mixmhcpred_results
from nose.tools import eq_

example_output =  """Peptide\tScore_bestAllele\tBestAllele\t%Rank_bestAllele\tScore_A0201\t%Rank_A0201
MLDDFSAGA\t0.182093\tA0201\t0.3\t0.182093\t0.3
SPEGEETII\t-0.655341\tA0201\t51.0\t-0.655341\t51.0
ILDRIITNA\t0.203906\tA0201\t0.3\t0.203906\t0.3"""

def test_parse_mixmhcpred_results():
    with NamedTemporaryFile(mode="r+") as f:
        f.write(example_output)
        f.flush()
        binding_results = parse_mixmhcpred_results(f.name)

        eq_(len(binding_results), 3)
        eq_(binding_results[0].peptide, "MLDDFSAGA")
        eq_(binding_results[1].peptide, "SPEGEETII")
        eq_(binding_results[2].peptide, "ILDRIITNA")
