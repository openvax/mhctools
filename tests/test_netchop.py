from numpy import testing

from mhctools import NetChop

# Peptides from http://tools.iedb.org/netchop/example/
peptides = """
MSLLTEVETPIRNEWGCRCNDSSDPLVVAASIIGIVHLILWIIDRLFSKSIYRIFKHGLKRGPSTEGVPESMREEYREEQQNAVDADDGHFVSIELE
MDSHTVSSFQVDCFLWHVRKQVADQDLGDAPFLDRLRRDQKSLKGRGSTLGLNIETATCVGKQIVERILKEESDEAFKMTMASALASRYLTDMTIEEMSRDWFMLMPKQKVAGPLCVRMDQAIMDKNIILKANFSVIFDRLENLTLLRAFTEEGAIVGEISPLPSLPGHTNEDVKNAIGVLIGGLEWNDNTVRVSETLQRFTWRSSNETGGPPFTPTQKRKMAGTIRSEV
MDSHTVSSFQDILMRMSKMQLGSSSGDLNGMITQFESLKLYRDSLGEAVMRLGDLHSLQHRNGKWREQLGQKFEEIRWLIEEVRHKLKTTENSFEQITFMQALQLLFEVEQEIRTFSFQLI
""".strip().split()


def test_simple():
    obj = NetChop()
    result = obj.predict(peptides)
    assert len(result) == 3
    assert len(result[0]) == len(peptides[0])
    assert len(result[1]) == len(peptides[1])
    assert len(result[2]) == len(peptides[2])

    # These numbers are from running http://tools.iedb.org/netchop
    # via the web interface on 12/19/2016.
    testing.assert_almost_equal(result[0][95], 0.976629)
    testing.assert_almost_equal(result[0][22], 0.022000)
    testing.assert_almost_equal(result[1][146], 0.977417)
    testing.assert_almost_equal(result[1][84], 0.285210)
    testing.assert_almost_equal(result[2][0], 0.547588)
    testing.assert_almost_equal(result[2][84], 0.104684)
