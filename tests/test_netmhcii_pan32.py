from mhctools import NetMHCIIpan43, NetMHCIIpan43_EL, NetMHCIIpan43_BA


def test_netmhciipan43():
    predictor = NetMHCIIpan43(alleles=['DRB1_0101'])
    predictions = predictor.predict_subsequences(
        ["AAGAARIIAVDINKD", "AAGIVESVGEGVTTV"]).to_dataframe()
    assert predictions.shape == (2, 9)
    assert all(pp is None for pp in predictions.affinity)


def test_netmhciipan43_ba():
    predictor_ba = NetMHCIIpan43_BA(
        alleles=['DRB1_0101', "HLA-DQA1*05:11-DQB1*03:02"])
    binding_predictions = predictor_ba.predict_subsequences(
        ["AAGAARIIAVDINKD", "AAGIVESVGEGVTTV"]).to_dataframe()
    assert binding_predictions.shape == (4, 9)
    assert not all(pp is None for pp in binding_predictions.affinity)


def test_netmhciipan43_el():
    predictor_el = NetMHCIIpan43_EL(
        alleles=['DRB1_0101', "HLA-DQA1*05:11-DQB1*03:02"])
    el_predictions = predictor_el.predict_subsequences(
        ["AAGAARIIAVDINKD", "AAGIVESVGEGVTTV"]).to_dataframe()
    assert el_predictions.shape == (4, 9)
    assert all(pp is None for pp in el_predictions.affinity)
