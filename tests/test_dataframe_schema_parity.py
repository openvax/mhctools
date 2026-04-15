# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#       http://www.apache.org/licenses/LICENSE-2.0

"""Regression test for openvax/mhctools#193.

predict_peptides_dataframe (legacy/deprecated) and predict_proteins_dataframe
(current) must emit the same canonical column schema so downstream code can
treat rows from either path uniformly.
"""

from mhctools import RandomBindingPredictor
from mhctools.pred import COLUMNS


def test_peptides_and_proteins_dataframes_share_schema():
    p = RandomBindingPredictor(
        alleles=["HLA-A*02:01"], default_peptide_lengths=[9],
    )
    df_peptides = p.predict_peptides_dataframe(["SIINFEKLA"])
    df_proteins = p.predict_proteins_dataframe({"src": "MASIINFEKLA"})

    assert list(df_peptides.columns) == list(COLUMNS), (
        "predict_peptides_dataframe must emit canonical COLUMNS schema "
        "(was previously missing predictor_version, kind, value and used "
        "prediction_method_name; see #193)"
    )
    assert list(df_proteins.columns) == list(COLUMNS)
    assert list(df_peptides.columns) == list(df_proteins.columns)


def test_peptides_dataframe_has_predictor_identity_columns():
    """Downstream (e.g. topiary CachedPredictor) needs a stable
    (predictor_name, predictor_version) identity on every row."""
    p = RandomBindingPredictor(alleles=["HLA-A*02:01"], default_peptide_lengths=[9])
    df = p.predict_peptides_dataframe(["SIINFEKLA"])
    assert "predictor_name" in df.columns
    assert "predictor_version" in df.columns
    assert "kind" in df.columns
    assert "value" in df.columns
    # legacy name must not reappear
    assert "prediction_method_name" not in df.columns
    # canonical name is populated (RandomBindingPredictor sets it)
    assert df["predictor_name"].iloc[0] != ""
