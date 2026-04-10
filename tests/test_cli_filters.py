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

"""
Tests for CLI --max-affinity, --max-percentile-rank, --min-score filters.

Uses the ``random`` predictor so no external tool is needed.
"""

import pandas as pd

from mhctools.cli.script import parse_args, run_predictor, apply_filters


PEPTIDES = ["SIINFEKL", "AAAAAAAAA"]


def _make_df():
    """Run the random predictor and get a DataFrame."""
    args = parse_args([
        "--mhc-predictor", "random",
        "--sequence", *PEPTIDES,
        "--mhc-alleles", "HLA-A*02:01",
    ])
    return run_predictor(args).to_dataframe()


def test_no_filters_returns_all():
    args = parse_args([
        "--mhc-predictor", "random",
        "--sequence", *PEPTIDES,
        "--mhc-alleles", "HLA-A*02:01",
    ])
    df = _make_df()
    filtered = apply_filters(df, args)
    assert len(filtered) == len(df)


def test_max_affinity_zero_excludes_all():
    args = parse_args([
        "--mhc-predictor", "random",
        "--sequence", *PEPTIDES,
        "--mhc-alleles", "HLA-A*02:01",
        "--max-affinity", "0",
    ])
    df = _make_df()
    # random produces affinity in (0, 10000), so 0 threshold should drop most/all
    filtered = apply_filters(df, args)
    assert len(filtered) <= len(df)
    assert all(filtered["affinity"] <= 0)


def test_max_affinity_large_keeps_all():
    args = parse_args([
        "--mhc-predictor", "random",
        "--sequence", *PEPTIDES,
        "--mhc-alleles", "HLA-A*02:01",
        "--max-affinity", "999999",
    ])
    df = _make_df()
    filtered = apply_filters(df, args)
    assert len(filtered) == len(df)


def test_max_percentile_rank():
    args = parse_args([
        "--mhc-predictor", "random",
        "--sequence", *PEPTIDES,
        "--mhc-alleles", "HLA-A*02:01",
        "--max-percentile-rank", "50",
    ])
    df = _make_df()
    filtered = apply_filters(df, args)
    assert all(filtered["percentile_rank"] <= 50)


def test_min_score():
    args = parse_args([
        "--mhc-predictor", "random",
        "--sequence", *PEPTIDES,
        "--mhc-alleles", "HLA-A*02:01",
        "--min-score", "0.5",
    ])
    df = _make_df()
    filtered = apply_filters(df, args)
    assert all(filtered["score"] >= 0.5)


def test_combined_filters():
    """Multiple filters should all apply (AND logic)."""
    args = parse_args([
        "--mhc-predictor", "random",
        "--sequence", *PEPTIDES,
        "--mhc-alleles", "HLA-A*02:01",
        "--max-affinity", "5000",
        "--max-percentile-rank", "50",
        "--min-score", "0.5",
    ])
    df = _make_df()
    filtered = apply_filters(df, args)
    if len(filtered) > 0:
        assert all(filtered["affinity"] <= 5000)
        assert all(filtered["percentile_rank"] <= 50)
        assert all(filtered["score"] >= 0.5)


def test_filter_with_nan_values():
    """Rows with NaN in a filtered column should be kept (not silently dropped)."""
    df = pd.DataFrame({
        "peptide": ["AAA", "BBB"],
        "affinity": [100.0, float("nan")],
        "percentile_rank": [1.0, float("nan")],
        "score": [0.9, float("nan")],
    })
    args = parse_args([
        "--mhc-predictor", "random",
        "--sequence", "SIINFEKL",
        "--mhc-alleles", "HLA-A*02:01",
        "--max-affinity", "500",
    ])
    filtered = apply_filters(df, args)
    # Both rows should survive: one passes the filter, the other has NaN
    assert len(filtered) == 2
