
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

from argparse import ArgumentTypeError, Namespace

import pytest

from mhctools.cli import mhc_predictors, make_mhc_arg_parser, predictors_from_args, _cls_accepts
from mhctools.cli.args import (
    _parse_predictor_names,
    _build_predictor,
    mhc_binding_predictor_from_args,
)
from mhctools.bigmhc import BigMHC, BigMHC_EL, BigMHC_IM
from mhctools.random_predictor import RandomBindingPredictor
from mhctools.pepsickle import Pepsickle


# ── Registry completeness ──────────────────────────────────────────

def test_all_models_in_registry():
    """Verify all major model classes are available via CLI."""
    expected = [
        "netmhcpan42", "netmhcpan42-ba", "netmhcpan42-el",
        "netmhciipan43", "netmhciipan43-ba", "netmhciipan43-el",
        "bigmhc", "bigmhc-el", "bigmhc-im",
        "netmhcstabpan", "netchop", "pepsickle",
        # legacy entries that should still be present
        "netmhcpan", "netmhcpan4", "netmhcpan41",
        "netmhciipan", "netmhciipan4",
        "mhcflurry", "mixmhcpred", "random",
    ]
    for name in expected:
        assert name in mhc_predictors, f"{name} missing from mhc_predictors registry"


def test_all_registry_values_are_callable():
    """Every registry value should be callable (class or factory function)."""
    for name, cls in mhc_predictors.items():
        assert callable(cls), (
            f"mhc_predictors[{name!r}] = {cls!r} is not callable"
        )


# ── BigMHC_EL / BigMHC_IM subclass registry entries ───────────────

def test_bigmhc_el_maps_to_subclass():
    assert mhc_predictors["bigmhc-el"] == BigMHC_EL


def test_bigmhc_im_maps_to_subclass():
    assert mhc_predictors["bigmhc-im"] == BigMHC_IM


def test_bigmhc_maps_to_base():
    assert mhc_predictors["bigmhc"] == BigMHC


def test_bigmhc_el_is_subclass():
    assert issubclass(BigMHC_EL, BigMHC)


def test_bigmhc_im_is_subclass():
    assert issubclass(BigMHC_IM, BigMHC)


def test_bigmhc_el_no_mode_arg():
    """BigMHC_EL should not accept a mode argument — it's fixed to 'el'."""
    import inspect
    sig = inspect.signature(BigMHC_EL.__init__)
    assert "mode" not in sig.parameters


def test_bigmhc_im_no_mode_arg():
    """BigMHC_IM should not accept a mode argument — it's fixed to 'im'."""
    import inspect
    sig = inspect.signature(BigMHC_IM.__init__)
    assert "mode" not in sig.parameters


# ── _cls_accepts introspection ─────────────────────────────────────

def test_cls_accepts_bigmhc_path():
    assert _cls_accepts(BigMHC, "bigmhc_path") is True
    assert _cls_accepts(BigMHC_EL, "bigmhc_path") is True
    assert _cls_accepts(BigMHC_IM, "bigmhc_path") is True


def test_cls_accepts_bigmhc_rejects_program_name():
    assert _cls_accepts(BigMHC, "program_name") is False


def test_cls_accepts_bigmhc_rejects_default_peptide_lengths():
    assert _cls_accepts(BigMHC, "default_peptide_lengths") is False
    assert _cls_accepts(BigMHC_EL, "default_peptide_lengths") is False


def test_cls_accepts_random_default_peptide_lengths():
    assert _cls_accepts(RandomBindingPredictor, "default_peptide_lengths") is True


def test_cls_accepts_random_rejects_bigmhc_path():
    assert _cls_accepts(RandomBindingPredictor, "bigmhc_path") is False


def test_cls_accepts_pepsickle_rejects_alleles():
    assert _cls_accepts(Pepsickle, "alleles") is False


def test_cls_accepts_pepsickle_accepts_dpl():
    assert _cls_accepts(Pepsickle, "default_peptide_lengths") is True


# ── _parse_predictor_names ─────────────────────────────────────────

def test_parse_predictor_names_single():
    assert _parse_predictor_names("random") == ["random"]


def test_parse_predictor_names_comma_separated():
    result = _parse_predictor_names("random,bigmhc,pepsickle")
    assert result == ["random", "bigmhc", "pepsickle"]


def test_parse_predictor_names_strips_whitespace():
    result = _parse_predictor_names("  random , bigmhc ")
    assert result == ["random", "bigmhc"]


def test_parse_predictor_names_case_insensitive():
    result = _parse_predictor_names("Random,BIGMHC")
    assert result == ["random", "bigmhc"]


def test_parse_predictor_names_unknown_raises():
    with pytest.raises(ArgumentTypeError, match="Unknown predictor"):
        _parse_predictor_names("nonexistent")


# ── Argparser construction ─────────────────────────────────────────

def test_parser_single_predictor():
    parser = make_mhc_arg_parser()
    args = parser.parse_args(["--mhc-predictor", "random", "--mhc-alleles", "HLA-A*02:01"])
    # nargs="+" with type that returns list → list of lists
    assert args.mhc_predictor == [["random"]]


def test_parser_multiple_space_separated():
    parser = make_mhc_arg_parser()
    args = parser.parse_args([
        "--mhc-predictor", "random", "bigmhc",
        "--mhc-alleles", "HLA-A*02:01",
    ])
    assert args.mhc_predictor == [["random"], ["bigmhc"]]


def test_parser_comma_separated():
    parser = make_mhc_arg_parser()
    args = parser.parse_args([
        "--mhc-predictor", "random,bigmhc",
        "--mhc-alleles", "HLA-A*02:01",
    ])
    assert args.mhc_predictor == [["random", "bigmhc"]]


def test_parser_mixed_space_and_comma():
    parser = make_mhc_arg_parser()
    args = parser.parse_args([
        "--mhc-predictor", "random,bigmhc", "pepsickle",
        "--mhc-alleles", "HLA-A*02:01",
    ])
    assert args.mhc_predictor == [["random", "bigmhc"], ["pepsickle"]]


# ── _build_predictor kwarg filtering ──────────────────────────────

def _stub_args(**overrides):
    defaults = dict(
        mhc_predictor=None,
        mhc_alleles="HLA-A*02:01",
        mhc_alleles_file=None,
        mhc_peptide_lengths=None,
        mhc_epitope_lengths=None,
        mhc_predictor_models_path=None,
        mhc_predictor_path=None,
        do_not_raise_on_error=False,
    )
    defaults.update(overrides)
    return Namespace(**defaults)


def test_build_predictor_skips_dpl_for_bigmhc():
    """BigMHC doesn't accept default_peptide_lengths; should not receive it."""
    args = _stub_args()
    pred = _build_predictor(
        RandomBindingPredictor, "random",
        alleles=["HLA-A*02:01"], peptide_lengths=[9, 10], args=args)
    assert pred.default_peptide_lengths == [9, 10]

    # BigMHC would TypeError if dpl were passed; FileNotFoundError means
    # it got through kwarg filtering and failed on path lookup
    with pytest.raises(FileNotFoundError):
        _build_predictor(
            BigMHC, "bigmhc",
            alleles=["HLA-A*02:01"], peptide_lengths=[9, 10], args=args)


def test_build_predictor_bigmhc_path_mapping():
    """--mhc-predictor-path should map to bigmhc_path for BigMHC."""
    args = _stub_args(mhc_predictor_path="/some/path")
    with pytest.raises(FileNotFoundError):
        _build_predictor(
            BigMHC_EL, "bigmhc-el",
            alleles=["HLA-A*02:01"], peptide_lengths=None, args=args)


def test_build_predictor_skips_alleles_for_pepsickle():
    """Pepsickle doesn't accept alleles; should not receive them."""
    args = _stub_args()
    pred = _build_predictor(
        Pepsickle, "pepsickle",
        alleles=["HLA-A*02:01"], peptide_lengths=[9], args=args)
    assert pred is not None


# ── predictors_from_args (full integration) ────────────────────────

def _make_args(predictor_tokens, **overrides):
    """Build a Namespace matching what argparse produces with nargs='+' type."""
    names = []
    for tok in predictor_tokens:
        names.append(_parse_predictor_names(tok))
    defaults = dict(
        mhc_predictor=names,
        mhc_alleles="HLA-A*02:01",
        mhc_alleles_file=None,
        mhc_peptide_lengths=None,
        mhc_epitope_lengths=None,
        mhc_predictor_models_path=None,
        mhc_predictor_path=None,
        do_not_raise_on_error=False,
    )
    defaults.update(overrides)
    return Namespace(**defaults)


def test_predictors_from_args_single():
    args = _make_args(["random"])
    result = predictors_from_args(args)
    assert len(result) == 1
    assert isinstance(result[0], RandomBindingPredictor)


def test_predictors_from_args_multiple():
    args = _make_args(["random", "random"])
    result = predictors_from_args(args)
    assert len(result) == 2
    assert all(isinstance(r, RandomBindingPredictor) for r in result)


def test_predictors_from_args_comma_separated():
    args = _make_args(["random,random"])
    result = predictors_from_args(args)
    assert len(result) == 2


def test_predictors_from_args_no_allele_predictor():
    """Processing predictors should not require alleles."""
    args = _make_args(["pepsickle"], mhc_alleles="")
    result = predictors_from_args(args)
    assert len(result) == 1
    assert isinstance(result[0], Pepsickle)


def test_predictors_from_args_mixed_allele_and_no_allele():
    """Alleles should be fetched when at least one predictor needs them."""
    args = _make_args(["random,pepsickle"])
    result = predictors_from_args(args)
    assert len(result) == 2
    assert isinstance(result[0], RandomBindingPredictor)
    assert isinstance(result[1], Pepsickle)


def test_predictors_from_args_peptide_lengths_only_where_accepted():
    """default_peptide_lengths passed to random, skipped for bigmhc."""
    args = _make_args(["random"], mhc_peptide_lengths=[8, 9])
    result = predictors_from_args(args)
    assert result[0].default_peptide_lengths == [8, 9]


def test_predictors_from_args_bigmhc_path():
    args = _make_args(["bigmhc-el"], mhc_predictor_path="/some/path")
    with pytest.raises(FileNotFoundError):
        predictors_from_args(args)


# ── mhc_binding_predictor_from_args (legacy compat) ───────────────

def test_legacy_binding_predictor_single():
    args = _make_args(["random"])
    result = mhc_binding_predictor_from_args(args)
    assert isinstance(result, RandomBindingPredictor)


def test_legacy_binding_predictor_rejects_multiple():
    args = _make_args(["random", "random"])
    with pytest.raises(ValueError, match="exactly one"):
        mhc_binding_predictor_from_args(args)
