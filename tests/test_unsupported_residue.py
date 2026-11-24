# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#       http://www.apache.org/licenses/LICENSE-2.0

"""Tests for typed UnsupportedResidueError across base, MHCflurry, and adapters."""

import pickle
import types

import pandas as pd
import pytest

import mhctools
from mhctools import (
    Calis,
    CapHLA,
    DeepTAP,
    ERAMER,
    MHCflurry,
    MHCflurry_Affinity,
    MixMHC2pred,
    MixMHCpred,
    PRIME,
    PeptideInput,
    PeptiVerse,
    PlifePred2,
    SMM,
    TCR,
    UnsupportedResidueError,
)
from mhctools.base_commandline_predictor import BaseCommandlinePredictor
from mhctools.base_predictor import BasePredictor
from mhctools.binding_prediction import BindingPrediction
from mhctools.binding_prediction_collection import BindingPredictionCollection
from mhctools.mhcflurry import translate_mhcflurry_residue_error
from mhctools.mixtcrpred import _validate_tcr
from mhctools.nettcr import _encode_feature


def _make_fake_mhcflurry_predictor(
    *,
    supported=("HLA-A*02:01",),
    affinity_error=None,
    presentation_error=None,
    reject_nonstandard=True,
    bad_presentation_rows=False,
    duplicate_presentation_row=False,
    unrequested_presentation_allele=None,
    missing_presentation_row=False,
):
    """Build a fake Class1PresentationPredictor for offline tests."""
    affinity_calls = []
    presentation_calls = []
    canonical = set("ACDEFGHIKLMNPQRSTVWY")

    def _check_native_peptides(peptides):
        out_of_range = [p for p in peptides if not (8 <= len(p) <= 15)]
        if out_of_range:
            raise ValueError(
                "%d peptides have lengths outside of supported range [8, 15]: %s"
                % (len(out_of_range), str(pd.Series(out_of_range).unique()))
            )
        if reject_nonstandard:
            nonstandard = [
                p for p in peptides if any(ch.upper() not in canonical for ch in p)
            ]
            if nonstandard:
                raise ValueError(
                    "%d peptides have nonstandard amino acids: %s"
                    % (len(nonstandard), str(pd.Series(nonstandard).unique()))
                )

    def predict_to_dataframe(peptides, alleles, include_percentile_ranks=True):
        affinity_calls.append(
            {
                "peptides": list(peptides),
                "alleles": list(alleles),
                "include_percentile_ranks": include_percentile_ranks,
            }
        )
        if affinity_error is not None:
            raise affinity_error
        _check_native_peptides(peptides)
        data = {
            "peptide": list(peptides),
            "allele": list(alleles),
            "prediction": [250.0] * len(peptides),
        }
        if include_percentile_ranks:
            data["prediction_percentile"] = [1.2] * len(peptides)
        return pd.DataFrame(data)

    affinity_predictor = types.SimpleNamespace(
        predict_to_dataframe=predict_to_dataframe,
        supported_alleles=list(supported),
        affinity_calls=affinity_calls,
    )

    def predict(
        peptides,
        alleles,
        sample_names=None,
        n_flanks=None,
        c_flanks=None,
        include_affinity_percentile=False,
        verbose=0,
        throw=True,
    ):
        presentation_calls.append(
            {
                "peptides": list(peptides),
                "alleles": alleles,
                "n_flanks": None if n_flanks is None else list(n_flanks),
                "c_flanks": None if c_flanks is None else list(c_flanks),
            }
        )
        if presentation_error is not None:
            raise presentation_error
        _check_native_peptides(peptides)
        rows = []
        best_allele = unrequested_presentation_allele or supported[0]
        for idx, peptide in enumerate(peptides):
            rows.append(
                {
                    "peptide": peptide,
                    "peptide_num": idx,
                    "best_allele": best_allele,
                    "presentation_score": 0.8,
                    "presentation_percentile": 1.0,
                    "processing_score": 0.6,
                }
            )
        if bad_presentation_rows:
            rows = rows[:-1]
        if duplicate_presentation_row and rows:
            rows[-1] = dict(rows[0])
        if missing_presentation_row and len(rows) >= 2:
            rows[1]["peptide_num"] = 0
        return pd.DataFrame(rows)

    return types.SimpleNamespace(
        affinity_predictor=affinity_predictor,
        predict=predict,
        supported_alleles=list(supported),
        affinity_calls=affinity_calls,
        presentation_calls=presentation_calls,
    )


# ---------------------------------------------------------------------------
# Public exception API
# ---------------------------------------------------------------------------


def test_unsupported_residue_error_public_export_and_attributes():
    assert "UnsupportedResidueError" in mhctools.__all__
    assert mhctools.UnsupportedResidueError is UnsupportedResidueError
    assert issubclass(UnsupportedResidueError, ValueError)

    error = UnsupportedResidueError(
        peptides=["SIINXEKL", "SIINXEKL", "AAAA*BBB"],
        residues=["X", "*", "X"],
        input_indices=[1, 3, 1],
        field="peptide",
        predictor_name="example",
    )
    assert isinstance(error, ValueError)
    assert error.peptide == "SIINXEKL"
    assert error.peptides == ("SIINXEKL", "AAAA*BBB")
    assert error.residues == ("*", "X")
    assert error.input_indices == (1, 3)
    assert error.field == "peptide"
    assert error.predictor_name == "example"
    assert "example" in str(error)
    assert "'X'" in str(error)

    unpickled = pickle.loads(pickle.dumps(error))
    assert type(unpickled) is UnsupportedResidueError
    assert str(unpickled) == str(error)
    assert unpickled.peptides == error.peptides
    assert unpickled.residues == error.residues
    assert unpickled.input_indices == error.input_indices
    assert unpickled.field == error.field
    assert unpickled.predictor_name == error.predictor_name


# ---------------------------------------------------------------------------
# BasePredictor._check_peptide_inputs
# ---------------------------------------------------------------------------


def test_base_predictor_raises_typed_error_for_X_and_invalid_chars():
    predictor = BasePredictor(min_peptide_length=8, max_peptide_length=15)

    with pytest.raises(UnsupportedResidueError, match="Invalid character 'X'") as exc_info:
        predictor._check_peptide_inputs(["SIINFEKL", "SIINXEKL"])
    error = exc_info.value
    assert isinstance(error, ValueError)
    assert error.peptide == "SIINXEKL"
    assert error.peptides == ("SIINXEKL",)
    assert error.residues == ("X",)
    assert error.input_indices == (1,)
    assert error.field == "peptide"
    assert error.predictor_name == "BasePredictor"

    with pytest.raises(UnsupportedResidueError, match="Invalid characters") as exc_info:
        predictor._check_peptide_inputs(["SIINFEKL", "SIIN*XKL"])
    error = exc_info.value
    assert error.peptides == ("SIIN*XKL",)
    assert error.residues == ("*", "X")
    assert error.input_indices == (1,)


def test_base_predictor_respects_allow_X_and_lowercase_flags():
    predictor_allowing_x = BasePredictor(
        min_peptide_length=8, max_peptide_length=15, allow_X_in_peptides=True
    )
    predictor_allowing_x._check_peptide_inputs(["SIINXEKL"])

    predictor_allowing_lower = BasePredictor(
        min_peptide_length=8,
        max_peptide_length=15,
        allow_X_in_peptides=False,
        allow_lowercase_in_peptides=True,
    )
    with pytest.raises(UnsupportedResidueError) as exc_info:
        predictor_allowing_lower._check_peptide_inputs(["siinxekl"])
    assert exc_info.value.residues == ("x",)
    assert exc_info.value.peptides == ("siinxekl",)


@pytest.mark.parametrize(
    "bad_peptide, match",
    [
        ("", "Invalid characters"),
        ("siinfekl", "Invalid lowercase"),
        ("SIIN", "too short"),
        ("A" * 16, "too long"),
    ],
)
def test_base_predictor_non_residue_failures_remain_plain_value_error(bad_peptide, match):
    predictor = BasePredictor(min_peptide_length=8, max_peptide_length=15)
    with pytest.raises(ValueError, match=match) as exc_info:
        predictor._check_peptide_inputs([bad_peptide])
    assert not isinstance(exc_info.value, UnsupportedResidueError)


def test_base_predictor_predict_remaps_indices_across_duplicate_peptides():
    class CheckedDummyPredictor(BasePredictor):
        def predict_peptides(self, peptides):
            self._check_peptide_inputs(peptides)
            return BindingPredictionCollection(
                [
                    BindingPrediction(
                        allele=allele,
                        peptide=peptide,
                        affinity=100.0,
                        prediction_method_name="dummy",
                    )
                    for allele in self.alleles
                    for peptide in peptides
                ]
            )

    predictor = CheckedDummyPredictor(alleles=["HLA-A*02:01"])
    with pytest.raises(UnsupportedResidueError) as exc_info:
        predictor.predict(["SIINFEKL", "SIINFEKL", "SIINXEKL", "SIINXEKL"])
    assert exc_info.value.peptides == ("SIINXEKL",)
    assert exc_info.value.residues == ("X",)
    assert exc_info.value.input_indices == (2, 3)


# ---------------------------------------------------------------------------
# MHCflurry & MHCflurry_Affinity native translation
# ---------------------------------------------------------------------------


def test_mhcflurry_predict_translates_native_residue_rejection_for_batch_recovery():
    fake = _make_fake_mhcflurry_predictor()
    predictor = MHCflurry(
        alleles=["HLA-A*02:01"],
        predictor=fake,
        include_affinity_percentile_ranks=False,
    )
    batch = ["SIINFEKL", "SIINXEKL", "GILGFVFTL", "AAAAUBBBB"]

    with pytest.raises(UnsupportedResidueError, match="nonstandard amino acids") as exc_info:
        predictor.predict(batch)
    error = exc_info.value
    assert isinstance(error, ValueError)
    assert isinstance(error.__cause__, ValueError)
    assert not isinstance(error.__cause__, UnsupportedResidueError)
    assert error.peptides == ("SIINXEKL", "AAAAUBBBB")
    assert error.peptide == "SIINXEKL"
    assert error.residues == ("B", "U", "X")
    assert error.input_indices == (1, 3)
    assert error.field == "peptide"
    assert error.predictor_name == "mhcflurry"

    # Safe batch recovery using offending indices
    recoverable = [p for i, p in enumerate(batch) if i not in set(error.input_indices)]
    results = predictor.predict(recoverable)
    assert [r.peptide for r in results] == ["SIINFEKL", "GILGFVFTL"]


def test_mhcflurry_predict_peptides_and_predict_proteins_translate_residue_error():
    fake = _make_fake_mhcflurry_predictor()
    predictor = MHCflurry(
        alleles=["HLA-A*02:01"],
        predictor=fake,
        include_affinity_percentile_ranks=False,
    )

    with pytest.raises(UnsupportedResidueError) as exc_info:
        predictor.predict_peptides(["SIINFEKL", "SIINXEKL"])
    assert exc_info.value.peptides == ("SIINXEKL",)
    assert exc_info.value.residues == ("X",)
    assert exc_info.value.input_indices == (1,)

    with pytest.raises(UnsupportedResidueError) as exc_info:
        predictor.predict_proteins({"prot": "AAASIINXEKLAAA"}, peptide_lengths=[8])
    assert exc_info.value.residues == ("X",)
    assert all("X" in p for p in exc_info.value.peptides)


def test_mhcflurry_presentation_call_residue_error_is_translated():
    fake = _make_fake_mhcflurry_predictor(
        reject_nonstandard=False,
        presentation_error=ValueError(
            "1 peptides have nonstandard amino acids: ['SIINXEKL']"
        ),
    )
    predictor = MHCflurry(
        alleles=["HLA-A*02:01"],
        predictor=fake,
        include_affinity_percentile_ranks=False,
    )
    with pytest.raises(UnsupportedResidueError) as exc_info:
        predictor.predict(["SIINFEKL", "SIINXEKL"])
    assert exc_info.value.peptides == ("SIINXEKL",)
    assert exc_info.value.residues == ("X",)
    assert exc_info.value.input_indices == (1,)


def test_mhcflurry_affinity_translates_and_remaps_duplicate_indices():
    fake = _make_fake_mhcflurry_predictor()
    predictor = MHCflurry_Affinity(
        alleles=["HLA-A*02:01"],
        predictor=fake.affinity_predictor,
        include_affinity_percentile_ranks=False,
    )
    with pytest.raises(UnsupportedResidueError) as exc_info:
        predictor.predict(["SIINFEKL", "SIINFEKL", "SIINXEKL"])
    assert exc_info.value.peptides == ("SIINXEKL",)
    assert exc_info.value.residues == ("X",)
    assert exc_info.value.input_indices == (2,)
    assert exc_info.value.predictor_name == "mhcflurry"


def test_mhcflurry_does_not_assume_fixed_alphabet_or_alter_sequences_and_flanks():
    fake = _make_fake_mhcflurry_predictor(reject_nonstandard=False)
    predictor = MHCflurry(
        alleles=["HLA-A*02:01"],
        predictor=fake,
        include_affinity_percentile_ranks=False,
    )
    results = predictor.predict(
        ["SIINXEKL"], n_flanks=["AAXAA"], c_flanks=["BBZBB"]
    )
    assert [r.peptide for r in results] == ["SIINXEKL"]
    assert fake.affinity_calls[-1]["peptides"] == ["SIINXEKL"]
    assert fake.presentation_calls[-1]["peptides"] == ["SIINXEKL"]
    assert fake.presentation_calls[-1]["n_flanks"] == ["AAXAA"]
    assert fake.presentation_calls[-1]["c_flanks"] == ["BBZBB"]


@pytest.mark.parametrize(
    "kwargs, match",
    [
        (
            {"affinity_error": ValueError("No affinity models available")},
            "No affinity models available",
        ),
        (
            {"affinity_error": ValueError("No sequences for allele(s): HLA-A*99:99.")},
            "No sequences for allele",
        ),
        (
            {
                "presentation_error": ValueError(
                    "Specify both or neither of n_flanks, c_flanks"
                )
            },
            "Specify both or neither",
        ),
        ({"bad_presentation_rows": True}, "presentation row"),
        ({"duplicate_presentation_row": True}, "duplicate presentation row"),
        (
            {"unrequested_presentation_allele": "HLA-B*07:02"},
            "unrequested model allele",
        ),
    ],
)
def test_mhcflurry_does_not_misclassify_config_length_or_output_errors(kwargs, match):
    fake = _make_fake_mhcflurry_predictor(**kwargs)
    predictor = MHCflurry(
        alleles=["HLA-A*02:01"],
        predictor=fake,
        include_affinity_percentile_ranks=False,
    )
    with pytest.raises(ValueError, match=match) as exc_info:
        predictor.predict(["SIINFEKL", "GILGFVFTL"])
    assert not isinstance(exc_info.value, UnsupportedResidueError)


def test_mhcflurry_length_error_remains_plain_value_error():
    fake = _make_fake_mhcflurry_predictor()
    predictor = MHCflurry(
        alleles=["HLA-A*02:01"],
        predictor=fake,
        include_affinity_percentile_ranks=False,
    )
    with pytest.raises(ValueError, match="lengths outside of supported range") as exc_info:
        predictor.predict(["SHORT"])
    assert not isinstance(exc_info.value, UnsupportedResidueError)


def test_translate_mhcflurry_residue_error_handles_custom_alphabet_fallback():
    native = ValueError("1 peptides have nonstandard amino acids: ['SIINFEKL']")
    translated = translate_mhcflurry_residue_error(native, ["GILGFVFTL", "SIINFEKL"])
    assert isinstance(translated, UnsupportedResidueError)
    assert translated.peptides == ("SIINFEKL",)
    assert translated.input_indices == (1,)
    assert translated.residues == ()

    already_typed = UnsupportedResidueError("already typed", peptides=("SIINXEKL",))
    assert translate_mhcflurry_residue_error(already_typed, ["SIINXEKL"]) is already_typed
    assert translate_mhcflurry_residue_error(ValueError("unrelated"), ["SIINXEKL"]) is None


# ---------------------------------------------------------------------------
# Adapter conformance tests
# ---------------------------------------------------------------------------


def test_adapter_conformance_commandline_and_base_validators(monkeypatch):
    monkeypatch.setattr(
        BaseCommandlinePredictor,
        "_determine_supported_alleles",
        staticmethod(lambda *args: {"HLA-A*02:01"}),
    )
    cmd = BaseCommandlinePredictor(
        program_name="fixture",
        alleles=["HLA-A*02:01"],
        parse_output_fn=lambda _: [],
        supported_alleles_flag="-list",
        input_file_flag="-f",
        length_flag="-l",
        allele_flag="-a",
    )
    with pytest.raises(UnsupportedResidueError) as exc_info:
        cmd.predict(["SIINXEKL"])
    assert exc_info.value.residues == ("X",)
    assert exc_info.value.peptides == ("SIINXEKL",)

    with pytest.raises(ValueError) as exc_info:
        cmd.predict(["SHORT"])
    assert not isinstance(exc_info.value, UnsupportedResidueError)


def test_adapter_conformance_across_sequence_predictors():
    # Calis (20 standard AAs)
    calis = Calis()
    with pytest.raises(UnsupportedResidueError) as exc_info:
        calis.predict(["SIINFEKL", "SIINXEKL"])
    assert exc_info.value.peptides == ("SIINXEKL",)
    assert exc_info.value.residues == ("X",)
    assert exc_info.value.input_indices == (1,)
    with pytest.raises(ValueError) as exc_info:
        calis.predict(["A"])
    assert not isinstance(exc_info.value, UnsupportedResidueError)

    # ERAMER (20 standard AAs)
    eramer = object.__new__(ERAMER)
    eramer.epitope_length = 8
    with pytest.raises(UnsupportedResidueError) as exc_info:
        eramer._check_peptides(["AAASIINXEKL"])
    assert exc_info.value.residues == ("X",)
    with pytest.raises(ValueError) as exc_info:
        eramer._check_peptides(["SHORT"])
    assert not isinstance(exc_info.value, UnsupportedResidueError)

    # CapHLA (20 standard AAs)
    caphla = object.__new__(CapHLA)
    with pytest.raises(UnsupportedResidueError) as exc_info:
        caphla._check_peptides(["SIINXEKL"])
    assert exc_info.value.residues == ("X",)
    with pytest.raises(ValueError) as exc_info:
        caphla._check_peptides(["SHORT"])
    assert not isinstance(exc_info.value, UnsupportedResidueError)

    # DeepTAP accepts X (21-AA alphabet) and rejects U with UnsupportedResidueError
    deeptap = object.__new__(DeepTAP)
    deeptap.max_peptide_length = 17
    deeptap._check_peptides(["SIINXEKL"])
    with pytest.raises(UnsupportedResidueError) as exc_info:
        deeptap._check_peptides(["SIINUEKL"])
    assert exc_info.value.residues == ("U",)
    with pytest.raises(ValueError) as exc_info:
        deeptap._check_peptides([""])
    assert not isinstance(exc_info.value, UnsupportedResidueError)

    # NetTCR (20 standard AAs)
    with pytest.raises(UnsupportedResidueError) as exc_info:
        _encode_feature(["SIINXEKL"], "pep")
    assert exc_info.value.residues == ("X",)
    assert exc_info.value.field == "pep"

    # MixTCRpred accepts X and '-' but rejects '*'
    _validate_tcr(TCR(cdr3a="CAVXG-F", cdr3b="CASSX-YF"))
    with pytest.raises(UnsupportedResidueError) as exc_info:
        _validate_tcr(TCR(cdr3a="CAV*GF", cdr3b="CASSYF"))
    assert exc_info.value.residues == ("*",)
    assert exc_info.value.field == "cdr3a"

    # PeptideInput / PeptiVerse / PlifePred2
    with pytest.raises(UnsupportedResidueError) as exc_info:
        PeptideInput("SIINXEKL")
    assert exc_info.value.residues == ("X",)
    assert exc_info.value.peptides == ("SIINXEKL",)
    with pytest.raises(ValueError) as exc_info:
        PeptideInput("siinfekl")
    assert not isinstance(exc_info.value, UnsupportedResidueError)

    peptiverse = object.__new__(PeptiVerse)
    peptiverse.max_peptide_length = 30
    with pytest.raises(UnsupportedResidueError) as exc_info:
        peptiverse._check_peptides(["SIINXEKL"])
    assert exc_info.value.residues == ("X",)

    plifepred2 = object.__new__(PlifePred2)
    with pytest.raises(UnsupportedResidueError) as exc_info:
        plifepred2._check_peptides(["SIINXEKLGGALQAKK"])
    assert exc_info.value.residues == ("X",)

    # BasePredictor subclasses using _check_peptide_inputs
    for cls in (MixMHCpred, MixMHC2pred, PRIME, SMM):
        assert cls._check_peptide_inputs is BasePredictor._check_peptide_inputs
