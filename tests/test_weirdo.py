"""Exercise the actual released WEIRDO model and its optional dependency boundary."""

import json
import subprocess
import sys

import numpy as np
import pandas as pd
import pytest

from mhctools import Kind, Prediction, WEIRDO
from mhctools.artifacts import artifact_status
from mhctools.cli.args import make_mhc_arg_parser, predictors_from_args
from mhctools.cli.script import main as cli_main
from mhctools.integrations import integration_status
from mhctools.pred import best_direction


@pytest.fixture(scope="module")
def owner():
    weirdo = pytest.importorskip("weirdo")
    return weirdo.load_pretrained_model()


def test_mhctools_import_does_not_import_optional_weirdo():
    code = """
import builtins
original = builtins.__import__
def guarded(name, globals=None, locals=None, fromlist=(), level=0):
    if level == 0 and (name == 'weirdo' or name.startswith('weirdo.')):
        raise ImportError('optional backend deliberately absent')
    return original(name, globals, locals, fromlist, level)
builtins.__import__ = guarded
import mhctools
try:
    mhctools.WEIRDO()
except ImportError as error:
    assert 'mhctools[weirdo]' in str(error)
else:
    raise AssertionError('missing optional dependency was not reported')
"""
    subprocess.run([sys.executable, "-c", code], check=True, capture_output=True, text=True)


def test_reference_output_duplicates_and_provenance(owner):
    peptides = ["acdefghi", "MTMDKSELV", "ACDEFGHI"]
    predictor = WEIRDO(batch_size=1)
    results = predictor.predict(iter(peptides))
    expected = owner.predict_dataframe(peptides, batch_size=1)
    assert [r.preds[0].peptide for r in results] == expected.peptide.tolist()
    np.testing.assert_array_equal([r.preds[0].score for r in results], expected.foreignness)
    for result in results:
        prediction = result.foreignness
        assert prediction.kind == Kind.foreignness
        assert prediction.allele == ""
        assert prediction.value is None
        assert prediction.predictor_name == "weirdo"
        assert owner.prediction_identity in prediction.predictor_version
        assert prediction.measurement_context.score_semantics == owner.score_policy
        assert Prediction.from_dict(json.loads(json.dumps(prediction.to_dict()))) == prediction
    frame = predictor.predict_dataframe(peptides)
    assert frame.allele.tolist() == ["", "", ""]
    assert set(frame.predictor_version) == {predictor.predictor_version}
    assert predictor.kind_support() == {
        Kind.foreignness: {"mhc_dependence": "none", "mhc_class": "none"},
    }
    assert best_direction(Kind.foreignness, "score") == "max"


def test_empty_batch(owner):
    predictor = WEIRDO()
    assert predictor.predict([]) == []
    assert predictor.predict_dataframe([]).empty


@pytest.mark.parametrize("peptides", [[""], ["ACDEFGHX"], ["ACDEFGHI "],
                                      ["A" * 12], [None], [123]])
@pytest.mark.parametrize("method", ["predict", "predict_dataframe"])
def test_invalid_inputs_match_owner_contract(owner, peptides, method):
    with pytest.raises(ValueError):
        owner.score(peptides)
    with pytest.raises(ValueError):
        getattr(WEIRDO(), method)(peptides)


@pytest.mark.parametrize("options", [{"host_taxon_id": 10090}, {"score_policy": "other"},
                                    {"model": "other"}, {"batch_size": 0},
                                    {"default_peptide_lengths": [7, 9]}])
def test_unsupported_settings_fail(owner, options):
    with pytest.raises(ValueError):
        WEIRDO(**options)


def test_foreignness_cannot_be_projected_onto_an_allele():
    with pytest.raises(ValueError, match="MHC-independent"):
        Prediction(kind=Kind.foreignness, score=0.5, peptide="ACDEFGHI", allele="HLA-A*02:01")


def test_cli_constructor_requires_no_alleles(owner):
    args = make_mhc_arg_parser().parse_args(["--mhc-predictor", "weirdo", "--mhc-peptide-lengths", "9"])
    predictor, = predictors_from_args(args)
    assert isinstance(predictor, WEIRDO)
    assert predictor.default_peptide_lengths == (9,)
    assert predictor.predict(["MTMDKSELV"])[0].foreignness is not None


def test_artifact_and_reproduced_capability(owner):
    status = artifact_status("weirdo")
    assert status.status == "ready"
    assert status.manager == "weirdo package"
    assert status.fetchable is False
    assert owner.manifest["artifact_id"] in status.version
    assert integration_status("weirdo", check="reproduced", timeout=30).reproduced is True


def test_documented_table_command_scores_real_model(owner, tmp_path):
    source = tmp_path / "peptides.csv"
    output = tmp_path / "scored.csv"
    peptides = ["ACDEFGHI", "MTMDKSELV", "ACDEFGHI"]
    pd.DataFrame({"peptide": peptides}).to_csv(source, index=False)
    cli_main(["predict-table", "--input", str(source), "--out", str(output),
              "--predictor", "weirdo:weirdo_foreignness:foreignness"])
    frame = pd.read_csv(output)
    assert frame.peptide.tolist() == peptides
    np.testing.assert_allclose(frame.weirdo_foreignness, owner.score(peptides), atol=1e-6)
