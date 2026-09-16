"""Model provenance follows the loaded weights through every prediction door."""

import sys
import types

import pandas as pd
import pytest

from mhctools import MHCflurry, MHCflurry_Affinity
from mhctools import mhcflurry as wrapper


MHCFLURRY_TWINS = (MHCflurry, MHCflurry_Affinity)
PREDICTION_DOORS = (
    "predict_dataframe", "predict_proteins_dataframe", "predict_peptides_dataframe",
)


def fake_backend():
    alleles = ["HLA-A*02:01"]

    def affinity(peptides, alleles, include_percentile_ranks=True):
        return pd.DataFrame(dict(
            peptide=peptides, allele=alleles, prediction=[12.5] * len(peptides),
            prediction_percentile=[1.5] * len(peptides),
        ))

    def presentation(peptides, alleles, **kwargs):
        return pd.DataFrame(dict(
            peptide=peptides, peptide_num=range(len(peptides)),
            best_allele=["HLA-A*02:01"] * len(peptides),
            presentation_score=[0.75] * len(peptides),
            presentation_percentile=[2.0] * len(peptides),
            processing_score=[0.25] * len(peptides),
        ))

    aff = types.SimpleNamespace(
        supported_alleles=alleles, predict_to_dataframe=affinity,
    )
    return types.SimpleNamespace(
        supported_alleles=alleles, affinity_predictor=aff, predict=presentation,
    )


@pytest.fixture
def installed_models(monkeypatch, tmp_path):
    package = types.ModuleType("mhcflurry")
    package.__version__ = "2.2.1"
    downloads = types.ModuleType("mhcflurry.downloads")
    state = {"release": "2.2.0", "loads": []}
    downloads.get_current_release = lambda: state["release"]
    downloads.get_path = lambda name, subdir, test_exists=True: str(
        tmp_path / str(state["release"]) / name / subdir)
    downloads.get_default_class1_presentation_models_dir = lambda test_exists=True: downloads.get_path(
        "models_class1_presentation", "models")
    downloads.get_default_class1_models_dir = lambda test_exists=True: downloads.get_path(
        "models_class1_pan", "models.combined")

    def load(kind, path):
        state["loads"].append((kind, path))
        backend = fake_backend()
        return backend if kind == "presentation" else backend.affinity_predictor

    package.Class1PresentationPredictor = types.SimpleNamespace(
        load=lambda path=None: load("presentation", path))
    package.Class1AffinityPredictor = types.SimpleNamespace(
        load=lambda path=None: load("affinity", path))
    package.downloads = downloads
    monkeypatch.setitem(sys.modules, "mhcflurry", package)
    monkeypatch.setitem(sys.modules, "mhcflurry.downloads", downloads)
    monkeypatch.setattr(wrapper, "_model_cache", {})
    return package, downloads, state


@pytest.mark.parametrize("cls", MHCFLURRY_TWINS)
def test_default_model_version_reaches_every_prediction_door(cls, installed_models):
    model = cls(alleles=["HLA-A*02:01"], default_peptide_lengths=[9])
    assert model.predictor_version == "2.2.1+release-2.2.0"
    frames = []
    for door in PREDICTION_DOORS:
        inputs = {"protein": "SIINFEKLA"} if "proteins" in door else ["SIINFEKLA"]
        if door == "predict_peptides_dataframe":
            with pytest.warns(DeprecationWarning):
                frame = getattr(model, door)(inputs)
        else:
            frame = getattr(model, door)(inputs)
        assert set(frame.predictor_version) == {model.predictor_version}
        frames.append(frame)
    for frame in frames[1:]:
        columns = ["peptide", "allele", "kind", "value", "score", "percentile_rank"]
        pd.testing.assert_frame_equal(frames[0][columns], frame[columns])


@pytest.mark.parametrize("cls", MHCFLURRY_TWINS)
@pytest.mark.parametrize("source", ["injected", "path", "environment", "no_release"])
@pytest.mark.parametrize("version", [None, "custom-weights-7"])
def test_custom_models_are_never_labeled_as_default_release(
    cls, source, version, installed_models, tmp_path,
):
    _, downloads, state = installed_models
    options = {"predictor_version": version}
    if source == "injected":
        backend = fake_backend()
        options["predictor"] = backend if cls is MHCflurry else backend.affinity_predictor
    elif source == "path":
        options["models_path"] = str(tmp_path / "custom")
    elif source == "environment":
        downloads.get_default_class1_presentation_models_dir = lambda: str(tmp_path / "custom")
        downloads.get_default_class1_models_dir = lambda: str(tmp_path / "custom")
    else:
        state["release"] = None
    model = cls(alleles=["HLA-A*02:01"], **options)
    assert model.predictor_version == version
    frame = model.predict_dataframe(["SIINFEKLA"])
    if version is None:
        assert frame.predictor_version.isna().all()
    else:
        assert set(frame.predictor_version) == {version}


@pytest.mark.parametrize("cls", MHCFLURRY_TWINS)
def test_cached_object_keeps_its_version_and_release_switch_loads_new_weights(cls, installed_models):
    package, _, state = installed_models
    first = cls(alleles=["HLA-A*02:01"])
    package.__version__ = "3.0.0"
    reused = cls(alleles=["HLA-A*02:01"])
    assert reused.predictor is first.predictor
    assert reused.predictor_version == first.predictor_version == "2.2.1+release-2.2.0"
    assert len(state["loads"]) == 1
    state["release"] = "3.0.0"
    switched = cls(alleles=["HLA-A*02:01"])
    assert switched.predictor is not first.predictor
    assert switched.predictor_version == "3.0.0+release-3.0.0"
    assert first.predictor_version == "2.2.1+release-2.2.0"
    assert len(state["loads"]) == 2


@pytest.mark.parametrize("cls", MHCFLURRY_TWINS)
@pytest.mark.parametrize("version", ["", " ", 12])
def test_explicit_version_must_be_a_nonempty_string(cls, version, installed_models):
    with pytest.raises(ValueError, match="predictor_version"):
        cls(alleles=["HLA-A*02:01"], predictor_version=version)


def test_public_composite_version_checks_actual_model_path(installed_models, tmp_path):
    from mhctools import mhcflurry_composite_version

    assert mhcflurry_composite_version() == "2.2.1+release-2.2.0"
    with pytest.raises(RuntimeError, match="custom"):
        mhcflurry_composite_version(str(tmp_path / "custom"))


def test_legacy_conversion_accepts_provenance_without_changing_legacy_schema():
    from mhctools import BindingPrediction, BindingPredictionCollection

    binding = BindingPrediction(peptide="SIINFEKLA", allele="HLA-A*02:01", affinity=12.5)
    original = binding.to_dict()
    collection = BindingPredictionCollection([binding])
    outputs = [
        binding.to_pred(predictor_version="experiment-1"),
        collection.to_preds(predictor_version="experiment-1")[0],
        collection.to_peptide_preds(predictor_version="experiment-1")[0].preds[0],
    ]
    assert all(pred.predictor_version == "experiment-1" for pred in outputs)
    assert binding.to_dict() == original
    assert binding.to_pred().predictor_version == ""
    assert collection.to_preds()[0].predictor_version == ""
    assert collection.to_peptide_preds()[0].preds[0].predictor_version == ""


def test_unversioned_legacy_predictor_retains_the_original_empty_string():
    from mhctools import RandomBindingPredictor

    model = RandomBindingPredictor(alleles=["HLA-A*02:01"])
    assert model.predict_dataframe(["SIINFEKLA"]).predictor_version.tolist() == [""]


@pytest.mark.parametrize("custom", [False, True])
def test_affinity_fallback_uses_the_presentation_objects_actual_provenance(
    installed_models, tmp_path, custom,
):
    _, downloads, state = installed_models

    def absent_standalone(**kwargs):
        raise RuntimeError("standalone affinity bundle missing")

    downloads.get_default_class1_models_dir = absent_standalone
    if custom:
        downloads.get_default_class1_presentation_models_dir = \
            lambda **kwargs: str(tmp_path / "custom")
    presentation = MHCflurry(alleles=["HLA-A*02:01"])
    affinity = MHCflurry_Affinity(alleles=["HLA-A*02:01"])
    assert affinity.predictor is presentation.predictor.affinity_predictor
    expected = None if custom else "2.2.1+release-2.2.0"
    assert affinity.predictor_version == presentation.predictor_version == expected
    assert len(state["loads"]) == 1
    assert affinity.predict_dataframe(["SIINFEKLA"]).value.tolist() == [12.5]
