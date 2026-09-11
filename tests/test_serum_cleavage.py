"""Source-observation, activation, and chemistry regressions for serum candidates."""

from dataclasses import replace
from importlib.resources import files
import json

import pytest

from mhctools import CleavageInput, get_cleavage_model, predict_cleavage
from mhctools.benchmark import AssayMeasurement, evaluate_benchmark, predict_cleavage_measurements
from mhctools.cli.script import main


def test_ace_dipeptide_and_exceptional_chemistry():
    ace = get_cleavage_model("ace-dipeptidyl")
    result = ace.predict(CleavageInput("DRVYIHPFHL", source_start=40))
    assert [(s.bond, s.status) for s in result.sites] == [(8, "matched")]
    assert result.to_dict()["sites"][0]["source_bond"] == 48
    assert ace.predict("DRVYIHPF").sites[0].status == "not_matched"
    assert ace.predict(CleavageInput("SDKP", n_term="acetylated")).sites[0].bond == 2
    assert ace.predict(CleavageInput("RPKPQQFFGLM", c_term="amidated")).unsupported_reason
    for seq in ("AAE", "AAD", "APA"):
        assert ace.predict(seq).sites[0].status == "not_matched"


def test_mme_site_context_and_input_scope():
    mme = get_cleavage_model("mme-hydrophobic")
    result = mme.predict("YGGFL")
    assert next(s for s in result.sites if s.bond == 3).status == "matched"
    result = mme.predict(CleavageInput("RPKPQQFFGLM", c_term="amidated"))
    assert {s.bond for s in result.sites if s.status == "matched"} == {6, 7, 9}
    assert mme.predict("A" * 31).unsupported_reason
    assert mme.predict(CleavageInput("YGGFL", c_term="unknown")).unsupported_reason


@pytest.mark.parametrize("state", [None, "unknown", "zymogen", "inactive"])
def test_cpb2_without_active_enzyme_abstains(state):
    result = get_cleavage_model("cpb2-basic", enzyme_state=state).predict("YFPGQFAFSK")
    assert result.unsupported_reason
    assert not result.sites
    assert dict(result.conditions)["enzyme_state"] == (state or "unknown")


def test_cpb2_active_is_separate_from_cpn():
    results = predict_cleavage("YFPGQFAFSK", models=["cpn-basic", "cpb2-basic"], enzyme_states={"CPB2": "active"})
    assert all([(s.bond, s.status) for s in r.sites] == [(9, "matched")] for r in results)
    assert results[0].conditions == ()
    assert dict(results[1].conditions) == {"enzyme_state": "active"}
    assert get_cleavage_model("cpb2-basic", enzyme_state="active").predict(
        CleavageInput("YFPGQFAFSK", c_term="amidated")).unsupported_reason
    with pytest.raises(ValueError, match="enzyme_state"):
        get_cleavage_model("cpb2-basic", enzyme_state="activated-ish")
    with pytest.raises(ValueError, match="only"):
        get_cleavage_model("cpn-basic", enzyme_state="active")
    with pytest.raises(ValueError, match="without selecting"):
        predict_cleavage("AK", models="cpn-basic", enzyme_states={"CPB2": "active"})
    with pytest.raises(ValueError, match="only"):
        predict_cleavage("AK", enzyme_states={"CPN1": "active"})


def test_serum_reference_observations_and_abstentions():
    data = json.loads(files("mhctools").joinpath("data/serum_cleavage_reference.json").read_text())
    measurements = [AssayMeasurement(**m) for m in data["measurements"]]
    predictions = predict_cleavage_measurements(measurements, data["models"])
    report = evaluate_benchmark(measurements, predictions, evaluation="reproduction")
    comparable = [r for r in report["records"] if r["exclusion"] is None]
    assert len(comparable) == 15
    assert sum(r["measurement"]["value"] == 0 for r in comparable) == 2
    assert all(r["measurement"]["value"] == r["prediction"]["value"] for r in comparable)
    cpb2 = next(m for m in measurements if m.enzyme == "CPB2")
    inactive = replace(cpb2, conditions={"enzyme_state": "zymogen"})
    prediction = predict_cleavage_measurements([inactive], ["cpb2-basic"])[0]
    assert prediction.value is None and "zymogen" in prediction.reason


def test_serum_cli_preserves_activation_and_reference_report(capsys):
    main(["cleavage", "--sequence", "YFPGQFAFSK", "--model", "cpb2-basic", "--enzyme-state", "CPB2=active"])
    result = json.loads(capsys.readouterr().out)["results"][0]
    assert result["conditions"] == [["enzyme_state", "active"]]
    assert result["sites"][0]["bond"] == 9
    main(["benchmark", "--reference-cleavage", "serum"])
    report = json.loads(capsys.readouterr().out)
    assert report["evaluation"] == "reproduction"
    assert all(d["measurement_count"] == 0 for d in report["requested_domains"])
    with pytest.raises(SystemExit):
        main(["cleavage", "--sequence", "AK", "--enzyme-state", "CPB2=active", "--enzyme-state", "CPB2=zymogen"])


def test_enzyme_state_gate_is_driven_by_requires_activation_not_a_name():
    # get_cleavage_model checks each rule's own requires_activation flag,
    # not a hardcoded "cpb2-basic" literal, so this generalizes correctly
    # if a second activation-gated enzyme is ever curated.
    from mhctools import get_cleavage_model
    with pytest.raises(ValueError, match="only supported for models that require activation"):
        get_cleavage_model("cpn-basic", enzyme_state="active")
    with pytest.raises(ValueError, match="not supported for dpp4-qpisa"):
        get_cleavage_model("dpp4-qpisa", enzyme_state="active")
    with pytest.raises(ValueError, match="not supported for source-reference models"):
        get_cleavage_model("thop1-observed", enzyme_state="active")
    # cpb2-basic itself is unaffected by the generalization.
    active = get_cleavage_model("cpb2-basic", enzyme_state="active")
    assert dict(active.predict("YFPGQFAFSK").conditions)["enzyme_state"] == "active"
