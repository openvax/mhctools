"""Single-step PWM mapping and provenance, independent of external weights."""

import hashlib

import pytest

from mhctools import CleavageInput, ERAMERCleavage


@pytest.fixture
def pwm(tmp_path):
    openpyxl = pytest.importorskip("openpyxl")
    workbook = openpyxl.Workbook()
    workbook.remove(workbook.active)
    for length in range(9, 17):
        sheet = workbook.create_sheet("Length %d" % length)
        sheet.append(["index", "aa"] + list(range(length)))
        for index, aa in enumerate("ACDEFGHIKLMNPQRSTVWY"):
            sheet.append([index, aa] + [length / 100] * length)
    path = tmp_path / "PWM.xlsx"
    workbook.save(path)
    return path


def test_one_step_score_and_loaded_artifact_identity(pwm):
    digest = hashlib.sha256(pwm.read_bytes()).hexdigest()
    predictor = ERAMERCleavage(pwm_path=pwm)
    # The in-memory snapshot retains its identity even if the file changes.
    pwm.write_bytes(b"changed after loading")
    result = predictor.predict(CleavageInput("LAAAFGAAA", source_start=40))
    assert result.model.version == "pwm-sha256:" + digest
    assert len(result.sites) == 1
    assert result.sites[0].bond == 1
    assert result.sites[0].score == pytest.approx(0.09)
    assert result.to_dict()["sites"][0]["source_bond"] == 41
    assert predictor.predict("A" * 16).sites[0].score == pytest.approx(0.16)


def test_erap_domain_and_terminal_chemistry(pwm):
    predictor = ERAMERCleavage(pwm_path=pwm)
    for peptide in ("A" * 8, "A" * 17, CleavageInput("A" * 9, n_term="acetylated")):
        result = predictor.predict(peptide)
        assert result.unsupported_reason
        assert not result.sites


def test_a_corrupt_pwm_file_fails_that_one_model_without_crashing_the_batch(tmp_path, monkeypatch):
    # Loading an external asset can raise almost any exception type (here,
    # openpyxl/zipfile's BadZipFile, which is neither ValueError, OSError
    # nor ImportError); this must still degrade to a "failed" prediction for
    # this one model, never abort predictions for every other model too.
    pytest.importorskip("openpyxl")
    from mhctools.benchmark import AssayMeasurement, predict_cleavage_measurements
    (tmp_path / "PWM.xlsx").write_bytes(b"not a real xlsx file at all")
    monkeypatch.setenv("ERAMER_HOME", str(tmp_path))
    m = AssayMeasurement(measurement_id="m1", source_measurement_id="m1", source="https://x",
        dataset="d", study="s", assay="a", sequence="RPPGFSPFR", chemistry="linear_L_free",
        endpoint="site_cleavage", units="binary", species="Homo sapiens", matrix="buffer",
        value=1, enzyme="XPNPEP1", bond=1)
    predictions = predict_cleavage_measurements([m], ["eramer-step", "app1-xp"])
    by_model = {p.model: p for p in predictions}
    assert by_model["eramer-step"].status == "failed"
    assert "zip file" in by_model["eramer-step"].reason
    assert by_model["app1-xp"].status == "scored"
    assert by_model["app1-xp"].value == 1


def test_cleavage_cli_reports_a_corrupt_asset_cleanly_not_a_traceback(tmp_path, capsys, monkeypatch):
    pytest.importorskip("openpyxl")
    from mhctools.cli.script import main
    (tmp_path / "PWM.xlsx").write_bytes(b"not a real xlsx file at all")
    monkeypatch.setenv("ERAMER_HOME", str(tmp_path))
    with pytest.raises(SystemExit) as excinfo:
        main(["cleavage", "--sequence", "RPPGFSPFR", "--model", "eramer-step"])
    assert excinfo.value.code == 2
    assert "zip file" in capsys.readouterr().err
