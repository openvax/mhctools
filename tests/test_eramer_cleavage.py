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
