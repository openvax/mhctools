"""Scientific topology, abstention, and public reporting regressions."""

import json

import pytest

from mhctools import (
    CleavageInput, cleavage_models, get_cleavage_model, predict_cleavage,
)
from mhctools.cli.script import main


@pytest.mark.parametrize("model,sequence,bond", [
    ("cpn-basic", "RPPGFSPFR", 8),  # C-terminal Arg of bradykinin
    ("app2-xp", "RPPGFSPFR", 1),  # first residue, NOT the Pro at position 2
    ("fap-dipeptidyl", "APGAA", 2),
    ("fap-endo-gp", "TSGPNQ", 4),  # alpha-2-antiplasmin recognition context
    ("enpep-acidic", "DRVYIHPF", 1),  # angiotensin II N-terminal Asp
    ("anpep-ala", "AAGF", 1),
    ("dpp8-xp-xa", "APGAA", 2),
    ("dpp9-xp-xa", "VPYGSFKHV", 2),  # RU1 epitope N-terminal dipeptide
    ("prep-pro", "RPPGFSPFR", 7),
    ("erap2-basic", "RSLYNTVATL", 1),
])
def test_recognition_bonds(model, sequence, bond):
    result = get_cleavage_model(model).predict(CleavageInput(sequence, source_start=20))
    assert result.unsupported_reason is None
    site = next(s for s in result.sites if s.bond == bond)
    assert site.status == "matched"
    assert site.score is None
    assert result.model.evidence == "motif_rule"
    assert result.model.score_units is None
    assert any(s["source_bond"] == 20 + bond for s in result.to_dict()["sites"])


@pytest.mark.parametrize("model,sequence", [
    ("cpn-basic", "RPPGFSPFA"),
    ("app2-xp", "AAPG"),
    ("fap-dipeptidyl", "APPG"),
    ("fap-endo-gp", "TSGPPQ"),
    ("enpep-acidic", "AAGF"),
    ("anpep-ala", "RAGF"),
    ("dpp8-xp-xa", "VPPGSFKHV"),
    ("dpp9-xp-xa", "VPPGSFKHV"),
    ("prep-pro", "AAAA"),
    ("erap2-basic", "ASLYNTVATL"),
])
def test_nonmatches_retain_rule_semantics(model, sequence):
    result = get_cleavage_model(model).predict(sequence)
    assert result.sites
    assert {s.status for s in result.sites} == {"not_matched"}
    assert all(s.score is None for s in result.sites)
    assert result.unsupported_reason is None


def test_blocked_termini_distinguish_endo_and_exo():
    peptide = CleavageInput("GPNQ", n_term="acetylated")
    assert get_cleavage_model("fap-dipeptidyl").predict(peptide).unsupported_reason
    result = get_cleavage_model("fap-endo-gp").predict(peptide)
    assert next(s for s in result.sites if s.bond == 2).status == "matched"
    assert get_cleavage_model("cpn-basic").predict(
        CleavageInput("AARK", c_term="amidated")).unsupported_reason


def test_all_models_abstain_for_unknown_chemistry_and_too_short_input():
    for model in cleavage_models():
        predictor = get_cleavage_model(model.name)
        for peptide in (CleavageInput("AAPRG", n_term="unknown"), CleavageInput("A")):
            result = predictor.predict(peptide)
            assert result.unsupported_reason, model.name
            assert not result.sites


def test_internal_dpp_motif_requires_conditional_fragment():
    peptide = CleavageInput("GGVPYGSFKHV", source_id="precursor", source_start=5)
    initial = get_cleavage_model("dpp9-xp-xa").predict(peptide)
    assert [(s.bond, s.status) for s in initial.sites] == [(2, "not_matched")]
    fragment = peptide.fragment(2, 11, n_term="free", c_term="free")
    followup = get_cleavage_model("dpp9-xp-xa").predict(fragment)
    assert followup.to_dict()["sites"][0]["source_bond"] == 9
    assert followup.sites[0].status == "matched"


def test_prep_length_scope():
    assert get_cleavage_model("prep-pro").predict("A" * 29 + "PA").unsupported_reason
    assert get_cleavage_model("prep-pro").predict("PAA").unsupported_reason


def test_compartment_filter_and_explicit_selection():
    results = predict_cleavage("VPYGSFKHV", compartment="cytosol")
    assert {r.model.enzyme for r in results} == {"DPP8", "DPP9", "PREP", "TPP2", "NPEPPS", "XPNPEP1", "THOP1", "NLN"}
    assert len(predict_cleavage("HAE", models=["dpp4-qpisa", "dpp4-qpisa"])) == 1
    with pytest.raises(ValueError, match="not annotated"):
        predict_cleavage("HAE", models="dpp4-qpisa", compartment="er")
    with pytest.raises(ValueError, match="Unknown compartment"):
        predict_cleavage("HAE", compartment="blood")
    with pytest.raises(ValueError, match="Unknown cleavage model"):
        predict_cleavage("HAE", models="dpp3")


def test_cli_json_preserves_scores_chemistry_and_provenance(capsys, tmp_path):
    main(["cleavage", "--sequence", "HAE", "--model", "dpp4-qpisa",
          "--source-start", "10", "--source-id", "construct"])
    data = json.loads(capsys.readouterr().out)
    result = data["results"][0]
    assert data["schema_version"] == 1
    assert result["sites"][0]["score"] == pytest.approx(2.1694)
    assert result["sites"][0]["source_bond"] == 12
    assert result["model"]["references"]
    path = tmp_path / "evidence.json"
    main(["cleavage", "--sequence", "HAE", "--model", "dpp4-qpisa",
          "--n-term", "acetylated", "--out", str(path)])
    result = json.loads(path.read_text())["results"][0]
    assert result["peptide"]["n_term"] == "acetylated"
    assert result["unsupported_reason"]
    assert result["sites"] == []


def test_cli_lists_optional_models_without_loading_assets(capsys, monkeypatch):
    monkeypatch.setenv("ERAMER_HOME", "/does-not-exist")
    main(["cleavage", "--list-models"])
    data = json.loads(capsys.readouterr().out)
    assert len(data["models"]) == 21
    assert any(m["name"] == "eramer-step" for m in data["models"])


@pytest.mark.parametrize("args", [[], ["--sequence", "HAX"],
    ["--sequence", "HAE", "--model", "unknown"],
    ["--list-models", "--sequence", "HAE"]])
def test_cli_invalid_requests_fail(args):
    with pytest.raises(SystemExit) as error:
        main(["cleavage"] + args)
    assert error.value.code == 2


def test_position_track_spans_internal_scan_and_fragment_cascade():
    """Locks in the coordinate contract a downstream tool (e.g. vaxrank)
    needs to overlay cleavage evidence from several models as one track
    indexed by position in a full parent sequence.

    Internal-topology models (here, neprilysin) assess every bond of
    whatever peptide they are given in a single predict() call. Terminal-
    topology models (here, CPN1) only ever assess the CURRENTLY exposed
    end of their input; to ask about an internal position of a longer
    precursor, the caller models that trimming step explicitly via
    ``fragment()`` and reads the result's ``source_bond``, which still maps
    back to the parent's absolute coordinates. Both modes must resolve into
    the same coordinate system so their evidence can share one track.
    """
    parent = CleavageInput("RPPGFSPFRSSRQ", source_id="precursor")
    track = {}

    def record(result):
        for site in result.to_dict()["sites"]:
            track.setdefault(site["source_bond"], []).append((result.model.name, site["status"]))

    # Mode 1: internal topology scans the whole input directly.
    internal = get_cleavage_model("mme-hydrophobic").predict(parent)
    assert internal.unsupported_reason is None
    assert len(internal.sites) == len(parent.sequence) - 1
    record(internal)
    assert track[4] == [("mme-hydrophobic", "matched")]
    assert track[7] == [("mme-hydrophobic", "matched")]
    assert track[1] == [("mme-hydrophobic", "not_matched")]

    # Mode 2: terminal topology requires an explicit fragment cascade step.
    # Here, a hypothesized prior trimming exposes residues 10-13 (SSRQ).
    fragment = parent.fragment(9, len(parent.sequence), n_term="free", c_term="free")
    trimmed = get_cleavage_model("cpn-basic").predict(fragment)
    assert trimmed.unsupported_reason is None
    record(trimmed)
    # The fragment's local bond 3 (its own last internal bond) lands on the
    # parent's absolute position 12 -- the same coordinate space as mode 1.
    assert trimmed.sites[0].bond == 3
    assert track[12] == [("mme-hydrophobic", "not_matched"), ("cpn-basic", "not_matched")]

    # One combined track, addressable by absolute parent position, carrying
    # evidence from models with entirely different assessment strategies.
    assert set(track) == set(range(1, len(parent.sequence)))


def test_compartment_help_text_matches_the_real_compartment_set(capsys):
    # The --compartment help text is a hand-maintained string, not derived
    # from cleavage_models() at runtime (that would force loading the whole
    # motif/reference-catalog panel on every CLI invocation just to render
    # --help). This test is the guardrail instead: it fails the moment a
    # compartment is added or removed without updating the help text.
    from mhctools.peptidases import cleavage_models
    real_compartments = {c for m in cleavage_models() for c in m.compartments}
    with pytest.raises(SystemExit):
        main(["cleavage", "--help"])
    help_text = " ".join(capsys.readouterr().out.split())
    marker = "Filter enzyme locations:"
    assert marker in help_text
    listed_text = help_text.split(marker, 1)[1].split(" --", 1)[0]
    listed = {token.strip(" ,.") for token in listed_text.split(",")}
    assert listed == real_compartments
