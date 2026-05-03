# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

import json
import subprocess

import pytest

from mhctools.pepsickle import Pepsickle
from mhctools.processing_predictor import (
    score_cterm,
    score_cterm_anti_max_internal,
    score_nterm_cterm_anti_max_internal,
)
from mhctools.proteasome_predictor import ProteasomePredictor
from mhctools.pred import Kind, PeptideResult, COLUMNS

pepsickle = pytest.importorskip("pepsickle")


@pytest.fixture(scope="module")
def predictor():
    return Pepsickle()


PROTEIN = "MFVFLVLLPLVSSQCVNLTTRTQLPPAYTNSFTRGVY"
PEPTIDES = ["SIINFEKL", "GILGFVFTL", "NLVPMVATV"]


# -- init --

def test_init_defaults():
    p = Pepsickle()
    assert p.threshold == 0.5
    assert p.human_only is False
    assert p.isolate_subprocess is False
    assert p.default_peptide_lengths == [9]
    assert p.scoring is score_cterm_anti_max_internal
    assert not hasattr(p, "alleles")


def test_is_proteasome_predictor():
    assert isinstance(Pepsickle(), ProteasomePredictor)


def test_init_custom_scoring_callable():
    p = Pepsickle(scoring=score_cterm)
    assert p.scoring is score_cterm


def test_init_custom_scoring_string():
    p = Pepsickle(scoring="cterm")
    assert p.scoring is score_cterm


def test_str():
    s = str(Pepsickle())
    assert "Pepsickle" in s


# -- cleavage_probs --

def test_cleavage_probs(predictor):
    probs = predictor.cleavage_probs(PROTEIN)
    assert len(probs) == len(PROTEIN)
    assert all(0.0 <= p <= 1.0 for p in probs)


def test_isolated_cleavage_probs_matches_in_process():
    in_process = Pepsickle().cleavage_probs(PROTEIN)
    isolated = Pepsickle(isolate_subprocess=True).cleavage_probs(PROTEIN)
    assert isolated == pytest.approx(in_process)


def test_isolated_cleavage_probs_batches_unique_sequences(monkeypatch):
    calls = []

    def fake_run(args, input, text, capture_output, timeout):
        calls.append({
            "args": args,
            "request": json.loads(input),
            "text": text,
            "capture_output": capture_output,
            "timeout": timeout,
        })
        results = {
            sequence: [0.25] * len(sequence)
            for sequence in calls[-1]["request"]["sequences"]
        }
        return subprocess.CompletedProcess(
            args=args,
            returncode=0,
            stdout=json.dumps({"results": results}),
            stderr="",
        )

    monkeypatch.setattr("mhctools.pepsickle.subprocess.run", fake_run)

    predictor = Pepsickle(
        threshold=0.25,
        human_only=True,
        isolate_subprocess=True,
        subprocess_timeout=7,
    )
    result = predictor.cleavage_probs_many([PROTEIN, PROTEIN, "AC"])

    assert result[PROTEIN] == [0.25] * len(PROTEIN)
    assert result["AC"] == [0.25, 0.25]
    assert len(calls) == 1
    assert calls[0]["request"] == {
        "human_only": True,
        "threshold": 0.25,
        "sequences": [PROTEIN, "AC"],
    }
    assert calls[0]["text"] is True
    assert calls[0]["capture_output"] is True
    assert calls[0]["timeout"] == 7


# -- predict --

def test_predict_returns_peptide_preds(predictor):
    results = predictor.predict(PEPTIDES)
    assert isinstance(results, list)
    assert len(results) == len(PEPTIDES)
    for pp in results:
        assert isinstance(pp, PeptideResult)
        assert len(pp.preds) == 1
        pred = pp.preds[0]
        assert pred.kind == Kind.proteasome_cleavage
        assert 0.0 <= pred.score <= 1.0
        assert pred.predictor_name == "pepsickle"


def test_predict_with_flanks(predictor):
    # With flanking context, C-terminal should get non-zero scores
    results = predictor.predict(
        ["SIINFEKL"],
        n_flanks=["MFVFLVLL"],
        c_flanks=["PLVSSQCV"],
    )
    pred = results[0].preds[0]
    assert pred.n_flank == "MFVFLVLL"
    assert pred.c_flank == "PLVSSQCV"
    assert pred.kind == Kind.proteasome_cleavage


def test_predict_dataframe(predictor):
    df = predictor.predict_dataframe(PEPTIDES)
    assert list(df.columns) == list(COLUMNS)
    assert len(df) == len(PEPTIDES)
    assert (df["kind"] == "proteasome_cleavage").all()
    assert (df["predictor_name"] == "pepsickle").all()


# -- predict_proteins --

def test_predict_proteins(predictor):
    result = predictor.predict_proteins({"spike": PROTEIN})
    assert "spike" in result
    preds_list = result["spike"]
    assert all(isinstance(pp, PeptideResult) for pp in preds_list)
    expected_count = len(PROTEIN) - 9 + 1
    assert len(preds_list) == expected_count
    for pp in preds_list:
        pred = pp.preds[0]
        assert pred.kind == Kind.proteasome_cleavage
        assert pred.source_sequence_name == "spike"
        assert 0.0 <= pred.score <= 1.0


def test_predict_proteins_scores_vary(predictor):
    result = predictor.predict_proteins({"spike": PROTEIN})
    scores = [pp.preds[0].score for pp in result["spike"]]
    assert len(set(round(s, 4) for s in scores)) > 1


def test_predict_proteins_with_flank_length(predictor):
    result = predictor.predict_proteins({"spike": PROTEIN}, flank_length=5)
    preds_list = result["spike"]
    # Check that flanks are recorded
    mid = preds_list[len(preds_list) // 2].preds[0]
    assert len(mid.n_flank) > 0
    assert len(mid.c_flank) > 0


def test_predict_proteins_dataframe(predictor):
    df = predictor.predict_proteins_dataframe(
        {"spike": PROTEIN}, sample_name="sample1")
    assert list(df.columns) == list(COLUMNS)
    assert (df["source_sequence_name"] == "spike").all()
    assert (df["sample_name"] == "sample1").all()


def test_predict_proteins_string_input(predictor):
    result = predictor.predict_proteins(PROTEIN)
    assert "seq" in result


def test_predict_proteins_multiple_lengths(predictor):
    p = Pepsickle(default_peptide_lengths=[8, 9, 10])
    result = p.predict_proteins({"s": PROTEIN})
    expected = sum(len(PROTEIN) - k + 1 for k in [8, 9, 10])
    assert len(result["s"]) == expected


# -- predict_cleavage_sites --

def test_predict_cleavage_sites(predictor):
    result = predictor.predict_cleavage_sites({"spike": PROTEIN})
    probs = result["spike"]
    assert len(probs) == len(PROTEIN)
    assert all(0.0 <= p <= 1.0 for p in probs)


# -- scoring --

def test_scoring_cterm_only():
    p = Pepsickle(scoring=score_cterm)
    result = p.predict_proteins({"s": PROTEIN})
    probs = p.cleavage_probs(PROTEIN)
    # First 9-mer: offset=0, c_term = probs[8]
    first = result["s"][0].preds[0]
    assert first.score == pytest.approx(probs[8])


def test_scoring_methods_produce_different_scores():
    fns = [score_cterm, score_cterm_anti_max_internal,
           score_nterm_cterm_anti_max_internal]
    scores_by_fn = {}
    for fn in fns:
        p = Pepsickle(scoring=fn)
        result = p.predict_proteins({"s": PROTEIN})
        scores_by_fn[fn.__name__] = [
            pp.preds[0].score for pp in result["s"]]
    unique = set(
        tuple(round(s, 6) for s in v) for v in scores_by_fn.values())
    assert len(unique) > 1

