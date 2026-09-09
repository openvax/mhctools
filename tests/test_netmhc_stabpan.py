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

import os

from numpy.testing import assert_allclose
from mhctools import NetMHCstabpan
from mhctools.base_commandline_predictor import BaseCommandlinePredictor
from mhctools.binding_prediction_collection import BindingPredictionCollection


DEFAULT_ALLELE = 'HLA-A*02:01'

# Protein sequences to test with against known values of netMHCstabpan web server.
protein_sequences = [
    "ALDKNLHQL",
    "ALHEEVVCV",
    "ALPPTVYEV",
    "AVLGSFSYV",
    "EMASVLFKA",
]

web_server_predictions = [
    4.6393,
    10.6011,
    11.8201,
    4.6722,
    1.8404,
]

# REMINDER: a program called "netMHCstabpan" MUST be installed and working for
# this test suite to succeeed. Also all peptides must be the same length.


def test_netmhc_stabpan_accuracy():    
    # Check that the netMHCstabpan program is working and returning th eexpected outputs.
    predictor = NetMHCstabpan(
        alleles=[DEFAULT_ALLELE], program_name='netMHCstabpan')

    binding_predictions = predictor.predict_peptides(protein_sequences)
    stability_predictions = [p.score for p in binding_predictions]
    rank_predictions = [p.percentile_rank for p in binding_predictions]

    assert len(web_server_predictions) == len(binding_predictions)
    assert len(stability_predictions) == len(binding_predictions)

    for prank in rank_predictions:
        # Make sure that correct mapping is done by checking percentiles aren't above 100.
        assert prank < 100
    
    for i, (expected, actual) in enumerate(zip(web_server_predictions, stability_predictions)):
        # Check to make sure that the stability predictions are within 0.01 of the webserver values.
        # This could be the result of different versions of dependencies or the nature of the ANN itself.
        assert_allclose(expected, actual, atol=0.01, err_msg="Peptide %d: expected %f but got %f" % (i, expected, actual))


def test_netmhc_stabpan_groups_mixed_length_peptides(monkeypatch):
    def fake_collect(self, commands, input_filenames, temp_dir_list,
                     sequence_key_mapping=None):
        seen_groups = []
        for path in input_filenames:
            with open(path) as fd:
                seen_groups.append([line.strip() for line in fd])
            os.remove(path)
        for output_file in commands:
            output_file.close()
            os.remove(output_file.name)
        self.seen_groups = seen_groups
        return BindingPredictionCollection([])

    monkeypatch.setattr(
        BaseCommandlinePredictor,
        "_determine_supported_alleles",
        staticmethod(lambda command, flag: {"HLA-A02:01"}))
    monkeypatch.setattr(
        NetMHCstabpan,
        "_run_commands_and_collect_predictions",
        fake_collect)
    monkeypatch.setattr(
        NetMHCstabpan,
        "_check_results",
        lambda self, binding_predictions, peptides, alleles: None)

    predictor = NetMHCstabpan(alleles=[DEFAULT_ALLELE])
    predictor.predict_peptides(["SIINFEKL", "SIINFEKLL", "SIINFEKLQY"])

    assert predictor.group_peptides_by_length is True
    assert sorted(predictor.seen_groups) == sorted([
        ["SIINFEKL"],
        ["SIINFEKLL"],
        ["SIINFEKLQY"],
    ])


_STABPAN_FIXTURE = "\n".join([
    "# NetMHCstabpan version 1.0",
    "-" * 100,
    " pos      HLA         peptide    Identity   Prediction  Thalf(h) %Rank_Stab",
    "-" * 100,
    "    0  HLA-A*02:01   AAAAAAAAAA   PEPLIST      0.075      0.27      19.00",
    "-" * 100,
])


def _zero_half_life_fixture():
    return "\n".join([
        "# NetMHCstabpan version 1.0",
        "-" * 100,
        " pos      HLA         peptide    Identity   Prediction  Thalf(h) %Rank_Stab",
        "-" * 100,
        "    0  HLA-A*02:01   AAAAAAAAAA   PEPLIST      0.000      0.00      99.00",
        "-" * 100,
    ])


def test_stability_value_is_half_life_in_hours():
    # Thalf(h) must reach `value`, the units-bearing field, because
    # Kind.pMHC_stability declares "hours" in VALUE_UNITS.
    from mhctools.parsing import parse_netmhcstabpan
    from mhctools.pred import Kind, value_unit

    predictions = parse_netmhcstabpan(_STABPAN_FIXTURE)
    assert len(predictions) == 1
    pred = predictions[0].to_pred(kind=Kind.pMHC_stability)
    assert pred.kind == Kind.pMHC_stability
    assert pred.value == 0.27
    assert pred.score == 0.27
    assert value_unit(pred.kind) == "hours"


def test_half_life_never_reaches_the_legacy_affinity_field():
    # Regression: `affinity` is documented as an IC50. A duration written
    # there silently changes units for every legacy consumer.
    from mhctools.binding_prediction_collection import BindingPredictionCollection
    from mhctools.parsing import parse_netmhcstabpan

    predictions = parse_netmhcstabpan(_STABPAN_FIXTURE)
    assert predictions[0].affinity is None

    frame = BindingPredictionCollection(predictions).to_dataframe()
    assert frame["affinity"].isna().all()
    # The half-life is still there, in the column that has always carried it.
    assert frame["score"].tolist() == [0.27]


def test_legacy_serialization_roundtrip_keeps_affinity_empty():
    from mhctools.binding_prediction import BindingPrediction
    from mhctools.parsing import parse_netmhcstabpan

    original = parse_netmhcstabpan(_STABPAN_FIXTURE)[0]
    restored = BindingPrediction.from_dict(original.to_dict())
    assert restored.affinity is None
    assert restored.score == 0.27


def test_best_by_value_finds_the_stability_half_life():
    from mhctools.parsing import parse_netmhcstabpan
    from mhctools.pred import Kind, PeptideResult

    pred = parse_netmhcstabpan(_STABPAN_FIXTURE)[0].to_pred(
        kind=Kind.pMHC_stability)
    result = PeptideResult(preds=(pred,))
    best = result.best_by_value(Kind.pMHC_stability)
    assert best is not None, "the registered max-value direction was unreachable"
    assert best.value == 0.27


def test_zero_half_life_is_kept_and_not_reconstructed_as_an_affinity():
    # A Thalf of 0 is a real reading of a very unstable complex, not a missing
    # IC50 to rebuild as 50000 ** (1 - score).
    from mhctools.parsing import parse_netmhcstabpan
    from mhctools.pred import Kind

    prediction = parse_netmhcstabpan(_zero_half_life_fixture())[0]
    assert prediction.affinity is None
    assert prediction.to_pred(kind=Kind.pMHC_stability).value == 0.0
