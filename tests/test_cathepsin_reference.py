"""Source-backed, chemistry-specific CatL evidence through the batch workflow."""

from copy import deepcopy
import json
from pathlib import Path

from mhctools import load_cleavage_batch, predict_cleavage_batch, write_cleavage_batch


FIXTURE = Path(__file__).parent / "data/tusar2023/catl-protected-peptides.json"


def test_catl_products_conditions_and_chemical_form_survive_batch_reload(tmp_path):
    request = json.loads(FIXTURE.read_text())
    inputs = deepcopy(request["inputs"])
    inputs.append(dict(inputs[0], id="unprotected", n_term="free", c_term="free"))
    inputs.append(dict(id="unseen", scope="construct", sequence="ADLLHPSP",
                       n_term="acetylated", c_term="amidated"))
    report = predict_cleavage_batch(inputs, request["scenarios"],
                                    reference_panels=request["reference_panels"], raise_on_error=True)
    p3, p4, unprotected, unseen = report["assessments"]
    assert [r["status"] for r in (p3, p4, unprotected, unseen)] == [
        "assessed", "assessed", "unsupported", "unsupported"]
    assert [s["bond"] for s in p3["result"]["sites"]] == [3, 4]
    assert [s["bond"] for s in p4["result"]["sites"]] == [4, 5, 8]
    for row in (p3, p4):
        assert all(s["status"] == "reported" and s["score"] is None for s in row["result"]["sites"])
        conditions = dict(row["result"]["conditions"])
        assert conditions["pH"] == "5.5" and conditions["incubation_hours"] == "2"
        assert conditions["source_measurement_id"].startswith("Supplementary Data 4 page 1")
        value = next(v for v in request["inputs"] if v["id"] == row["input_id"])
        assert all(value["sequence"][p["start"]:p["end"]] == p["sequence"]
                   for p in value["observed_products"])
    assert p3["overlays"][0]["c_boundary"]["status"] == "sequence_endpoint"
    assert p4["overlays"][0]["internal"][0]["status"] == "unassessed"
    path, html = tmp_path / "catl.json", tmp_path / "catl.html"
    write_cleavage_batch(report, path, html_path=html)
    assert load_cleavage_batch(path) == report
    assert report["reference_panels"] == request["reference_panels"]
    assert "MALDI-TOF" in html.read_text() and "acetylated" in html.read_text()
