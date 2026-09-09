"""Offline HTTP timeout and error-policy regressions for IEDB predictors."""

import io
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer
from threading import Event, Thread

import pytest

from mhctools import iedb


RESPONSE = (
    b"allele seq_num start end length peptide ic50 percentile_rank\n"
    b"HLA-A*02:01 1 1 9 9 SIINFEKLA 123.0 1.5\n")


@pytest.mark.parametrize("predictor_class", [
    iedb.IedbNetMHCcons, iedb.IedbNetMHCpan, iedb.IedbSMM,
    iedb.IedbSMM_PMBEC, iedb.IedbNetMHCIIpan,
])
def test_public_predictors_forward_timeout_and_close_response(
        predictor_class, monkeypatch):
    handles = []
    timeouts = []

    def open_response(request, timeout):
        assert request.get_method() == "POST"
        timeouts.append(timeout)
        handle = io.BytesIO(RESPONSE)
        handles.append(handle)
        return handle

    monkeypatch.setattr(iedb, "urlopen", open_response)
    for configured in (None, 2.5):
        options = {} if configured is None else {"request_timeout": configured}
        predictor = predictor_class(alleles=["HLA-A*02:01"], **options)
        predictions = predictor.predict_subsequences(
            {"input": "SIINFEKLA"}, peptide_lengths=[9])
        assert predictions[0].affinity == 123.0
        assert handles[-1].closed
    assert timeouts == [iedb.DEFAULT_REQUEST_TIMEOUT, 2.5]


@pytest.mark.parametrize("timeout", [0, -1, float("inf"), float("nan"), None, True])
def test_invalid_timeout_is_rejected(timeout):
    with pytest.raises(ValueError, match="finite positive"):
        iedb.IedbNetMHCpan(alleles=["HLA-A*02:01"], request_timeout=timeout)


def test_read_timeout_closes_response_and_respects_error_policy(monkeypatch):
    class StalledResponse(io.BytesIO):
        def read(self, *args):
            raise TimeoutError("IEDB response stalled")

    handles = []

    def open_response(request, timeout):
        handle = StalledResponse()
        handles.append(handle)
        return handle

    monkeypatch.setattr(iedb, "urlopen", open_response)
    predictor = iedb.IedbNetMHCpan(alleles=["HLA-A*02:01"])
    with pytest.raises(TimeoutError, match="stalled"):
        predictor.predict_subsequences({"input": "SIINFEKLA"}, [9])
    assert handles[-1].closed
    predictor.raise_on_error = False
    assert len(predictor.predict_subsequences({"input": "SIINFEKLA"}, [9])) == 0
    assert handles[-1].closed


def test_http_body_stall_obeys_socket_timeout():
    """Exercise urllib's actual transport against a local stalled response."""
    release = Event()

    class Handler(BaseHTTPRequestHandler):
        def do_POST(self):
            self.send_response(200)
            self.send_header("Content-Length", "100")
            self.end_headers()
            self.wfile.flush()
            # Bound the fixture even if the client timeout regresses.
            release.wait(1)

        def log_message(self, *args):
            pass

    server = ThreadingHTTPServer(("127.0.0.1", 0), Handler)
    thread = Thread(target=lambda: server.serve_forever(poll_interval=0.01))
    thread.start()
    try:
        url = "http://127.0.0.1:%d/" % server.server_port
        with pytest.raises(TimeoutError):
            iedb._query_iedb({"method": "fixture"}, url, timeout=0.05)
    finally:
        release.set()
        server.shutdown()
        server.server_close()
        thread.join()
