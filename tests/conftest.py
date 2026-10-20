"""Headless test setup and the optional fully provisioned release gate."""

import pytest


def pytest_addoption(parser):
    parser.addoption(
        "--require-all", action="store_true",
        help="Fail the run if any test is skipped or marked as an expected failure.")


def pytest_configure(config):
    import matplotlib

    # Static plot tests must not open a GUI on macOS framework Python (#607).
    matplotlib.use("Agg")
    if config.getoption("--require-all"):
        config.pluginmanager.register(_RequireAllTests(), "mhctools-require-all")


@pytest.fixture(scope="session")
def netmhc3_runtime():
    from tests.common import require_executable
    require_executable("netMHC-3.4", args=("-h",))


@pytest.fixture(scope="session")
def netmhccons_runtime(netmhc3_runtime):
    from tests.common import require_executable
    require_executable("netMHCcons", args=("-h",))


class _RequireAllTests:
    def __init__(self):
        self.incomplete = set()

    def pytest_collectreport(self, report):
        if report.skipped:
            self.incomplete.add(report.nodeid)

    def pytest_runtest_logreport(self, report):
        if report.skipped or hasattr(report, "wasxfail"):
            self.incomplete.add(report.nodeid)

    def pytest_sessionfinish(self, session, exitstatus):
        if self.incomplete and exitstatus == pytest.ExitCode.OK:
            session.exitstatus = pytest.ExitCode.TESTS_FAILED

    def pytest_terminal_summary(self, terminalreporter):
        if self.incomplete:
            terminalreporter.section("--require-all: incomplete test coverage", red=True)
            for nodeid in sorted(self.incomplete):
                terminalreporter.write_line(nodeid)
