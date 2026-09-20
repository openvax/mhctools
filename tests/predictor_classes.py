import mhctools
import pytest

from mhctools.optional_backend import probe_executable

mhc1_predictor_classes = [
    mhctools.NetMHCpan,
    mhctools.NetMHC,
    mhctools.NetMHC4
]

# Availability, not host architecture, determines whether these integrations
# can run: the test runtime supports the legacy Linux tools on Apple Silicon.
netmhc3 = probe_executable("netMHC-3.4", args=("-h",))
netmhccons = probe_executable("netMHCcons", args=("-h",))
for predictor, capabilities in (
        (mhctools.NetMHC3, (netmhc3,)),
        (mhctools.NetMHCcons, (netmhc3, netmhccons))):
    reasons = "; ".join(c.reason for c in capabilities if not c.runnable)
    mhc1_predictor_classes.append(pytest.param(
        predictor, marks=pytest.mark.skipif(bool(reasons), reason=reasons)))

mhc2_predictor_classes = [
    mhctools.NetMHCIIpan,
]
