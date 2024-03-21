import mhctools
from .arch import apple_silicon

mhc1_predictor_classes = [
    mhctools.NetMHCpan,
    mhctools.NetMHC,
    mhctools.NetMHC4,
    mhctools.IedbNetMHCcons,
    mhctools.IedbNetMHCpan,
    mhctools.IedbSMM,
    mhctools.IedbSMM_PMBEC,
]

if not apple_silicon: 
    mhc1_predictor_classes += [
        mhctools.NetMHC3,
        mhctools.NetMHCcons,
    ]

mhc2_predictor_classes = [
    mhctools.NetMHCIIpan,
    mhctools.IedbNetMHCIIpan,
]