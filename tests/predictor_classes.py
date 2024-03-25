import mhctools
from .arch import apple_silicon

mhc1_predictor_classes = [
    mhctools.NetMHCpan,
    mhctools.NetMHC,
    mhctools.NetMHC4
]

# for now excluding because of 403 errors
mhc1_iedb_predictors = [
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
]

# for now excluding because of 403 errors
mhc2_iedb_predictors = [
    mhctools.IedbNetMHCIIpan,
]
