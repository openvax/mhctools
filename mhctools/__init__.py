from .binding_prediction import BindingPrediction
from .binding_prediction_collection import BindingPredictionCollection
from .pred import Prediction, Pred, PeptideResult, PeptidePreds, Kind, preds_from_rows
from .sample import MultiSample
from .iedb import (
    IedbNetMHCcons,
    IedbNetMHCpan,
    IedbSMM,
    IedbSMM_PMBEC,
    IedbNetMHCIIpan,
)
from .mixmhcpred import MixMHCpred
from .processing_predictor import (
    ProcessingPredictor,
    SCORING_MODES,
    resolve_scoring,
    score_cterm,
    score_nterm_cterm,
    score_cterm_anti_max_internal,
    score_cterm_anti_mean_internal,
    score_nterm_cterm_anti_max_internal,
    score_nterm_cterm_anti_mean_internal,
)
from .proteasome_predictor import ProteasomePredictor
from .netchop import NetChop
from .pepsickle import Pepsickle
from .netmhc import NetMHC
from .netmhc3 import NetMHC3
from .netmhc4 import NetMHC4
from .netmhc_cons import NetMHCcons
from .netmhc_pan import NetMHCpan
from .netmhc_pan28 import NetMHCpan28
from .netmhc_pan3 import NetMHCpan3
from .netmhc_pan4 import NetMHCpan4, NetMHCpan4_BA, NetMHCpan4_EL
from .netmhc_pan41 import NetMHCpan41, NetMHCpan41_BA, NetMHCpan41_EL
from .netmhc_pan42 import NetMHCpan42, NetMHCpan42_BA, NetMHCpan42_EL
from .netmhcii_pan import NetMHCIIpan, NetMHCIIpan3, NetMHCIIpan4, NetMHCIIpan4_BA, NetMHCIIpan4_EL, NetMHCIIpan43, NetMHCIIpan43_BA, NetMHCIIpan43_EL
from .random_predictor import RandomBindingPredictor
from .netmhcstabpan import NetMHCstabpan
from .unsupported_allele import UnsupportedAllele


# Lazy-load predictors that pull in heavy optional dependencies (torch, Keras/TF).
# Accessing mhctools.BigMHC or mhctools.MHCflurry triggers the import on first use.
_LAZY_IMPORTS = {
    "BigMHC": (".bigmhc", "BigMHC"),
    "BigMHC_EL": (".bigmhc", "BigMHC_EL"),
    "BigMHC_IM": (".bigmhc", "BigMHC_IM"),
    "MHCflurry": (".mhcflurry", "MHCflurry"),
    "MHCflurry_Affinity": (".mhcflurry", "MHCflurry_Affinity"),
}


def __getattr__(name):
    """PEP 562: lazy-load heavy predictor modules on first access."""
    if name in _LAZY_IMPORTS:
        from importlib import import_module
        module_name, attr = _LAZY_IMPORTS[name]
        module = import_module(module_name, package=__name__)
        value = getattr(module, attr)
        globals()[name] = value
        return value
    raise AttributeError(
        "module %r has no attribute %r" % (__name__, name))

__version__ = "3.13.0"

__all__ = [
    "Prediction",
    "Pred",  # backward compat alias
    "PeptideResult",
    "PeptidePreds",  # backward compat alias
    "Kind",
    "preds_from_rows",
    "MultiSample",
    "BindingPrediction",
    "BindingPredictionCollection",
    "IedbNetMHCcons",
    "IedbNetMHCpan",
    "IedbSMM",
    "IedbSMM_PMBEC",
    "IedbNetMHCIIpan",
    "MixMHCpred",
    "MHCflurry",
    "MHCflurry_Affinity",
    "ProcessingPredictor",
    "ProteasomePredictor",
    "SCORING_MODES",
    "resolve_scoring",
    "score_cterm",
    "score_nterm_cterm",
    "score_cterm_anti_max_internal",
    "score_cterm_anti_mean_internal",
    "score_nterm_cterm_anti_max_internal",
    "score_nterm_cterm_anti_mean_internal",
    "NetChop",
    "Pepsickle",
    "NetMHC",
    "NetMHC3",
    "NetMHC4",
    "NetMHCcons",
    "NetMHCpan",
    "NetMHCpan28",
    "NetMHCpan3",
    "NetMHCpan4",
    "NetMHCpan41",
    "NetMHCpan4_BA",
    "NetMHCpan4_EL",
    "NetMHCpan41_BA",
    "NetMHCpan41_EL",
    "NetMHCpan42",
    "NetMHCpan42_BA",
    "NetMHCpan42_EL",
    "NetMHCIIpan",
    "NetMHCIIpan3",
    "NetMHCIIpan4",
    "NetMHCIIpan4_BA",
    "NetMHCIIpan4_EL",
    "NetMHCIIpan43",
    "NetMHCIIpan43_BA",
    "NetMHCIIpan43_EL",
    "NetMHCstabpan",
    "BigMHC",
    "BigMHC_EL",
    "BigMHC_IM",
    "RandomBindingPredictor",
    "UnsupportedAllele",
]
