
# Copyright (c) 2016. Mount Sinai School of Medicine
#
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

"""
Commandline options for MHC Binding Prediction
"""

from argparse import ArgumentParser
import inspect
import logging

from ..allele_normalization import normalize_allele_name

from .parsing_helpers import parse_int_list
from .. import (
    NetMHC,
    NetMHC3,
    NetMHC4,
    NetMHCpan,
    NetMHCpan28,
    NetMHCpan3,
    NetMHCpan4,
    NetMHCpan4_EL,
    NetMHCpan4_BA,
    NetMHCpan41,
    NetMHCpan41_EL,
    NetMHCpan41_BA,
    NetMHCpan42,
    NetMHCpan42_EL,
    NetMHCpan42_BA,
    NetMHCIIpan,
    NetMHCIIpan3,
    NetMHCIIpan4,
    NetMHCIIpan4_EL,
    NetMHCIIpan4_BA,
    NetMHCIIpan43,
    NetMHCIIpan43_EL,
    NetMHCIIpan43_BA,
    NetMHCcons,
    NetMHCstabpan,
    NetChop,
    Pepsickle,
    RandomBindingPredictor,
    IedbNetMHCpan,
    IedbNetMHCcons,
    IedbSMM,
    IedbSMM_PMBEC,
    IedbNetMHCIIpan,
    MixMHCpred,
)


class _LazyPredictor:
    """Lazy-import wrapper for heavy predictors (BigMHC, MHCflurry).

    Behaves like the class itself: calling it instantiates, subclass
    checks work, and __name__ is preserved for introspection.
    """
    def __init__(self, module_name, class_name):
        self._module_name = module_name
        self._class_name = class_name
        self._cls = None
        self.__name__ = class_name

    def _resolve(self):
        if self._cls is None:
            from importlib import import_module
            module = import_module(self._module_name, package="mhctools")
            self._cls = getattr(module, self._class_name)
        return self._cls

    def __call__(self, *args, **kwargs):
        return self._resolve()(*args, **kwargs)

    def __getattr__(self, item):
        return getattr(self._resolve(), item)

    def __eq__(self, other):
        if isinstance(other, _LazyPredictor):
            return (self._module_name, self._class_name) == (
                other._module_name, other._class_name)
        return self._resolve() is other

    def __hash__(self):
        return hash((self._module_name, self._class_name))


BigMHC = _LazyPredictor(".bigmhc", "BigMHC")
BigMHC_EL = _LazyPredictor(".bigmhc", "BigMHC_EL")
BigMHC_IM = _LazyPredictor(".bigmhc", "BigMHC_IM")
MHCflurry = _LazyPredictor(".mhcflurry", "MHCflurry")
MHCflurry_Affinity = _LazyPredictor(".mhcflurry", "MHCflurry_Affinity")


logger = logging.getLogger(__name__)

mhc_predictors = {
    "netmhc": NetMHC,
    "netmhc3": NetMHC3,
    "netmhc4": NetMHC4,
    "netmhcpan": NetMHCpan,
    "netmhcpan4": NetMHCpan4,
    "netmhcpan4-ba": NetMHCpan4_BA,
    "netmhcpan4-el": NetMHCpan4_EL,
    "netmhcpan41": NetMHCpan41,
    "netmhcpan41-ba": NetMHCpan41_BA,
    "netmhcpan41-el": NetMHCpan41_EL,
    "netmhcpan42": NetMHCpan42,
    "netmhcpan42-ba": NetMHCpan42_BA,
    "netmhcpan42-el": NetMHCpan42_EL,
    "netmhcpan28": NetMHCpan28,
    "netmhcpan3": NetMHCpan3,
    "netmhciipan": NetMHCIIpan,
    "netmhciipan3": NetMHCIIpan3,
    "netmhciipan4": NetMHCIIpan4,
    "netmhciipan4-ba": NetMHCIIpan4_BA,
    "netmhciipan4-el": NetMHCIIpan4_EL,
    "netmhciipan43": NetMHCIIpan43,
    "netmhciipan43-ba": NetMHCIIpan43_BA,
    "netmhciipan43-el": NetMHCIIpan43_EL,
    "netmhccons": NetMHCcons,
    "netmhcstabpan": NetMHCstabpan,
    "bigmhc": BigMHC,
    "bigmhc-el": BigMHC_EL,
    "bigmhc-im": BigMHC_IM,
    "netchop": NetChop,
    "pepsickle": Pepsickle,
    "random": RandomBindingPredictor,
    # use NetMHCpan via IEDB's web API
    "netmhcpan-iedb": IedbNetMHCpan,
    # use NetMHCcons via IEDB's web API
    "netmhccons-iedb": IedbNetMHCcons,
    # use SMM via IEDB's web API
    "smm-iedb": IedbSMM,
    # use SMM-PMBEC via IEDB's web API
    "smm-pmbec-iedb": IedbSMM_PMBEC,
    # Class II MHC binding prediction using NetMHCIIpan via IEDB
    "netmhciipan-iedb": IedbNetMHCIIpan,
    "mhcflurry": MHCflurry,
    "mhcflurry-affinity": MHCflurry_Affinity,
    "mixmhcpred": MixMHCpred,
}


def _cls_accepts(cls, param_name):
    """Check whether a predictor class/factory accepts a given __init__ parameter."""
    # Unwrap lazy predictors so inspect sees the real class.
    if isinstance(cls, _LazyPredictor):
        cls = cls._resolve()
    sig = inspect.signature(cls)
    if any(p.kind == inspect.Parameter.VAR_KEYWORD for p in sig.parameters.values()):
        return True
    return param_name in sig.parameters


def _parse_predictor_names(value):
    """Normalize a single CLI token: lowercase, split on commas, validate."""
    names = [s.strip().lower() for s in value.split(",") if s.strip()]
    for name in names:
        if name not in mhc_predictors:
            from argparse import ArgumentTypeError
            raise ArgumentTypeError(
                "Unknown predictor %r. Available: %s"
                % (name, ", ".join(sorted(mhc_predictors.keys()))))
    return names


def add_mhc_args(arg_parser):
    mhc_options_arg_group = arg_parser.add_argument_group(
        title="Prediction Options",
        description="Model selection and prediction options")

    mhc_options_arg_group.add_argument(
        "--mhc-predictor",
        nargs="+",
        type=_parse_predictor_names,
        required=True,
        help=(
            "One or more predictor names, space or comma separated. "
            "E.g. '--mhc-predictor netmhcpan42' or "
            "'--mhc-predictor netmhcpan42 bigmhc-el pepsickle' or "
            "'--mhc-predictor netmhcpan42,bigmhc-el,pepsickle'. "
            "Available: %s" % ", ".join(sorted(mhc_predictors.keys()))
        ))

    mhc_options_arg_group.add_argument(
        "--mhc-predictor-path",
        help="Path to executable program or model directory",
        default=None,
        required=False)

    mhc_options_arg_group.add_argument(
        "--mhc-predictor-models-path",
        help="Path to models to use with predictor (e.g. mhcflurry).")

    mhc_options_arg_group.add_argument(
        "--mhc-peptide-lengths",
        type=parse_int_list,
        help="Lengths of epitopes to consider for binding prediction")

    mhc_options_arg_group.add_argument(
        "--mhc-epitope-lengths",
        type=parse_int_list,
        help="Deprecated name for --mhc-peptide-lengths")

    mhc_options_arg_group.add_argument(
        "--mhc-alleles-file",
        help="File with one HLA allele per line")

    mhc_options_arg_group.add_argument(
        "--mhc-alleles",
        default="",
        help="Comma or space separated list of allele (default HLA-A*02:01)")

    mhc_options_arg_group.add_argument(
        "--do-not-raise-on-error",
        action="store_true", default=False,
        help="Only applies to IEDB predictors: if this arg is present, will not crash on any "
        "errors, which can result from connection issues or supplying unsupported MHC alleles. "
        "In such cases, some predictions may get dropped from the returned result set.")

    return mhc_options_arg_group


def make_mhc_arg_parser(**kwargs):
    parser = ArgumentParser(**kwargs)
    add_mhc_args(parser)
    return parser

def mhc_alleles_from_args(args):
    alleles = [
        normalize_allele_name(allele.strip())
        for comma_group in args.mhc_alleles.split(",")
        for allele in comma_group.split(" ")
        if allele.strip()
    ]
    if args.mhc_alleles_file:
        with open(args.mhc_alleles_file, 'r') as f:
            for line in f:
                line = line.strip()
                if line:
                    alleles.append(normalize_allele_name(line))
    if len(alleles) == 0:
        raise ValueError(
            "MHC alleles required (use --mhc-alleles or --mhc-alleles-file)")
    return alleles

def _flatten_predictor_names(args):
    """Flatten the nested list from nargs="+" + comma-splitting type function."""
    raw = args.mhc_predictor  # list of lists
    return [name for group in raw for name in group]

def _build_predictor(cls, name, alleles, peptide_lengths, args):
    """Build a single predictor instance, passing only kwargs it accepts."""
    kwargs = {}
    if alleles is not None and _cls_accepts(cls, "alleles"):
        kwargs["alleles"] = alleles
    if peptide_lengths is not None and _cls_accepts(cls, "default_peptide_lengths"):
        kwargs["default_peptide_lengths"] = peptide_lengths
    if getattr(args, "mhc_predictor_models_path", None) and _cls_accepts(cls, "models_path"):
        kwargs["models_path"] = args.mhc_predictor_models_path
    if getattr(args, "mhc_predictor_path", None):
        if _cls_accepts(cls, "bigmhc_path"):
            kwargs["bigmhc_path"] = args.mhc_predictor_path
        elif _cls_accepts(cls, "program_name"):
            kwargs["program_name"] = args.mhc_predictor_path
    if getattr(args, "do_not_raise_on_error", False):
        if _cls_accepts(cls, "raise_on_error"):
            kwargs["raise_on_error"] = False
        else:
            logger.warning(
                "--do-not-raise-on-error ignored for predictor %s", name)
    logger.info("Building predictor %s(%s)",
                getattr(cls, "__name__", name), kwargs)
    return cls(**kwargs)


def predictors_from_args(args):
    """Build predictor instances from --mhc-predictor.

    Returns a list of instantiated predictor objects.
    """
    names = _flatten_predictor_names(args)

    # Fetch alleles only if at least one predictor needs them
    needs_alleles = any(_cls_accepts(mhc_predictors[n], "alleles") for n in names)
    alleles = None
    if needs_alleles:
        alleles = mhc_alleles_from_args(args)

    peptide_lengths = getattr(args, "mhc_peptide_lengths", None)
    if not peptide_lengths:
        peptide_lengths = getattr(args, "mhc_epitope_lengths", None)

    return [
        _build_predictor(mhc_predictors[name], name, alleles, peptide_lengths, args)
        for name in names
    ]


def mhc_binding_predictor_from_args(args):
    """Build a single predictor from --mhc-predictor (legacy entry point)."""
    names = _flatten_predictor_names(args)
    if len(names) != 1:
        raise ValueError(
            "mhc_binding_predictor_from_args expects exactly one predictor, "
            "got %d. Use predictors_from_args for multiple." % len(names))
    return predictors_from_args(args)[0]
