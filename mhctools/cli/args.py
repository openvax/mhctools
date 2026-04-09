
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
    BigMHC,
    NetChop,
    Pepsickle,
    RandomBindingPredictor,
    IedbNetMHCpan,
    IedbNetMHCcons,
    IedbSMM,
    IedbSMM_PMBEC,
    IedbNetMHCIIpan,
    MHCflurry,
    MixMHCpred,
)


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
    "mixmhcpred": MixMHCpred,
}

# Predictors that don't take alleles (processing/cleavage predictors)
_no_allele_predictors = {"netchop", "pepsickle"}

def _parse_predictor_list(value):
    """Parse a comma-separated list of predictor names."""
    names = [s.strip().lower() for s in value.split(",")]
    for name in names:
        if name not in mhc_predictors:
            from argparse import ArgumentTypeError
            available = sorted(mhc_predictors.keys())
            raise ArgumentTypeError(
                "Unknown predictor %r. Available: %s" % (name, available))
    return names


def add_mhc_args(arg_parser):
    mhc_options_arg_group = arg_parser.add_argument_group(
        title="Prediction Options",
        description="Model selection and prediction options")

    mhc_options_arg_group.add_argument(
        "--predictors",
        type=_parse_predictor_list,
        default=None,
        help=(
            "Comma-separated list of predictor names to run. "
            "Can specify one or more: e.g. 'netmhcpan42' or "
            "'netmhcpan42,bigmhc,pepsickle'. "
            "Available: %s" % ", ".join(sorted(mhc_predictors.keys()))
        ))

    mhc_options_arg_group.add_argument(
        "--mhc-predictor",
        choices=list(sorted(mhc_predictors.keys())),
        type=lambda s: s.lower().strip(),
        default=None,
        help="Deprecated: use --predictors instead. "
             "Single prediction method to use.")

    mhc_options_arg_group.add_argument(
        "--mhc-predictor-path",
        help="Path to executable program used for command-line predictors",
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

def predictors_from_args(args):
    """Build a list of predictor instances from --predictors (or --mhc-predictor).

    Returns a list of instantiated predictor objects.
    """
    predictor_names = getattr(args, "predictors", None)
    mhc_predictor = getattr(args, "mhc_predictor", None)

    if predictor_names and mhc_predictor:
        raise ValueError(
            "--predictors and --mhc-predictor are mutually exclusive. "
            "Use --predictors (which supports multiple models)."
        )

    if not predictor_names and not mhc_predictor:
        raise ValueError(
            "Either --predictors or --mhc-predictor is required."
        )

    if mhc_predictor:
        logger.warning(
            "--mhc-predictor is deprecated; use --predictors instead"
        )
        predictor_names = [mhc_predictor]

    # Determine if any model needs alleles
    needs_alleles = any(n not in _no_allele_predictors for n in predictor_names)
    alleles = None
    if needs_alleles:
        alleles = mhc_alleles_from_args(args)
    elif args.mhc_alleles or getattr(args, "mhc_alleles_file", None):
        logger.info("Alleles ignored — no selected predictor requires them")

    peptide_lengths = getattr(args, "mhc_peptide_lengths", None)
    if not peptide_lengths:
        peptide_lengths = getattr(args, "mhc_epitope_lengths", None)

    models = []
    for name in predictor_names:
        cls = mhc_predictors[name]
        kwargs = {}
        if name not in _no_allele_predictors:
            kwargs["alleles"] = alleles
        if peptide_lengths is not None:
            kwargs["default_peptide_lengths"] = peptide_lengths
        if getattr(args, "mhc_predictor_models_path", None):
            kwargs["models_path"] = args.mhc_predictor_models_path
        if getattr(args, "mhc_predictor_path", None):
            kwargs["program_name"] = args.mhc_predictor_path
        if getattr(args, "do_not_raise_on_error", False):
            if "iedb" in name:
                kwargs["raise_on_error"] = False
            else:
                logger.warning(
                    "--do-not-raise-on-error ignored for non-IEDB predictor %s",
                    name,
                )
        logger.info("Building predictor %s(%s)", cls.__name__, kwargs)
        models.append(cls(**kwargs))
    return models


def mhc_binding_predictor_from_args(args):
    mhc_class = mhc_predictors.get(args.mhc_predictor)
    if mhc_class is None:
        raise ValueError(
            "Invalid MHC prediction method: %s" % (args.mhc_predictor,))

    needs_alleles = args.mhc_predictor not in _no_allele_predictors

    if needs_alleles:
        alleles = mhc_alleles_from_args(args)
    else:
        alleles = None
        # Still allow alleles to be empty for processing predictors
        if args.mhc_alleles or args.mhc_alleles_file:
            logger.info(
                "Alleles ignored for processing predictor %s" %
                args.mhc_predictor)

    peptide_lengths = args.mhc_peptide_lengths
    if not peptide_lengths:
        peptide_lengths = args.mhc_epitope_lengths
    logger.info(
        ("Building MHC binding prediction %s"
         " for alleles %s"
         " and epitope lengths %s") % (
            mhc_class.__name__,
            alleles,
            peptide_lengths))

    kwargs = {}
    if needs_alleles:
        kwargs["alleles"] = alleles
    if peptide_lengths is not None:
        kwargs["default_peptide_lengths"] = peptide_lengths

    if args.mhc_predictor_models_path:
        kwargs["models_path"] = args.mhc_predictor_models_path

    if args.mhc_predictor_path:
        kwargs["program_name"] = args.mhc_predictor_path

    if args.do_not_raise_on_error:
        if 'iedb' in args.mhc_predictor.lower():
            kwargs["raise_on_error"] = False
        else:
            logger.warning('--do-not-raise-on-error ignored for non-IEDB predictor')

    return mhc_class(**kwargs)
