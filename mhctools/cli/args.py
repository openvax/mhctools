
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

from __future__ import print_function, division, absolute_import
from argparse import ArgumentParser
import logging

from mhcnames import normalize_allele_name

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
    NetMHCIIpan,
    NetMHCIIpan3,
    NetMHCIIpan4,
    NetMHCIIpan4_EL,
    NetMHCIIpan4_BA,
    NetMHCcons,
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
    "netmhcpan28": NetMHCpan28,
    "netmhcpan3": NetMHCpan3,
    "netmhciipan": NetMHCIIpan,
    "netmhciipan3": NetMHCIIpan3,
    "netmhciipan4": NetMHCIIpan4,
    "netmhciipan4-ba": NetMHCIIpan4_BA,
    "netmhciipan4-el": NetMHCIIpan4_EL,
    "netmhccons": NetMHCcons,
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
    # TODO implement SMM predictors in mhctools
    # "smm": None,
    # "smm-pmbec": None,
    "mhcflurry": MHCflurry,
    "mixmhcpred": MixMHCpred,
}

def add_mhc_args(arg_parser):
    mhc_options_arg_group = arg_parser.add_argument_group(
        title="MHC Prediction Options",
        description="MHC Binding Prediction Options")

    mhc_options_arg_group.add_argument(
        "--mhc-predictor",
        choices=list(sorted(mhc_predictors.keys())),
        type=lambda s: s.lower().strip(),
        required=True)

    mhc_options_arg_group.add_argument(
        "--mhc-predictor-path",
        help="Path to executable program used for MHC binding predictor",
        default=None,
        required=False)

    mhc_options_arg_group.add_argument(
        "--mhc-predictor-models-path",
        help="Path to models to use with predictor. Currently supported only "
        "by mhcflurry predictor.")

    mhc_options_arg_group.add_argument(
        "--mhc-peptide-lengths",
        type=parse_int_list,
        help="Lengths of epitopes to consider for MHC binding prediction")

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

def mhc_binding_predictor_from_args(args):
    mhc_class = mhc_predictors.get(args.mhc_predictor)
    if mhc_class is None:
        raise ValueError(
            "Invalid MHC prediction method: %s" % (args.mhc_predictor,))
    alleles = mhc_alleles_from_args(args)
    peptide_lengths = args.mhc_peptide_lengths
    if not peptide_lengths:
        peptide_lengths = args.mhc_epitope_lengths
    logger.info(
        ("Building MHC binding prediction %s"
         " for alleles %s"
         " and epitope lengths %s") % (
            mhc_class.__class__.__name__,
            alleles,
            peptide_lengths))
    kwargs = dict(alleles=alleles)
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
