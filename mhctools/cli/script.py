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

from __future__ import print_function, division, absolute_import
import logging
import logging.config
import pkg_resources
import sys

from pyensembl.fasta import parse_fasta_dictionary

from .args import make_mhc_arg_parser, mhc_binding_predictor_from_args

logging.config.fileConfig(pkg_resources.resource_filename(__name__, 'logging.conf'))
logger = logging.getLogger(__name__)


arg_parser = make_mhc_arg_parser(
    prog="mhctools",
    description=("Predict MHC ligands from protein sequences."))

def add_input_args(arg_parser):
    input_group = arg_parser.add_argument_group("Inputs")
    input_group.add_argument(
        "--sequence",
        nargs="*",
        help=(
            "Protein sequences from which MHC binding predictor will extract "
            "potentially shorter peptides."))
    input_group.add_argument(
        "--input-fasta-file",
        help="Path to FASTA file which contains protein sequences")
    return input_group

def add_output_args(parser):
    output_group = arg_parser.add_argument_group("Outputs")
    output_group.add_argument("--output-csv", default=None)
    return output_group

add_input_args(arg_parser)
add_output_args(arg_parser)

def main(args_list=None):
    """
    Script to make pMHC binding predictions from amino acid sequences.

    Usage example:
        mhctools
            --sequence SFFPIQQQQQAAALLLI \
            --sequence SILQQQAQAQQAQAASSSC \
            --mhc-predictor netmhc \
            --mhc-alleles HLA-A0201 H2-Db \
            --output-csv epitope.csv
    """
    if args_list is None:
        args_list = sys.argv[1:]
    args = arg_parser.parse_args(args_list)
    predictor = mhc_binding_predictor_from_args(args)
    if args.sequence:
        input_dictionary = {
            ("input-sequence-%d" % i): seq
            for (i, seq)
            in enumerate(args.sequence)
        }
    else:
        input_dictionary = {}

    if args.input_fasta_file:
        input_dictionary.update(parse_fasta_dictionary(args.input_fasta_file))

    if len(input_dictionary) == 0:
        raise ValueError(
            ("No input sequences provided, "
             "use either --sequence or --input-fasta-file"))
    epitope_collection = predictor.predict(input_dictionary)
    df = epitope_collection.to_dataframe()
    logger.info('\n%s', df)
    if args.output_csv:
        df.to_csv(args.output_csv)
