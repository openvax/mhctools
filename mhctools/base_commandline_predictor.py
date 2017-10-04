# Copyright (c) 2014-2017. Mount Sinai School of Medicine
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#       http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

from __future__ import print_function, division, absolute_import
from collections import defaultdict
import logging
from subprocess import check_output
import tempfile

from six import string_types
from typechecks import require_string, require_integer, require_iterable_of
from mhcnames import normalize_allele_name, AlleleParseError

from .base_predictor import BasePredictor
from .unsupported_allele import UnsupportedAllele
from .process_helpers import run_command
from .cleanup_context import CleanupFiles
from .input_file_formats import create_input_peptides_files
from .process_helpers import run_multiple_commands_redirect_stdout
from .binding_prediction_collection import BindingPredictionCollection

logger = logging.getLogger(__name__)

class BaseCommandlinePredictor(BasePredictor):
    """
    Base class for MHC binding predictors that run a local external
    program and write their output to a local file.
    """
    def __init__(
            self,
            program_name,
            alleles,
            parse_output_fn,
            supported_alleles_flag,
            input_file_flag,
            length_flag,
            allele_flag,
            peptide_mode_flags=["-p"],
            tempdir_flag=None,
            extra_flags=[],
            max_peptides_per_file=10 ** 4,
            process_limit=-1,
            default_peptide_lengths=[9],
            group_peptides_by_length=False,
            min_peptide_length=8,
            max_peptide_length=None,):
        """
        Parameters
        ----------
        program_name : str
            Name of prediction program to run
            (e.g. "netMHCcons" or "netMHCIIpan")

        alleles : list of str
            MHC alleles

        supported_alleles_flag : str
            Flag to pass to the predictor to get a list of supported alleles
            (e.g. "-A", "-list", "-listMHC")

        parse_output_fn : fn
            Takes the stdout string from the predictor and returns a collection
            of BindingPrediction objects

        input_file_flag : str
            How to specify the input FASTA file of source sequences (e.g. "-f")

        length_flag : str
            How to specify the desired predicted peptide length (e.g. "-length")

        allele_flag : str
            How to specify the allele we want predictions for (e.g. "-a")

        peptide_mode_flags : list of str
            How to switch from the default FASTA subsequences input mode to
            where peptides are explicitly given one per line of a text file.

        tempdir_flag : str, optional
            How to specify the predictor's temporary directory (e.g. "-tdir")

        extra_flags : list of str
            Extra flags to pass to the predictor

        max_peptides_per_file : int, optional
            Maximum number of lines per file when predicting peptides directly.

        process_limit : int, optional
            Maximum number of parallel processes to start
            (0 for no limit, -1 for use all available processors)

        default_peptide_lengths : list of int, optional
            When making predictions across subsequences of protein sequences,
            what peptide lengths to predict for.

        group_peptides_by_length : bool
            Run commandline predictor on groups of peptides of equal length

        min_peptide_length : int
            Shortest peptide this predictor can handle

        max_peptide_length : int
            Longest peptide this predictor can handle
        """
        require_string(program_name, "Predictor program name")
        self.program_name = program_name

        if supported_alleles_flag is not None:
            require_string(supported_alleles_flag, "Supported alleles flag")
        self.supported_alleles_flag = supported_alleles_flag

        require_string(input_file_flag, "Input file flag")
        self.input_file_flag = input_file_flag

        require_string(length_flag, "Peptide length flag")
        self.length_flag = length_flag

        require_string(allele_flag, "Allele flag")
        self.allele_flag = allele_flag

        require_iterable_of(peptide_mode_flags, string_types)
        self.peptide_mode_flags = peptide_mode_flags

        if tempdir_flag is not None:
            require_string(tempdir_flag, "Temporary directory flag")
        self.tempdir_flag = tempdir_flag

        require_iterable_of(extra_flags, string_types)
        self.extra_flags = extra_flags

        require_integer(
            max_peptides_per_file,
            "Maximum number of lines in a peptides input file")
        self.max_peptides_per_file = max_peptides_per_file

        require_integer(process_limit, "Maximum number of processes")
        self.process_limit = process_limit

        self.parse_output_fn = parse_output_fn

        if isinstance(default_peptide_lengths, int):
            default_peptide_lengths = [default_peptide_lengths]

        self.group_peptides_by_length = group_peptides_by_length

        if self.supported_alleles_flag:
            valid_alleles = self._determine_supported_alleles(
                self.program_name,
                self.supported_alleles_flag)
        else:
            # if we're not running the tool to determine supported alleles
            # then at least try running it by itself to determine if it's
            # it's present
            try:
                run_command([self.program_name])
            except:
                raise SystemError("Failed to run %s" % self.program_name)
            valid_alleles = None

        try:
            BasePredictor.__init__(
                self,
                alleles=alleles,
                valid_alleles=valid_alleles,
                default_peptide_lengths=default_peptide_lengths,
                min_peptide_length=min_peptide_length,
                max_peptide_length=max_peptide_length)
        except UnsupportedAllele as e:
            if self.supported_alleles_flag:
                additional_message = (
                    "\nRun command %s %s to see a list of valid alleles" % (
                        self.program_name,
                        self.supported_alleles_flag))
            else:
                additional_message = ""
            raise UnsupportedAllele(str(e) + additional_message)

    @staticmethod
    def _determine_supported_alleles(command, supported_allele_flag):
        """
        Try asking the commandline predictor (e.g. netMHCpan)
        which alleles it supports.
        """
        try:
            # convert to str since Python3 returns a `bytes` object
            supported_alleles_output = check_output([
                command, supported_allele_flag
            ])
            supported_alleles_str = supported_alleles_output.decode("ascii", "ignore")
            assert len(supported_alleles_str) > 0, \
                '%s returned empty allele list' % command
            supported_alleles = set([])
            for line in supported_alleles_str.split("\n"):
                line = line.strip()
                if not line.startswith('#') and len(line) > 0:
                    try:
                        # We need to normalize these alleles (the output of the predictor
                        # when it lists its supported alleles) so that they are comparable with
                        # our own alleles.
                        supported_alleles.add(normalize_allele_name(line))
                    except AlleleParseError as error:
                        logger.info("Skipping allele %s: %s", line, error)
                        continue
            if len(supported_alleles) == 0:
                raise ValueError("Unable to determine supported alleles")
            return supported_alleles
        except Exception as e:
            logger.exception(e)
            raise SystemError("Failed to run %s %s. Possibly an incorrect executable version?" % (
                command,
                supported_allele_flag))

    def prepare_allele_name(self, allele_name):
        """
        How does the predictor expect to see allele names?
        """
        return allele_name.replace("*", "")

    def _build_command(
            self,
            input_filename,
            allele,
            length=None,
            temp_dirname=None,
            peptide_mode=False):
        args = [self.program_name]
        if peptide_mode:
            args.extend(self.peptide_mode_flags)
        args.extend([self.allele_flag, self.prepare_allele_name(allele)])
        if length:
            args.extend([self.length_flag, str(length)])
        if self.tempdir_flag and temp_dirname:
            args.extend([self.tempdir_flag, temp_dirname])
        args.extend(self.extra_flags)
        if self.input_file_flag:
            args.extend([self.input_file_flag, input_filename])
        else:
            args.append(input_filename)
        return args

    def _run_commands_and_collect_predictions(
            self,
            commands,
            input_filenames,
            temp_dir_list,
            sequence_key_mapping=None):
        if sequence_key_mapping is None:
            sequence_key_mapping = defaultdict(lambda: "seq")
        binding_predictions = []

        # Cleanup either when finished or if an exception gets raised by
        # deleting the input and output files
        filenames_to_delete = list(input_filenames) + [
            f.name for f in commands.keys()]
        with CleanupFiles(
                filenames=filenames_to_delete,
                directories=temp_dir_list):
            run_multiple_commands_redirect_stdout(
                commands,
                print_commands=True,
                process_limit=self.process_limit)
            for output_file, command in commands.items():
                # closing/opening looks insane
                # but I was getting empty files otherwise
                output_file.close()
                with open(output_file.name, 'r') as f:
                    file_contents = f.read()
                    binding_predictions.extend(
                        self.parse_output_fn(
                            stdout=file_contents,
                            sequence_key_mapping=sequence_key_mapping,
                            prediction_method_name=self.program_name))

        if len(binding_predictions) == 0:
            logger.warn("No binding predictions from %s" % self.program_name)
        return BindingPredictionCollection(binding_predictions)

    def predict_peptides(self, peptides):
        self._check_peptide_inputs(peptides)
        input_filenames = create_input_peptides_files(
            peptides,
            max_peptides_per_file=self.max_peptides_per_file,
            group_by_length=self.group_peptides_by_length)
        logger.debug("Created %d input files" % len(input_filenames))
        commands = {}
        dirs = []

        for i, input_filename in enumerate(input_filenames):
            for j, allele in enumerate(self.alleles):
                if self.tempdir_flag:
                    temp_dirname = tempfile.mkdtemp(
                        prefix="tmp_%d_%d_%s" % (
                            i,
                            j,
                            self.program_name))
                    logger.debug(
                        "Created temporary directory %s for allele %s",
                        temp_dirname,
                        allele)
                    dirs.append(temp_dirname)
                else:
                    temp_dirname = None
                output_file = tempfile.NamedTemporaryFile(
                    "w",
                    prefix="%s_output_length_%d_%d" % (
                        self.program_name, i, j),
                    delete=False)
                commands[output_file] = self._build_command(
                    input_filename=input_filename,
                    allele=allele,
                    peptide_mode=True,
                    temp_dirname=temp_dirname)
        results = self._run_commands_and_collect_predictions(
            commands=commands,
            input_filenames=input_filenames,
            temp_dir_list=dirs)
        self._check_results(
            results,
            peptides=peptides,
            alleles=self.alleles)
        return results
