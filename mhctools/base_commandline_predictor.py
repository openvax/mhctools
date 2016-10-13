# Copyright (c) 2014. Mount Sinai School of Medicine
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
import logging
from subprocess import check_output
import tempfile

from typechecks import require_string, require_integer

from .common import check_sequence_dictionary
from .alleles import normalize_allele_name, AlleleParseError
from .base_predictor import BasePredictor, UnsupportedAllele
from .process_helpers import run_command
from .cleanup_context import CleanupFiles
from .epitope_collection import EpitopeCollection
from .file_formats import create_input_fasta_files
from .process_helpers import run_multiple_commands_redirect_stdout


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
            epitope_lengths,
            parse_output_fn,
            supported_alleles_flag,
            input_fasta_flag,
            length_flag,
            allele_flag,
            tempdir_flag=None,
            extra_flags=[],
            max_file_records=None,
            process_limit=0):
        """
        Parameters
        ----------
        program_name : str
            Name of prediction program to run
            (e.g. "netMHCcons" or "netMHCIIpan")

        alleles : list of str
            MHC alleles

        epitope_lengths : list of int

        supported_alleles_flag : str
            Flag to pass to the predictor to get a list of supported alleles
            (e.g. "-A", "-list", "-listMHC")

        parse_output_fn : fn
            Takes the stdout string from the predictor and returns a collection
            of BindingPrediction objects

        input_fasta_flag : str
            How to specify the input FASTA file of source sequences (e.g. "-f")

        length_flag : str
            How to specify the desired predicted epitope length (e.g. "-length")

        allele_flag : str
            How to specify the allele we want predictions for (e.g. "-a")

        tempdir_flag : str, optional
            How to specify the predictor's temporary directory (e.g. "-tdir")

        extra_flags : list of str
            Extra flags to pass to the predictor

        max_file_records : int, optional
            Maximum number of sequences per input FASTA file

        process_limit : int, optional
            Maximum number of parallel processes to start
        """
        require_string(program_name, "Predictor program name")
        self.program_name = program_name

        if supported_alleles_flag is not None:
            require_string(supported_alleles_flag, "Supported alleles flag")
        self.supported_alleles_flag = supported_alleles_flag

        require_string(input_fasta_flag, "Input FASTA file flag")
        self.input_fasta_flag = input_fasta_flag

        require_string(allele_flag, "Allele flag")
        self.allele_flag = allele_flag

        require_string(length_flag, "Peptide length flag")
        self.length_flag = length_flag

        if tempdir_flag is not None:
            require_string(tempdir_flag, "Temporary directory flag")
        self.tempdir_flag = tempdir_flag

        self.extra_flags = extra_flags

        if max_file_records is not None:
            require_integer(
                    max_file_records,
                    "Maximum number of sequences per input files")
        self.max_file_records = max_file_records

        require_integer(process_limit, "Maximum number of processes")
        self.process_limit = process_limit

        self.parse_output_fn = parse_output_fn

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
                alleles,
                epitope_lengths,
                valid_alleles=valid_alleles)
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
            return supported_alleles
        except:
            raise SystemError("Failed to run %s %s. Possibly an incorrect executable version?" % (
                command,
                supported_allele_flag))

    def prepare_allele_name(self, allele_name):
        """
        How does the predictor expect to see allele names?
        """
        return allele_name.replace("*", "")

    def _build_command(self, input_filename, allele, length, temp_dirname=None):
        args = [self.program_name]
        if self.input_fasta_flag:
            args.extend([self.input_fasta_flag, input_filename])
        else:
            args.append(input_filename)
        args.extend([
            self.allele_flag, allele,
            self.length_flag, str(length),
        ])
        if self.tempdir_flag and temp_dirname:
            args.extend([self.tempdir_flag, temp_dirname])
        args.extend(self.extra_flags)
        return args

    def predict(self, fasta_dictionary):
        """
        Run multiple predictors simultaneously, split across:
            1) multiple input files, each containing self.max_file_records
            2) every allele + peptide length combination
        """
        fasta_dictionary = check_sequence_dictionary(fasta_dictionary)
        input_filenames, sequence_key_mapping = create_input_fasta_files(
            fasta_dictionary,
            max_file_records=self.max_file_records)
        output_files = {}
        commands = {}
        dirs = []
        alleles = [
            self.prepare_allele_name(allele)
            for allele in self.alleles
        ]
        for i, input_filename in enumerate(input_filenames):
            for j, allele in enumerate(alleles):
                for length in self.epitope_lengths:
                    if self.tempdir_flag:
                        temp_dirname = tempfile.mkdtemp(
                            prefix="tmp_%s_length_%d" % (
                                self.program_name, length))
                        logger.debug(
                            "Created temporary directory %s for allele %s, length %d",
                            temp_dirname,
                            allele,
                            length)
                        dirs.append(temp_dirname)
                    else:
                        temp_dirname = None
                    output_file = tempfile.NamedTemporaryFile(
                        "w",
                        prefix="%s_output_%d_%d" % (self.program_name, i, j),
                        delete=False)
                    commands[output_file] = self._build_command(
                        input_filename=input_filename,
                        allele=allele,
                        length=length,
                        temp_dirname=temp_dirname)

        epitope_collections = []

        # Cleanup either when finished or if an exception gets raised by
        # deleting the input and output files
        filenames_to_delete = input_filenames
        for f in output_files.keys():
            filenames_to_delete.append(f.name)

        with CleanupFiles(
                filenames=filenames_to_delete,
                directories=dirs):
            run_multiple_commands_redirect_stdout(
                commands,
                print_commands=True,
                process_limit=self.process_limit)
            for output_file, command in commands.items():
                # closing/opening looks insane
                # but I was getting empty files otherwise
                output_file.close()
                with open(output_file.name, 'r') as f:
                    epitope_collection = self.parse_output_fn(
                        stdout=f.read(),
                        fasta_dictionary=fasta_dictionary,
                        sequence_key_mapping=sequence_key_mapping,
                        prediction_method_name=self.program_name)
                    epitope_collections.append(epitope_collection)

        if len(epitope_collections) != len(commands):
            raise ValueError(
                ("Expected an epitope collection for each "
                 "command (%d), but instead there are %d") % (
                    len(commands),
                    len(epitope_collections)))

        if len(epitope_collections) == 0:
            raise ValueError("No epitopes from %s" % self.program_name)

        # flatten all epitope collections into a single object
        return EpitopeCollection([
            binding_prediction
            for sublist in epitope_collections
            for binding_prediction in sublist])
