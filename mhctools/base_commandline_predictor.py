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

from .alleles import normalize_allele_name, AlleleParseError
from .base_predictor import BasePredictor, UnsupportedAllele
from .process_helpers import run_command


class BaseCommandlinePredictor(BasePredictor):
    """
    Base class for MHC binding predictors that run a local external
    program and write their output to a local file.
    """
    def __init__(
            self,
            name,
            command,
            alleles,
            epitope_lengths,
            supported_allele_flag='-listMHC'):
        self.name = name
        self.command = command
        self.supported_allele_flag = supported_allele_flag

        if not isinstance(command, str):
            raise TypeError(
                "Expected %s command to be string, got %s : %s" % (
                    name, command, type(command)))

        if supported_allele_flag:
            valid_alleles = self._determine_valid_alleles(
                command, supported_allele_flag)
        else:
            # if we're not running the tool to determine supported alleles
            # then at least try running it by itself to determine if it's
            # it's present
            try:
                run_command([command])
            except:
                raise SystemError("Failed to run %s" % command)
            valid_alleles = None

        try:
            BasePredictor.__init__(
                self,
                alleles,
                epitope_lengths,
                valid_alleles=valid_alleles)
        except UnsupportedAllele as e:
            additional_message = ("\nRun command %s %s to see a list of valid alleles" % (
                command, supported_allele_flag) if supported_allele_flag else "")
            raise UnsupportedAllele(e.message + additional_message)

    @staticmethod
    def _determine_valid_alleles(command, supported_allele_flag):
        """
        Try asking the commandline predictor (e.g. netMHCpan)
        which alleles it supports.
        """
        try:
            # convert to str since Python3 returns a `bytes` object
            valid_alleles = check_output([
                command, supported_allele_flag
            ])
            valid_alleles_str = valid_alleles.decode("ascii", "ignore")
            assert len(valid_alleles_str) > 0, \
                '%s returned empty allele list' % command
            valid_alleles = set([])
            for line in valid_alleles_str.split("\n"):
                line = line.strip()
                if not line.startswith('#') and len(line) > 0:
                    try:
                        # We need to normalize these alleles (the output of the predictor
                        # when it lists its supported alleles) so that they are comparable with
                        # our own alleles.
                        valid_alleles.add(normalize_allele_name(line))
                    except AlleleParseError as error:
                        logging.info("Skipping allele %s: %s" % (
                            line, error))
                        continue
            return valid_alleles
        except:
            logging.warning(
                "Failed to run %s %s",
                command,
                supported_allele_flag)
            raise
