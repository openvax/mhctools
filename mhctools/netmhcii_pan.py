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

from functools import partial
import logging
import os
from subprocess import check_output

from mhcnames import parse_classi_or_classii_allele_name

from .base_commandline_predictor import BaseCommandlinePredictor
from .parsing import parse_netmhciipan_stdout, parse_netmhciipan4_stdout

logger = logging.getLogger(__name__)


class NetMHCIIpanBase(BaseCommandlinePredictor):
    def __init__(
            self,
            alleles,
            parse_output_fn,
            default_peptide_lengths=[15, 16, 17, 18, 19, 20],
            program_name="netMHCIIpan",
            process_limit=-1,
            extra_flags=[]):
        BaseCommandlinePredictor.__init__(
            self,
            program_name=program_name,
            alleles=alleles,
            default_peptide_lengths=default_peptide_lengths,
            parse_output_fn=parse_output_fn,
            supported_alleles_flag="-list",
            input_file_flag="-f",
            allele_flag="-a",
            peptide_mode_flags=["-inptype", "1"],
            length_flag="-length",
            tempdir_flag="-tdir",
            process_limit=process_limit,
            extra_flags=extra_flags,
            min_peptide_length=9)

    def _prepare_drb_allele_name(self, parsed_beta_allele):
        """
        Assume that we're dealing with a human DRB allele
        which NetMHCIIpan treats differently because there is
        little population diversity in the DR-alpha gene
        """
        if "DRB" not in parsed_beta_allele.gene:
            raise ValueError("Unexpected allele %s" % parsed_beta_allele)
        return "%s_%s%s" % (
            parsed_beta_allele.gene,
            parsed_beta_allele.allele_family,
            parsed_beta_allele.allele_code)

    def prepare_allele_name(self, allele_name):
        """
        netMHCIIpan has some unique requirements for allele formats,
        expecting the following forms:
         - DRB1_0101 (for non-alpha/beta pairs)
         - HLA-DQA10501-DQB10636 (for alpha and beta pairs)

        Other than human class II alleles, the only other alleles that
        netMHCIIpan accepts are the following mouse alleles:
         - H-2-IAb
         - H-2-IAd
        """
        parsed_alleles = parse_classi_or_classii_allele_name(allele_name)
        if len(parsed_alleles) == 1:
            allele = parsed_alleles[0]
            if allele.species == "H-2":
                return "%s-%s%s" % (
                    allele.species,
                    allele.gene,
                    allele.allele_code)
            return self._prepare_drb_allele_name(allele)

        else:
            alpha, beta = parsed_alleles
            if "DRA" in alpha.gene:
                return self._prepare_drb_allele_name(beta)
            return "HLA-%s%s%s-%s%s%s" % (
                alpha.gene,
                alpha.allele_family,
                alpha.allele_code,
                beta.gene,
                beta.allele_family,
                beta.allele_code)


class NetMHCIIpan3(NetMHCIIpanBase):
    def __init__(
            self,
            alleles,
            default_peptide_lengths=[15, 16, 17, 18, 19, 20],
            program_name="netMHCIIpan",
            process_limit=-1,
            extra_flags=[]):
        NetMHCIIpanBase.__init__(
            self,
            alleles=alleles,
            parse_output_fn=parse_netmhciipan_stdout,
            default_peptide_lengths=default_peptide_lengths,
            program_name=program_name,
            process_limit=process_limit,
            extra_flags=extra_flags)


class NetMHCIIpan4(NetMHCIIpanBase):
    def __init__(
            self,
            alleles,
            default_peptide_lengths=[15, 16, 17, 18, 19, 20],
            program_name="netMHCIIpan",
            process_limit=-1,
            mode="elution_score",
            extra_flags=[]):
        """
        Wrapper for NetMHCIIpan 4.0, using a different parser.
        """

        if mode not in ['binding_affinity', 'elution_score', 'all_output']:
            raise ValueError("Unsupported mode", mode)

        # Default to including binding affinity data (-BA flag), though the score and %rank will
        # still be EL-based
        NetMHCIIpanBase.__init__(
            self,
            alleles=alleles,
            program_name=program_name,
            process_limit=process_limit,
            parse_output_fn=partial(parse_netmhciipan4_stdout, mode=mode),
            default_peptide_lengths=default_peptide_lengths,
            extra_flags=['-BA'] + extra_flags)

        self.mode = mode


class NetMHCIIpan4_EL(NetMHCIIpan4):
    """
    Wrapper for NetMHCIIpan4 when the preferred mode is elution score
    """
    def __init__(
            self,
            alleles,
            default_peptide_lengths=[15, 16, 17, 18, 19, 20],
            program_name="netMHCIIpan",
            process_limit=-1,
            extra_flags=[]):
        NetMHCIIpan4.__init__(
            self,
            alleles=alleles,
            default_peptide_lengths=default_peptide_lengths,
            program_name=program_name,
            process_limit=process_limit,
            mode="elution_score",
            extra_flags=extra_flags)


class NetMHCIIpan4_BA(NetMHCIIpan4):
    """
    Wrapper for NetMHCIIpan4 when the preferred mode is binding affinity
    """
    def __init__(
            self,
            alleles,
            default_peptide_lengths=[15, 16, 17, 18, 19, 20],
            program_name="netMHCIIpan",
            process_limit=-1,
            extra_flags=[]):
        NetMHCIIpan4.__init__(
            self,
            alleles=alleles,
            default_peptide_lengths=default_peptide_lengths,
            program_name=program_name,
            process_limit=process_limit,
            mode="binding_affinity",
            extra_flags=extra_flags)


def NetMHCIIpan(
        alleles,
        process_limit=-1,
        program_name="netMHCIIpan",
        default_peptide_lengths=[15, 16, 17, 18, 19, 20],
        extra_flags=[]):
    """
    This function wraps NetMHCIIpan4 and NetMHCIIpan3 to automatically detect which class to use.
    """
    # convert to str since Python3 returns a "bytes" object.
    with open(os.devnull, 'w') as devnull:
        output = check_output([program_name, "-h"], stderr=devnull)
    output_str = output.decode("ascii", "ignore")
    kwargs = {
        "alleles": alleles,
        "default_peptide_lengths": default_peptide_lengths,
        "program_name": program_name,
        "process_limit": process_limit,
        "extra_flags": extra_flags,
    }
    if "NetMHCIIpan-4.0" in output_str:
        logger.info("Using NetMHCIIpan 4.0")
        return NetMHCIIpan4(**kwargs)
    elif "NetMHCIIpan-3" in output_str:
        logger.info("Using NetMHCIIpan 3.x")
        return NetMHCIIpan3(**kwargs)
    else:
        raise ValueError("This software expects NetMHCIIpan version 3.x or 4.0")
