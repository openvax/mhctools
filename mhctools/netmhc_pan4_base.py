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

"""Shared implementation for the NetMHCpan 4.x wrappers (4.0, 4.1, 4.2)."""

from functools import partial

from .base_commandline_predictor import BaseCommandlinePredictor
from .parsing import parse_netmhcpan_to_preds
from .pred import Kind


class NetMHCpan4Base(BaseCommandlinePredictor):
    """
    Common command construction for NetMHCpan 4.x.

    Subclasses set ``parse_stdout_fn`` (a parser taking ``mode=``) and
    ``affinity_mode_kinds`` (the prediction kinds reported in binding
    affinity mode), and define ``__init__`` to give their own default
    ``program_name``.
    """

    parse_stdout_fn = None
    affinity_mode_kinds = (Kind.pMHC_affinity,)

    def __init__(
            self,
            alleles,
            default_peptide_lengths,
            program_name,
            process_limit,
            mode,
            extra_flags,
            max_peptides_per_file,
            max_alleles_per_command):
        if mode == "binding_affinity":
            # The -BA flag is required to predict binding affinity
            flags = ["-BA"]
        elif mode == "elution_score":
            flags = []
        else:
            raise ValueError("Unsupported mode", mode)
        self.mode = mode

        BaseCommandlinePredictor.__init__(
            self,
            program_name=program_name,
            alleles=alleles,
            default_peptide_lengths=default_peptide_lengths,
            parse_output_fn=partial(type(self).parse_stdout_fn, mode=mode),
            parse_to_preds_fn=parse_netmhcpan_to_preds,
            supported_alleles_flag="-listMHC",
            input_file_flag="-f",
            length_flag="-l",
            allele_flag="-a",
            extra_flags=flags + extra_flags,
            max_peptides_per_file=max_peptides_per_file,
            max_alleles_per_command=max_alleles_per_command,
            process_limit=process_limit)

    def kind_support(self):
        if self.mode == "binding_affinity":
            kinds = self.affinity_mode_kinds
        else:
            kinds = (Kind.pMHC_presentation,)
        return {
            kind: {
                "mhc_dependence": "single_allele",
                "mhc_class": "I",
            }
            for kind in kinds
        }
