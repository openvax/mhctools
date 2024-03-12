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

from .base_commandline_predictor import BaseCommandlinePredictor
from .parsing import parse_netmhcstabpan

class NetMHCstabpan(BaseCommandlinePredictor):
    def __init__(
            self,
            alleles,
            default_peptide_lengths=[],
            program_name="netMHCstabpan",
            process_limit=-1,
            flags=[]):
        """
        Wrapper for NetMHCstabpan.
        """
        BaseCommandlinePredictor.__init__(
            self,
            program_name=program_name,
            alleles=alleles,
            default_peptide_lengths=default_peptide_lengths,
            parse_output_fn=parse_netmhcstabpan,
            supported_alleles_flag="-listMHC",
            input_file_flag="-f",
            length_flag="-l",
            allele_flag="-a",
            extra_flags=flags,
            process_limit=process_limit)
        
    def predict_peptides(self, peptides):
        peptide_lengths = set(len(p) for p in peptides)
        if len(peptide_lengths) > 1:
            raise ValueError("All peptides must be the same length")
        return super().predict_peptides(peptides)
    
