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

import pandas as pd

from .pred import COLUMNS


class MultiSample:
    """
    Run a predictor across multiple samples, each with its own alleles.

    Parameters
    ----------
    samples : dict
        Mapping of sample_name -> list of allele strings.
        Example: {"pat001": ["HLA-A*02:01", "HLA-B*07:02"],
                  "pat002": ["HLA-A*01:01", "HLA-B*08:01"]}

    predictor_class : class
        A predictor class (e.g. NetMHCpan41) that accepts an ``alleles``
        keyword argument.

    **predictor_kwargs
        Additional keyword arguments forwarded to the predictor constructor.
    """

    def __init__(self, samples, predictor_class, **predictor_kwargs):
        self.samples = samples
        self.predictor_class = predictor_class
        self.predictor_kwargs = predictor_kwargs

    def _make_predictor(self, alleles):
        return self.predictor_class(alleles=alleles, **self.predictor_kwargs)

    # --- peptide predictions ---

    def predict(self, peptides):
        """
        Returns
        -------
        dict mapping sample_name -> list of PeptidePreds
        """
        results = {}
        for sample_name, alleles in self.samples.items():
            predictor = self._make_predictor(alleles)
            results[sample_name] = predictor.predict(peptides)
        return results

    def predict_dataframe(self, peptides):
        """predict() flattened to a DataFrame with sample_name column."""
        dfs = []
        for sample_name, pp_list in self.predict(peptides).items():
            for pp in pp_list:
                dfs.append(pp.to_dataframe(sample_name))
        if not dfs:
            return pd.DataFrame(columns=COLUMNS)
        return pd.concat(dfs, ignore_index=True)

    # --- protein scanning ---

    def predict_proteins(self, sequence_dict, peptide_lengths=None):
        """
        Returns
        -------
        dict mapping sample_name -> {sequence_name: list of PeptidePreds}
        """
        results = {}
        for sample_name, alleles in self.samples.items():
            predictor = self._make_predictor(alleles)
            results[sample_name] = predictor.predict_proteins(
                sequence_dict, peptide_lengths=peptide_lengths)
        return results

    def predict_proteins_dataframe(self, sequence_dict, peptide_lengths=None):
        """predict_proteins() flattened to a DataFrame with sample_name column."""
        dfs = []
        for sample_name, seq_dict in self.predict_proteins(
                sequence_dict, peptide_lengths).items():
            for seq_name, pp_list in seq_dict.items():
                for pp in pp_list:
                    dfs.append(pp.to_dataframe(sample_name))
        if not dfs:
            return pd.DataFrame(columns=COLUMNS)
        return pd.concat(dfs, ignore_index=True)
