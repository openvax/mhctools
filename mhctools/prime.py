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

"""Wrapper for PRIME — PRedictor of Immunogenic Epitopes (Gfeller lab).

PRIME predicts CD8+ T-cell immunogenicity of class-I peptides by combining
MHC-I binding (via MixMHCpred) with a TCR-recognition propensity model. This
wrapper shells out to a user-provided PRIME install (PRIME is academic /
non-commercial licensed, so mhctools does not vendor it) and emits one
``Kind.immunogenicity`` prediction per (peptide, allele).

PRIME itself calls MixMHCpred (v3.0+ recommended); point at it with
``mixmhcpred_path`` if it is not on ``PATH``.

Upstream: https://github.com/GfellerLab/PRIME
Cite: Gfeller et al., Cell Systems 2023 — "Improved predictions of antigen
presentation and TCR recognition with MixMHCpred2.2 and PRIME2.0".

Like every current CD8 immunogenicity predictor, PRIME ranks better than it
generalizes to novel neoepitopes (independent benchmarks put the field near
AUC 0.5-0.65) — a prioritization aid, not ground truth.
"""

from os import remove
from os.path import exists, join
from tempfile import NamedTemporaryFile, mkdtemp

import pandas as pd

from .allele_normalization import normalize_allele_name
from .base_predictor import BasePredictor, _check_flank_inputs
from .cleanup_context import CleanupFiles
from .pred import Kind, PeptideResult, Prediction
from .process_helpers import run_command


class PRIME(BasePredictor):
    """Wrapper for the PRIME immunogenicity predictor.

    Parameters
    ----------
    alleles : list of str
        Class-I alleles.
    default_peptide_lengths : list of int
    program_name : str
        PRIME executable (name on ``PATH`` or an absolute path).
    mixmhcpred_path : str, optional
        Absolute path to the MixMHCpred executable PRIME should use, forwarded
        as PRIME's ``-mix`` option. If omitted, PRIME finds MixMHCpred on
        ``PATH``.
    """

    def __init__(
            self,
            alleles,
            default_peptide_lengths=[9],
            program_name="PRIME",
            mixmhcpred_path=None):
        BasePredictor.__init__(
            self,
            alleles=alleles,
            default_peptide_lengths=default_peptide_lengths,
            min_peptide_length=8,
            max_peptide_length=14,
            allow_X_in_peptides=False,
            allow_lowercase_in_peptides=False)
        self.program_name = program_name
        self.mixmhcpred_path = mixmhcpred_path

    def predict(self, peptides, n_flanks=None, c_flanks=None):
        """Predict immunogenicity for a list of peptides.

        Flanks are accepted for a uniform API but ignored (PRIME does not use
        flanking context).

        Returns
        -------
        list of PeptideResult
            One entry per input peptide; each holds one
            ``Kind.immunogenicity`` prediction per allele (``score`` = PRIME
            score, higher = more immunogenic; ``percentile_rank`` = PRIME
            %Rank, lower = better).
        """
        peptide_list, _, _ = _check_flank_inputs(peptides, n_flanks, c_flanks)
        self._check_peptide_inputs(peptide_list)
        if not peptide_list:
            return []

        # PRIME accepts HLA-A*02:01 etc.; normalize for a stable identity that
        # matches the allele strings we hand back in each Prediction.
        alleles = [normalize_allele_name(a) for a in self.alleles]

        temp_dir = mkdtemp(prefix="mhctools", suffix="prime")
        input_file_path = join(temp_dir, "prime_input.txt")
        output_file_path = join(temp_dir, "prime_output.txt")
        with open(input_file_path, "w") as f:
            f.write("\n".join(peptide_list) + "\n")

        args = [
            self.program_name,
            "-i", input_file_path,
            "-o", output_file_path,
            "-a", ",".join(alleles),
        ]
        if self.mixmhcpred_path:
            args += ["-mix", self.mixmhcpred_path]

        with CleanupFiles(
                filenames=[input_file_path, output_file_path],
                directories=[temp_dir]):
            with NamedTemporaryFile(
                    prefix="PRIME_stdout", mode="w", delete=False) as stdout_file:
                stdout_file_name = stdout_file.name
                run_command(
                    args,
                    suppress_stderr=False,
                    redirect_stdout_file=stdout_file)
            if exists(output_file_path):
                preds_by_peptide = parse_prime_results(output_file_path, alleles)
            else:
                with open(stdout_file_name) as f:
                    stdout = f.read().strip()
                raise ValueError(
                    "PRIME produced no output for alleles %s; stdout: %s"
                    % (alleles, stdout))
            remove(stdout_file_name)

        # Preserve input peptide order; a peptide PRIME dropped becomes empty.
        return [
            PeptideResult(preds=tuple(preds_by_peptide.get(peptide, ())))
            for peptide in peptide_list
        ]

    def _default_pred_kind(self):
        return Kind.immunogenicity

    def kind_support(self):
        return {
            Kind.immunogenicity: {
                "mhc_dependence": "single_allele",
                "mhc_class": "I",
            },
        }


def parse_prime_results(filename, alleles):
    """Parse a PRIME output file into ``{peptide: [Prediction, ...]}``.

    PRIME's output has a ``#``-comment header block, then columns::

        Peptide  %Rank_bestAllele  Score_bestAllele  %RankBinding_bestAllele
        BestAllele  %Rank_<A1>  Score_<A1>  %RankBinding_<A1>  %Rank_<A2> ...

    The per-allele column triples (``%Rank`` / ``Score`` / ``%RankBinding``)
    follow the order the alleles were passed to PRIME. PRIME rewrites allele
    names into its own short form in the column headers, so we map the triples
    to ``alleles`` **by position** rather than by header name.

    Parameters
    ----------
    filename : str
    alleles : list of str
        The (normalized) alleles handed to PRIME, in order. Their identity is
        carried through to the returned predictions.

    Returns
    -------
    dict mapping peptide -> list of Prediction (one per allele)
    """
    df = pd.read_csv(filename, comment="#", sep="\t")
    columns = list(df.columns)
    if "Peptide" not in columns or "BestAllele" not in columns:
        raise ValueError(
            "Unexpected PRIME output columns: %s" % columns)

    # Per-allele columns start right after the fixed best-allele block.
    base = columns.index("BestAllele") + 1
    per_allele_columns = columns[base:]
    if len(per_allele_columns) != 3 * len(alleles):
        raise ValueError(
            "PRIME returned %d per-allele columns for %d alleles (%s)"
            % (len(per_allele_columns), len(alleles), per_allele_columns))

    peptide_idx = columns.index("Peptide")
    results = {}
    for row in df.itertuples(index=False):
        peptide = row[peptide_idx]
        preds = []
        for i, allele in enumerate(alleles):
            rank = float(row[base + 3 * i])       # %Rank_<allele>
            score = float(row[base + 3 * i + 1])  # Score_<allele>
            # row[base + 3 * i + 2] is %RankBinding (MixMHCpred %Rank); we
            # surface only PRIME's own immunogenicity signal here.
            preds.append(Prediction(
                kind=Kind.immunogenicity,
                score=score,
                peptide=peptide,
                allele=allele,
                percentile_rank=rank,
                predictor_name="prime"))
        results[peptide] = preds
    return results
