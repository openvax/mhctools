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

"""Wrapper for MixMHC2pred — pan-allele MHC class-II presentation (Gfeller lab).

MixMHC2pred predicts MHC class-II ligands / eluted-ligand likelihood. This
wrapper shells out to a user-provided MixMHC2pred install (it is academic /
non-commercial licensed, so mhctools does not vendor it) and emits one
``Kind.pMHC_presentation`` prediction per (peptide, allele), with ``score`` =
MixMHC2pred's raw presentation score (higher = better) and ``percentile_rank``
= its %Rank (lower = better).

Class-II alleles use MixMHC2pred's own nomenclature (``HLA-`` dropped, ``_``
separators, paired α/β chains joined by ``__``): e.g. ``DRB1_15_01``,
``DQA1_01_02__DQB1_06_02``. :func:`to_mixmhc2pred_allele` maps the usual
spellings (``HLA-DRB1*15:01``, the mhcgnomes-canonical
``HLA-DRA1*01:01-DRB1*15:01``, or the native ``DRB1_15_01``) onto that form.

Upstream: https://github.com/GfellerLab/MixMHC2pred
Cite: Racle et al., *Nat. Biotechnol.* 2019; Racle et al., bioRxiv 2022.
"""

from os import remove
from os.path import exists, join
from tempfile import NamedTemporaryFile, mkdtemp

import pandas as pd

from .base_predictor import BasePredictor, _check_flank_inputs
from .cleanup_context import CleanupFiles
from .pred import Kind, PeptideResult, Prediction
from .process_helpers import run_command


def to_mixmhc2pred_allele(allele):
    """Convert a class-II allele name to MixMHC2pred's nomenclature.

    Handles ``HLA-DRB1*15:01``, the mhcgnomes-canonical
    ``HLA-DRA1*01:01-DRB1*15:01`` (the explicit DRA1 α-chain is dropped —
    MixMHC2pred names DR alleles by their DRB chain only), paired DP/DQ chains,
    and names already in MixMHC2pred form. Examples::

        HLA-DRB1*15:01              -> DRB1_15_01
        HLA-DRA1*01:01-DRB1*15:01   -> DRB1_15_01
        HLA-DQA1*01:02-DQB1*06:02   -> DQA1_01_02__DQB1_06_02
        DRB1_15_01                  -> DRB1_15_01
    """
    text = str(allele).strip()
    if text.upper().startswith("HLA-"):
        text = text[4:]
    # mhcgnomes joins paired α/β chains with '-'; users may use '/'.
    text = text.replace("/", "-")
    chains = [
        chain.replace("*", "_").replace(":", "_")
        for chain in text.split("-") if chain
    ]
    # The DRA1 α-chain is implicit in MixMHC2pred's DR naming.
    chains = [chain for chain in chains if not chain.upper().startswith("DRA")]
    return "__".join(chains)


class MixMHC2pred(BasePredictor):
    """Wrapper for the MixMHC2pred class-II presentation predictor.

    Parameters
    ----------
    alleles : list of str
        Class-II alleles (any spelling :func:`to_mixmhc2pred_allele` accepts).
    default_peptide_lengths : list of int
    program_name : str
        MixMHC2pred executable (name on ``PATH`` or an absolute path — use
        ``MixMHC2pred_unix`` on Linux, ``MixMHC2pred`` on macOS). The allele
        definition folder (``PWMdef``) is resolved relative to the executable.
    """

    mhc_class = "II"

    def __init__(
            self,
            alleles,
            default_peptide_lengths=[15],
            program_name="MixMHC2pred"):
        BasePredictor.__init__(
            self,
            alleles=alleles,
            default_peptide_lengths=default_peptide_lengths,
            min_peptide_length=9,
            max_peptide_length=25,
            allow_X_in_peptides=False,
            allow_lowercase_in_peptides=False,
            # Keep MixMHC2pred-native names (DRB1_15_01, ...) verbatim; HLA
            # spellings are still canonicalized. to_mixmhc2pred_allele maps
            # either onto the CLI form.
            keep_unparseable_alleles=True)
        self.program_name = program_name

    def predict(self, peptides, n_flanks=None, c_flanks=None):
        """Predict class-II presentation for a list of peptides.

        Flanks are accepted for a uniform API but ignored (this wrapper runs
        MixMHC2pred in ``--no_context`` mode).

        Returns
        -------
        list of PeptideResult
            One entry per input peptide; each holds one
            ``Kind.pMHC_presentation`` prediction per allele.
        """
        peptide_list, _, _ = _check_flank_inputs(peptides, n_flanks, c_flanks)
        self._check_peptide_inputs(peptide_list)
        if not peptide_list:
            return []

        alleles = list(self.alleles)
        cli_names = [to_mixmhc2pred_allele(a) for a in alleles]

        temp_dir = mkdtemp(prefix="mhctools", suffix="mixmhc2pred")
        input_file_path = join(temp_dir, "mixmhc2pred_input.txt")
        output_file_path = join(temp_dir, "mixmhc2pred_output.txt")
        with open(input_file_path, "w") as f:
            f.write("\n".join(peptide_list) + "\n")

        args = [
            self.program_name,
            "-i", input_file_path,
            "-o", output_file_path,
            "-a", *cli_names,
            "--no_context",
            "--extra_out",
        ]

        with CleanupFiles(
                filenames=[input_file_path, output_file_path],
                directories=[temp_dir]):
            with NamedTemporaryFile(
                    prefix="MixMHC2pred_stdout", mode="w",
                    delete=False) as stdout_file:
                stdout_file_name = stdout_file.name
                run_command(
                    args,
                    suppress_stderr=False,
                    redirect_stdout_file=stdout_file)
            if exists(output_file_path):
                preds_by_peptide = parse_mixmhc2pred_results(
                    output_file_path, alleles, cli_names)
            else:
                with open(stdout_file_name) as f:
                    stdout = f.read().strip()
                raise ValueError(
                    "MixMHC2pred produced no output for alleles %s; stdout: %s"
                    % (cli_names, stdout))
            remove(stdout_file_name)

        return [
            PeptideResult(preds=tuple(preds_by_peptide.get(peptide, ())))
            for peptide in peptide_list
        ]

    def _default_pred_kind(self):
        return Kind.pMHC_presentation

    def kind_support(self):
        return {
            Kind.pMHC_presentation: {
                "mhc_dependence": "single_allele",
                "mhc_class": "II",
            },
        }


def parse_mixmhc2pred_results(filename, alleles, cli_names):
    """Parse a MixMHC2pred (``--extra_out``) output into ``{peptide: [preds]}``.

    MixMHC2pred's output has a ``#``-comment header, then columns including,
    per allele ``M`` (its MixMHC2pred name), ``%Rank_M`` and ``Score_M``. We
    look those up by name and attach the caller's original allele identity.

    Parameters
    ----------
    filename : str
    alleles : list of str
        The original allele identities (carried through to the predictions).
    cli_names : list of str
        The MixMHC2pred names passed for each allele, in the same order — the
        column-header suffixes to read.

    Returns
    -------
    dict mapping peptide -> list of Prediction (one per allele)
    """
    df = pd.read_csv(filename, comment="#", sep="\t")
    if "Peptide" not in df.columns:
        raise ValueError("Unexpected MixMHC2pred output columns: %s"
                         % list(df.columns))
    for cli_name in cli_names:
        for prefix in ("%Rank_", "Score_"):
            column = prefix + cli_name
            if column not in df.columns:
                raise ValueError(
                    "MixMHC2pred output missing column %r (columns: %s)"
                    % (column, list(df.columns)))

    # itertuples mangles '%Rank_...' names, so index by position. Resolve the
    # column positions once, not per row.
    peptide_idx = df.columns.get_loc("Peptide")
    rank_idx = [df.columns.get_loc("%Rank_" + name) for name in cli_names]
    score_idx = [df.columns.get_loc("Score_" + name) for name in cli_names]

    results = {}
    for row in df.itertuples(index=False):
        peptide = row[peptide_idx]
        preds = [
            Prediction(
                kind=Kind.pMHC_presentation,
                score=float(row[score_idx[i]]),
                peptide=peptide,
                allele=allele,
                percentile_rank=float(row[rank_idx[i]]),
                predictor_name="mixmhc2pred")
            for i, allele in enumerate(alleles)]
        results[peptide] = preds
    return results
