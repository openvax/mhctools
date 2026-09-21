"""Local SMM 1.0 and SMM-PMBEC 1.0 from IEDB's standalone MHC-I tools.

The upstream program owns the models and scoring algorithm. Install it with
``python scripts/setup_test_backends.py smm --accept-license`` or point
``IEDB_MHCI_EXECUTABLE`` at a launcher for its ``src/predict_binding.py``.
"""

from collections import defaultdict
import csv
import io
import math
import os
from pathlib import Path
import shutil
import subprocess
from tempfile import TemporaryDirectory

from .allele_normalization import normalize_allele_name
from .base_predictor import BasePredictor
from .binding_prediction import BindingPrediction
from .binding_prediction_collection import BindingPredictionCollection


def resolve_smm_executable(program_name=None):
    """Resolve exactly the executable that inference will use."""
    program = program_name or os.environ.get("IEDB_MHCI_EXECUTABLE") or "iedb-mhci"
    path = shutil.which(os.path.expanduser(str(program)))
    if path is None:
        raise FileNotFoundError(
            "Local IEDB SMM tools not found: %s. Set IEDB_MHCI_EXECUTABLE "
            "to a launcher for predict_binding.py; see docs/testing.md." % program)
    return os.path.abspath(path)


def parse_smm_output(text, peptides, method):
    """Validate native TSV output against the submitted FASTA records."""
    reader = csv.DictReader(io.StringIO(text), delimiter="\t")
    required = {"allele", "seq_num", "start", "end", "length", "peptide", "ic50", "rank"}
    if not required.issubset(reader.fieldnames or ()):
        raise ValueError("Invalid local IEDB output (expected prediction table): %s" % text)
    predictions = []
    seen = set()
    for row in reader:
        try:
            index = int(row["seq_num"]) - 1
            peptide = row["peptide"]
            allele = normalize_allele_name(row["allele"])
            affinity, rank = float(row["ic50"]), float(row["rank"])
            if (not 0 <= index < len(peptides) or peptides[index] != peptide
                    or int(row["start"]) != 1 or int(row["end"]) != len(peptide)
                    or int(row["length"]) != len(peptide)
                    or not math.isfinite(affinity) or affinity <= 0
                    or not math.isfinite(rank) or not 0 <= rank <= 100
                    or (index, allele) in seen):
                raise ValueError("Invalid prediction identity, coordinates or values")
            seen.add((index, allele))
        except (TypeError, ValueError, KeyError) as error:
            raise ValueError("Invalid local IEDB prediction: %r" % row) from error
        predictions.append(BindingPrediction(
            peptide=peptide, allele=allele, affinity=affinity,
            percentile_rank=rank, source_sequence_name="seq%d" % (index + 1),
            prediction_method_name=method))
    return predictions


class SMM(BasePredictor):
    """Run the official standalone SMM predictor locally, with IC50 in nM.

    Parameters
    ----------
    alleles : list of str
        MHC alleles. Unsupported allele/length pairs fail explicitly.
    default_peptide_lengths : list of int
        Window lengths used by protein prediction (default: nine residues).
    program_name : str, optional
        Executable launcher for IEDB's ``predict_binding.py``. Defaults to
        ``IEDB_MHCI_EXECUTABLE`` or ``iedb-mhci`` on PATH.
    timeout : float
        Maximum seconds for each local prediction batch.
    """

    prediction_method = "smm"
    predictor_version = "1.0"

    def __init__(self, alleles=None, default_peptide_lengths=None, program_name=None, timeout=120):
        super().__init__(alleles=alleles, default_peptide_lengths=default_peptide_lengths or [9])
        if isinstance(timeout, bool) or not isinstance(timeout, (int, float)) or not math.isfinite(timeout) or timeout <= 0:
            raise ValueError("timeout must be a finite positive number")
        self.program_name = resolve_smm_executable(program_name)
        self.timeout = timeout

    def predict_peptides(self, peptides):
        peptides = list(peptides)
        if not self.alleles:
            raise ValueError("SMM requires at least one MHC allele")
        self._check_peptide_inputs(peptides)
        groups = defaultdict(list)
        for peptide in dict.fromkeys(peptides):
            groups[len(peptide)].append(peptide)
        predictions = []
        with TemporaryDirectory(prefix="mhctools-smm-") as directory:
            fasta = Path(directory) / "peptides.fasta"
            for length, batch in sorted(groups.items()):
                fasta.write_text("".join(
                    ">seq%d\n%s\n" % (i + 1, peptide) for i, peptide in enumerate(batch)))
                alleles = sorted(self.alleles)
                command = [self.program_name, self.prediction_method, ",".join(alleles),
                           ",".join([str(length)] * len(alleles)), str(fasta)]
                result = subprocess.run(command, capture_output=True, text=True, timeout=self.timeout)
                if result.returncode:
                    raise RuntimeError("Local IEDB %s failed (%d):\n%s\n%s" % (
                        self.prediction_method, result.returncode, result.stdout, result.stderr))
                try:
                    batch_predictions = parse_smm_output(result.stdout, batch, self.prediction_method)
                    self._check_results(batch_predictions, batch, self.alleles)
                except ValueError as error:
                    raise ValueError("%s\nLocal IEDB stderr: %s" % (error, result.stderr)) from error
                predictions.extend(batch_predictions)
        by_pair = {(p.peptide, p.allele): p for p in predictions}
        return BindingPredictionCollection([
            by_pair[peptide, allele].clone_with_updates(source_sequence_name="seq%d" % (i + 1))
            for i, peptide in enumerate(peptides) for allele in sorted(self.alleles)
        ])


class SMMPMBEC(SMM):
    """Run IEDB's official local SMM-PMBEC 1.0 model."""

    prediction_method = "smmpmbec"
