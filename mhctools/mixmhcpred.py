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

"""Wrapper for MixMHCpred 3.0 class-I ligand presentation predictions.

MixMHCpred 3.0 is a pan-allele predictor which can also align user-provided
MHC-I sequences and render binding-motif artifacts. Its academic license does
not permit redistribution, so mhctools shells out to a user-provided official
checkout and does not vendor MixMHCpred code, models, or reference data.

The canonical :meth:`MixMHCpred.predict` API returns one
``Kind.pMHC_presentation`` prediction per peptide and allele. The
:meth:`MixMHCpred.predict_detailed` API additionally preserves MixMHCpred's
complete table, including its independently-computed best score, best allele,
best rank, closest training allele, distance score, and pan-allele status.
Sequence alignment, sequence-driven prediction, and motif/HTML generation are
available through :meth:`MixMHCpred.predict_allele_sequences`.

Upstream: https://github.com/GfellerLab/MixMHCpred
"""

from __future__ import annotations

import os
import re
import shutil
from dataclasses import dataclass
from pathlib import Path
from subprocess import STDOUT, CalledProcessError, check_output
from tempfile import TemporaryDirectory

import pandas as pd

from .allele_normalization import normalize_allele_name_or_raw
from .base_predictor import BasePredictor, _check_flank_inputs
from .binding_prediction import BindingPrediction
from .binding_prediction_collection import BindingPredictionCollection
from .pred import Kind, PeptideResult, Prediction
from .process_helpers import run_command

_VERSION_RE = re.compile(r"MixMHCpred\s*(?:\(v)?([0-9]+(?:\.[0-9]+)+)")
_CLOSEST_RE = re.compile(r"(.+?)\s*\(([-+0-9.eE]+)\)\s*$")
_SEQUENCE_QUALITY_RE = re.compile(
    r"<h2>(.*?)</h2>\s*"
    r"<h2>closest Allele = (.*?), distance = (.*?)</h2>\s*"
    r"<h2>closest Allele from database = (.*?), distance = (.*?)</h2>",
    re.IGNORECASE | re.DOTALL,
)


@dataclass(frozen=True)
class MixMHCpredAlleleInfo:
    """Quality/provenance information for one MixMHCpred output allele."""

    allele: str
    native_allele: str
    closest_training_allele: str = ""
    distance: float | None = None
    pan_allele: bool = False
    closest_database_allele: str = ""
    database_distance: float | None = None


@dataclass(frozen=True)
class MixMHCpredArtifacts:
    """Persistent files produced by a MixMHCpred 3.0 artifact run."""

    output_dir: str
    files: tuple
    binding_predictions: str = ""
    alignment: str = ""
    overview_html: str = ""


@dataclass
class MixMHCpredResult:
    """Complete parsed output from one MixMHCpred invocation.

    ``table`` is the exact upstream tabular output, retaining
    ``Score_bestAllele``, ``BestAllele``, ``%Rank_bestAllele``, and every
    allele-specific score/rank column. ``predictions`` maps those per-allele
    columns into mhctools' canonical API.
    """

    table: pd.DataFrame
    allele_info: tuple
    version: str = ""
    comments: tuple = ()
    stdout: str = ""
    aligned_sequences: tuple = ()
    artifacts: MixMHCpredArtifacts | None = None

    @property
    def predictions(self):
        """Return the detailed table as canonical ``PeptideResult`` objects."""
        if self.table.empty:
            return []

        peptide_index = self.table.columns.get_loc("Peptide")
        allele_columns = tuple(
            (
                info,
                self.table.columns.get_loc(f"Score_{info.native_allele}"),
                self.table.columns.get_loc(f"%Rank_{info.native_allele}"),
            )
            for info in self.allele_info
        )
        results = []
        for row in self.table.itertuples(index=False, name=None):
            peptide = str(row[peptide_index])
            preds = []
            for info, score_index, rank_index in allele_columns:
                preds.append(Prediction(
                    kind=Kind.pMHC_presentation,
                    score=float(row[score_index]),
                    peptide=peptide,
                    allele=info.allele,
                    percentile_rank=float(row[rank_index]),
                    predictor_name="mixmhcpred",
                    predictor_version=self.version,
                ))
            results.append(PeptideResult(preds=tuple(preds)))
        return results


def resolve_mixmhcpred_path(program_name=None):
    """Resolve a MixMHCpred executable.

    Resolution order is an explicit ``program_name``, ``MIXMHCPRED_PATH``,
    then ``MixMHCpred`` on ``PATH``. A path may name either the executable or
    the root of an official checkout.
    """
    candidate = program_name or os.environ.get("MIXMHCPRED_PATH")
    if not candidate:
        candidate = shutil.which("MixMHCpred")
    if not candidate:
        raise FileNotFoundError(
            "MixMHCpred was not found. Download the official v3.0 release "
            "from https://github.com/GfellerLab/MixMHCpred, accept its "
            "academic/non-commercial license, and set MIXMHCPRED_PATH to "
            "the checkout or executable.")

    path = Path(candidate).expanduser()
    if path.is_dir():
        path = path / "MixMHCpred"
    if not path.exists() and os.sep not in str(candidate):
        resolved = shutil.which(str(candidate))
        if resolved:
            path = Path(resolved)
    if not path.is_file():
        raise FileNotFoundError(
            f"MixMHCpred executable does not exist: {path}")
    if " " in str(path.resolve()):
        raise ValueError(
            "MixMHCpred 3.0 does not support spaces in its installation path: "
            f"{path}")
    return str(path.resolve())


def mixmhcpred_version(program_name):
    """Return the version reported by ``MixMHCpred --help``."""
    program = resolve_mixmhcpred_path(program_name)
    try:
        output = check_output(
            [program, "--help"], stderr=STDOUT, text=True)
    except (CalledProcessError, OSError) as e:
        raise RuntimeError(
            f"Could not query MixMHCpred version from {program}: {e}")
    match = _VERSION_RE.search(output)
    return match.group(1) if match else ""


def _native_alleles_from_table(table):
    return tuple(
        column[len("Score_"):]
        for column in table.columns
        if column.startswith("Score_") and column != "Score_bestAllele"
    )


def _normalize_output_allele(allele):
    try:
        return normalize_allele_name_or_raw(allele)
    except (TypeError, ValueError):
        return allele


def _parse_comments(filename):
    comments = []
    with open(filename) as f:
        for line in f:
            if line.startswith("#"):
                comments.append(line.rstrip("\n"))
    return tuple(comments)


def _comment_value(comments, labels):
    for comment in comments:
        text = comment.lstrip("#").strip()
        for label in labels:
            if text.lower().startswith(label.lower()):
                return text[len(label):].strip()
    return ""


def _sequence_quality_from_html(path):
    if not path or not Path(path).is_file():
        return {}
    text = Path(path).read_text()
    result = {}
    for match in _SEQUENCE_QUALITY_RE.finditer(text):
        allele, closest, distance, database, database_distance = match.groups()
        result[allele] = (
            closest,
            float(distance),
            database,
            float(database_distance),
        )
    return result


def parse_mixmhcpred_output(
        filename,
        alleles=None,
        artifacts=None,
        stdout="",
        aligned_sequences=(),
        sequence_mode=False):
    """Parse the complete output of MixMHCpred 2.x or 3.0.

    The returned :class:`MixMHCpredResult` preserves the raw table and header
    metadata. If ``alleles`` are supplied, their caller-facing spellings are
    retained while native names are taken from the output columns. Set
    ``sequence_mode=True`` for predictions generated from user-provided MHC-I
    sequences so their pan-allele provenance is retained.
    """
    comments = _parse_comments(filename)
    table = pd.read_csv(filename, comment="#", sep="\t")
    required = {
        "Peptide", "Score_bestAllele", "BestAllele", "%Rank_bestAllele",
    }
    missing = required - set(table.columns)
    if missing:
        raise ValueError(
            "Unexpected MixMHCpred output; missing columns "
            f"{sorted(missing)} (columns: {list(table.columns)})")

    native_alleles = _native_alleles_from_table(table)
    for native in native_alleles:
        rank_column = f"%Rank_{native}"
        if rank_column not in table.columns:
            raise ValueError(
                f"MixMHCpred output missing column {rank_column!r} "
                f"(columns: {list(table.columns)})")

    if alleles is None:
        public_alleles = tuple(
            _normalize_output_allele(name) for name in native_alleles)
    else:
        public_alleles = tuple(alleles)
        if len(public_alleles) != len(native_alleles):
            raise ValueError(
                f"MixMHCpred returned {len(native_alleles)} allele(s), "
                f"expected {len(public_alleles)}: {native_alleles}")

    version = ""
    for comment in comments:
        match = _VERSION_RE.search(comment)
        if match:
            version = match.group(1)
            break

    allele_comment = _comment_value(
        comments, ("Alleles:", "Predictions for Alleles :"))
    pan_native = set()
    if " - predicted motif for " in allele_comment:
        _, pan_text = allele_comment.split(" - predicted motif for ", 1)
        pan_native = {x.strip() for x in pan_text.split(",") if x.strip()}

    closest_text = _comment_value(
        comments, ("Closest Allele (distance score):",))
    closest = []
    if closest_text:
        for item in closest_text.split(" -- "):
            match = _CLOSEST_RE.match(item)
            if not match:
                raise ValueError(
                    "Could not parse MixMHCpred closest-allele metadata: "
                    f"{item!r}")
            closest.append((match.group(1), float(match.group(2))))
        if len(closest) != len(native_alleles):
            raise ValueError(
                f"MixMHCpred returned {len(closest)} closest-allele entries "
                f"for {len(native_alleles)} alleles")

    sequence_quality = _sequence_quality_from_html(
        artifacts.overview_html if artifacts else "")
    allele_info = []
    for i, (public, native) in enumerate(zip(public_alleles, native_alleles)):
        closest_name = ""
        distance = None
        database_name = ""
        database_distance = None
        if closest:
            closest_name, distance = closest[i]
        elif native in sequence_quality:
            (closest_name, distance,
             database_name, database_distance) = sequence_quality[native]
        allele_info.append(MixMHCpredAlleleInfo(
            allele=public,
            native_allele=native,
            closest_training_allele=_normalize_output_allele(closest_name),
            distance=distance,
            pan_allele=sequence_mode or native in pan_native,
            closest_database_allele=_normalize_output_allele(database_name),
            database_distance=database_distance,
        ))

    return MixMHCpredResult(
        table=table,
        allele_info=tuple(allele_info),
        version=version,
        comments=comments,
        stdout=stdout,
        aligned_sequences=tuple(aligned_sequences),
        artifacts=artifacts,
    )


def parse_mixmhcpred_results(filename):
    """Parse MixMHCpred output into legacy ``BindingPrediction`` objects."""
    result = parse_mixmhcpred_output(filename)
    return [
        BindingPrediction(
            peptide=str(peptide),
            allele=_normalize_output_allele(allele),
            score=float(score),
            percentile_rank=float(percentile_rank),
            prediction_method_name="mixmhcpred",
        )
        for peptide, allele, score, percentile_rank in zip(
            result.table["Peptide"],
            result.table["BestAllele"],
            result.table["Score_bestAllele"],
            result.table["%Rank_bestAllele"],
        )
    ]


def _validate_prediction_rows(result, peptides):
    observed = result.table["Peptide"].astype(str).tolist()
    if observed != list(peptides):
        raise RuntimeError(
            "MixMHCpred output peptides do not match its inputs. "
            f"Expected {list(peptides)}, got {observed}")
    return result


def _sequence_allele_info(names, artifacts=None):
    quality = _sequence_quality_from_html(
        artifacts.overview_html if artifacts else "")
    result = []
    for name in names:
        closest, distance, database, database_distance = quality.get(
            name, ("", None, "", None))
        result.append(MixMHCpredAlleleInfo(
            allele=name,
            native_allele=name,
            closest_training_allele=_normalize_output_allele(closest),
            distance=distance,
            pan_allele=True,
            closest_database_allele=_normalize_output_allele(database),
            database_distance=database_distance,
        ))
    return tuple(result)


def _read_fasta(path):
    records = []
    name = None
    sequence = []
    with open(path) as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if name is not None:
                    records.append((name, "".join(sequence)))
                name = line[1:]
                sequence = []
            elif name is None:
                raise ValueError("Sequence data appears before the first FASTA header")
            else:
                sequence.append(line)
    if name is not None:
        records.append((name, "".join(sequence)))
    return tuple(records)


def _write_fasta(path, records):
    with open(path, "w") as f:
        f.write("".join(f">{name}\n{sequence}\n" for name, sequence in records))


def _collect_artifacts(output_dir):
    root = Path(output_dir).resolve()
    files = tuple(sorted(
        str(path.resolve()) for path in root.rglob("*") if path.is_file()))

    def existing(name):
        path = root / name
        return str(path) if path.is_file() else ""

    return MixMHCpredArtifacts(
        output_dir=str(root),
        files=files,
        binding_predictions=existing("Binding_predictions.txt"),
        alignment=existing("final_alignment.fasta"),
        overview_html=existing("Data_overview.html"),
    )


class MixMHCpred(BasePredictor):
    """MixMHCpred class-I ligand-presentation predictor.

    Parameters
    ----------
    alleles : sequence of str, optional
        MHC-I alleles. Version 3.0 accepts its training alleles plus MHC-I
        alleles present in its sequence database for pan-allele inference.
    default_peptide_lengths : sequence of int
    program_name : str, optional
        Executable or official checkout directory. If omitted, resolve
        ``MIXMHCPRED_PATH`` and then ``MixMHCpred`` on ``PATH``.
    exclude_peptides_with_cysteine : bool
        Exclude cysteine-containing peptides in mhctools before invocation.
        This is version-independent; MixMHCpred 3.0 removed the old ``-c``
        option.
    """

    mhc_class = "I"

    def __init__(
            self,
            alleles=None,
            default_peptide_lengths=(9,),
            program_name=None,
            exclude_peptides_with_cysteine=False):
        if isinstance(alleles, str):
            alleles = alleles.split(",")
        normalized = []
        for allele in alleles or ():
            value = normalize_allele_name_or_raw(allele)
            if value not in normalized:
                normalized.append(value)
        BasePredictor.__init__(
            self,
            alleles=normalized,
            default_peptide_lengths=default_peptide_lengths,
            min_peptide_length=8,
            max_peptide_length=14,
            allow_X_in_peptides=False,
            allow_lowercase_in_peptides=False,
            keep_unparseable_alleles=True,
        )
        # BasePredictor historically de-duplicates through a set. Restore the
        # caller's stable order so native output columns have deterministic
        # identities.
        self.alleles = normalized
        self.program_name = program_name
        self.exclude_peptides_with_cysteine = exclude_peptides_with_cysteine
        self._version = None

    @property
    def executable(self):
        return resolve_mixmhcpred_path(self.program_name)

    @property
    def version(self):
        if self._version is None:
            self._version = mixmhcpred_version(self.executable)
        return self._version

    def _require_v3(self, feature):
        version = self.version
        if not version or int(version.split(".", 1)[0]) < 3:
            raise RuntimeError(
                f"{feature} requires MixMHCpred 3.0 or newer; detected "
                f"{version or 'an unknown version'} at {self.executable}")

    @staticmethod
    def _validate_output_dir(output_dir):
        path = Path(output_dir).expanduser().resolve()
        if path.exists():
            raise FileExistsError(
                "Refusing to let MixMHCpred replace existing output path: "
                f"{path}")
        if " " in str(path):
            raise ValueError(
                "MixMHCpred 3.0 does not support spaces in output paths: "
                f"{path}")
        return str(path)

    def _run_command(self, args, stdout_path):
        try:
            with open(stdout_path, "w") as stdout_file:
                run_command(
                    args,
                    suppress_stderr=False,
                    redirect_stdout_file=stdout_file,
                )
        except (CalledProcessError, OSError) as e:
            stdout = (
                Path(stdout_path).read_text()
                if Path(stdout_path).exists() else "")
            raise RuntimeError(
                f"MixMHCpred failed: {e}\n{stdout.strip()}") from e
        return Path(stdout_path).read_text()

    def predict_detailed(self, peptides, output_dir=None, output_motifs=False):
        """Run allele-based prediction and retain every MixMHCpred output.

        Set ``output_motifs=True`` and provide a new ``output_dir`` to retain
        v3.0's motif matrices/images, peptide-length distributions, HTML, and
        binding-prediction table. The method refuses an existing output path
        because upstream deletes its output directory.
        """
        peptide_list, _, _ = _check_flank_inputs(peptides)
        self._check_peptide_inputs(peptide_list)
        if self.exclude_peptides_with_cysteine:
            peptide_list = [p for p in peptide_list if "C" not in p]
        if not peptide_list or not self.alleles:
            return MixMHCpredResult(
                table=pd.DataFrame(columns=(
                    "Peptide", "Score_bestAllele", "BestAllele",
                    "%Rank_bestAllele")),
                allele_info=(),
            )
        if output_motifs and output_dir is None:
            raise ValueError("output_dir is required when output_motifs=True")
        if output_motifs:
            self._require_v3("Motif generation")

        persistent_dir = (
            self._validate_output_dir(output_dir) if output_motifs else None)
        with TemporaryDirectory(prefix="mhctools-mixmhcpred-") as temp_dir:
            input_path = Path(temp_dir) / "peptides.txt"
            stdout_path = Path(temp_dir) / "stdout.txt"
            input_path.write_text("\n".join(peptide_list) + "\n")
            output_path = (
                Path(persistent_dir) if output_motifs
                else Path(temp_dir) / "predictions.txt")
            args = [
                self.executable,
                "-i", str(input_path),
                "-o", str(output_path),
                "-a", ",".join(self.alleles),
            ]
            if output_motifs:
                args.extend(("-m", "1"))
            stdout = self._run_command(args, stdout_path)
            prediction_path = (
                output_path / "Binding_predictions.txt"
                if output_motifs else output_path)
            if not prediction_path.is_file():
                raise RuntimeError(
                    "MixMHCpred exited without creating predictions at "
                    f"{prediction_path}. stdout: {stdout.strip()}")
            artifacts = (
                _collect_artifacts(persistent_dir) if output_motifs else None)
            result = parse_mixmhcpred_output(
                prediction_path,
                alleles=self.alleles,
                artifacts=artifacts,
                stdout=stdout,
            )
            if not result.version:
                result.version = self.version
            return _validate_prediction_rows(result, peptide_list)

    def predict(self, peptides, n_flanks=None, c_flanks=None):
        """Predict class-I presentation, preserving one result per input."""
        peptide_list, _, _ = _check_flank_inputs(
            peptides, n_flanks, c_flanks)
        self._check_peptide_inputs(peptide_list)
        if not peptide_list:
            return []
        if not self.alleles:
            return [PeptideResult() for _ in peptide_list]

        if self.exclude_peptides_with_cysteine:
            allowed = [p for p in peptide_list if "C" not in p]
        else:
            allowed = peptide_list
        predicted = iter(self.predict_detailed(allowed).predictions)
        results = []
        for peptide in peptide_list:
            if self.exclude_peptides_with_cysteine and "C" in peptide:
                results.append(PeptideResult())
            else:
                results.append(next(predicted))
        return results

    def predict_peptides(self, peptides):
        """Legacy ``BindingPredictionCollection`` prediction API."""
        return BindingPredictionCollection([
            BindingPrediction.from_pred(pred)
            for peptide_result in self.predict(peptides)
            for pred in peptide_result.preds
        ])

    def predict_allele_sequences(
            self,
            allele_sequences,
            peptides=None,
            output_dir=None,
            output_motifs=False):
        """Align MHC-I sequences and optionally predict peptides/motifs.

        Parameters
        ----------
        allele_sequences : mapping or str/path-like
            Mapping of allele identifiers to unaligned MHC-I sequences, or a
            FASTA filename.
        peptides : sequence of str, optional
            Peptides to score against sequence-derived motifs. Omit for
            alignment and/or motif generation only.
        output_dir : str/path-like, optional
            New directory in which to retain all upstream artifacts. Required
            for ``output_motifs=True``. If omitted for alignment/prediction,
            files are parsed from a temporary directory and only structured
            Python output is retained.
        output_motifs : bool
            Generate PWM tables/images, peptide-length distributions, and the
            HTML overview.
        """
        if hasattr(allele_sequences, "items"):
            records = tuple(
                (str(name), str(sequence))
                for name, sequence in allele_sequences.items())
        else:
            records = _read_fasta(allele_sequences)
        if not records:
            raise ValueError("At least one MHC-I allele sequence is required")
        names = [name for name, _ in records]
        if any(not name for name in names) or len(set(names)) != len(names):
            raise ValueError("MHC-I FASTA identifiers must be non-empty and unique")
        for name, sequence in records:
            if not sequence or not sequence.isalpha() or not sequence.isupper():
                raise ValueError(
                    f"Invalid MHC-I sequence for {name!r}: expected uppercase "
                    "amino acids")

        peptide_list = None
        if peptides is not None:
            peptide_list, _, _ = _check_flank_inputs(peptides)
            self._check_peptide_inputs(peptide_list)
            if self.exclude_peptides_with_cysteine:
                peptide_list = [p for p in peptide_list if "C" not in p]
            if not peptide_list:
                peptide_list = None
        if output_motifs and output_dir is None:
            raise ValueError("output_dir is required when output_motifs=True")
        self._require_v3("Allele-sequence alignment and prediction")

        persistent_dir = (
            self._validate_output_dir(output_dir)
            if output_dir is not None else None)
        with TemporaryDirectory(
                prefix="mhctools-mixmhcpred-sequences-") as temp_dir:
            sequence_path = Path(temp_dir) / "alleles.fasta"
            peptide_path = Path(temp_dir) / "peptides.txt"
            stdout_path = Path(temp_dir) / "stdout.txt"
            _write_fasta(sequence_path, records)
            if peptide_list is not None:
                peptide_path.write_text("\n".join(peptide_list) + "\n")
            run_output_dir = Path(persistent_dir or (Path(temp_dir) / "output"))
            args = [
                self.executable,
                "-s", str(sequence_path),
                "-o", str(run_output_dir),
            ]
            if peptide_list is not None:
                args.extend(("-i", str(peptide_path), "-p", "1"))
            if output_motifs:
                args.extend(("-m", "1"))
            stdout = self._run_command(args, stdout_path)

            alignment_path = run_output_dir / "final_alignment.fasta"
            if not alignment_path.is_file() or alignment_path.stat().st_size == 0:
                raise RuntimeError(
                    "MixMHCpred exited without a non-empty sequence alignment. "
                    "Ensure MAFFT is installed and runnable. stdout: "
                    f"{stdout.strip()}")
            aligned = _read_fasta(alignment_path)
            artifacts = (
                _collect_artifacts(persistent_dir) if persistent_dir else None)

            if peptide_list is None:
                return MixMHCpredResult(
                    table=pd.DataFrame(columns=(
                        "Peptide", "Score_bestAllele", "BestAllele",
                        "%Rank_bestAllele")),
                    allele_info=_sequence_allele_info(names, artifacts),
                    version=self.version,
                    stdout=stdout,
                    aligned_sequences=aligned,
                    artifacts=artifacts,
                )

            prediction_path = run_output_dir / "Binding_predictions.txt"
            if not prediction_path.is_file():
                raise RuntimeError(
                    "MixMHCpred exited without sequence-driven predictions at "
                    f"{prediction_path}. stdout: {stdout.strip()}")
            result = parse_mixmhcpred_output(
                prediction_path,
                alleles=names,
                artifacts=artifacts,
                stdout=stdout,
                aligned_sequences=aligned,
                sequence_mode=True,
            )
            if not result.version:
                result.version = self.version
            return _validate_prediction_rows(result, peptide_list)

    def _default_pred_kind(self):
        return Kind.pMHC_presentation

    def kind_support(self):
        return {
            Kind.pMHC_presentation: {
                "mhc_dependence": "single_allele",
                "mhc_class": "I",
            },
        }
