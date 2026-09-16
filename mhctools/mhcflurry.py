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

import logging
import math
import os

from .base_predictor import BasePredictor
from .base_predictor import _check_flank_inputs
from .binding_prediction import BindingPrediction
from .binding_prediction_collection import BindingPredictionCollection
from .pred import Kind, Prediction
from .unsupported_allele import UnsupportedAllele

logger = logging.getLogger(__name__)

# Loaded model and inferred provenance, keyed by (kind, resolved models path).
_model_cache = {}
_PERCENT_RANK_SUPPORT_UNKNOWN = object()


def _normalize_models_path(models_path):
    """Normalize a models directory path for cache-key deduplication.

    Different strings pointing at the same directory (relative vs.
    absolute, with ``~``, via a symlink) should share a single cache
    entry. ``None`` stays ``None`` — mhcflurry resolves that to the
    package-default path internally, and we treat it as its own key.
    """
    if models_path is None:
        return None
    return os.path.realpath(os.path.expanduser(models_path))


def mhcflurry_composite_version(models_path=None):
    """Identify an official MHCflurry model bundle and its Python package.

    Parameters
    ----------
    models_path : str or path-like, optional
        Directory actually loaded. If absent, use the configured default
        presentation directory, including environment overrides. Official
        presentation and affinity bundle paths share the release identity.

    Returns
    -------
    str
        Package and release identity, such as ``2.2.1+release-2.2.0``.
        This identifies the configured official release, not a content hash
        or verification that its files have never been modified.

    Raises
    ------
    RuntimeError
        If the package/release is unknown or the selected path is custom.
        Custom and injected models need caller-supplied provenance; labeling
        them with the active default release would describe different weights.
    """
    try:
        import mhcflurry
        from mhcflurry import downloads
    except ImportError as error:
        raise RuntimeError(
            "mhcflurry is not installed; supply predictor_version explicitly."
        ) from error
    package_version = getattr(mhcflurry, "__version__", None)
    if not package_version:
        raise RuntimeError("mhcflurry exposes no __version__; supply predictor_version explicitly.")
    try:
        release = downloads.get_current_release()
    except Exception as error:
        raise RuntimeError(
            "Could not read mhcflurry's current model release; "
            "supply predictor_version explicitly."
        ) from error
    if not release:
        raise RuntimeError(
            "mhcflurry has no active model release. Run mhcflurry-downloads fetch "
            "or supply predictor_version explicitly for custom models.")
    presentation_path = downloads.get_path(
        "models_class1_presentation", "models", test_exists=False)
    official_paths = {
        _normalize_models_path(presentation_path),
        _normalize_models_path(os.path.join(presentation_path, "affinity_predictor")),
        _normalize_models_path(downloads.get_path(
            "models_class1_pan", "models.combined", test_exists=False)),
    }
    selected_path = (
        downloads.get_default_class1_presentation_models_dir(test_exists=False)
        if models_path is None else models_path)
    if _normalize_models_path(selected_path) not in official_paths:
        raise RuntimeError(
            "Cannot infer provenance for custom MHCflurry models at %r; "
            "supply predictor_version explicitly." % str(selected_path))
    return "%s+release-%s" % (package_version, release)


def _load_predictor(model_kind, predictor, models_path, predictor_version):
    """Keep loaded weights and the provenance captured at load time together."""
    if predictor_version is not None and (
            not isinstance(predictor_version, str) or not predictor_version.strip()):
        raise ValueError("predictor_version must be a non-empty string or None")
    if predictor is not None:
        return predictor, predictor_version

    from mhcflurry import Class1AffinityPredictor, Class1PresentationPredictor, downloads
    if model_kind == "presentation":
        loader = Class1PresentationPredictor.load
        default_path = downloads.get_default_class1_presentation_models_dir
    else:
        loader = Class1AffinityPredictor.load
        default_path = downloads.get_default_class1_models_dir
    if models_path is None:
        try:
            models_path = default_path()
        except RuntimeError as missing_affinity:
            if model_kind != "affinity":
                raise
            # Match Class1AffinityPredictor.load(): when the standalone bundle
            # is absent, use the presentation bundle's actual affinity model.
            # Its captured provenance follows that object, including overrides.
            try:
                presentation, version = _load_predictor(
                    "presentation", None, None, predictor_version)
            except RuntimeError:
                raise missing_affinity
            return presentation.affinity_predictor, version
    path = _normalize_models_path(models_path)
    cache_key = (model_kind, path)
    if cache_key not in _model_cache:
        try:
            inferred_version = mhcflurry_composite_version(path)
        except RuntimeError:
            # Custom/default-overridden weights remain usable but unversioned.
            # A cache consumer must require the caller to identify those weights.
            inferred_version = None
        logger.info("Loading MHCflurry %s models from %s", model_kind, path)
        _model_cache[cache_key] = (loader(path), inferred_version)
    loaded, inferred_version = _model_cache[cache_key]
    return loaded, predictor_version if predictor_version is not None else inferred_version


def _affinity_percent_rank_calibrated_allele(affinity_predictor, allele):
    """Return the allele whose affinity percentile calibration can be used."""
    helper = getattr(
        affinity_predictor, "percent_rank_calibrated_allele", None)
    if callable(helper):
        return helper(allele)

    transforms = getattr(
        affinity_predictor, "allele_to_percent_rank_transform", None)
    if transforms is None:
        return _PERCENT_RANK_SUPPORT_UNKNOWN

    canonicalize = getattr(
        affinity_predictor, "canonicalize_allele_name", None)
    normalized_allele = (
        canonicalize(allele) if callable(canonicalize) else allele)
    if normalized_allele in transforms:
        return normalized_allele

    allele_to_sequence = getattr(affinity_predictor, "allele_to_sequence", None)
    if (
            not allele_to_sequence
            or normalized_allele not in allele_to_sequence):
        return None

    sequence = allele_to_sequence[normalized_allele]
    for other_allele in sorted(allele_to_sequence):
        if (
                allele_to_sequence[other_allele] == sequence
                and other_allele in transforms):
            return other_allele
    return None


def _check_affinity_percent_rank_support(affinity_predictor, alleles):
    """Raise if requested alleles cannot get affinity percentile ranks."""
    missing_alleles = []
    for allele in alleles:
        calibrated = _affinity_percent_rank_calibrated_allele(
            affinity_predictor, allele)
        if calibrated is _PERCENT_RANK_SUPPORT_UNKNOWN:
            return
        if calibrated is None:
            missing_alleles.append(allele)
    if missing_alleles:
        raise ValueError(
            "MHCflurry affinity percentile ranks are unavailable for "
            "allele(s): %s. Raw affinity prediction may still be supported. "
            "Pass include_affinity_percentile_ranks=False to omit affinity "
            "percentile ranks, or calibrate MHCflurry percentile ranks for "
            "these alleles."
            % ", ".join(sorted(missing_alleles)))


class MHCflurry(BasePredictor):
    """
    MHCflurry predictor using the modern Class1PresentationPredictor API.

    Produces per-allele ``pMHC_affinity`` predictions. For presentation,
    ``presentation_allele_mode`` controls whether mhctools treats the allele
    set as one class-I haplotype or as a panel of independent one-allele
    samples. It also surfaces MHCflurry's ``antigen_processing`` (cleavage)
    score — computed by the presentation model from the peptide + flanks and
    allele-independent — as one allele-less prediction per peptide. The legacy
    ``predict_peptides`` method returns BindingPrediction objects based on
    affinity values for backward compat.

    See https://github.com/openvax/mhcflurry
    """
    uses_flanking_sequences = True
    flank_length = 15
    max_haplotype_alleles = 6
    presentation_allele_modes = frozenset((
        "auto",
        "haplotype",
        "per_allele",
    ))

    @classmethod
    def fetch(cls, version=None):
        """Fetch MHCflurry presentation models using its native manager."""
        from .artifacts import fetch
        return fetch("mhcflurry", version=version)

    def __init__(
            self,
            alleles,
            default_peptide_lengths=[9],
            predictor=None,
            models_path=None,
            include_affinity_percentile_ranks=True,
            presentation_allele_mode="auto",
            predictor_version=None):
        """
        Parameters
        -----------
        alleles : list of str

        default_peptide_lengths : list of int

        predictor : mhcflurry.Class1PresentationPredictor (optional)
            MHCflurry presentation predictor to use

        models_path : string
            Models dir to use if predictor argument is None

        predictor_version : str, optional
            Explicit model provenance, required to cache custom or injected
            weights. Official default models infer package plus model release
            at load time. Custom/injected weights otherwise remain unversioned.

        include_affinity_percentile_ranks : bool
            Whether to request affinity percentile ranks. Enabled by default.
            If enabled, requested alleles must have MHCflurry affinity
            percentile-rank calibration, either directly or through an allele
            with the same pseudosequence.

        presentation_allele_mode : {"auto", "haplotype", "per_allele"}
            How to interpret the requested alleles for presentation scoring.
            ``"haplotype"`` treats the alleles as one sample genotype and
            emits one presentation prediction per peptide. ``"per_allele"``
            treats each allele as a separate one-allele synthetic sample and
            emits one presentation prediction per peptide/allele pair.
            ``"auto"`` uses haplotype mode for up to six alleles and
            per-allele mode for larger allele panels.
        """
        BasePredictor.__init__(
            self,
            alleles=alleles,
            default_peptide_lengths=default_peptide_lengths,
            min_peptide_length=8,
            max_peptide_length=15)
        self.predictor, self.predictor_version = _load_predictor(
            "presentation", predictor, models_path, predictor_version)

        self.include_affinity_percentile_ranks = \
            include_affinity_percentile_ranks
        self.presentation_allele_mode = self._resolve_presentation_allele_mode(
            presentation_allele_mode)

        for allele in self.alleles:
            if allele not in self.predictor.supported_alleles:
                raise UnsupportedAllele(allele)
        if self.include_affinity_percentile_ranks:
            _check_affinity_percent_rank_support(
                self.predictor.affinity_predictor, self.alleles)

    def _resolve_presentation_allele_mode(self, presentation_allele_mode):
        if presentation_allele_mode not in self.presentation_allele_modes:
            raise ValueError(
                "presentation_allele_mode must be one of %s, got %r" % (
                    sorted(self.presentation_allele_modes),
                    presentation_allele_mode))
        if presentation_allele_mode == "auto":
            if len(self.alleles) <= self.max_haplotype_alleles:
                return "haplotype"
            return "per_allele"
        if (
                presentation_allele_mode == "haplotype" and
                len(self.alleles) > self.max_haplotype_alleles):
            raise ValueError(
                "MHCflurry presentation haplotype mode accepts at most %d "
                "alleles, got %d. Use presentation_allele_mode='per_allele' "
                "for allele panels." % (
                    self.max_haplotype_alleles,
                    len(self.alleles)))
        return presentation_allele_mode

    def _predict_protein_flank_lengths(self):
        processing_predictor = getattr(
            self.predictor, "processing_predictor_with_flanks", None)
        sequence_lengths = getattr(
            processing_predictor, "sequence_lengths", None)
        if sequence_lengths:
            return (
                int(sequence_lengths.get("n_flank", self.flank_length)),
                int(sequence_lengths.get("c_flank", self.flank_length)),
            )
        return super()._predict_protein_flank_lengths()

    def predict_peptides(self, peptides):
        """
        Predict MHC binding affinity and presentation for peptides.

        Returns a BindingPredictionCollection (legacy API) using affinity
        values for backward compatibility.
        """
        peptide_list = list(peptides)
        allele_list = list(self.alleles)

        # Build cross product for batch prediction
        batch_peptides = peptide_list * len(allele_list)
        batch_alleles = [a for a in allele_list for _ in peptide_list]

        df = self.predictor.affinity_predictor.predict_to_dataframe(
            peptides=batch_peptides,
            alleles=batch_alleles,
            include_percentile_ranks=self.include_affinity_percentile_ranks,
        )
        binding_predictions = []
        for row in df.itertuples(index=False):
            binding_predictions.append(BindingPrediction(
                allele=row.allele,
                peptide=row.peptide,
                affinity=row.prediction,
                percentile_rank=(
                    row.prediction_percentile
                    if hasattr(row, 'prediction_percentile') else None),
                prediction_method_name="mhcflurry",
            ))
        return BindingPredictionCollection(binding_predictions)

    def predict(self, peptides, n_flanks=None, c_flanks=None):
        """
        Predict for a list of peptide sequences.

        Returns a list of PeptideResult, each containing one pMHC_affinity
        Prediction per allele. Presentation predictions are haplotype-level
        when ``presentation_allele_mode`` is ``"haplotype"`` and per-allele
        when it is ``"per_allele"``.

        Uses batch prediction across alleles for affinity. For presentation,
        haplotype mode passes the allele list as one MHCflurry sample genotype;
        per-allele mode passes a sample-to-one-allele dict.
        """
        from .pred import PeptideResult

        peptide_list, n_flank_list, c_flank_list = _check_flank_inputs(
            peptides, n_flanks, c_flanks)
        if n_flank_list is not None or c_flank_list is not None:
            if n_flank_list is None:
                n_flank_list = [""] * len(peptide_list)
            if c_flank_list is None:
                c_flank_list = [""] * len(peptide_list)
        allele_list = list(self.alleles)

        # Build cross product
        batch_peptides = peptide_list * len(allele_list)
        batch_alleles = [a for a in allele_list for _ in peptide_list]
        batch_indices = list(range(len(peptide_list))) * len(allele_list)

        # Single batched call for affinity
        aff_df = self.predictor.affinity_predictor.predict_to_dataframe(
            peptides=batch_peptides,
            alleles=batch_alleles,
            include_percentile_ranks=self.include_affinity_percentile_ranks,
        )

        if self.presentation_allele_mode == "haplotype":
            presentation_alleles = allele_list
        else:
            presentation_alleles = {
                allele: [allele]
                for allele in allele_list
            }

        kwargs = {
            "peptides": peptide_list,
            "alleles": presentation_alleles,
            "include_affinity_percentile": False,
            "verbose": 0,
        }
        if n_flank_list is not None:
            kwargs["n_flanks"] = n_flank_list
        if c_flank_list is not None:
            kwargs["c_flanks"] = c_flank_list
        pres_df = self.predictor.predict(**kwargs)
        expected_presentation_rows = len(peptide_list)
        if self.presentation_allele_mode == "per_allele":
            expected_presentation_rows *= len(allele_list)
        if len(pres_df) != expected_presentation_rows:
            raise ValueError(
                "MHCflurry returned %d presentation row(s) for %d "
                "expected peptide/context row(s)" % (
                    len(pres_df),
                    expected_presentation_rows))

        pres_by_peptide_index = {i: [] for i in range(len(peptide_list))}
        # MHCflurry's presentation predict() also returns a processing_score
        # (its antigen-processing / cleavage head). It depends only on the
        # peptide + flanks, not the allele, so it's identical across the
        # per-allele rows of a peptide; we keep the first seen per peptide.
        processing_by_peptide_index = {}
        seen_presentation_keys = set()
        for row_position, row in enumerate(pres_df.itertuples(index=False)):
            row_index = int(getattr(row, "peptide_num", row_position))
            if self.presentation_allele_mode == "haplotype":
                allele = getattr(row, "best_allele", getattr(row, "allele", ""))
                key = (row_index, "")
            else:
                allele = getattr(
                    row,
                    "sample_name",
                    getattr(row, "best_allele", getattr(row, "allele", "")))
                key = (row_index, allele)
            if key in seen_presentation_keys:
                raise ValueError(
                    "MHCflurry returned duplicate presentation row for "
                    "peptide index %d and allele '%s'" % (row_index, allele))
            seen_presentation_keys.add(key)
            pres_by_peptide_index[row_index].append((
                row.presentation_score,
                row.presentation_percentile,
                allele,
            ))
            processing_by_peptide_index.setdefault(
                row_index, getattr(row, "processing_score", None))

        groups = [list() for _ in peptide_list]
        for row_index, row in zip(batch_indices, aff_df.itertuples(index=False)):
            pep = row.peptide
            allele = row.allele
            affinity_nM = row.prediction
            affinity_pct = (
                row.prediction_percentile
                if hasattr(row, 'prediction_percentile') else None)

            aff_score = max(0.0, min(1.0,
                1.0 - math.log(max(affinity_nM, 1e-6)) / math.log(50000)))

            n_flank = (
                n_flank_list[row_index] if n_flank_list is not None else "")
            c_flank = (
                c_flank_list[row_index] if c_flank_list is not None else "")

            groups[row_index].append(Prediction(
                kind=Kind.pMHC_affinity,
                score=aff_score,
                peptide=pep,
                allele=allele,
                n_flank=n_flank,
                c_flank=c_flank,
                value=affinity_nM,
                percentile_rank=affinity_pct,
                predictor_name="mhcflurry",
                predictor_version=self.predictor_version,
            ))

        for row_index, pep in enumerate(peptide_list):
            if not pres_by_peptide_index[row_index]:
                raise ValueError(
                    "MHCflurry: missing presentation score for "
                    "peptide index %d, peptide='%s'" % (row_index, pep))
            n_flank = (
                n_flank_list[row_index] if n_flank_list is not None else "")
            c_flank = (
                c_flank_list[row_index] if c_flank_list is not None else "")
            for pres_score, pres_pct, presentation_allele in (
                    pres_by_peptide_index[row_index]):
                groups[row_index].append(Prediction(
                    kind=Kind.pMHC_presentation,
                    score=pres_score,
                    peptide=pep,
                    allele=presentation_allele or "",
                    n_flank=n_flank,
                    c_flank=c_flank,
                    percentile_rank=pres_pct,
                    predictor_name="mhcflurry",
                    predictor_version=self.predictor_version,
                ))

            # Surface MHCflurry's antigen-processing (cleavage) score, which
            # its presentation predictor already computes from the peptide +
            # flanks. Allele-independent, so emit once per peptide, allele-less.
            processing_score = processing_by_peptide_index.get(row_index)
            if processing_score is not None:
                groups[row_index].append(Prediction(
                    kind=Kind.antigen_processing,
                    score=processing_score,
                    peptide=pep,
                    allele="",
                    n_flank=n_flank,
                    c_flank=c_flank,
                    predictor_name="mhcflurry",
                    predictor_version=self.predictor_version,
                ))

        return [PeptideResult(preds=tuple(preds)) for preds in groups]

    def predict_with_flanks(self, peptides, n_flanks, c_flanks):
        return self.predict(
            peptides,
            n_flanks=n_flanks,
            c_flanks=c_flanks)

    def _default_pred_kind(self):
        return Kind.pMHC_affinity

    def kind_support(self):
        presentation_dependence = (
            "haplotype"
            if self.presentation_allele_mode == "haplotype"
            else "single_allele")
        return {
            Kind.pMHC_affinity: {
                "mhc_dependence": "single_allele",
                "mhc_class": "I",
            },
            Kind.pMHC_presentation: {
                "mhc_dependence": presentation_dependence,
                "mhc_class": "I",
            },
            Kind.antigen_processing: {
                "mhc_dependence": "none",
                "mhc_class": "none",
            },
        }


class MHCflurry_Affinity(BasePredictor):
    """
    MHCflurry predictor using the older Class1AffinityPredictor API.

    Only produces ``pMHC_affinity`` predictions. Use this if you only need
    affinity scores and don't want the presentation model overhead.

    See https://github.com/openvax/mhcflurry
    """

    @classmethod
    def fetch(cls, version=None):
        """Fetch MHCflurry affinity models using its native manager."""
        from .artifacts import fetch
        return fetch("mhcflurry-affinity", version=version)

    def __init__(
            self,
            alleles,
            default_peptide_lengths=[9],
            predictor=None,
            models_path=None,
            include_affinity_percentile_ranks=True,
            predictor_version=None):
        """
        Parameters
        -----------
        alleles : list of str

        default_peptide_lengths : list of int

        predictor : mhcflurry.Class1AffinityPredictor (optional)
            MHCflurry affinity predictor to use

        models_path : string
            Models dir to use if predictor argument is None

        predictor_version : str, optional
            Explicit model provenance. Official default models infer package
            plus model release at load time. Custom/injected weights otherwise
            remain unversioned and require a label before caching predictions.

        include_affinity_percentile_ranks : bool
            Whether to request affinity percentile ranks. Enabled by default.
            If enabled, requested alleles must have MHCflurry affinity
            percentile-rank calibration, either directly or through an allele
            with the same pseudosequence.
        """
        BasePredictor.__init__(
            self,
            alleles=alleles,
            default_peptide_lengths=default_peptide_lengths,
            min_peptide_length=8,
            max_peptide_length=15)
        self.predictor, self.predictor_version = _load_predictor(
            "affinity", predictor, models_path, predictor_version)

        self.include_affinity_percentile_ranks = \
            include_affinity_percentile_ranks

        for allele in self.alleles:
            if allele not in self.predictor.supported_alleles:
                raise UnsupportedAllele(allele)
        if self.include_affinity_percentile_ranks:
            _check_affinity_percent_rank_support(self.predictor, self.alleles)

    def predict_peptides(self, peptides):
        """
        Predict MHC affinity for peptides.
        """
        peptide_list = list(peptides)
        allele_list = list(self.alleles)

        batch_peptides = peptide_list * len(allele_list)
        batch_alleles = [a for a in allele_list for _ in peptide_list]

        df = self.predictor.predict_to_dataframe(
            peptides=batch_peptides,
            alleles=batch_alleles,
            include_percentile_ranks=self.include_affinity_percentile_ranks,
        )
        binding_predictions = []
        for row in df.itertuples(index=False):
            binding_predictions.append(BindingPrediction(
                allele=row.allele,
                peptide=row.peptide,
                affinity=row.prediction,
                percentile_rank=(
                    row.prediction_percentile
                    if hasattr(row, 'prediction_percentile') else None),
                prediction_method_name="mhcflurry",
            ))
        return BindingPredictionCollection(binding_predictions)
