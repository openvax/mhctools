"""Optional WEIRDO category scoring with an explicit, allele-free kind."""

from .pred import Kind, MeasurementContext, PeptideResult, Prediction
from .wrapper_base import AlleleFreePredictor


class WEIRDO(AlleleFreePredictor):
    """Score peptides with WEIRDO's verified, bundled category model.

    Scores are descriptive human-category complements. They have no confirmed
    immune-response ranking benefit and must be explicitly opted into.
    WEIRDO owns validation, normalization, score policy and artifact identity.
    Importing mhctools does not import the optional weirdo package.

    Parameters
    ----------
    model : str, optional
        Immutable WEIRDO model identifier; None selects its released default.
    models_path : path-like, optional
        Explicit identical copy of the released artifact. No downloads occur.
    score_policy : str, optional
        Named WEIRDO policy; None selects the released human-reference policy.
    host_taxon_id : int, optional
        Reference host. None uses the policy default; currently only human is
        supported. Unsupported hosts fail rather than reusing human scores.
    batch_size : int
        Maximum distinct peptides evaluated in one WEIRDO batch.
    default_peptide_lengths : int or iterable of int, optional
        Default lengths within the model's evaluated 8–11-residue range.
    """

    prediction_method_name = "weirdo"

    def __init__(self, model=None, models_path=None, score_policy=None,
                 host_taxon_id=None, batch_size=1024,
                 default_peptide_lengths=None):
        try:
            from weirdo import load_pretrained_model
        except ImportError as error:
            raise ImportError(
                "WEIRDO requires its released model API; install mhctools[weirdo]"
            ) from error
        options = {"model_dir": models_path}
        for name, value in (("model", model), ("score_policy", score_policy),
                            ("host_taxon_id", host_taxon_id)):
            if value is not None:
                options[name] = value
        self._model = load_pretrained_model(**options)
        self._model.predict_dataframe([], batch_size=batch_size)
        self.batch_size = batch_size
        contract = self._model.manifest["input"]
        self.min_peptide_length = contract["min_length"]
        self.max_peptide_length = contract["max_length"]
        supported = range(self.min_peptide_length, self.max_peptide_length + 1)
        if default_peptide_lengths is None:
            default_peptide_lengths = supported
        elif type(default_peptide_lengths) is int:
            default_peptide_lengths = [default_peptide_lengths]
        self.default_peptide_lengths = tuple(default_peptide_lengths)
        if not self.default_peptide_lengths or any(
            type(length) is not int or length not in supported
            for length in self.default_peptide_lengths
        ):
            raise ValueError("default_peptide_lengths must be within the released model's supported range")

    @property
    def predictor_version(self):
        """Exact adapter, package, artifact, policy and host identity for caches."""
        from . import __version__
        return "mhctools:%s:%s" % (__version__, self._model.prediction_identity)

    def __str__(self):
        return "WEIRDO(%s)" % self.predictor_version

    def _default_pred_kind(self):
        return Kind.foreignness

    def predict(self, peptides):
        """Return one allele-free result per input, including duplicates.

        Parameters
        ----------
        peptides : str or iterable of str
            Canonical 8–11-mers; lowercase is normalized by WEIRDO. Invalid
            input raises ValueError. Empty input returns an empty list.

        Returns
        -------
        list of PeptideResult
            ``Kind.foreignness`` scores with empty allele, complete predictor
            identity and explicit policy semantics. Higher scores mean a lower
            human-category output, not a higher immune-response probability.
        """
        frame = self._model.predict_dataframe(peptides, batch_size=self.batch_size)
        context = MeasurementContext(
            estimate_type="ml_predicted",
            score_semantics=self._model.score_policy,
            detail="Experimental category summary; reference host taxon %s; no confirmed immune-response ranking benefit"
            % self._model.host_taxon_id,
        )
        return [PeptideResult(preds=(Prediction(
            kind=Kind.foreignness, score=float(row.foreignness), peptide=row.peptide,
            predictor_name=self.prediction_method_name,
            predictor_version=self.predictor_version,
            measurement_context=context,
        ),)) for row in frame.itertuples(index=False)]
