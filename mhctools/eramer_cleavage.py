"""Expose one ERAMER PWM trimming assessment with its exact artifact identity."""

from dataclasses import replace
import hashlib
from io import BytesIO
from pathlib import Path

from .cleavage import CleavageInput, CleavageModel, CleavageResult, CleavageSite
from .eramer import _find_pwm_path, _specificity, load_pwm


class ERAMERCleavage:
    """Score the initial ERAP1 trimming step of a 9-16mer with ERAMER's PWM.

    This is the existing ERAMER intermediate-specificity score, not the
    cascade average. The external GPL-licensed PWM is loaded, never vendored.
    """

    model = CleavageModel(
        name="eramer-step", version="2024-external-PWM", enzyme="ERAP1",
        uniprot="Q9NZ08", species="Homo sapiens", compartments=("er",),
        evidence="quantitative_model",
        references=("https://pubmed.ncbi.nlm.nih.gov/38925438/",
                    "https://github.com/aalokaily/ERAMER"),
        assay="ERAMER length-specific position weight matrices for ERAP1 specificity",
        limitations=("Requires separately fetched ERAMER PWM and openpyxl. "
                     "Single intermediate score; not a probability, kinetic rate or cascade average. "
                     "No allotype, enzyme exposure or MHC protection model."),
        score_name="intermediate_pwm_specificity", score_units="native ERAMER specificity")

    def __init__(self, eramer_home=None, pwm_path=None):
        path = Path(_find_pwm_path(eramer_home, pwm_path)).expanduser().resolve()
        data = path.read_bytes()
        self._weights = load_pwm(BytesIO(data))
        self.model = replace(self.model, version="pwm-sha256:" + hashlib.sha256(data).hexdigest())

    def predict(self, peptide):
        """Assess only the initial N-terminal bond; fragments require new input."""
        if isinstance(peptide, str):
            peptide = CleavageInput(peptide)
        if not isinstance(peptide, CleavageInput):
            raise TypeError("Expected CleavageInput or canonical peptide string")
        if peptide.n_term != "free" or peptide.c_term != "free":
            return CleavageResult(peptide, self.model, unsupported_reason="ERAMER requires unmodified free termini")
        length = len(peptide.sequence)
        if not 9 <= length <= 16:
            return CleavageResult(peptide, self.model, unsupported_reason="ERAMER has PWMs only for 9-16 residues")
        score = _specificity(peptide.sequence, self._weights[length])
        return CleavageResult(peptide, self.model, (CleavageSite(
            1, "scored", "Mean length-specific PWM weight for this exposed precursor", score),))
