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

"""The Calis (IEDB) class-I T-cell immunogenicity model.

Calis et al. 2013 is the classic sequence-only MHC-I immunogenicity predictor:
a fixed per-amino-acid log-enrichment scale, weighted by per-position
importance, with the anchor positions (P1, P2, and the C-terminus) masked out
because their identity is driven by MHC binding rather than TCR recognition.

The ~30 parameters come from the open-access paper (Calis et al.,
PLoS Comput. Biol. 2013, 10.1371/journal.pcbi.1003266) and are reproduced here,
so the predictor is self-contained. It is allele-independent: it always masks
P1/P2/C-terminus, the default the IEDB tool uses when no allele is given. A
score > 0 leans immunogenic, < 0 non-immunogenic.

Like every current CD8 immunogenicity predictor it ranks better than it
generalizes (independent benchmarks put the field near AUC 0.5-0.65 on unseen
neoepitopes) — a prioritization aid, not ground truth.
"""

from .pred import Kind, PeptideResult, Prediction
from .wrapper_base import AlleleFreePredictor

# Per-amino-acid log-enrichment score, log(freq(aa | immunogenic) /
# freq(aa | non-immunogenic)) from the Calis 2013 training set. Higher = more
# associated with immunogenicity (W most, K/M least).
IMMUNOSCALE = {
    "A": 0.127, "C": -0.175, "D": 0.072, "E": 0.325, "F": 0.380,
    "G": 0.110, "H": 0.105, "I": 0.432, "K": -0.700, "L": -0.036,
    "M": -0.570, "N": -0.021, "P": -0.036, "Q": -0.376, "R": 0.168,
    "S": -0.537, "T": 0.126, "V": 0.134, "W": 0.719, "Y": -0.012,
}

# Per-position importance weights for a 9-mer (KL divergence between the
# immunogenic and non-immunogenic training peptides at each position). Index 0
# is P1. The interior positions (P3-P8) carry the signal; P1/P2/P9 are ~0 and
# are also masked below.
IMMUNOWEIGHT = [0.00, 0.00, 0.10, 0.31, 0.30, 0.29, 0.26, 0.18, 0.00]

# Amino acids the scale is defined for.
_VALID_AMINO_ACIDS = frozenset(IMMUNOSCALE)


def position_weights(peptide_length):
    """Per-position importance weights for a peptide of the given length.

    For 9-mers this is :data:`IMMUNOWEIGHT` verbatim. For longer peptides the
    Calis rule pads the interior (after P5) with ``0.30`` weights, matching the
    IEDB tool. For shorter peptides it returns ``IMMUNOWEIGHT`` unchanged and the
    caller uses only the leading ``peptide_length`` entries.
    """
    if peptide_length > 9:
        return (
            IMMUNOWEIGHT[:5]
            + [0.30] * (peptide_length - 9)
            + IMMUNOWEIGHT[5:])
    return IMMUNOWEIGHT


def immunogenicity_score(peptide):
    """Calis immunogenicity score for one peptide (default P1/P2/C-term mask).

    Parameters
    ----------
    peptide : str
        Amino-acid sequence (upper-case, standard residues).

    Returns
    -------
    float
        Sum over unmasked positions of ``position_weight * immunoscale``,
        rounded to 5 decimals (as the IEDB tool does). > 0 leans immunogenic.
    """
    peptide = peptide.upper()
    n = len(peptide)
    # Default anchor mask (0-indexed): P1, P2, and the C-terminal residue.
    masked = {0, 1, n - 1}
    weights = position_weights(n)
    score = 0.0
    for i, aa in enumerate(peptide):
        if i in masked:
            continue
        score += weights[i] * IMMUNOSCALE[aa]
    return round(score, 5)


class Calis(AlleleFreePredictor):
    """The Calis (IEDB) class-I immunogenicity predictor.

    Allele-independent: ``predict()`` returns one
    :class:`~mhctools.pred.Prediction` per peptide with an empty ``allele`` and
    ``kind == Kind.immunogenicity``.

    Parameters
    ----------
    min_peptide_length : int
        Shortest peptide accepted. Default 8 (class-I minimum).
    max_peptide_length : int
        Longest peptide accepted. Default 11 (the class-I range the model was
        designed for; longer peptides are handled by the interior-padding rule
        but are outside its intended scope).
    """

    mhc_class = "I"

    def __init__(self, min_peptide_length=8, max_peptide_length=11):
        if min_peptide_length < 3:
            # With P1/P2/C-term masked, peptides shorter than this have no
            # scoring positions left.
            raise ValueError("min_peptide_length must be >= 3")
        if max_peptide_length < min_peptide_length:
            raise ValueError(
                "max_peptide_length (%d) < min_peptide_length (%d)"
                % (max_peptide_length, min_peptide_length))
        self.min_peptide_length = min_peptide_length
        self.max_peptide_length = max_peptide_length

    def __str__(self):
        return "Calis(min_peptide_length=%d, max_peptide_length=%d)" % (
            self.min_peptide_length, self.max_peptide_length)

    def _default_pred_kind(self):
        return Kind.immunogenicity

    def _check_peptides(self, peptides):
        for peptide in peptides:
            n = len(peptide)
            if n < self.min_peptide_length or n > self.max_peptide_length:
                raise ValueError(
                    "Calis supports peptides of length %d-%d; got %r (length %d)"
                    % (self.min_peptide_length, self.max_peptide_length,
                       peptide, n))
            invalid = set(peptide) - _VALID_AMINO_ACIDS
            if invalid:
                raise ValueError(
                    "Peptide %r contains non-standard amino acid(s): %s"
                    % (peptide, "".join(sorted(invalid))))

    def predict(self, peptides):
        """Predict immunogenicity for a list of peptides.

        Returns
        -------
        list of PeptideResult
            One entry per input peptide, each holding a single
            ``Kind.immunogenicity`` prediction (empty ``allele``; ``score`` is
            the Calis score, higher = more immunogenic).
        """
        peptide_list = self._normalize_peptides(peptides)
        self._check_peptides(peptide_list)
        return [
            PeptideResult(preds=(Prediction(
                kind=Kind.immunogenicity,
                score=immunogenicity_score(peptide),
                peptide=peptide,
                predictor_name="calis"),))
            for peptide in peptide_list
        ]
