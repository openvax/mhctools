"""Native ITCell cathepsin specificity profiles, without structural inference.

The independently implemented scoring formula follows the released author
script. Bundled, unmodified count matrices and background frequencies are
LGPL-2.1 data; attribution, source hashes and license accompany the data.
"""

from functools import lru_cache
import hashlib
import json
import math
from pathlib import Path

from .cleavage import CleavageModel, CleavageResult, CleavageSite, coerce_peptide


ITCELL_MODELS = {
    "itcell-cat" + enzyme.lower() + "-" + str(minutes): dict(enzyme=enzyme, minutes=minutes)
    for enzyme in ("B", "S", "H") for minutes in (15, 60, 240)
}
UNIPROT = {"B": "P07858", "S": "P25774", "H": "P09668"}
DATA_SHA256 = "da0ec86118f2039eaed4e6023ec0c320c0712a059f15304c472a8b830154aa32"


@lru_cache(maxsize=1)
def _data():
    raw = (Path(__file__).parent / "data/itcell_profiles.json").read_bytes()
    if hashlib.sha256(raw).hexdigest() != DATA_SHA256:
        raise ValueError("ITCell profile data checksum mismatch")
    return json.loads(raw)


def _weights(profile, background):
    """Use the author's shared maximum-column denominator and +1 counts."""
    raw = profile["raw_counts"]
    if hashlib.sha256(raw.encode()).hexdigest() != profile["sha256"]:
        raise ValueError("ITCell raw matrix checksum mismatch")
    counts = {fields[0]: fields[1:] for fields in
              (line.split() for line in raw.splitlines())}
    if set(counts) != set(background) or any(len(v) != 8 for v in counts.values()):
        raise ValueError("Expected all 20 residues and eight positions")
    totals = [sum(float(v[pos]) for v in counts.values() if v[pos] != "X")
              for pos in range(8)]
    denominator = max(totals) + 20
    return {aa: tuple(0.0 if value == "X" else
                     math.log2(((float(value) + 1) / denominator) / background[aa])
                     for value in values) for aa, values in counts.items()}


class ITCellCleavage:
    """Score one released cathepsin B/S/H profile on canonical sequence.

    Parameters
    ----------
    enzyme : str
        ``B`` and ``S`` assess endoprotease sites. ``H`` assesses only the
        initial N-terminal mono-aminopeptidase step.
    minutes : int
        Select the source MSP-MS count profile (15, 60 or 240 minutes).
        This selects training observations, not a requested incubation time.

    Notes
    -----
    Native log2 enrichment scores are not cleavage probabilities, depletion
    percentages or serum half-life. Unavailable flanks contribute zero as in
    the author script. Cathepsin B carboxydipeptidase activity is not modeled.
    No subsequent fragments or repeated trimming are inferred.
    """

    def __init__(self, enzyme="S", minutes=240):
        self.model = self.catalog_model(enzyme, minutes)
        self.enzyme, self.minutes = enzyme, minutes
        self.topology = "n_terminal" if enzyme == "H" else "internal"
        self.removed = 1 if enzyme == "H" else None
        data = _data()
        self.profile = data["profiles"]["cat%s_%d_count.txt" % (enzyme, minutes)]
        self.weights = _weights(self.profile, data["background_frequencies"])
        self.threshold = 2.0 if enzyme == "H" else 3.0

    @staticmethod
    def catalog_model(enzyme="S", minutes=240):
        """Return source-scoped metadata without invoking an external tool."""
        if enzyme not in UNIPROT or type(minutes) is not int or minutes not in (15, 60, 240):
            raise ValueError("ITCell requires enzyme B/S/H and profile 15/60/240")
        profile = _data()["profiles"]["cat%s_%d_count.txt" % (enzyme, minutes)]
        return CleavageModel(
            name="itcell-cat%s-%d" % (enzyme.lower(), minutes),
            version="matrix-sha256:" + profile["sha256"],
            enzyme="Cathepsin " + enzyme, uniprot=UNIPROT[enzyme], species="Homo sapiens",
            compartments=("endosome",), evidence="quantitative_model",
            references=("https://pmc.ncbi.nlm.nih.gov/articles/PMC6219782/",
                        "https://github.com/salilab/itcell-lib", "https://doi.org/10.5281/zenodo.3227044"),
            assay=("ITCell recombinant human cathepsin MSP-MS specificity; 228 tetradecapeptides, "
                   "0.2 ug/mL enzyme, pH 6.5, 1 mM TCEP; %d-minute count profile" % minutes),
            limitations=(
                "Sequence specificity only; not calibrated probability, degradation rate or serum survival. "
                "Source uses norleucine in place of methionine; cysteine is unrepresented. "
                "Missing sequence flanks contribute zero. Human background frequencies retained as published. "
                "Profile time is training-assay scope, not a user incubation-time prediction. "
                + ("Only initial N-terminal aminopeptidase step; no internal cuts or cascade." if enzyme == "H"
                   else "Endoprotease sites only; cathepsin B C-terminal trimming is outside scope.")),
            score_name="itcell_log2_enrichment", score_units="sum of log2 profile/background ratios",
            scored_endpoint="site_cleavage")

    def predict(self, peptide):
        """Retain every assessed native score, including below-threshold sites."""
        peptide = coerce_peptide(peptide)
        if peptide.n_term != "free" or peptide.c_term != "free":
            return CleavageResult(peptide, self.model,
                                  unsupported_reason="ITCell requires canonical sequence with assumed free termini")
        if len(peptide.sequence) < 2:
            return CleavageResult(peptide, self.model, unsupported_reason="No internal peptide bond")
        padded = "XXX" + peptide.sequence + "XXX"
        bonds = (1,) if self.enzyme == "H" else range(1, len(peptide.sequence))
        sites = []
        for bond in bonds:
            context = padded[bond - 1:bond + 7]
            score = sum(self.weights[aa][pos] for pos, aa in enumerate(context) if aa != "X")
            reason = "Native ITCell P4-P4prime profile; context " + context
            if "M" in context:
                reason += "; M uses source norleucine surrogate"
            if "C" in context:
                reason += "; C unrepresented in source assay"
            if "X" in context:
                reason += "; absent flanks contribute zero"
            sites.append(CleavageSite(bond, "scored", reason, score))
        return CleavageResult(peptide, self.model, tuple(sites), conditions=(
            ("matrix_sha256", self.profile["sha256"]), ("bundled_data_sha256", DATA_SHA256),
            ("profile_minutes", str(self.minutes)), ("pH", "6.5"),
            ("candidate_threshold_strict_gt", str(self.threshold)),
            ("topology", self.topology),
            ("input_scope", "sequence specificity; no structural or exposure inference")))
