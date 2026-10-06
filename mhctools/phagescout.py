"""Independent native PhageScout sequence scoring with licensed source data.

No author notebook code is included. Unmodified deposited matrices and
peptide profiles are CC BY 4.0; their attribution and license accompany them.
"""

import csv
from functools import lru_cache
import hashlib
from io import StringIO
import json
import math
from pathlib import Path
from types import MappingProxyType

from .cleavage import AMINO_ACIDS, CleavageModel, CleavageResult, CleavageSite, coerce_peptide


DATA_SHA256 = "d9feb51e0ee12dafaad92aaea2d4c2c28a474102129b9276348e14b1a6160f9b"
ENZYMES = {
    "ELANE": ("elastase", "Neutrophil elastase", "P08246", "5 nM, 1.5 hours"),
    "CTSG": ("cathepsin G", "Cathepsin G", "P08311", "10 nM, 6 hours"),
}
PROFILES = (
    "pwm-deseq2", "pwm-relaxed-unaligned", "pwm-relaxed-aligned",
    "peptide-relaxed-unaligned", "peptide-relaxed-aligned",
)
OPTIONAL_PROFILES = ("peptide-deseq2",)
PHAGESCOUT_MODELS = {
    "phagescout-%s-%s" % (enzyme.lower(), profile): dict(enzyme=enzyme, profile=profile)
    for enzyme in ENZYMES for profile in PROFILES
}
PHAGESCOUT_OPTIONAL_MODELS = {
    "phagescout-%s-%s" % (enzyme.lower(), profile): dict(enzyme=enzyme, profile=profile)
    for enzyme in ENZYMES for profile in OPTIONAL_PROFILES
}


@lru_cache(maxsize=1)
def _data():
    raw = (Path(__file__).parent / "data/phagescout_profiles.json").read_bytes()
    if hashlib.sha256(raw).hexdigest() != DATA_SHA256:
        raise ValueError("PhageScout bundled data checksum mismatch")
    return json.loads(raw)


def _asset(enzyme, profile):
    source = ENZYMES[enzyme][0]
    suffix = profile.split("-", 1)[1].replace("-", "_")
    if profile.startswith("pwm-"):
        suffix = "relaxed" if suffix == "relaxed_unaligned" else suffix
        name = "%s_%s_pwm.txt" % (source, suffix)
    else:
        name = "peptide_profile_%s_%s.txt" % (source, suffix)
    asset = _data()["profiles"][name]
    if hashlib.sha256(asset["raw"].encode()).hexdigest() != asset["sha256"]:
        raise ValueError("PhageScout source asset checksum mismatch: " + name)
    return name, asset


@lru_cache(maxsize=10)
def _weights(enzyme, profile):
    _, asset = _asset(enzyme, profile)
    if profile.startswith("pwm-"):
        rows = list(csv.reader(StringIO(asset["raw"]), delimiter="\t"))
        residues = rows.pop(0)
        width = 9 if profile.endswith("-aligned") else 5
        if (len(rows) != width or not AMINO_ACIDS <= set(residues) or
                len(set(residues)) != len(residues) or
                any(len(row) != len(residues) for row in rows)):
            raise ValueError("Invalid PhageScout matrix shape")
        values = [[float(value) for value in row] for row in rows]
        if not all(math.isfinite(value) for row in values for value in row):
            raise ValueError("Nonfinite PhageScout matrix weight")
        return MappingProxyType({aa: tuple(row[i] for row in values) for i, aa in enumerate(residues)})
    return _peptide_weights(StringIO(asset["raw"]), 9 if profile.endswith("-aligned") else 5)


def _peptide_weights(source, width):
    rows = csv.DictReader(source, delimiter="\t")
    if rows.fieldnames != ["peptide", "log2FoldChange"]:
        raise ValueError("Invalid PhageScout peptide profile header")
    groups = {}
    for order, row in enumerate(rows):
        pattern, score = row["peptide"], float(row["log2FoldChange"])
        if (len(pattern) != width or not set(pattern) <= AMINO_ACIDS | {"-"} or
                not math.isfinite(score)):
            raise ValueError("Invalid PhageScout peptide profile")
        positions = tuple(i for i, aa in enumerate(pattern) if aa != "-")
        key = "".join(pattern[i] for i in positions)
        # First released row wins, including ties across wildcard masks.
        groups.setdefault(positions, {}).setdefault(key, (order, score))
    return MappingProxyType({positions: MappingProxyType(lookup) for positions, lookup in groups.items()})


@lru_cache(maxsize=2)
def _full_weights(path):
    # Only explicit optional model construction reaches these large tables.
    with path.open() as source:
        return _peptide_weights(source, 5)


class PhageScout:
    """Score released human ELANE/CTSG sequence-specificity profiles.

    Parameters
    ----------
    enzyme : str
        Human ``ELANE`` (neutrophil elastase) or ``CTSG`` (cathepsin G).
    profile : str
        Exact entry in :data:`PROFILES` or :data:`OPTIONAL_PROFILES`.
        ``peptide-deseq2`` requires ``mhctools fetch phagescout``.
        PWM models assess every bond with
        available context. Peptide models assess only released-profile
        matches; unmatched bonds are unassessed, never assigned zero.
    profile_dir : str or pathlib.Path, optional
        Full DESeq2 table directory; overrides ``PHAGESCOUT_HOME`` and the
        managed data directory. Applies only to ``peptide-deseq2``.

    Notes
    -----
    Scores prioritize sequence contexts under active-enzyme exposure. The
    aligned P1 anchor is inferred from phage alignment, not a measured cut.
    No structural classifier, normalized score, threshold or survival model
    is substituted for the native sequence score.
    """

    def __init__(self, enzyme="ELANE", profile="pwm-relaxed-aligned", profile_dir=None):
        self.model = self.catalog_model(enzyme, profile)
        self.enzyme, self.profile = enzyme, profile
        self.aligned = profile.endswith("-aligned")
        self.is_pwm = profile.startswith("pwm-")
        self.bundled = profile not in OPTIONAL_PROFILES
        if self.bundled:
            if profile_dir is not None:
                raise ValueError("profile_dir applies only to the optional peptide-deseq2 profile")
            self.weights = _weights(enzyme, profile)
            self.asset_name, asset = _asset(enzyme, profile)
            self.asset_sha256 = asset["sha256"]
        else:
            from .phagescout_artifacts import selected_profile
            path, self.asset_sha256 = selected_profile(enzyme, profile_dir)
            self.asset_name = path.name
            self.weights = _full_weights(path)

    @staticmethod
    def catalog_model(enzyme="ELANE", profile="pwm-relaxed-aligned"):
        """Describe an exact source model without external inference."""
        if enzyme not in ENZYMES or profile not in PROFILES + OPTIONAL_PROFILES:
            raise ValueError("PhageScout requires enzyme ELANE/CTSG and an exact released profile")
        source, label, uniprot, assay = ENZYMES[enzyme]
        if profile in OPTIONAL_PROFILES:
            from .phagescout_artifacts import asset as optional_asset
            digest = optional_asset(enzyme)[2]
        else:
            _, asset = _asset(enzyme, profile)
            digest = asset["sha256"]
        aligned, pwm = profile.endswith("-aligned"), profile.startswith("pwm-")
        units = ("sum of positional log2 enrichment weights" if aligned else
                 "mean of five-mer sums of positional log2 enrichment weights")
        if profile == "pwm-deseq2":
            units = "mean of five-mer sums of DESeq2-derived positional weights"
        elif not pwm:
            units = ("first matching aligned phage log2 fold change" if aligned else
                     "mean matched five-mer phage log2 fold change")
        return CleavageModel(
            name="phagescout-%s-%s" % (enzyme.lower(), profile),
            version="asset-sha256:" + digest, enzyme=label, uniprot=uniprot,
            species="Homo sapiens", compartments=("extracellular",),
            evidence="quantitative_model",
            references=("https://doi.org/10.3390/ijms27177593",
                        "https://zenodo.org/records/21387981"),
            assay="PhageScout randomized five-mer phage-display library; active %s, %s; "
                  "FLAG-eluted control; source pH not established here" % (source, assay),
            limitations=(
                "Phage-derived recognition feature, not calibrated probability, loss, rate or half-life. "
                "Enzyme exposure, activation, inhibitors and substrate accessibility are not inferred. "
                "Nine-mer P5-P4prime alignment assigns an inferred P1 anchor; phage enrichment does not "
                "experimentally identify the scissile bond. Canonical sequence only; terminal chemistry "
                "effects and phage versus free-peptide transfer are unvalidated. "
                + ("Missing PWM flanks are omitted as in source scoring." if pwm and aligned else
                   "Complete five-mer windows spanning each bond are averaged." if not aligned else
                   "Requires full nine-mer context; first source-ordered wildcard match wins.")
                + (" Unmatched profiles are unassessed, not resistant." if not pwm else "")),
            score_name="phagescout_" + profile.replace("-", "_"), score_units=units,
            scored_endpoint="site_cleavage")

    def _score(self, sequence, bond):
        if self.aligned:
            start, end = max(0, bond - 5), min(len(sequence), bond + 4)
            if self.is_pwm:
                return sum(self.weights[sequence[i]][i - bond + 5] for i in range(start, end))
            if end - start != 9:
                return None
            context = sequence[start:end]
            matches = [match for positions, lookup in self.weights.items()
                       if (match := lookup.get("".join(context[i] for i in positions))) is not None]
            return min(matches)[1] if matches else None
        scores = []
        for start in range(max(0, bond - 4), min(bond, len(sequence) - 4)):
            context = sequence[start:start + 5]
            if self.is_pwm:
                scores.append(sum(self.weights[aa][i] for i, aa in enumerate(context)))
            else:
                match = self.weights[(0, 1, 2, 3, 4)].get(context)
                if match is not None:
                    scores.append(match[1])
        return sum(scores) / len(scores) if scores else None

    def predict(self, peptide):
        """Preserve native scores, canonical bonds and explicit coverage gaps."""
        peptide = coerce_peptide(peptide)
        if peptide.n_term != "free" or peptide.c_term != "free":
            return CleavageResult(peptide, self.model, unsupported_reason=
                                  "PhageScout terminal chemistry effects are unestablished; assumed free termini required")
        if len(peptide.sequence) < 2:
            return CleavageResult(peptide, self.model, unsupported_reason="No internal peptide bond")
        sites, missing = [], {}
        for bond in range(1, len(peptide.sequence)):
            score = self._score(peptide.sequence, bond)
            if score is None:
                missing[str(bond)] = ("Requires complete nine-mer context" if self.aligned and
                                      not 5 <= bond <= len(peptide.sequence) - 4 else
                                      "No matching released peptide profile" if not self.is_pwm and
                                      len(peptide.sequence) >= 5 else
                                      "Requires a complete five-mer window")
            else:
                reason = ("P5-P4prime context; inferred P1 anchor" if self.aligned else
                          "Mean of available five-mer windows spanning the bond")
                if self.is_pwm and self.aligned and not 5 <= bond <= len(peptide.sequence) - 4:
                    reason += "; missing flanks omitted"
                sites.append(CleavageSite(bond, "scored", reason, score))
        conditions = (("asset", self.asset_name), ("asset_sha256", self.asset_sha256),
                      ("normalization", "none"),
                      ("exposure", "conditional active enzyme; not inferred"),
                      ("unassessed_bonds", json.dumps(missing, sort_keys=True)))
        if self.bundled:
            conditions += (("bundled_data_sha256", DATA_SHA256),)
        return CleavageResult(peptide, self.model, tuple(sites),
                              unsupported_reason=None if sites else "No assessable source-profile context; see unassessed_bonds",
                              conditions=conditions)
