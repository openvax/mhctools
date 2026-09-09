"""Curated human peptidase recognition rules and model selection.

Rules flag partial substrate-recognition patterns, not enzyme activity,
probabilities or whole-peptide stability. Every rule retains its scope.
"""

from dataclasses import dataclass
import re
from typing import Optional, Tuple

from .cleavage import CleavageInput, CleavageModel, CleavageResult, CleavageSite
from .dpp4 import DPP4qPISA


@dataclass(frozen=True)
class PeptidaseMotif:
    """A partial specificity rule with explicit terminal or internal topology."""

    model: CleavageModel
    topology: str
    removed: int
    left: str
    right: str
    description: str
    min_length: int = 2
    max_length: Optional[int] = None
    n_chemistry: Tuple[str, ...] = ("free",)

    def __post_init__(self):
        if self.topology not in ("n_terminal", "c_terminal", "internal"):
            raise ValueError("Unknown cleavage topology")
        if (self.topology != "internal" and
                (self.removed < 1 or self.min_length <= self.removed)):
            raise ValueError("Terminal rule must leave a peptide bond to assess")
        re.compile(self.left)
        re.compile(self.right)

    def predict(self, peptide):
        """Assess eligible bonds; a non-match is not evidence of resistance."""
        if isinstance(peptide, str):
            peptide = CleavageInput(peptide)
        if not isinstance(peptide, CleavageInput):
            raise TypeError("Expected CleavageInput or canonical peptide string")
        seq = peptide.sequence
        reason = None
        if peptide.n_term not in self.n_chemistry or peptide.c_term != "free":
            reason = "Terminal chemistry is outside this rule's documented input domain"
        elif len(seq) < self.min_length:
            reason = "Rule requires at least %d residues" % self.min_length
        elif self.max_length is not None and len(seq) > self.max_length:
            reason = "Rule is limited to peptides of at most %d residues" % self.max_length
        if reason:
            return CleavageResult(peptide, self.model, unsupported_reason=reason)
        if self.topology == "n_terminal":
            bonds = (self.removed,)
        elif self.topology == "c_terminal":
            bonds = (len(seq) - self.removed,)
        else:
            bonds = range(1, len(seq))
        sites = []
        for bond in bonds:
            match = (re.search(self.left + "$", seq[:bond]) is not None and
                     re.match(self.right, seq[bond:]) is not None)
            sites.append(CleavageSite(
                bond, "matched" if match else "not_matched", self.description))
        return CleavageResult(peptide, self.model, tuple(sites))


def _rule(name, enzyme, accession, compartments, papers, topology, removed,
          left, right, description, assay, limitations, **kwargs):
    metadata = CleavageModel(
        name=name, version="1", enzyme=enzyme, uniprot=accession,
        species="Homo sapiens", compartments=compartments, evidence="motif_rule",
        references=tuple(papers) + (
            "https://www.uniprot.org/uniprotkb/%s/entry" % accession,),
        assay=assay,
        limitations=limitations + " Partial specificity rule; no kinetic or matrix calibration.")
    return PeptidaseMotif(metadata, topology, removed, left, right, description, **kwargs)


_RULES = (
    _rule("cpn-basic", "CPN1", "P15169", ("plasma", "serum", "extracellular"),
          ("https://doi.org/10.1016/0003-9861(75)90104-6",),
          "c_terminal", 1, ".", "[KR]", "Free C-terminal Lys/Arg removal",
          "Human plasma CPN peptide-substrate assays",
          "Penultimate residues affect rates; Lys and Arg are not kinetically equivalent."),
    _rule("app2-xp", "XPNPEP2", "O43895", ("plasma", "extracellular"),
          ("https://pubmed.ncbi.nlm.nih.gov/15361070/",),
          "n_terminal", 1, ".", "P", "Exposed N-terminal X|Pro",
          "Recombinant human membrane aminopeptidase P and kinin substrates",
          "Removes the FIRST residue; this is not DPP-like dipeptide removal."),
    _rule("fap-dipeptidyl", "FAP", "Q12884", ("extracellular", "plasma"),
          ("https://pubmed.ncbi.nlm.nih.gov/16410248/",),
          "n_terminal", 2, ".P", "[^P]", "Exposed N-terminal X-Pro|non-Pro",
          "Purified human FAP dipeptide substrate profiling",
          "Dipeptidyl activity only; P2 preferences and whole-peptide context omitted.",
          min_length=3),
    _rule("fap-endo-gp", "FAP", "Q12884", ("extracellular", "plasma"),
          ("https://pubmed.ncbi.nlm.nih.gov/16480718/",
           "https://pubmed.ncbi.nlm.nih.gov/16410248/",
           "https://pubmed.ncbi.nlm.nih.gov/22750443/"),
          "internal", 0, "GP", "[^P]", "Gly-Pro|non-Pro endopeptidase recognition",
          "Human FAP substrate profiling around alpha-2-antiplasmin cleavage site",
          "P3 preferences and accessibility omitted; independent of dipeptidyl activity.",
          min_length=3, n_chemistry=("free", "acetylated")),
    _rule("enpep-acidic", "ENPEP", "Q07075", ("extracellular",),
          ("https://pubmed.ncbi.nlm.nih.gov/23888046/",),
          "n_terminal", 1, "[DE]", ".", "Exposed N-terminal Asp/Glu|X",
          "Human aminopeptidase A structural/substrate specificity experiments",
          "Calcium-dependent preference; enzyme exposure and calcium are not modeled."),
    _rule("anpep-ala", "ANPEP", "P15144", ("extracellular",),
          ("https://pubmed.ncbi.nlm.nih.gov/22932899/",),
          "n_terminal", 1, "A", ".", "Preferred N-terminal Ala|X",
          "Human aminopeptidase N structural/biochemical specificity",
          "Only an Ala preference flag. Broad alternative substrates and X-Pro dipeptide removal omitted."),
    _rule("dpp8-xp-xa", "DPP8", "Q6V1X1", ("cytosol",),
          ("https://pubmed.ncbi.nlm.nih.gov/11012666/",
           "https://pmc.ncbi.nlm.nih.gov/articles/PMC3656252/"),
          "n_terminal", 2, ".[PA]", "[^P]", "Exposed X-(Pro/Ala)|non-Pro",
          "Human DPP8 biochemical characterization and substrate degradomics",
          "Qualitative rule only; DPP4 and C. elegans DPF-3 coefficients are inapplicable.",
          min_length=3),
    _rule("dpp9-xp-xa", "DPP9", "Q86TI2", ("cytosol",),
          ("https://pubmed.ncbi.nlm.nih.gov/19667070/",
           "https://pmc.ncbi.nlm.nih.gov/articles/PMC3656252/"),
          "n_terminal", 2, ".[PA]", "[^P]", "Exposed X-(Pro/Ala)|non-Pro",
          "Human DPP9 peptide/antigen processing assays and substrate degradomics",
          "Shared recognition pattern with DPP8 does not imply identical substrate rates.",
          min_length=3),
    _rule("prep-pro", "PREP", "P48147", ("cytosol", "extracellular", "serum"),
          ("https://pubmed.ncbi.nlm.nih.gov/22750443/",),
          "internal", 0, ".P", ".", "Internal post-proline X-Pro|X recognition",
          "Recombinant human POP/FAP peptide substrate profiling",
          "Conservative 4-30-residue model domain; not an absolute enzyme size cutoff. "
          "Weaker non-Pro cleavage and flanking charge preferences omitted; abundance unmodeled.",
          min_length=4, max_length=30),
    _rule("erap2-basic", "ERAP2", "Q6P179", ("er",),
          ("https://pubmed.ncbi.nlm.nih.gov/12799365/",
           "https://pubmed.ncbi.nlm.nih.gov/26381406/"),
          "n_terminal", 1, "[RK]", ".", "Preferred exposed N-terminal Arg/Lys|X",
          "Human ERAP2 biochemical specificity and peptide-complex structures",
          "A basic-residue preference flag only; length, allotype and internal sequence affect trimming."),
)


def cleavage_models(include_optional=False):
    """List built-in model metadata without loading optional runtimes."""
    models = (DPP4qPISA.model,) + tuple(rule.model for rule in _RULES)
    if include_optional:
        from .eramer_cleavage import ERAMERCleavage
        models += (ERAMERCleavage.model,)
    return models


def get_cleavage_model(name):
    """Select an exact model name; never silently substitute another enzyme."""
    if name == "dpp4-qpisa":
        return DPP4qPISA()
    if name == "eramer-step":
        from .eramer_cleavage import ERAMERCleavage
        return ERAMERCleavage()
    for rule in _RULES:
        if rule.model.name == name:
            return rule
    raise ValueError("Unknown cleavage model %r; choices: %s" % (
        name, ", ".join(model.name for model in cleavage_models(include_optional=True))))


def predict_cleavage(peptide, models=None, compartment=None):
    """Return distinct model results without aggregation or enzyme ranking.

    ``compartment`` filters annotated enzyme locations, not assay validation.
    Explicit selections incompatible with that filter are rejected.
    """
    if compartment is not None and compartment not in {
            c for model in cleavage_models() for c in model.compartments}:
        raise ValueError("Unknown compartment %r" % compartment)
    if models is None:
        models = [m.name for m in cleavage_models()
                  if compartment is None or compartment in m.compartments]
    elif isinstance(models, str):
        models = [models]
    results = []
    for name in dict.fromkeys(models):
        predictor = get_cleavage_model(name)
        if compartment is not None and compartment not in predictor.model.compartments:
            raise ValueError("Model %s is not annotated for compartment %s" % (name, compartment))
        results.append(predictor.predict(peptide))
    return tuple(results)
