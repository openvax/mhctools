"""Curated human peptidase recognition rules and model selection.

Rules flag partial substrate-recognition patterns, not enzyme activity,
probabilities or whole-peptide stability. Every rule retains its scope.
"""

from dataclasses import dataclass, replace
import re
from typing import Optional, Tuple

from .cleavage import CleavageInput, CleavageModel, CleavageResult, CleavageSite
from .dpp4 import DPP4qPISA
from .substrate_reference import substrate_references


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
    c_chemistry: Tuple[str, ...] = ("free",)
    requires_activation: bool = False
    enzyme_state: str = "unknown"

    def __post_init__(self):
        if self.topology not in ("n_terminal", "c_terminal", "internal"):
            raise ValueError("Unknown cleavage topology")
        if (self.topology != "internal" and
                (self.removed < 1 or self.min_length <= self.removed)):
            raise ValueError("Terminal rule must leave a peptide bond to assess")
        re.compile(self.left)
        re.compile(self.right)
        if self.enzyme_state not in ("unknown", "active", "zymogen", "inactive"):
            raise ValueError("enzyme_state must be unknown, active, zymogen or inactive")

    def predict(self, peptide):
        """Assess eligible bonds; a non-match is not evidence of resistance."""
        if isinstance(peptide, str):
            peptide = CleavageInput(peptide)
        if not isinstance(peptide, CleavageInput):
            raise TypeError("Expected CleavageInput or canonical peptide string")
        seq = peptide.sequence
        reason = None
        conditions = (("enzyme_state", self.enzyme_state),) if self.requires_activation else ()
        if self.requires_activation and self.enzyme_state != "active":
            reason = "This rule requires explicitly active %s; supplied state is %s" % (
                self.model.enzyme, self.enzyme_state)
        elif peptide.n_term not in self.n_chemistry or peptide.c_term not in self.c_chemistry:
            reason = "Terminal chemistry is outside this rule's documented input domain"
        elif len(seq) < self.min_length:
            reason = "Rule requires at least %d residues" % self.min_length
        elif self.max_length is not None and len(seq) > self.max_length:
            reason = "Rule is limited to peptides of at most %d residues" % self.max_length
        if reason:
            return CleavageResult(peptide, self.model, unsupported_reason=reason, conditions=conditions)
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
        return CleavageResult(peptide, self.model, tuple(sites), conditions=conditions)


def _rule(name, enzyme, accession, compartments, papers, topology, removed,
          left, right, description, assay, limitations, *, strictness, basis, **kwargs):
    metadata = CleavageModel(
        name=name, version="1", enzyme=enzyme, uniprot=accession,
        species="Homo sapiens", compartments=compartments, evidence="motif_rule",
        references=tuple(papers) + (
            "https://www.uniprot.org/uniprotkb/%s/entry" % accession,),
        assay=assay,
        limitations=limitations + " Partial specificity rule; no kinetic or matrix calibration.",
        motif_strictness=strictness, strictness_basis=basis)
    return PeptidaseMotif(metadata, topology, removed, left, right, description, **kwargs)


_RULES = (
    _rule("app1-xp", "XPNPEP1", "Q9NQW7", ("cytosol",),
          ("https://pubmed.ncbi.nlm.nih.gov/11106490/", "https://doi.org/10.1074/jbc.M710274200"),
          "n_terminal", 1, ".", "P", "Exposed N-terminal X|Pro",
          "Recombinant human cytosolic aminopeptidase P; RPP and bradykinin hydrolysis",
          "Removes the first residue. Manganese dependence and full substrate context affect activity; "
          "distinct from extracellular XPNPEP2 and dipeptidyl processing.",
          strictness="required",
          basis="Aminopeptidase P is defined by hydrolysis of the X-Pro bond; the source reports "
                "the recombinant human cytosolic enzyme hydrolysing the X-Pro bond of bradykinin "
                "and substance P."),
    _rule("tpp2-tripeptidyl", "TPP2", "P29144", ("cytosol",),
          ("https://doi.org/10.4049/jimmunol.169.8.4161",),
          "n_terminal", 3, "..[^P]", "[^P]", "Exposed tripeptide with non-Pro P1 and P1-prime",
          "Immunopurified human TPP2 processing of RU1 antigen precursors",
          "Partial proline constraint, not proof of turnover for every matching peptide. "
          "Assembly, other residues and longer substrate context are unmodeled. "
          "Not bacterial Xaa-Xaa-Pro tripeptidase specificity and not an endopeptidase model.",
          min_length=4,
          strictness="permissive",
          basis="Removing an N-terminal tripeptide is this enzyme's topology, not a selective "
                "substrate rule. The proline constraints are partial and the source does not "
                "establish turnover of every matching peptide."),
    _rule("npepps-n-terminal", "NPEPPS", "P55786", ("cytosol",),
          ("https://doi.org/10.4049/jimmunol.169.8.4161",),
          "n_terminal", 1, "[^GP]", "[^GP]", "Exposed first bond lacking poor Gly/Pro contexts",
          "Human puromycin-sensitive aminopeptidase processing of RU1 precursors",
          "Broad, low-specificity context flag only. Gly-containing bonds can be slowly cleaved; "
          "non-matches do not establish resistance. Other residues, competing aminopeptidases "
          "and full precursor sequence affect processing.",
          strictness="permissive",
          basis="Puromycin-sensitive aminopeptidase acts broadly; the Gly/Pro exclusion is a weak "
                "context flag and Gly-containing bonds can still be cleaved slowly."),
    _rule("ace-dipeptidyl", "ACE", "P12821", ("plasma", "serum", "extracellular"),
          ("https://doi.org/10.1042/BJ20040634", "https://doi.org/10.1038/srep13742"),
          "c_terminal", 2, ".", "[^P][^DE]", "Free C-terminal dipeptide: non-Pro followed by non-Asp/Glu",
          "Human ACE dipeptide hydrolysis; angiotensins and N-acetyl-SDKP",
          "Ordinary dipeptide recognition only. Domain, chloride, concentration and context affect activity. "
          "Exceptional endopeptide/tripeptide cleavage (including amidated substance P) is omitted.",
          min_length=3, n_chemistry=("free", "acetylated"),
          strictness="required",
          basis="Within the ordinary dipeptidyl-carboxypeptidase route the source's angiotensin "
                "series requires a free C-terminal dipeptide with non-Pro at P1-prime. "
                "Exceptional endopeptidase cleavages are outside this rule."),
    _rule("mme-hydrophobic", "MME", "P08473", ("extracellular",),
          ("https://pubmed.ncbi.nlm.nih.gov/6349683/", "https://pubmed.ncbi.nlm.nih.gov/2417254/"),
          "internal", 0, ".", "[FILY]", "Selected hydrophobic P1-prime preference: Phe/Ile/Leu/Tyr",
          "Purified human kidney neprilysin cleavage of enkephalin, kinins and angiotensins",
          "Incomplete preference flag; whole sequence affects site selection. Other residues can be cleaved. "
          "The 2-30-residue input scope is conservative, not an absolute enzyme size cutoff. "
          "Membrane exposure and circulating activity are unmodeled.",
          max_length=30, c_chemistry=("free", "amidated"),
          strictness="preferred",
          basis="The source reports a preference for hydrophobic residues after the scissile "
                "bond, and other residues are also cleaved, so a non-match is weak evidence."),
    _rule("cpb2-basic", "CPB2", "Q96IY4", ("plasma", "serum", "extracellular"),
          ("https://doi.org/10.1074/jbc.274.49.35046", "https://pmc.ncbi.nlm.nih.gov/articles/PMC2613638/"),
          "c_terminal", 1, ".", "[KR]", "Activated CPB2 removal of free C-terminal Lys/Arg",
          "Human plasma-derived active TAFI/CPB2 peptide hydrolysis",
          "Requires explicit enzyme_state=active. Zymogen abundance does not establish active enzyme. "
          "Activation, spontaneous inactivation, inhibitors and serum clotting effects are unmodeled.",
          requires_activation=True,
          strictness="required",
          basis="Basic carboxypeptidase activity requires a free C-terminal Lys or Arg. The grade "
                "covers recognition only and still presumes an explicitly activated enzyme."),
    _rule("cpn-basic", "CPN1", "P15169", ("plasma", "serum", "extracellular"),
          ("https://doi.org/10.1016/0003-9861(75)90104-6",),
          "c_terminal", 1, ".", "[KR]", "Free C-terminal Lys/Arg removal",
          "Human plasma CPN peptide-substrate assays",
          "Penultimate residues affect rates; Lys and Arg are not kinetically equivalent.",
          strictness="required",
          basis="Carboxypeptidase N removes a free C-terminal Lys or Arg; penultimate residues "
                "change the rate but not the requirement."),
    _rule("app2-xp", "XPNPEP2", "O43895", ("plasma", "extracellular"),
          ("https://pubmed.ncbi.nlm.nih.gov/15361070/",),
          "n_terminal", 1, ".", "P", "Exposed N-terminal X|Pro",
          "Recombinant human membrane aminopeptidase P and kinin substrates",
          "Removes the FIRST residue; this is not DPP-like dipeptide removal.",
          strictness="required",
          basis="Aminopeptidase P is defined by hydrolysis of the X-Pro bond; the source's kinin "
                "substrates all carry proline in the second position."),
    _rule("fap-dipeptidyl", "FAP", "Q12884", ("extracellular", "plasma"),
          ("https://pubmed.ncbi.nlm.nih.gov/16410248/",),
          "n_terminal", 2, ".P", "[^P]", "Exposed N-terminal X-Pro|non-Pro",
          "Purified human FAP dipeptide substrate profiling",
          "Dipeptidyl activity only; P2 preferences and whole-peptide context omitted.",
          min_length=3,
          strictness="required",
          basis="FAP dipeptidyl activity requires proline in the second position; the source's "
                "profiling is built on X-Pro dipeptide substrates."),
    _rule("fap-endo-gp", "FAP", "Q12884", ("extracellular", "plasma"),
          ("https://pubmed.ncbi.nlm.nih.gov/16480718/",
           "https://pubmed.ncbi.nlm.nih.gov/16410248/",
           "https://pubmed.ncbi.nlm.nih.gov/22750443/"),
          "internal", 0, "GP", "[^P]", "Gly-Pro|non-Pro endopeptidase recognition",
          "Human FAP substrate profiling around alpha-2-antiplasmin cleavage site",
          "P3 preferences and accessibility omitted; independent of dipeptidyl activity.",
          min_length=3, n_chemistry=("free", "acetylated"),
          strictness="required",
          basis="The source's endopeptidase profiling establishes Gly-Pro before the scissile "
                "bond and a non-proline after it; P3 preferences and accessibility remain "
                "unmodeled."),
    _rule("enpep-acidic", "ENPEP", "Q07075", ("extracellular",),
          ("https://pubmed.ncbi.nlm.nih.gov/23888046/",),
          "n_terminal", 1, "[DE]", ".", "Exposed N-terminal Asp/Glu|X",
          "Human aminopeptidase A structural/substrate specificity experiments",
          "Calcium-dependent preference; enzyme exposure and calcium are not modeled.",
          strictness="preferred",
          basis="Aminopeptidase A prefers an acidic N-terminal residue and the preference is "
                "calcium dependent, so a non-match does not exclude slower cleavage."),
    _rule("anpep-ala", "ANPEP", "P15144", ("extracellular",),
          ("https://pubmed.ncbi.nlm.nih.gov/22932899/",),
          "n_terminal", 1, "A", ".", "Preferred N-terminal Ala|X",
          "Human aminopeptidase N structural/biochemical specificity",
          "Only an Ala preference flag. Broad alternative substrates and X-Pro dipeptide removal omitted.",
          strictness="permissive",
          basis="Aminopeptidase N acts broadly; alanine is a favoured first residue but many "
                "other N-termini are cleaved, so a match adds little."),
    _rule("dpp8-xp-xa", "DPP8", "Q6V1X1", ("cytosol",),
          ("https://pubmed.ncbi.nlm.nih.gov/11012666/",
           "https://pmc.ncbi.nlm.nih.gov/articles/PMC3656252/"),
          "n_terminal", 2, ".[PA]", "[^P]", "Exposed X-(Pro/Ala)|non-Pro",
          "Human DPP8 biochemical characterization and substrate degradomics",
          "Qualitative rule only; DPP4 and C. elegans DPF-3 coefficients are inapplicable.",
          min_length=3,
          strictness="required",
          basis="Dipeptidyl peptidase activity requires proline or alanine in the second "
                "position; the source's degradomics substrates share that constraint."),
    _rule("dpp9-xp-xa", "DPP9", "Q86TI2", ("cytosol",),
          ("https://pubmed.ncbi.nlm.nih.gov/19667070/",
           "https://pmc.ncbi.nlm.nih.gov/articles/PMC3656252/"),
          "n_terminal", 2, ".[PA]", "[^P]", "Exposed X-(Pro/Ala)|non-Pro",
          "Human DPP9 peptide/antigen processing assays and substrate degradomics",
          "Shared recognition pattern with DPP8 does not imply identical substrate rates.",
          min_length=3,
          strictness="required",
          basis="Dipeptidyl peptidase activity requires proline or alanine in the second "
                "position; the source's antigen-processing and degradomics substrates share that "
                "constraint."),
    _rule("prep-pro", "PREP", "P48147", ("cytosol", "extracellular", "serum"),
          ("https://pubmed.ncbi.nlm.nih.gov/22750443/",),
          "internal", 0, ".P", ".", "Internal post-proline X-Pro|X recognition",
          "Recombinant human POP/FAP peptide substrate profiling",
          "Conservative 4-30-residue model domain; not an absolute enzyme size cutoff. "
          "Weaker non-Pro cleavage and flanking charge preferences omitted; abundance unmodeled.",
          min_length=4, max_length=30,
          strictness="required",
          basis="Prolyl oligopeptidase requires proline before the scissile bond; documented "
                "weaker non-proline cleavage is outside this rule's scope."),
    _rule("erap2-basic", "ERAP2", "Q6P179", ("er",),
          ("https://pubmed.ncbi.nlm.nih.gov/12799365/",
           "https://pubmed.ncbi.nlm.nih.gov/26381406/"),
          "n_terminal", 1, "[RK]", ".", "Preferred exposed N-terminal Arg/Lys|X",
          "Human ERAP2 biochemical specificity and peptide-complex structures",
          "A basic-residue preference flag only; length, allotype and internal sequence affect trimming.",
          strictness="preferred",
          basis="The source reports a basic-residue preference at the N terminus rather than a "
                "requirement; ERAP2 trims many other residues more slowly."),
)


def cleavage_models(include_optional=False):
    """List built-in model metadata without loading optional runtimes."""
    models = (DPP4qPISA.model,) + tuple(rule.model for rule in _RULES) + tuple(
        reference.model for reference in substrate_references())
    names = [m.name for m in models]
    if len(set(names)) != len(names):
        raise ValueError("Duplicate cleavage model name in the built-in panel: %r" % (
            sorted({n for n in names if names.count(n) > 1}),))
    if include_optional:
        from .eramer_cleavage import ERAMERCleavage
        models += (ERAMERCleavage.model,)
    return models


def get_cleavage_model(name, *, enzyme_state=None):
    """Select an exact model name; never silently substitute another enzyme."""
    if enzyme_state is not None and name != "cpb2-basic":
        raise ValueError("Explicit enzyme state is currently supported only for cpb2-basic")
    if name == "dpp4-qpisa":
        return DPP4qPISA()
    if name == "eramer-step":
        from .eramer_cleavage import ERAMERCleavage
        return ERAMERCleavage()
    for rule in _RULES:
        if rule.model.name == name:
            return replace(rule, enzyme_state=enzyme_state) if enzyme_state is not None else rule
    for reference in substrate_references():
        if reference.model.name == name:
            return reference
    raise ValueError("Unknown cleavage model %r; choices: %s" % (
        name, ", ".join(model.name for model in cleavage_models(include_optional=True))))


def predict_cleavage(peptide, models=None, compartment=None, *, enzyme_states=None):
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
    models = tuple(dict.fromkeys(models))
    enzyme_states = {} if enzyme_states is None else dict(enzyme_states)
    if set(enzyme_states) - {"CPB2"}:
        raise ValueError("Explicit enzyme state is currently supported only for CPB2")
    if enzyme_states and "cpb2-basic" not in models:
        raise ValueError("CPB2 state supplied without selecting cpb2-basic")
    results = []
    for name in models:
        predictor = get_cleavage_model(name, enzyme_state=enzyme_states.get("CPB2") if name == "cpb2-basic" else None)
        if compartment is not None and compartment not in predictor.model.compartments:
            raise ValueError("Model %s is not annotated for compartment %s" % (name, compartment))
        results.append(predictor.predict(peptide))
    return tuple(results)
