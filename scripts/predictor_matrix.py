"""Source of truth for ``docs/predictor-matrix.md``.

The matrix is curated here (install route, license tier, input shape) and
checked against the code by ``tests/test_docs_predictor_matrix.py``, which
fails when an exported predictor class or CLI name has no row, or when a row
disagrees with what an instance reports. Regenerate the page with::

    python scripts/predictor_matrix.py            # rewrite docs/predictor-matrix.md
    python scripts/predictor_matrix.py --check    # exit 1 if the page is stale
"""

import argparse
import sys
from pathlib import Path

PAGE = Path(__file__).resolve().parent.parent / "docs" / "predictor-matrix.md"

# License tiers, defined once in docs/licensing.md.
LICENSE_TIERS = {
    "open": "Open source; code and weights can be fetched or installed freely",
    "academic": "Academic / non-commercial; review the upstream license first",
    "dtu": "DTU academic license, identity-bound; you request it from DTU",
    "unlicensed": "Upstream publishes no license; `--accept-license` records your own authorization",
    "builtin": "Built into mhctools; nothing to download",
}

FAMILIES = [
    ("binding", "MHC binding and presentation", "predictors/binding.md"),
    ("processing", "Antigen processing", "predictors/processing.md"),
    ("immunogenicity", "Immunogenicity", "predictors/immunogenicity.md"),
    ("tcr", "TCR specificity", "predictors/tcr.md"),
    ("peptide-pk", "Peptide half-life", "predictors/peptide-pk.md"),
]

# Fields: family, name, classes, cli, kinds, mhc_class, inputs, lengths,
# artifact (`mhctools fetch` / `ls` name(s), or None), install, license, page.
ROWS = [
    dict(family="processing", name="PhageScout sequence profiles",
         classes=["PhageScout"], cli=[], kinds=[], mhc_class="none",
         inputs="canonical sequence; human ELANE/CTSG, inferred site contexts",
         lengths="2+ for aligned PWMs; 5+ for five-mer scoring; 9+ for aligned peptide profiles",
         artifact="phagescout", install="bundled profiles; optional `mhctools fetch phagescout`; no external runtime",
         license="open", page="cleavage/phagescout.md"),
    dict(family="processing", name="ITCell cathepsin profiles",
         classes=["ITCellCleavage"], cli=[], kinds=["endolysosomal_cleavage"], mhc_class="none",
         inputs="canonical sequence; B/S internal, H initial N-terminal trimming",
         lengths="2+; missing flanks contribute zero as in author code", artifact=None,
         install="bundled human B/S/H profiles; no external runtime",
         license="open", page="cleavage/itcell.md"),
    dict(family="processing", name="CleaveNet",
         classes=["CleaveNet"], cli=[], kinds=["substrate_cleavage"], mhc_class="none",
         inputs="whole substrates or ten-residue windows; dedicated native result",
         lengths="1-10; centered padding for shorter inputs", artifact="cleavenet",
         install="`mhctools fetch cleavenet` + isolated TensorFlow 2.18.0",
         license="open", page="cleavage/cleavenet.md"),
    dict(family="binding", name="NetMHCpan 4.1 / 4.2",
         classes=["NetMHCpan", "NetMHCpan41", "NetMHCpan41_BA",
                  "NetMHCpan41_EL", "NetMHCpan42", "NetMHCpan42_BA",
                  "NetMHCpan42_EL"],
         cli=["netmhcpan", "netmhcpan41", "netmhcpan41-ba", "netmhcpan41-el",
              "netmhcpan42", "netmhcpan42-ba", "netmhcpan42-el"],
         kinds=["pMHC_affinity", "pMHC_presentation"], mhc_class="I",
         inputs="peptides + alleles", lengths="9",
         artifact="netmhcpan", install="your licensed DTU install",
         license="dtu", page="predictors/binding.md#netmhcpan"),
    dict(family="binding", name="NetMHCpan 4.0",
         classes=["NetMHCpan4", "NetMHCpan4_BA", "NetMHCpan4_EL"],
         cli=["netmhcpan4", "netmhcpan4-ba", "netmhcpan4-el"],
         kinds=["pMHC_affinity", "pMHC_presentation"], mhc_class="I",
         inputs="peptides + alleles", lengths="9",
         artifact="netmhcpan", install="your licensed DTU install",
         license="dtu", page="predictors/binding.md#netmhcpan"),
    dict(family="binding", name="NetMHCpan 3.0 / 2.8",
         classes=["NetMHCpan3", "NetMHCpan28"],
         cli=["netmhcpan3", "netmhcpan28"],
         kinds=["pMHC_affinity"], mhc_class="I",
         inputs="peptides + alleles", lengths="9",
         artifact="netmhcpan", install="your licensed DTU install",
         license="dtu", page="predictors/binding.md#netmhcpan"),
    dict(family="binding", name="NetMHC 3.4 / 4.0",
         classes=["NetMHC", "NetMHC3", "NetMHC4"],
         cli=["netmhc", "netmhc3", "netmhc4"],
         kinds=["pMHC_affinity"], mhc_class="I",
         inputs="peptides + alleles", lengths="9",
         artifact="netmhc", install="your licensed DTU install",
         license="dtu", page="predictors/binding.md#netmhc"),
    dict(family="binding", name="NetMHCcons",
         classes=["NetMHCcons"], cli=["netmhccons"],
         kinds=["pMHC_affinity"], mhc_class="I",
         inputs="peptides + alleles", lengths="9",
         artifact="netmhccons", install="your licensed DTU install",
         license="dtu", page="predictors/binding.md#netmhccons"),
    dict(family="binding", name="NetMHCIIpan 4.3 / 4.1 / 4.0 / 3.x",
         classes=["NetMHCIIpan", "NetMHCIIpan3", "NetMHCIIpan4",
                  "NetMHCIIpan4_BA", "NetMHCIIpan4_EL", "NetMHCIIpan43",
                  "NetMHCIIpan43_BA", "NetMHCIIpan43_EL"],
         cli=["netmhciipan", "netmhciipan3", "netmhciipan4",
              "netmhciipan4-ba", "netmhciipan4-el", "netmhciipan43",
              "netmhciipan43-ba", "netmhciipan43-el"],
         kinds=["pMHC_affinity", "pMHC_presentation"], mhc_class="II",
         inputs="peptides + alleles", lengths="15-20",
         artifact="netmhciipan", install="your licensed DTU install",
         license="dtu", page="predictors/binding.md#netmhciipan"),
    dict(family="binding", name="NetMHCstabpan",
         classes=["NetMHCstabpan"], cli=["netmhcstabpan"],
         kinds=["pMHC_stability"], mhc_class="I",
         inputs="peptides + alleles", lengths="(none; pass lengths)",
         artifact="netmhcstabpan", install="your licensed DTU install",
         license="dtu", page="predictors/binding.md#netmhcstabpan"),
    dict(family="binding", name="MHCflurry",
         classes=["MHCflurry", "MHCflurry_Affinity"],
         cli=["mhcflurry", "mhcflurry-affinity"],
         kinds=["pMHC_affinity", "pMHC_presentation", "antigen_processing"],
         mhc_class="I", inputs="peptides + alleles (flanks optional)",
         lengths="9",
         artifact=["mhcflurry", "mhcflurry-affinity"],
         install="`mhctools fetch mhcflurry`",
         license="open", page="predictors/binding.md#mhcflurry"),
    dict(family="binding", name="BigMHC",
         classes=["BigMHC", "BigMHC_EL", "BigMHC_IM"],
         cli=["bigmhc", "bigmhc-el", "bigmhc-im"],
         kinds=["pMHC_presentation", "immunogenicity"], mhc_class="I",
         inputs="peptides + alleles", lengths="(set per call)",
         artifact="bigmhc",
         install="`mhctools fetch bigmhc --accept-license` + PyTorch",
         license="academic", page="predictors/binding.md#bigmhc"),
    dict(family="binding", name="CapHLA",
         classes=["CapHLA", "CapHLA_BA", "CapHLA_EL"],
         cli=["caphla", "caphla-ba", "caphla-el"],
         kinds=["pMHC_affinity", "pMHC_presentation"], mhc_class="I and II",
         inputs="peptides + alleles", lengths="(7-25 supported)",
         artifact="caphla",
         install="`pip install \"mhctools[caphla]\"` + `mhctools fetch caphla`",
         license="open", page="predictors/binding.md#caphla"),
    dict(family="binding", name="MixMHCpred",
         classes=["MixMHCpred"], cli=["mixmhcpred"],
         kinds=["pMHC_presentation"], mhc_class="I",
         inputs="peptides + alleles", lengths="9 (8-14 supported)",
         artifact="mixmhcpred",
         install="upstream release + `MIXMHCPRED_PATH`",
         license="academic", page="predictors/binding.md#mixmhcpred"),
    dict(family="binding", name="MixMHC2pred",
         classes=["MixMHC2pred"], cli=["mixmhc2pred"],
         kinds=["pMHC_presentation"], mhc_class="II",
         inputs="peptides + alleles", lengths="15",
         artifact="mixmhc2pred", install="upstream release (with `PWMdef/`)",
         license="academic", page="predictors/binding.md#mixmhc2pred"),
    dict(family="binding", name="SMM / SMM-PMBEC",
         classes=["SMM", "SMMPMBEC"], cli=["smm", "smm-pmbec"],
         kinds=["pMHC_affinity"], mhc_class="I",
         inputs="peptides + alleles", lengths="9",
         artifact="smm",
         install="`python scripts/setup_test_backends.py smm --accept-license`"
                 " or `IEDB_MHCI_EXECUTABLE`",
         license="open", page="predictors/binding.md#smm-and-smm-pmbec"),
    dict(family="binding", name="Historical IEDB names",
         classes=["IedbNetMHCpan", "IedbNetMHCcons", "IedbNetMHCIIpan",
                  "IedbSMM", "IedbSMM_PMBEC"],
         cli=["netmhcpan-iedb", "netmhccons-iedb", "netmhciipan-iedb",
              "smm-iedb", "smm-pmbec-iedb"],
         kinds=["pMHC_affinity"], mhc_class="I and II",
         inputs="peptides + alleles", lengths="8-11 (class I), 15-20 (II)",
         artifact=None, install="the local predictor each one maps to",
         license="dtu",
         page="predictors/binding.md#compatibility-names-for-the-old-iedb-predictors"),
    dict(family="binding", name="RandomBindingPredictor",
         classes=["RandomBindingPredictor"], cli=["random"],
         kinds=["pMHC_affinity"], mhc_class="I",
         inputs="peptides + alleles", lengths="9",
         artifact=None, install="built in", license="builtin",
         page="predictors/binding.md#randombindingpredictor"),

    dict(family="processing", name="Pepsickle",
         classes=["Pepsickle"], cli=["pepsickle"],
         kinds=["proteasome_cleavage"], mhc_class="none",
         inputs="peptides (flanks recommended)", lengths="9",
         artifact="pepsickle", install="`pip install pepsickle`",
         license="open", page="predictors/processing.md#pepsickle"),
    dict(family="processing", name="NetChop",
         classes=["NetChop"], cli=["netchop"],
         kinds=["proteasome_cleavage"], mhc_class="none",
         inputs="peptides (flanks recommended)", lengths="9",
         artifact="netchop",
         install="your licensed DTU install (`NETCHOP_HOME`)",
         license="dtu", page="predictors/processing.md#netchop"),
    dict(family="processing", name="NetCleave",
         classes=["NetCleave", "NetCleave_I", "NetCleave_II"],
         cli=["netcleave", "netcleave-i", "netcleave-ii"],
         kinds=["proteasome_cleavage", "endolysosomal_cleavage"],
         mhc_class="I or II", inputs="peptides + C-terminal flank (>= 3 residues)",
         lengths="9 (I), 15 (II)", artifact="netcleave",
         install="`mhctools fetch netcleave --accept-license`",
         license="unlicensed", page="predictors/processing.md#netcleave"),
    dict(family="processing", name="DeepTAP",
         classes=["DeepTAP"], cli=["deeptap"],
         kinds=["tap_transport"], mhc_class="none",
         inputs="peptides only", lengths="(any)",
         artifact="deeptap", install="`mhctools fetch deeptap` + torch Python",
         license="open", page="predictors/processing.md#deeptap"),
    dict(family="processing", name="ERAMER",
         classes=["ERAMER"], cli=["eramer"],
         kinds=["erap_trimming"], mhc_class="none (class I context)",
         inputs="peptides only (9-16mer precursors)", lengths="(9-16 supported)",
         artifact="eramer", install="`mhctools fetch eramer` + `openpyxl`",
         license="open", page="predictors/processing.md#eramer"),

    dict(family="immunogenicity", name="Calis",
         classes=["Calis"], cli=["calis"],
         kinds=["immunogenicity"], mhc_class="I",
         inputs="peptides only", lengths="(any)",
         artifact="calis", install="built in", license="builtin",
         page="predictors/immunogenicity.md#calis"),
    dict(family="immunogenicity", name="PRIME",
         classes=["PRIME"], cli=["prime"],
         kinds=["immunogenicity"], mhc_class="I",
         inputs="peptides + alleles", lengths="9",
         artifact="prime", install="upstream clone + MixMHCpred 3.0+",
         license="academic", page="predictors/immunogenicity.md#prime"),
    dict(family="immunogenicity", name="DeepImmuno",
         classes=["DeepImmuno"], cli=["deepimmuno"],
         kinds=["immunogenicity"], mhc_class="I",
         inputs="peptides + alleles", lengths="9-10 only",
         artifact="deepimmuno",
         install="`mhctools fetch deepimmuno` + Keras-2-capable Python",
         license="open", page="predictors/immunogenicity.md#deepimmuno"),
    dict(family="immunogenicity", name="TLimmuno2",
         classes=["TLimmuno2"], cli=["tlimmuno2"],
         kinds=["immunogenicity"], mhc_class="II",
         inputs="peptides + class II alleles", lengths="(set per call)",
         artifact="tlimmuno2",
         install="`mhctools fetch tlimmuno2 --accept-license`",
         license="unlicensed", page="predictors/immunogenicity.md#tlimmuno2"),

    dict(family="tcr", name="NetTCR",
         classes=["NetTCR"], cli=[],
         kinds=["pMHC_TCR_binding"], mhc_class="I",
         inputs="(peptide, TCR) pairs", lengths="(n/a)",
         artifact="nettcr",
         install="`mhctools fetch nettcr --accept-license` + `mhctools[nettcr]`",
         license="academic", page="predictors/tcr.md#nettcr"),
    dict(family="tcr", name="Tulip",
         classes=["Tulip"], cli=[],
         kinds=["pMHC_TCR_binding"], mhc_class="I",
         inputs="(peptide, TCR) pairs + `mhc=`", lengths="(n/a)",
         artifact="tulip", install="`mhctools fetch tulip` + isolated Python 3.11",
         license="open", page="predictors/tcr.md#tulip"),
    dict(family="tcr", name="MixTCRpred",
         classes=["MixTCRpred"], cli=[],
         kinds=["pMHC_TCR_binding"], mhc_class="I or II (per model)",
         inputs="TCRs against one fixed pMHC target per model",
         lengths="(n/a)", artifact="mixtcrpred",
         install="`mhctools[mixtcrpred]` + `mhctools fetch mixtcrpred --accept-license`",
         license="academic", page="predictors/tcr.md#mixtcrpred"),

    dict(family="peptide-pk", name="PeptiVerse",
         classes=["PeptiVerse"], cli=["peptiverse"],
         kinds=["peptide_half_life"], mhc_class="none",
         inputs="peptides or `PeptideInput`", lengths="(sequence only)",
         artifact=None,
         install="pinned snapshots (`PEPTIVERSE_HOME`, `PEPTIVERSE_ESM_HOME`)",
         license="open", page="predictors/peptide-pk.md#peptiverse"),
    dict(family="peptide-pk", name="PlifePred2",
         classes=["PlifePred2"], cli=["plifepred2"],
         kinds=["peptide_half_life"], mhc_class="none",
         inputs="peptides only (12-100 residues, natural)", lengths="(12-100)",
         artifact=None,
         install="`plifepred2==1.0` + Pfeature (`PLIFEPRED2_HOME`, `PFEATURE_HOME`)",
         license="open", page="predictors/peptide-pk.md#plifepred2"),
]

# Names that are not predictors even though they are exported and have
# `predict`/`kind_support` (shared base classes).
NOT_PREDICTORS = {"ProcessingPredictor", "ProteasomePredictor"}


def covered_classes():
    return {c for row in ROWS for c in row["classes"]}


def covered_cli_names():
    return {c for row in ROWS for c in row["cli"]}


def covered_artifacts():
    names = set()
    for row in ROWS:
        artifact = row["artifact"]
        if isinstance(artifact, str):
            names.add(artifact)
        elif artifact:
            names.update(artifact)
    return names


def _code(items):
    return ", ".join("`%s`" % i for i in items) if items else "(none)"


def render():
    out = [
        "# Predictor matrix",
        "",
        "<!-- Generated by scripts/predictor_matrix.py; do not edit by hand. -->",
        "",
        "Reference for supported predictors, Python classes, command-line "
        "names, inputs, and installation routes. For model selection, see "
        "[choosing a predictor](choosing.md).",
        "",
        "- Kinds are defined in [prediction kinds](kinds.md). A predictor "
        "emits only the kinds that apply (for example `NetMHCpan41_EL` emits "
        "presentation only).",
        "- **MHC class** is the class the model scores. None means the "
        "prediction is MHC-independent and `Prediction.allele` is empty.",
        "- **Default lengths** apply to `predict_proteins()` and the "
        "command line when you pass no lengths; they are narrower than "
        "what the model supports. See [peptide lengths](predictors/index.md#peptide-lengths).",
        "- **License tiers** are defined in [licensing](licensing.md); "
        "`mhctools fetch` behavior is in [getting models](artifacts.md).",
        "",
    ]
    for key, title, page in FAMILIES:
        rows = [r for r in ROWS if r["family"] == key]
        out += ["## %s" % title, "",
                "See the [%s guide](%s) for examples and model notes." % (
                    title.lower(), page), "",
                "| Predictor | MHC class | Input |",
                "|---|---|---|"]
        for r in rows:
            out.append("| [%s](%s) | %s | %s |" % (
                r["name"], r["page"], r["mhc_class"], r["inputs"]))
        out.append("")
        for r in rows:
            out += [
                "### %s" % r["name"], "",
                "- Python classes: %s" % _code(r["classes"]),
                "- CLI names: %s" % _code(r["cli"]),
                "- Prediction kinds: %s" % _code(r["kinds"]),
                "- Default scanning lengths: %s" % r["lengths"],
                "- Installation: %s" % r["install"],
                "- License: %s" % r["license"],
                "",
            ]
    out += ["## License tiers", "", "| Tier | Meaning |", "|---|---|"]
    for tier, meaning in LICENSE_TIERS.items():
        out.append("| %s | %s |" % (tier, meaning))
    out += ["",
            "TCR predictors have no `mhctools` command-line prediction name: "
            "their input is a peptide plus a TCR, not an allele. "
            "`MixTCRpred` has its own `mhctools mixtcrpred` subcommand; see "
            "[the command line](cli.md).",
            "",
            "Cleavage models (`mhctools cleavage`) are a separate panel; see "
            "[cleavage models](cleavage/models.md).",
            ""]
    return "\n".join(out)


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true",
                        help="exit 1 if docs/predictor-matrix.md is stale")
    args = parser.parse_args(argv)
    text = render()
    if args.check:
        current = PAGE.read_text() if PAGE.exists() else ""
        if current != text:
            print("docs/predictor-matrix.md is stale; run "
                  "python scripts/predictor_matrix.py", file=sys.stderr)
            return 1
        return 0
    PAGE.write_text(text)
    return 0


if __name__ == "__main__":
    sys.exit(main())
