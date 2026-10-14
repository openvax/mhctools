# Vaccine construct processing report

Use `mhctools vaccine-report` to inspect where proposed cuts fall relative to an
intended epitope in a synthetic long peptide (SLP) or an RNA-encoded protein.
It produces a PDF sequence map and files you can audit. A cut beside a target
may help release it; a cut inside the target may destroy it. The report keeps
those observations separate by model and processing route.

This report **does not estimate** the fraction of an injected dose that reaches
lymph-node APCs, the half-life of an intact core, or the fraction of APCs that
present it. For a clearly labeled, uncalibrated dose sensitivity analysis, use
[candidate window exposure](vaccine-windows.md). For the report's per-target
evidence summary, see [vaccine target summary](vaccine-target-summary.md).

## Run the example

PDF generation is included in a standard install. This example first runs a
built-in cleavage rule, then gives its actual JSON output to the PDF report.
Copy and run it from one directory:

```sh
pip install mhctools
mhctools cleavage \
  --sequence AAAKIILNFTTLAAA \
  --model mme-hydrophobic \
  --out cleavage.json
cat > vaccine-report.json <<'JSON'
{
  "schema_version": 1,
  "constructs": [{
    "id": "slp-demo",
    "sequence": "AAAKIILNFTTLAAA",
    "delivery": "synthetic_long_peptide",
    "intended_epitopes": [
      {"label": "KIILNFTTL", "start": 4, "end": 12}
    ],
    "cleavage_tracks": [{
      "context": "extracellular_interstitial",
      "cleavage_json": "cleavage.json"
    }]
  }]
}
JSON
mhctools vaccine-report \
  --input vaccine-report.json \
  --output-dir reports
```

The command prints a timestamped directory inside `reports/`. Open
`vaccine-processing-report.pdf` there. The manifest is also available as a
[downloadable JSON file](vaccine-report-simple.json); run the `cleavage`
command above beside it first. The report reads the CLI's schema-v2 output and
copies its model identity, evidence type, site decisions, and provenance into
`normalized-input.json`. You need only declare `context`: the biological route
in which you want to inspect that evidence. A model's enzyme-location
annotation alone does not establish exposure after vaccination. `name` is
optional and defaults to the model name. If the cleavage JSON contains several
results for the same sequence, add `"model": "mme-hydrophobic"` to select one.

For this input, the tool reports `motif_rule` evidence with matched bonds 4, 5,
6, 8, and 11. The rule flags a partial extracellular MME preference for a
hydrophobic residue after a bond; it is **not** a measured cut, cleavage rate,
or probability. No MHC predictions are supplied, so the report does not draw
MHC-window bars.

## Inspect `cleavage.json`

The `mhctools cleavage` command above writes this complete JSON. `models`
holds the model identity, evidence type, source references, and limits;
`results` holds the assessed sequence and each internal bond. Here, every
`score` is `null` because this is categorical motif evidence; `matched` and
`not_matched` are the actual site decisions. `source_bond` equals `bond`
because this input has `source_start: 0`.

```json
{
  "schema_version": 2,
  "models": {
    "mme-hydrophobic": {
      "name": "mme-hydrophobic",
      "version": "1",
      "enzyme": "MME",
      "uniprot": "P08473",
      "species": "Homo sapiens",
      "compartments": [
        "extracellular"
      ],
      "evidence": "motif_rule",
      "references": [
        "https://pubmed.ncbi.nlm.nih.gov/6349683/",
        "https://pubmed.ncbi.nlm.nih.gov/2417254/",
        "https://www.uniprot.org/uniprotkb/P08473/entry"
      ],
      "assay": "Purified human kidney neprilysin cleavage of enkephalin, kinins and angiotensins",
      "limitations": "Incomplete preference flag; whole sequence affects site selection. Other residues can be cleaved. The 2-30-residue input scope is conservative, not an absolute enzyme size cutoff. Membrane exposure and circulating activity are unmodeled. Partial specificity rule; no kinetic or matrix calibration.",
      "score_name": null,
      "score_units": null,
      "scored_endpoint": null,
      "motif_strictness": "preferred",
      "strictness_basis": "The source reports a preference for hydrophobic residues after the scissile bond, and other residues are also cleaved, so a non-match is weak evidence."
    }
  },
  "results": [
    {
      "peptide": {
        "sequence": "AAAKIILNFTTLAAA",
        "n_term": "free",
        "c_term": "free",
        "source_id": null,
        "source_start": 0
      },
      "model": "mme-hydrophobic",
      "sites": [
        {
          "bond": 1,
          "status": "not_matched",
          "reason": "Selected hydrophobic P1-prime preference: Phe/Ile/Leu/Tyr",
          "score": null,
          "source_bond": 1
        },
        {
          "bond": 2,
          "status": "not_matched",
          "reason": "Selected hydrophobic P1-prime preference: Phe/Ile/Leu/Tyr",
          "score": null,
          "source_bond": 2
        },
        {
          "bond": 3,
          "status": "not_matched",
          "reason": "Selected hydrophobic P1-prime preference: Phe/Ile/Leu/Tyr",
          "score": null,
          "source_bond": 3
        },
        {
          "bond": 4,
          "status": "matched",
          "reason": "Selected hydrophobic P1-prime preference: Phe/Ile/Leu/Tyr",
          "score": null,
          "source_bond": 4
        },
        {
          "bond": 5,
          "status": "matched",
          "reason": "Selected hydrophobic P1-prime preference: Phe/Ile/Leu/Tyr",
          "score": null,
          "source_bond": 5
        },
        {
          "bond": 6,
          "status": "matched",
          "reason": "Selected hydrophobic P1-prime preference: Phe/Ile/Leu/Tyr",
          "score": null,
          "source_bond": 6
        },
        {
          "bond": 7,
          "status": "not_matched",
          "reason": "Selected hydrophobic P1-prime preference: Phe/Ile/Leu/Tyr",
          "score": null,
          "source_bond": 7
        },
        {
          "bond": 8,
          "status": "matched",
          "reason": "Selected hydrophobic P1-prime preference: Phe/Ile/Leu/Tyr",
          "score": null,
          "source_bond": 8
        },
        {
          "bond": 9,
          "status": "not_matched",
          "reason": "Selected hydrophobic P1-prime preference: Phe/Ile/Leu/Tyr",
          "score": null,
          "source_bond": 9
        },
        {
          "bond": 10,
          "status": "not_matched",
          "reason": "Selected hydrophobic P1-prime preference: Phe/Ile/Leu/Tyr",
          "score": null,
          "source_bond": 10
        },
        {
          "bond": 11,
          "status": "matched",
          "reason": "Selected hydrophobic P1-prime preference: Phe/Ile/Leu/Tyr",
          "score": null,
          "source_bond": 11
        },
        {
          "bond": 12,
          "status": "not_matched",
          "reason": "Selected hydrophobic P1-prime preference: Phe/Ile/Leu/Tyr",
          "score": null,
          "source_bond": 12
        },
        {
          "bond": 13,
          "status": "not_matched",
          "reason": "Selected hydrophobic P1-prime preference: Phe/Ile/Leu/Tyr",
          "score": null,
          "source_bond": 13
        },
        {
          "bond": 14,
          "status": "not_matched",
          "reason": "Selected hydrophobic P1-prime preference: Phe/Ile/Leu/Tyr",
          "score": null,
          "source_bond": 14
        }
      ],
      "unsupported_reason": null,
      "conditions": [],
      "substrate_observation": null
    }
  ]
}
```

## Choose an extracellular context

The cleavage tool labels this MME model `extracellular` in its
`compartments` metadata. That is an enzyme-location annotation. The PDF
manifest uses a more specific **vaccine exposure route**:

| Report `context` | What it represents for an SLP |
|---|---|
| `extracellular_interstitial` | Free peptide around the injection site and other interstitial spaces before uptake; the simple IM SLP example uses this route. |
| `tumor_stroma` | Extracellular exposure in tumor tissue. For an SLP, declare `"exposures": ["tumor_stroma"]` only if that exposure is part of the scenario. |
| `circulation` | Free peptide in blood/plasma. For an SLP, declare `"exposures": ["circulation"]` only if systemic exposure is part of the scenario; IM injection alone does not establish it. |

The other report contexts—`endolysosomal`, `cytosolic_proteasome`, and
`er_trimming`—are intracellular. There is currently no separate
extracellular lymphatic-transit or draining-node context. Neither
`extracellular_interstitial` nor `circulation` measures peptide arrival at
lymph-node APCs; see the [delivery evidence review](research/vaccine-node-delivery.md).
A model compatible with an extracellular location still does not establish
that its enzyme is present or active along a particular vaccine route.

## Read the result

In this SLP example, the target spans residues 4–12. All five MME motif
matches are *inside* the target, so they flag possible damage during
extracellular exposure. The rule did not match at boundary bonds 3 and 12;
those non-matches do not prove the target is protected or cannot be released.
The report cannot tell us how much intact target survives, whether it enters
an APC, or whether it loads onto MHC. Coordinates refer to the complete
sequence supplied for that construct.

| Output file | What to inspect |
|---|---|
| `vaccine-processing-report.pdf` | Sequence map, intended target, MHC windows, and cleavage tracks grouped by route. |
| `target-summary.json` | Per-target flags, assessed and unassessed bonds, and missing expression, uptake, loading, or serum evidence. |
| `placement-assessments.json` | Cuts at target boundaries versus cuts within the target, separately for each track. |
| `processing-routes.csv` | Which processing contexts are primary, conditional, or not applicable for the declared delivery. These labels are not probabilities. |
| `selected-mhc-windows.csv` | Supplied MHC windows drawn in the PDF, if any. `display_rank` identifies each bar; a percentile rank is a predictor output, not a presentation probability. |
| `normalized-input.json`, `SHA256SUMS.json` | Exact inputs used and checksums for the run. |

A numerical site flag means a score crossed *that model's* declared threshold;
a categorical flag records a motif match or source-reported site. Neither
establishes that a cut occurs in vivo. Native scores from different models
cannot be multiplied or compared as calibrated probabilities. An unassessed
bond or missing route is not evidence of survival.

## Prepare a manifest

The [larger SLP and RNA example](vaccine-report-example.json) shows both
delivery types, intended epitopes, MHC windows, and per-bond tracks. In that
example, the RNA construct has an initiator methionine, so its target is at
residues 10–18 rather than the SLP's residues 9–17. Supply the
full delivered peptide sequence for an SLP. For `rna_encoded`, supply the
**amino-acid sequence of the translated construct**, including any relevant
signal peptide, linkers, and junctions; this field does not accept nucleotide
sequence. Declare `routing` as `cytosolic`, `secreted`, or
`lysosomal_targeted`. The example's added initial methionine is a schematic
translation product, not a validated vaccine ORF.

Target and MHC-window `start`/`end` coordinates are **one-based, inclusive
residue positions**. A hand-written numerical track has one value per internal
bond, so a sequence of length 23 needs 22 values. Bond 8 is between residues
8 and 9; `null` means that bond was not numerically assessed. Do not add a
terminal-bond value. Imported categorical results keep their `matched`,
`not_matched`, or `reported` site statuses in canonical evidence even though
their numerical scores are all `null`.
For a hand-written numerical track, supply its model `name`, route `context`,
`evidence_type`, native score `units`, provenance, and an optional threshold.
For a `cleavage_json` reference, the tool output supplies the model fields and
per-bond evidence; only `context` must be declared in the manifest. The
referenced JSON must come from `mhctools cleavage --out` for the complete
construct sequence, in the same directory as the manifest or at an absolute
path. Native scores retain their own scale and version rather than becoming
presentation probabilities.

For an injected SLP, extracellular exposure and uptake are relevant;
endolysosomal processing is a primary display context. Cytosolic
cross-presentation is conditional. For a cytosolic RNA construct, endogenous
proteasome/ER processing is primary; endolysosomal access is conditional.
Secreted or lysosome-targeted constructs need their declared route. The
report does not assume that IM injection puts free peptide in blood: declare
`"exposures": ["circulation"]` only when appropriate. These are route
assumptions for displaying evidence, not measured route weights. See the
[primary evidence on mRNA antigen presentation](evidence/mrna-antigen-presentation.md)
and [SLP cross-presentation in human dendritic cells](https://pmc.ncbi.nlm.nih.gov/articles/PMC3937416/).

## Python and Vaxrank integration

```python
from mhctools import VaccineReportInput, generate_vaccine_report

request = VaccineReportInput.from_json("vaccine-report.json")
run_directory = generate_vaccine_report(request, "reports")
```

When assembling a manifest from Vaxrank, use each final translated construct,
including linkers and junctions, so target and cut coordinates refer to the
sequence actually delivered. `placement_assessments` exposes boundary and
internal cuts separately. It does not produce an overall placement score or
change Vaxrank ranking on its own. Vaxrank's own PDF reports use its separately
installed PDF backends; this mhctools report needs no extra install group.
