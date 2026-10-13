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

Install the PDF dependency and download the [example manifest](vaccine-report-example.json),
then run:

```sh
pip install "mhctools[vaccine-report]"
curl -L -o vaccine-report-example.json \
  https://openvax.github.io/mhctools/vaccine-report-example.json
mhctools vaccine-report \
  --input vaccine-report-example.json \
  --output-dir reports
```

The command prints a timestamped directory inside `reports/`. Open
`vaccine-processing-report.pdf` there. The example compares an injected SLP
with a cytosolic RNA construct containing the same nine-residue target,
`KIILNFTTL`. **Its cleavage scores and MHC ranks are invented to demonstrate
report behavior; they are not predictions for that sequence.**

## Read the result

In the SLP example, the target spans residues 9–17. The endolysosomal track
flags a possible cut after residue 13, *inside* the target. It also flags the
bond after residue 8, immediately *before* the target. The first flag suggests
possible destruction; the second suggests possible release. Neither tells us
which event wins, how much intact target survives, or whether it loads onto MHC.
The RNA construct has an initiator methionine, so the same target is at
residues 10–18. Coordinates always refer to the complete sequence supplied
for that construct.

| Output file | What to inspect |
|---|---|
| `vaccine-processing-report.pdf` | Sequence map, intended target, MHC windows, and cleavage tracks grouped by route. |
| `target-summary.json` | Per-target flags, assessed and unassessed bonds, and missing expression, uptake, loading, or serum evidence. |
| `placement-assessments.json` | Cuts at target boundaries versus cuts within the target, separately for each track. |
| `processing-routes.csv` | Which processing contexts are primary, conditional, or not applicable for the declared delivery. These labels are not probabilities. |
| `selected-mhc-windows.csv` | MHC windows drawn in the PDF. `display_rank` identifies each bar; the supplied percentile rank is a predictor output, not a presentation probability. |
| `normalized-input.json`, `SHA256SUMS.json` | Exact inputs used and checksums for the run. |

A supported cut means a track crossed *its own* declared threshold. Native
scores from different models cannot be multiplied or compared as calibrated
probabilities. An unassessed bond or missing route is not evidence of survival.

## Prepare a manifest

The [complete example manifest](vaccine-report-example.json) shows both
delivery types, intended epitopes, MHC windows, and per-bond tracks. Supply the
full delivered peptide sequence for an SLP. For `rna_encoded`, supply the
**amino-acid sequence of the translated construct**, including any relevant
signal peptide, linkers, and junctions; this field does not accept nucleotide
sequence. Declare `routing` as `cytosolic`, `secreted`, or
`lysosomal_targeted`. The example's added initial methionine is a schematic
translation product, not a validated vaccine ORF.

Target and MHC-window `start`/`end` coordinates are **one-based, inclusive
residue positions**. A cleavage track has one value per internal bond, so a
sequence of length 23 needs 22 values. Bond 8 is between residues 8 and 9;
`null` means that bond was not assessed. Do not add a terminal-bond value.
For each track, name the model, route `context`, evidence type, units,
provenance, and threshold if it has one. Retain the model's native scale and
version instead of relabeling its score as a probability.

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

request = VaccineReportInput.from_dict(manifest)
run_directory = generate_vaccine_report(request, "reports")
```

When assembling a manifest from Vaxrank, use each final translated construct,
including linkers and junctions, so target and cut coordinates refer to the
sequence actually delivered. `placement_assessments` exposes boundary and
internal cuts separately. It does not produce an overall placement score or
change Vaxrank ranking on its own.
