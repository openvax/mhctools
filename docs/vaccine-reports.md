# Route-aware vaccine processing reports

`mhctools vaccine-report` turns a structured vaccine-construct manifest into a
timestamped PDF and machine-readable audit files. The same entry point is
available in Python as `generate_vaccine_report`.

The report distinguishes two delivery modalities:

- `synthetic_long_peptide`: free peptide is first exposed extracellularly and
  during uptake. Endolysosomal processing is a primary display context.
  Cytosolic proteasome and ER trimming remain a separate, conditional route
  because human and mouse dendritic-cell studies found proteasome- and
  TAP-dependent SLP cross-presentation after uptake and cytosolic export.
- `rna_encoded`: `routing` is required semantically and defaults to
  `cytosolic`. Cytosolic translation makes proteasome/ER/class-I processing
  primary. Endolysosomal/class-II processing is conditional for a cytosolic
  construct because autophagy can deliver endogenous antigen to MHC-II loading
  compartments. `secreted` and `lysosomal_targeted` make those routes explicit.

Serum, plasma, and free-peptide exopeptidase tracks are not shown for a
cytosolic RNA construct. They become relevant only when secretion or another
extracellular exposure is declared. Likewise, an injected SLP does not imply
circulation; add `"exposures": ["circulation"]` only when that route is part of
the intended experiment or formulation.

Primary route evidence:

- SLP uptake and cross-presentation in human dendritic cells:
  <https://pmc.ncbi.nlm.nih.gov/articles/PMC3937416/>
- SLP class-I and class-II processing in mouse and human dendritic cells:
  <https://pubmed.ncbi.nlm.nih.gov/23836147/>
- Endogenous antigen delivery to MHC-II compartments by autophagy:
  <https://pubmed.ncbi.nlm.nih.gov/17182262/>

These papers establish possible processing routes, not that every construct
uses each route or that a native predictor score is an in-vivo probability.

## Input schema

Coordinates are one-based and inclusive for intended epitopes and MHC windows.
Cleavage arrays contain exactly `length - 1` internal-bond values; array index
zero is bond 1 (`sequence[:1] | sequence[1:]`). `null` means unassessed. An
endpoint sentinel must never be supplied.

```json
{
  "schema_version": 1,
  "constructs": [
    {
      "id": "example-slp",
      "sequence": "RISVTPGEKIILNFTTLDLYRSR",
      "delivery": "synthetic_long_peptide",
      "exposures": ["tumor_stroma"],
      "intended_epitopes": [
        {"label": "mutant target", "start": 9, "end": 17}
      ],
      "mhc_windows": [
        {
          "mhc_class": "I",
          "allele": "HLA-A*02:01",
          "start": 9,
          "end": 17,
          "percentile_rank": 0.4,
          "predictor": "example predictor with exact version"
        }
      ],
      "cleavage_tracks": [
        {
          "name": "example proteasome model",
          "context": "cytosolic_proteasome",
          "evidence_type": "quantitative_model",
          "scores": [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 0.8,
                     0.7, 0.6, 0.5, 0.4, 0.3, 0.2, 0.1, 0.2, 0.3, 0.4,
                     0.5, 0.6],
          "threshold": 0.5,
          "units": "dimensionless native model output",
          "provenance": "package/model/weights identity"
        }
      ]
    },
    {
      "id": "example-rna",
      "sequence": "MALWMRLLPLLALLALWGPDPAAA",
      "delivery": "rna_encoded",
      "routing": "cytosolic"
    }
  ]
}
```

```sh
mhctools vaccine-report \
  --input docs/vaccine-report-example.json \
  --output-dir reports \
  --max-mhc-windows 10
```

Each run contains `vaccine-processing-report.pdf`, normalized input,
`processing-routes.csv`, `selected-mhc-windows.csv`, placement assessments,
and checksums. Selection retains the strongest intended-target overlap, then
allele diversity, then native percentile rank. All supplied windows remain in
the normalized input. Every selected window is drawn: the figure adds ligand
lanes as needed instead of dropping one, and the number on each bar is that
window's `display_rank` in `selected-mhc-windows.csv`.

## Python and Vaxrank integration

```python
from mhctools import VaccineReportInput, generate_vaccine_report

request = VaccineReportInput.from_dict(manifest)
run_directory = generate_vaccine_report(request, "reports")
```

Vaxrank can build one `VaccineConstruct` for each final translated construct,
including linkers and junctions, and call `placement_assessments` before any
ranking change. The result reports internal cuts (candidate target disruption)
and boundary cuts (candidate liberation) separately for every model and route.
It deliberately has no aggregate “good placement” number. A future Vaxrank
policy can compare placements only after choosing which routes, models, and
failure semantics are admissible; unassessed bonds must remain unassessed.
