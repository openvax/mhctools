# Assay-aware benchmarks

```sh
mhctools benchmark --reference-cleavage --out reference-report.json
mhctools benchmark --reference-cleavage serum --out serum-reference-report.json
mhctools benchmark --lineage-inventory --out lineage.json
mhctools benchmark --input measurements-and-predictions.json --out evaluation.json
mhctools benchmark --input site-measurements.json --model cpn-basic --model dpp9-xp-xa
```

The benchmark accepts JSON records represented by `AssayMeasurement`,
`BenchmarkPrediction` and `ModelLineage` in `mhctools.benchmark`. Each
measurement retains its own ID, original source measurement ID, source URL,
study, assay, exact sequence and chemistry, endpoint, native units, species,
matrix, optional cell type, censoring and experimental conditions. Site
records also require an enzyme and a peptide bond. No missing bond label
is converted into an experimental non-cleavage.

The input object has `measurements`, optional `predictions`, optional
`lineages`, and optional `requested_domains` lists. To run a cleavage model
use `--model` instead of supplying predictions. A minimal supplied evaluation:

```json
{
  "measurements": [{
    "measurement_id": "experiment-1",
    "source_measurement_id": "table-2-row-7",
    "source": "https://example.org/replace-with-primary-source",
    "dataset": "my-data", "study": "study-1", "assay": "incubation-1",
    "sequence": "HAEGT", "chemistry": "linear_L_free",
    "endpoint": "serum_half_life", "units": "hours",
    "species": "Homo sapiens", "matrix": "serum", "value": 2.0,
    "family": "GLP1", "conditions": {"temperature_C": "37"}
  }],
  "predictions": [{
    "measurement_id": "experiment-1", "model": "my-model",
    "endpoint": "serum_half_life", "units": "hours", "value": 3.0
  }],
  "requested_domains": [{
    "name": "primary human DC cytosolic delivery",
    "species": "Homo sapiens", "cell_type": "primary dendritic cell",
    "endpoint": "cytosolic_delivery"
  }]
}
```

These numbers and URL are illustrative fixtures, not experimental data.
Use explicit `unknown` descriptions for unspecified conditions and `null`
for an unknown measurement value. Censored observations use `less_than`,
`greater_than` or `unknown`; they remain visible but are excluded from ordinary
point-error metrics. Chemistry descriptions are preserved verbatim. The
cleavage runner can currently represent only `linear_L_free`,
`linear_L_N_acetylated` and `linear_L_C_amidated`; other descriptions abstain.

## Metrics and interpretation

Reports stratify by model, endpoint, native units/output scale, study, assay,
conditions, species, matrix, cell type and enzyme. Serum, plasma, whole blood,
systemic half-life/clearance, fluorescence uptake, CPP classification,
cytosolic delivery and antigen presentation remain distinct endpoints.
There is no implicit unit conversion or transformed-score calibration.

Compatible native regression outputs receive MAE, RMSE and signed bias.
Binary decisions receive confusion counts; explicit probabilities receive a
Brier score. Motif decisions are limited recognition rules, so their confusion
counts measure that rule's behavior on supplied labels. Native qPISA/PWM
scores cannot acquire probability metrics merely because the observed label
is binary. Intervals are evaluated only when explicitly identified as
individual-measurement prediction intervals in the native scale, with a
stated nominal coverage. Model disagreement is not an uncertainty interval.

Counts distinguish measurements, source measurements, unique sequences and
chemical forms. Missing predictions, unsupported chemistry, runtime failures,
unassessed sites and incompatible/censored records retain their reasons.
The optional `serum` panel adds 17 source observations, including two experimental
non-cleavages and two unsupported ACE/amide examples; its activation conditions
are applied separately to each CPB2 measurement.
The source-linked starter set reproduces four positive observations (two
CPN-like plasma observations for bradykinin, aminopeptidase P on RPP and DPP9
on the RU1 epitope). Eight model/observation combinations are off-enzyme and
unassessed. This small, positive-only set cannot estimate specificity,
calibration or independent predictive performance. CPN attribution is to
CPN-like inhibitor-sensitive plasma activity, not an isolated purified enzyme.

## Training overlap and applicability

`ModelLineage` records dataset, study, assay, sequence, exact chemical-form,
family and cell-type inventories. Provenance is `unknown`, `partial` or
`complete`; complete provenance requires actual sequence, study and family
inventories. Every evaluation record is checked against these inventories
and against any supplied `split="train"` partition. Related chemistry does
not erase sequence overlap. Family assignments are supplied explicitly; no
arbitrary sequence-similarity threshold or confidence score is invented.

The report labels an external evaluation unverified when training provenance,
family assignments or independence cannot be established. `split="reference"`
records are used only for `--evaluation reproduction`; the built-in reference
command always selects reproduction. Metrics remain descriptive, and model
agreement or reproduction is never represented as independent validation.

Requested domains report available measurements, comparable predictions and
missing evidence. A length restriction is a user-selected target population,
not a learned applicability threshold. The supplied reference set reports no
human serum half-life or primary-DC cytosolic-delivery observations.

## Dataset/model lineage inventory

The installed [inventory](../mhctools/data/model_lineage.json) records sources,
known relationships and unresolved questions rather than assuming models
have independent training data:

- [Tan 2024](https://doi.org/10.1093/bib/bbae350) draws on PEPlife, PepTherDia,
  THPdb and literature. Its species/organ/chemistry subsets require original
  assay review before distinguishing incubation half-life from systemic PK.
  Complete raw datasets are described as available on request.
- [PEPlife2](https://www.mdpi.com/2673-5601/6/2/26) updates half-life curation
  and retains repeated sequences measured under different conditions.
- [pepADMET](https://pubmed.ncbi.nlm.nih.gov/41512092/) spans multiple ADMET
  endpoints. Its exact half-life row overlap with Tan/PEPlife is unresolved;
  shared authorship does not establish either identity or independence.
- [POSEIDON](https://doi.org/10.1186/s13321-024-00810-7) revisits CPPsite2.0
  primary studies and retains repeated peptides under different conditions.
  Its regression endpoint is fluorescence uptake.
- [PerseuCPP](https://doi.org/10.1093/bioadv/vbaf213) reuses MLCPP2.0 training
  sets and adds CPPsite2.0 data. Some negative CPP labels are generated rather
  than experimentally observed. This shared ancestry with POSEIDON requires
  row-level auditing before claiming an independent comparison.
- qPISA coefficient reproduction and ERAMER adapter agreement are documented
  in the [cleavage guide](cleavage.md); their reproduction is separate from
  new assay validation.

No third-party raw exposure/uptake datasets or model weights are redistributed.
The included small factual cleavage curation is offered under CC0 with source
measurement identifiers and explicit caveats. Broad external validation on
long vaccine peptides or primary human DCs remains an empirical data task;
this report makes missing evidence visible and supplies the evaluation path.
