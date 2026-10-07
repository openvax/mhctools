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

## Released CPP split audit and reproduction

The PeptiVerse command consumes a **local** copy of the exact released CPP
metadata. It verifies SHA-256 before any model is loaded. Obtain it separately;
the raw dataset is not bundled with mhctools:

```sh
curl -fL -o permeability_meta_with_split.csv \
  https://huggingface.co/ChatterjeeLab/PeptiVerse/resolve/8cf0b21dae356278ae96b414a088e4360357d16c/training_data_cleaned/permeability_penetrance/permeability_meta_with_split.csv
mhctools benchmark --peptiverse-cpp-metadata permeability_meta_with_split.csv \
  --out cpp-audit.json
```

The expected file hash is
`c924f1b92fac14f2007afb1b4b3641047896219a17a5780beb865ee0f4b35ec8`.
Audit-only mode needs no ML runtime and runs no predictions. For optional
offline source-label reproduction, provision the
[CPP runtime and weights](predictors/uptake.md#peptiverse-cpp), then run:

```sh
mhctools benchmark --peptiverse-cpp-metadata permeability_meta_with_split.csv \
  --predict-cpp --out cpp-source-reproduction.json
# A small smoke cohort, explicitly reported as partial:
mhctools benchmark --peptiverse-cpp-metadata permeability_meta_with_split.csv \
  --predict-cpp --source-id seq_6 --source-id seq_7 --out cpp-partial.json
```

All 1,859 training and 465 validation records are audited. The splits have
no duplicate sequences, conflicting labels or exact sequence overlap.
Training lengths span 3–61 residues; validation lengths span 3–52, with 56
validation records at least 30 residues long. Original assay, study, chemical
form, species, matrix, cell type and family assignments are absent from this
four-column file. The report leaves them unknown. A generic released split
script demonstrates non-fouling preprocessing; it does not establish the
actual CPP cluster assignments. Exact nonoverlap alone cannot certify
study/family independence.

Every validation ID remains in the report, including rows outside a selected
cohort, unsupported inputs and runtime failures. Probability metrics reproduce
the source CPP labels; confusion counts use the native **0.5493** threshold.
The classifier's canonical/free-terminus input assumption is reported
separately from unknown experimental chemistry. No score is interpreted as a
physical uptake fraction, and predictive entropy is not an accuracy interval.
The output includes source/evaluator hashes, mhctools version and, after
successful prediction, the adapter's asset inventory.

The [recorded full-split CPU run](https://github.com/openvax/mhctools/blob/master/tests/data/peptiverse_cpp_validation_summary.json)
scored **465/465** validation records with no failures or unsupported inputs.
At the native threshold, it produced 203 true positives, 231 true negatives,
13 false positives and 18 false negatives against the released labels;
descriptive Brier score was **0.06989**. This aggregate record includes source,
evaluator and asset hashes and runtime versions, without redistributing raw
sequences or labels. Routine CI audits the full metadata and scores three
validation records; it does not repeat the full CPU run. Small cross-platform
embedding differences can affect scores, so this record is run evidence,
not a demand for bitwise equality on another CPU stack.

This is **source-label reproduction**, always reported as such, even if
`--evaluation external_validation` is supplied. The
[publication](https://doi.org/10.1038/s41467-026-74167-w) describes selection on
validation performance; this is not a newly untouched experimental test set.
Requested-domain reports find no verified primary human DC cytosolic-delivery,
antigen-presentation or beyond-training-length (62+ residue) evidence in these
metadata. That means evidence is missing, not that these biological processes
cannot occur. Broad [#291](https://github.com/openvax/mhctools/issues/291) and
[#302](https://github.com/openvax/mhctools/issues/302) validation work remains open.

## Dataset/model lineage inventory

The installed [inventory](https://github.com/openvax/mhctools/blob/master/mhctools/data/model_lineage.json) records sources,
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
- [PeptiVerse](https://doi.org/10.1038/s41467-026-74167-w) releases CPP training
  and validation sequence membership. Its missing assay/chemical-form/family
  metadata limits the audit above to source-label reproduction.
- qPISA coefficient reproduction and [ERAMER](predictors/processing.md#eramer) adapter agreement are documented
  in the [peptidase activity guide](cleavage/index.md); their reproduction is separate from
  new assay validation.

No third-party raw exposure/uptake datasets or model weights are redistributed.
The included small factual cleavage curation is offered under CC0 with source
measurement identifiers and explicit caveats. Broad external validation on
long vaccine peptides or primary human DCs remains an empirical data task;
this report makes missing evidence visible and supplies the evaluation path.
