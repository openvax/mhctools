# API reference

Use this reference to look up constructors, fields, and methods. For a complete
first prediction, start with [getting started](getting-started.md).

| Task | Reference | Explanation and examples |
|---|---|---|
| Read predictions and scores | [Results](#results) | [Results and DataFrames](results.md) |
| Describe a peptide or TCR | [Inputs](#inputs) | [Input shapes](predictors/index.md#input-shapes) |
| Add predictions to a table | [Table annotation](#annotating-tables) | [Annotation recipe](recipes.md) |
| Download or inspect models | [Models and artifacts](#models-and-artifacts) | [Getting models](artifacts.md) |
| Assess individual peptide bonds | [Peptidase activity](#peptidase-activity) | [Activity guide](cleavage/index.md) |

Predictor-specific setup and examples are in the [family guides](predictors/index.md).
The [predictor matrix](predictor-matrix.md) lists Python classes and command-line
names for every supported model.

## Results

`predict()` returns one [`PeptideResult`][mhctools.PeptideResult] per peptide.
Its predictions are [`Prediction`][mhctools.Prediction] objects, each describing
one endpoint and MHC or TCR context. [`Kind`][mhctools.Kind] names the endpoint;
[`MeasurementContext`][mhctools.MeasurementContext] records its measurement
semantics. [`MultiSample`][mhctools.MultiSample] runs a predictor for several
named MHC allele sets.

Read [results and DataFrames](results.md) for examples and
[prediction kinds](kinds.md) for score meanings and units.
Source: [prediction types](https://github.com/openvax/mhctools/blob/master/mhctools/pred.py)
and [sample helper](https://github.com/openvax/mhctools/blob/master/mhctools/sample.py).

::: mhctools.PeptideResult

::: mhctools.Prediction

::: mhctools.Kind

::: mhctools.MeasurementContext

::: mhctools.MultiSample

## Inputs

[`TCR`][mhctools.TCR] carries receptor sequences and gene identifiers.
[`PeptideInput`][mhctools.PeptideInput] records a peptide's sequence, terminal
chemistry, attachments, and source occurrence.
[`PeptideContext`][mhctools.PeptideContext] adds study, formulation, and assay
context when those details are known.

See [input shapes](predictors/index.md#input-shapes), the
[TCR guide](predictors/tcr.md), and [peptide exposure inputs](exposure-results.md).
Source: [TCR types](https://github.com/openvax/mhctools/blob/master/mhctools/tcr.py)
and [peptide types](https://github.com/openvax/mhctools/blob/master/mhctools/peptide_input.py).

::: mhctools.TCR

::: mhctools.PeptideInput

::: mhctools.PeptideContext

## Annotating tables

[`annotate_table()`][mhctools.annotate_table] adds predictions to a DataFrame.
An [`AnnotationSpec`][mhctools.AnnotationSpec] selects the model, endpoint, and
filter to apply. Start with the [annotation recipe](recipes.md) for a complete
example.

Source: [table annotation](https://github.com/openvax/mhctools/blob/master/mhctools/annotate.py).

::: mhctools.annotate_table

::: mhctools.AnnotationSpec

## Models and artifacts

[`fetch()`][mhctools.fetch] installs managed model artifacts.
[`list_artifacts()`][mhctools.list_artifacts] reports their installation status;
[`integration_status()`][mhctools.integration_status] checks whether a predictor
can run. See [getting models](artifacts.md) for the installation workflow and
[licensing](licensing.md) for upstream terms.

Source: [artifact management](https://github.com/openvax/mhctools/blob/master/mhctools/artifacts.py)
and [integration checks](https://github.com/openvax/mhctools/blob/master/mhctools/integrations.py).

::: mhctools.fetch

::: mhctools.list_artifacts

::: mhctools.integration_status

<a id="cleavage"></a>

## Peptidase activity

[`predict_cleavage()`][mhctools.predict_cleavage] assesses individual bonds for
selected peptidases. [`predict_cleavage_batch()`][mhctools.predict_cleavage_batch]
assesses a set of inputs and processing scenarios.
[`CleavageInput`][mhctools.CleavageInput] records sequence and terminal chemistry.
The Python names retain “cleavage”; the biological workflow is explained in the
[peptidase activity guide](cleavage/index.md).

See [choosing processing models](cleavage/choosing.md),
[reading the evidence](cleavage/evidence.md), and [batch assessments](cleavage/batch.md).
Source: [model dispatch](https://github.com/openvax/mhctools/blob/master/mhctools/peptidases.py),
[input and result types](https://github.com/openvax/mhctools/blob/master/mhctools/cleavage.py),
and [batch assessment](https://github.com/openvax/mhctools/blob/master/mhctools/cleavage_batch.py).

::: mhctools.predict_cleavage

::: mhctools.predict_cleavage_batch

::: mhctools.CleavageInput
