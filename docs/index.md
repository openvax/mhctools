# mhctools

mhctools is a Python library for running MHC binding, presentation,
immunogenicity, and antigen-processing predictors. It provides a common
interface to tools such as [NetMHCpan](predictors/binding.md#netmhcpan) and [MHCflurry](predictors/binding.md#mhcflurry), with results you can inspect
in Python or export as a pandas DataFrame.

## Quickstart

Install mhctools and download the [MHCflurry](predictors/binding.md#mhcflurry) model weights:

```sh
pip install mhctools
mhctools fetch mhcflurry
```

```python
from mhctools import MHCflurry

predictor = MHCflurry(alleles=["HLA-A*02:01"])
results = predictor.predict(["SIINFEKL", "GILGFVFTL"])
df = predictor.predict_dataframe(["SIINFEKL", "GILGFVFTL"])
```

The [getting started guide](getting-started.md) explains the results and shows
how to use another predictor.

## Choose what to predict

| Question | Predictors |
|---|---|
| Class I binding affinity | [NetMHCpan](predictors/binding.md#netmhcpan) or [MHCflurry](predictors/binding.md#mhcflurry) |
| Class I presentation | [NetMHCpan](predictors/binding.md#netmhcpan), [MHCflurry](predictors/binding.md#mhcflurry), [MixMHCpred](predictors/binding.md#mixmhcpred), [BigMHC (EL)](predictors/binding.md#bigmhc), [CapHLA](predictors/binding.md#caphla) |
| Class II binding or presentation | [NetMHCIIpan](predictors/binding.md#netmhciipan); [MixMHC2pred](predictors/binding.md#mixmhc2pred) for presentation |
| Peptide-MHC complex stability | [NetMHCstabpan](predictors/binding.md#netmhcstabpan) |
| Proteasomal cleavage | [Pepsickle](predictors/processing.md#pepsickle) or [NetChop](predictors/processing.md#netchop) |
| Class II cleavage | [NetCleave (class II)](predictors/processing.md#netcleave) |
| TAP transport | [DeepTAP](predictors/processing.md#deeptap) |
| ERAP1 trimming | [ERAMER](predictors/processing.md#eramer) |
| T-cell immunogenicity | [Calis](predictors/immunogenicity.md#calis), [PRIME](predictors/immunogenicity.md#prime), [BigMHC (IM)](predictors/binding.md#bigmhc), [DeepImmuno](predictors/immunogenicity.md#deepimmuno); [TLimmuno2](predictors/immunogenicity.md#tlimmuno2) for class II |
| Recognition by a specific TCR | [NetTCR](predictors/tcr.md#nettcr), [Tulip](predictors/tcr.md#tulip), [MixTCRpred](predictors/tcr.md#mixtcrpred) |
| Free-peptide half-life | [PeptiVerse](predictors/peptide-pk.md#peptiverse), [PlifePred2](predictors/peptide-pk.md#plifepred2) |
| Per-bond peptidase evidence | [Peptidase activity](cleavage/index.md) |

The [selection guide](choosing.md) explains input, installation, and license
constraints. Read each model's guide for its output and validation limits.

<a id="documentation"></a>

## Start here

Use the [recipes](recipes.md) for protein scans, multiple samples and table
annotation, or the [API reference](api.md) for class and function details.
For installation help, see [troubleshooting](troubleshooting.md).

<a id="choose-and-evaluate-models"></a>

Read [known limits](limitations.md) and [benchmarks](benchmarks.md) when
interpreting predictions. [Research and evidence](research/index.md) traces
primary studies, validation audits and open calibration questions.

<a id="processing-and-vaccine-analysis"></a>

For processing workflows, start with the
[antigen-processing guide](predictors/processing.md), then use
[batch assessments](cleavage/batch.md) or [vaccine reports](vaccine-reports.md).
