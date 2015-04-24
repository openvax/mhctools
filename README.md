# mhctools
Python interface to running command-line and web-based MHC binding predictors

## Example

```
from mhctools import NetMHCpan
# Run NetMHCpan for alleles HLA-A*01:01 and HLA-A*02:01
predictor = NetMHCpan(alleles=["A*02:01", "hla-a0101"])
# scan the short proteins 1L2Y and 1L3Y for epitopes
protein_sequences = (
  {"1L2Y": "NLYIQWLKDGGPSSGRPPPS",
  "1L3Y": "ECDTINCERYNGQVCGGPGRGLCFCGKCRCHPGFEGSACQA"
}
epitope_collection = predictor.predict(protein_sequences)
# flatten binding predictions into a Pandas DataFrame
df = epitope_collection.dataframe()
# epitope collection is sorted by percentile rank
# of binding predictions
strongest_predicted_binder = epitope_collection[0]
```
