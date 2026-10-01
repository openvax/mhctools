# Recipes

Short answers to common tasks. Each assumes you have a predictor built as in the
[quickstart](index.md#quickstart).

## Scan proteins instead of peptides

`predict_proteins()` takes a dictionary of sequences and returns
`{sequence_name: list[PeptideResult]}`, with each result's `offset` set:

```python
proteins = predictor.predict_proteins(
    {"TP53": "MEEPQSDPSVEPPLSQETFS...", "KRAS": "MTEYKLVVVGAGGVGKS..."},
    peptide_lengths=[9, 10],
)

for r in proteins["TP53"]:
    if r.affinity and r.affinity.value < 500:
        print(f"  offset={r.offset} {r.peptide} IC50={r.affinity.value:.0f}")
```

## Run many samples with different genotypes

```python
from mhctools import MultiSample, MHCflurry

ms = MultiSample(
    samples={
        "pat001": ["HLA-A*02:01", "HLA-B*07:02"],
        "pat002": ["HLA-A*01:01", "HLA-B*08:01"],
    },
    predictor_class=MHCflurry,
)

results = ms.predict(["SIINFEKL", "GILGFVFTL"])       # {sample: [PeptideResult]}
protein_results = ms.predict_proteins({"TP53": "MEEPQ..."})  # {sample: {seq: [...]}}

df = ms.predict_dataframe(["SIINFEKL"])               # flat, with sample_name
df = ms.predict_proteins_dataframe({"TP53": "MEEPQ..."})
```

## Add predictor scores to an existing table

Evaluation workflows often start from an annotated benchmark table — columns
like `sample_id`, `hit`, `peptide`, and a per-row genotype — and just need
scores appended. `annotate_table` is I/O-free and works on any `DataFrame`:

```python
from mhctools import annotate_table, AnnotationSpec, NetMHCpan42_BA

annotated = annotate_table(
    df,
    [AnnotationSpec(
        predictor=lambda alleles: NetMHCpan42_BA(alleles=alleles),
        output_column="netmhcpan4.2.ba",
        field="affinity")],
    peptide_column="peptide",
    allele_column="hla")
```

There's a [CLI equivalent](cli.md#annotate-a-table-predict-table) that reads and
writes CSV.
