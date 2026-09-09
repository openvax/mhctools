# Peptidase cleavage evidence

```python
from mhctools import CleavageInput, DPP4qPISA

peptide = CleavageInput("HAEGTFTSD", source_id="GLP-1 fragment")
result = DPP4qPISA().predict(peptide)
print(result.sites[0].bond)   # 2: HA | EGTFTSD
print(result.sites[0].score)  # 2.1694, native qPISA score
print(result.to_dict())      # input, model, assay, limitations and source bonds
```

This API describes individual peptide bonds. Bond `b` splits
`sequence[:b] | sequence[b:]`; `source_bond` adds the zero-based
`source_start` offset. Existing proteasomal scoring APIs remain available.

Inputs describe canonical, linear L-peptides. Strings assume free N and C
termini. Use `CleavageInput` to explicitly record `n_term="acetylated"` or
`c_term="amidated"`, or `unknown`. These forms are recorded but the initial
qPISA implementation abstains because they are outside its supported domain.
Other modifications, D-residues, cyclization and conjugates are unsupported;
do not strip their chemistry to obtain a score.

Exopeptidases assess exposed termini. An internal `HAE` is not an immediate
DPP4 site. To ask about a hypothesized product, use
`parent.fragment(start, end, n_term="free", c_term="free")` and predict that
fragment. The source offset is retained. This conditional analysis does not
predict formation of the fragment, cleavage order, rates or competition.

## Human DPP4 qPISA

The [Gudipati et al. 2024 paper](https://doi.org/10.1038/s44320-024-00071-4)
reports a model fitted to substrate depletion by purified human DPP4. The
source assay used tryptic HeLa peptides in HEPES pH 7.4, at 21 C for 4 hours.
The implementation independently evaluates the three rearranged terms in
Dataset EV2: `P1 + P2:P1 + P1:P1-prime`, where the first three peptide
residues are P2, P1 and P1-prime. Only bond 2 is assessed.

Higher scores indicate greater predicted log2 depletion relative to buffer
control in that assay. Negative values are retained. Scores are **not
cleavage probabilities, serum half-lives or stability ranks across enzymes**.
Compartment metadata records where the enzyme can act, not where the model
has been calibrated. Physiological exposure, structure, concentration and
competing enzymes are not modeled.

All triplets with complete coefficients can be evaluated, including those
without Pro or Ala at P1. Of 8,000 canonical triplets, 6,420 have complete
coefficients; the other 1,580 return an explicit missing-coefficient reason.
Evaluability does not establish that an individual triplet occurred in the
training data. The related C. elegans DPF-3 model is not used for human DPP8/9.

### Parameter provenance

`mhctools/data/dpp4_qpisa.json` contains the numeric cells from
`44320_2024_71_MOESM3_ESM.xlsx`, sheet `dpp4_modelParams`. Published `NA`
cells become JSON `null`; no numerical imputation or refitting is performed.
The source workbook SHA-256 is
`ee449da13b5ec66fd6fb08c16203da8c44f2e9c10572a355abdcb985694c6ddc`.
It was retrieved via the [Europe PMC supplementary archive](https://www.ebi.ac.uk/europepmc/webservices/rest/PMC11612144/supplementaryFiles).
The article assigns associated data to [CC0](https://creativecommons.org/publicdomain/zero/1.0/),
unless otherwise credited; Dataset EV2 has no separate credit restriction.
Attribution: Rajani Kanth Gudipati and colleagues, 2024, DOI above. No figures
or upstream R source code are redistributed.

## Interpreting rule-based evidence

The shared result contract also supports `motif_rule` evidence, reported as
`matched` or `not_matched` with no numerical score. A rule match indicates a
known recognition pattern; a non-match does not establish resistance.
`unsupported_reason` means the input could not be assessed and carries no
score. These states must remain separate in downstream displays and ranking.
