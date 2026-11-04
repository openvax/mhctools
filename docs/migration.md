# Migration guide

Old spellings that still work, and what replaces them.

## MixMHCpred FASTA handling (4.6.1)

`MixMHCpred.predict_allele_sequences()` now uses Biotite's FASTA reader and
writer for allele inputs and native alignments. Full headers, record order,
wrapped sequence lines and alignment gaps are preserved. Semicolon comment
lines are accepted and ignored. Sequence data before the first header,
empty inputs, duplicate/empty identifiers and empty sequences are rejected
before running MixMHCpred. The private `_read_fasta` and `_write_fasta`
helpers were removed.

## Sequence utilities and normal installation (4.6.0)

Default SLP mass-to-molecule conversion now uses Biotite's
`ProteinSequence.get_molecular_weight(monoisotopic=False)`. It still assumes
a free, linear, unmodified canonical L-peptide, includes one terminal water,
and converts 1 microgram to `6.02214076e17 / mass_g_per_mol` molecules.
Unknown/noncanonical residues are rejected, modified chemistry remains
unassessed, and an explicitly supplied mass is preserved. The dose `source`
records the Biotite version; RNA still needs its own supplied mass or dose.

This is an average isotopic mass, not a monoisotopic mass. Biotite uses
[ExPASy residue masses](https://web.expasy.org/findmod/findmod_masses.html)
and [adds terminal water](https://www.biotite-python.org/latest/apidoc/biotite.sequence.ProteinSequence.html#biotite.sequence.ProteinSequence.get_molecular_weight).
Its water constant is rounded to 18.015 Da (ExPASy: 18.01524 Da).
The residue tables differ slightly from Biopython 1.88: `ACDEFGHIKLMN`
changes from 1377.588 to 1377.5956 g/mol, and the sequence containing each
canonical residue once changes from 2395.7134 to 2395.7356 g/mol.
The direction depends on composition: a 100-cysteine chain changes by
about -40 ppm, whereas a 100-phenylalanine chain changes by about +19 ppm.
Derived molecule counts change inversely; previously saved results are
not rewritten. Independent regression values use the ExPASy table plus
water with a 0.0003 Da tolerance for that terminal-water rounding.

Biotite, ReportLab, WEIRDO, Pepsickle and PyTorch are normal requirements.
Biotite 2.0 requires Python 3.12 or later. The requirement `biotite>=1.6,<3`
lets Python 3.11 use the newest compatible stable release (1.6.0 at this
release); newer Python versions use 2.0. Both expose the same mass and FASTA
APIs used here. The Python 3.11–3.14 CI matrix validates this boundary.
The `weirdo`, `pepsickle`, `vaccine-report` and `caphla` extra names remain
accepted compatibility aliases and add no dependencies. No model download
or license acceptance happens during installation.

| Retained extra/runtime | Reason |
|---|---|
| `nettcr` | Its native LiteRT interpreter has a separate platform/Python wheel support boundary; TensorFlow is an alternative runtime. Licensed model files still require explicit fetching. |
| `mixtcrpred` | Torchvision and Lightning are specific to this separately licensed model runtime; inference uses reviewed checkpoint compatibility handling. |
| Pepsickle gradient-boosted models | Historical scikit-learn 0.23.2 artifacts require the isolated legacy runtime; the current host neural models install normally. |
| NetCleave, TULIP and other sidecars | Their pinned framework/artifact environments remain isolated as documented under backend setup. |
| `docs`, `dev` | Documentation and development tools are not prediction dependencies. |

mhctools no longer directly imports or requires Biopython. Pepsickle still
requires it, and the isolated NetCleave runtime uses it: a normal environment
therefore still contains Biopython. This is not a claim of removing it from
the installed dependency graph.

## NetMHCpan 4.1 coordinates (4.5.4)

`NetMHCpan41` and its BA/EL variants now return zero-based offsets from the
legacy `predict_subsequences()` API, matching the 4.0/4.2 wrappers and Python
sequence slicing. Previously these offsets were one residue too large.
Scores and allele identities are unchanged; the modern `predict_proteins()`
path already used zero-based coordinates. The duplicate, misnamed internal
`parse_netmhc41_stdout` function was removed; use `parse_netmhcpan41_stdout`
for the NetMHCpan 4.1 output format.

## Prediction types

| Old | New |
|---|---|
| `BindingPrediction`, `BindingPredictionCollection` | `Prediction` and `PeptideResult` (see [results](results.md)) |
| `Pred`, `PeptidePreds` | Aliases of `Prediction` and `PeptideResult`; use the new names |
| `predict_peptides(peptides)` | `predict(peptides)` |
| `predict_subsequences(sequences, peptide_lengths)` | `predict_proteins(sequences, peptide_lengths=...)` |
| `predict_peptides_dataframe()` | `predict_dataframe()`; now returns the canonical `mhctools.pred.COLUMNS` schema rather than the legacy columns (`prediction_method_name`, `length`, ...) |

Convert legacy results with `collection.to_preds()` (a list of `Prediction`) or
`collection.to_peptide_preds()` (a list of `PeptideResult`). The standalone
`to_peptide_preds()` conversion groups by `(peptide, offset,
source_sequence_name)` because it has no input list, so unlike `predict()` it
does not preserve repeated input peptides.

```python
collection = predictor.predict_subsequences({"1L2Y": "NLYIQWLKDGGPSSGRPPPS"},
                                            peptide_lengths=[9])
df = collection.to_dataframe()
preds = collection.to_preds()
```

## Predictor names

| Old | New |
|---|---|
| `IedbNetMHCpan`, `IedbNetMHCcons`, `IedbNetMHCIIpan`, `IedbSMM`, `IedbSMM_PMBEC` and their `*-iedb` CLI names | Still work as local wrappers. Prefer the explicit local names (`NetMHCpan41_BA`, `NetMHCcons`, `NetMHCIIpan43_BA`, `SMM`, `SMMPMBEC`) when recording provenance. See [compatibility names](predictors/binding.md#compatibility-names-for-the-old-iedb-predictors). |
| IEDB HTTP constructor arguments `url`, `request_timeout`, `raise_on_error` and the CLI option `--do-not-raise-on-error` | Removed. There is no HTTP fallback; missing installations and unsupported inputs raise errors. |

## Commands and options

| Old | New |
|---|---|
| `mhctools integrations` | `mhctools predictors` (the old spelling is an alias) |
| `--mhc-epitope-lengths` | `--mhc-peptide-lengths` (deprecated name still accepted) |

## Kinds

`serum_half_life`, `plasma_half_life`, `blood_half_life` and
`systemic_elimination_half_life` are accepted as input and canonicalized to
`peptide_half_life`, with the matrix kept in the [measurement
context](kinds.md#measurement-context). `canonical_kind()` performs the
mapping.
