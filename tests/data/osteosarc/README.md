# Osteosarc vaccine regression fixtures

Source: [Sid Sijbrandij's osteosarcoma dataset](https://registry.opendata.aws/sid-osteosarc/),
accessed September 21, 2026 (CC0 1.0). Exported through
[`osteosarc.Dataset.vaccine_peptides()`](https://github.com/iskandr/osteosarc)
with osteosarc 0.1.2. These are software regression fixtures, not a vaccine
recommendation or a predictor-accuracy benchmark.

`inputs/source.json` retains all 94 API records, their variant IDs, peptide
ordinals, source links, exact sequences, vaccine membership, and original
experimental annotations. Of these, 78 belong to a disclosed vaccine and 16
are associated experimental peptides without vaccine membership. The 33
variant–vaccine memberships with no disclosed sequence are listed separately.
An identical sequence attributed to both EPG5 and EXOC4 remains two source
records; sequence deduplication retains both references.

The immutable snapshot ID is
`9c8938b30873f6dbaf3a712e806fd1b813ee94dac9364ddfa49bb494748e745d`.
Its source receipts include retrieval times, URLs, sizes and SHA-256 hashes.
The full curation report is retained: 22 corrections applied, two already
fixed upstream, and eight stale corrections left unapplied by osteosarc.
Those statuses describe this snapshot and client version, not current data.

`validation` is a **variant-level** annotation; `experiments` retains the
separate peptide-level records. Neither becomes a predicted score or an
assumed label for a derived window. In particular, an untested ELISPOT entry
is not a negative result; see the [source's explanation](https://osteosarc.com/vaccines/).

## Predictor inputs and outputs

`inputs/panels.json` links every predictor input to its source record(s) with
zero-based offsets. The matching FASTA and plain-text files are the exact
inputs used for recording. Alleles are regression-test parameters; no allele
restriction or clinical suitability is inferred for any peptide.

| Capture | Input panel | Native output |
| --- | --- | --- |
| NetMHCpan 4.1b | 19 disclosed 9/10-mers; scan lengths 9 and 10 | 46 rows, each with BA and EL endpoints |
| NetMHCIIpan 4.3e | 55 distinct centered 15-residue windows | 110 rows, each with BA and EL endpoints |
| NetMHCcons 1.1 | 17 disclosed 9-mers | 17 rows per allele, two captures |
| MixMHCpred 3.0 | 19 disclosed 9/10-mers | 19 rows, two alleles |
| PRIME 2.1 with MixMHCpred 3.0 | 19 disclosed 9/10-mers | 19 rows, two alleles |
| MixMHC2pred 2.1-beta1 | 55 distinct centered 15-residue windows, no context | 55 rows, two alleles |
| NetChop 3.1 distribution, C-term 3.0 / 20S 3.0 models | All 77 distinct disclosed vaccine sequences | 1,705 residue scores per model, two captures |

Class I uses HLA-A*01:01 and HLA-B*08:01; class II uses HLA-DRB1*03:01 and
HLA-DRB1*08:01. The class II windows are deterministic test inputs, **not
experimentally established epitopes**. NetMHCpan's scan also produces two
9-mer windows from each disclosed 10-mer. Tests check those offsets explicitly.

`outputs/` contains native stdout/stderr and output tables from real local
executions. No prediction values were invented or copied from the source's
experimental labels. `outputs/manifest.json` records commands, input/output
checksums, runtime details, and hashes of the installed code/model trees.
Tree hashes cover sorted relative paths and file checksums, excluding scratch
and hidden files; see `tree_identity` in the generator for the exact algorithm.
Licensed binaries and model weights are not redistributed. Upstream citation
and licensing notices remain in the native outputs.

The tests replay native files through mhctools parsers and wrapper execution
boundaries. They check numeric anchors, both score endpoints, every Gfeller
table cell used by the adapters, allele identity, positions, and duplicate
input ordering. Subprocess and socket calls are blocked during these tests.
Existing live integration tests continue to exercise inference separately.

```sh
python -m pytest tests/test_osteosarc_fixtures.py --require-all
```

## Regeneration

Only regeneration needs osteosarc (Python 3.10+) and installed predictors.
Use an isolated environment with `osteosarc==0.1.2`. Open the existing immutable
cache snapshot; the exporter verifies its ID and makes no network requests:

```sh
env/osteosarc-fixtures/bin/python scripts/osteosarc_fixtures.py inputs \
  --cache env/osteosarc-cache \
  --snapshot mhctools-vaccine-fixtures-2026-09-21 \
  --snapshot-id 9c8938b30873f6dbaf3a712e806fd1b813ee94dac9364ddfa49bb494748e745d \
  --output /tmp/osteosarc-inputs-new

source env/test-backends/activate.sh
python scripts/osteosarc_fixtures.py record \
  --inputs /tmp/osteosarc-inputs-new --output /tmp/osteosarc-outputs-new
```

Both commands refuse an existing output directory. The recorder requires the
licensed NetMHC bundle, configured Gfeller executables, Docker, the legacy
NetMHC compatibility image, and the pinned NetChop image. The repository's
`scripts/setup_test_backends.py` provisions the optional test runtimes.
Native headers include run times and local paths, so a fresh inference capture
is not byte-identical. Review output and provenance changes before replacing
fixtures or updating numeric anchors.

If the old cache is unavailable, create a **new** named snapshot explicitly
with `Dataset.sync(name, cache=Cache(path))`, record its new ID, and review
the changed source data. Mutable upstream URLs cannot recreate old bytes by
themselves. Offline CI always uses the checked-in, checksummed export.
