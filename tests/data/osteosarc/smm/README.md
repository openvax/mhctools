# Local SMM vaccine fixtures

These are native outputs from SMM 1.0 and SMM-PMBEC 1.0 in the official
IEDB MHC-I 3.1.7 standalone release, recorded on 2026-09-21. No scores are
fabricated or taken from a web prediction service. The 19 sequences are the
unchanged [osteosarc class-I panel](../inputs/class_i.txt), grouped by length
(17 nine-residue peptides, two ten-residue peptides). The original
[source records and caveats](../README.md) apply. Alleles HLA-A*01:01 and
HLA-B*08:01 match the other recorded class-I predictors. These are model
outputs, not experimental validation.

Each method produces 38 peptide/allele predictions. TSV `seq_num` identifies
the input FASTA row, **not** output order: upstream sorts rows by affinity.
`ic50` is nM and `rank` is the upstream percentile. SMM may return affinities
above 50,000 nM; the wrapper preserves them. The accompanying stderr files
are empty on this capture. `manifest.json` records the release checksum,
installed code/model tree hash, exact commands, runtime, input-panel checksum,
and all captured file hashes. Upstream code and model files are not vendored.

Regenerate explicitly into a new directory after installing the pinned runtime:

```sh
python scripts/setup_test_backends.py smm --accept-license
source env/test-backends/activate.sh
python scripts/record_smm_fixtures.py \
  --archive env/test-backends/IEDB_MHC_I-3.1.7.tar.gz \
  --installation env/test-backends/iedb-3.1.7/mhc_i \
  --output /tmp/smm-capture-new
```

`tests/test_smm.py` replays these files without processes or network access.
`tests/test_smm_integration.py` runs both real local predictors and compares
all 76 values and percentiles with the recordings.
