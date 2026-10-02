# Maintaining the IEDB SMM subset

This page explains how to rebuild the pinned subset used by the SMM backends.
For ordinary installation, use the [SMM setup guide](../backends.md#local-smm-and-smm-pmbec).

## Rebuild the subset

Run these commands from the repository root:

```sh
mkdir -p env/test-backends
curl --fail --location \
    --output env/test-backends/IEDB_MHC_I-3.1.7.tar.gz \
    https://downloads.iedb.org/tools/mhci/3.1.7/IEDB_MHC_I-3.1.7.tar.gz
python scripts/build_iedb_smm_subset.py \
    --archive env/test-backends/IEDB_MHC_I-3.1.7.tar.gz \
    --output dist/IEDB_MHC_I-3.1.7-smm-subset.tar.gz
```

Both output paths are git-ignored. The
[builder](https://github.com/openvax/mhctools/blob/master/scripts/build_iedb_smm_subset.py)
verifies the official archive's SHA-256, copies the allowed files without
rewriting them, and normalizes archive order and metadata.

The command succeeds only when the resulting digest matches
`EXPECTED_SUBSET_SHA256`. A failed build is not an artifact to publish.

## Check a failed build

| Failure | Meaning and next step |
|---|---|
| Source checksum mismatch | The input is not the pinned official archive. Check the download before rebuilding. |
| Output checksum mismatch | The compressed bytes differ from the hosted subset. Check the compression runtime and extracted contents. |
| Empty allowed prefix | An expected part of the archive is missing. Inspect the upstream layout. |
| Unsupported member type | A selected entry is not a regular file or directory. Review it before changing the builder. |

The same compression runtime reproduces the subset checksum. Different zlib
versions can produce different compressed bytes even when the extracted files
match. The independent content audit is tracked in
[#465](https://github.com/openvax/mhctools/issues/465); keep the checksum gate in
place while that work is pending.

## Files in the subset

The allowlist is defined by `INCLUDED_PREFIXES` in the builder.

| Files | Purpose |
|---|---|
| `LIAI_license.txt`, `Copenhagen_license.txt`, upstream `README` | Preserve the upstream terms and release notes. |
| `src/`, `method/allele-info/`, `method/iedbtools-utilities/` | Provide the code and allele metadata needed by the SMM backends. |
| SMM and SMM-PMBEC training data | Supply the two supported methods. |
| Consensus training data | Supply `distribution_consensus_bin.cpickle` for percentile ranks, even though the consensus method is not used. |

Bundled DTU executables are excluded. The NetMHC family uses its own licensed
installation; see [predictor licenses](../licensing.md#predictor-licenses).

## Hosted artifact

The subset is served from the
[IEDB SMM subset 3.1.7 release](https://github.com/openvax/mhctools/releases/tag/iedb-smm-subset-3.1.7).

| Unpacked artifact | Size |
|---|---|
| Full official release | 1,031 MB across 38,236 members |
| Paths needed for SMM and SMM-PMBEC | 9.4 MB |

The full release is dominated by bundled DTU executables and per-method
training data. Hosting the small subset also avoids making every CI run depend
on the upstream download server. On 2026-09-28, that server timed out on four
attempts from GitHub runners after the Actions cache was evicted, blocking
unrelated PRs.
