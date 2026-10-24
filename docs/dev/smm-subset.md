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

The command succeeds only when the SHA-256 of the uncompressed tar stream
matches `EXPECTED_SUBSET_TAR_SHA256`. This covers file names, contents,
normalized metadata, order and tar padding. Different gzip encodings can pass
the same audit. A failed build is not an artifact to publish.

The output reports two checksums:

| Field | Meaning |
|---|---|
| `sha256` | Exact bytes of the generated `.tar.gz`; may vary across compression runtimes. |
| `tar_sha256` | Normalized tar stream after gzip decompression; must match the content pin. |

To check the tar stream independently, without extracting files:

```sh
gzip -dc dist/IEDB_MHC_I-3.1.7-smm-subset.tar.gz | shasum -a 256
```

The expected value is
`f1e3d6168644a33e64a8be00d77b4e56808db1bd9b72bcc704d4954828ad494e`.

## Check a failed build

| Failure | Meaning and next step |
|---|---|
| Source checksum mismatch | The input is not the pinned official archive. Check the download before rebuilding. |
| Tar checksum mismatch | The normalized tar stream differs from the pinned content. Inspect the selected files and tar metadata before publishing. |
| Empty allowed prefix | An expected part of the archive is missing. Inspect the upstream layout. |
| Unsupported member type | A selected entry is not a regular file or directory. Review it before changing the builder. |

The installer still requires the exact hosted `.tar.gz` checksum, recorded as
`SMM_SUBSET_SHA256` in `scripts/setup_test_backends.py` and
`EXPECTED_SUBSET_SHA256` in the builder. Passing the content audit does not make
a differently compressed rebuild interchangeable with the published download.
Publishing a replacement artifact requires explicitly updating its download pin.

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
