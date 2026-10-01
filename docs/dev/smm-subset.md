# Maintaining the IEDB SMM subset

Maintainer notes for the pinned IEDB subset that `mhctools` serves from its own
GitHub release. Users do not need this page; see
[installing SMM](../backends.md#local-smm-and-smm-pmbec).

The subset is served from this repo's
[`iedb-smm-subset-3.1.7`](https://github.com/openvax/mhctools/releases/tag/iedb-smm-subset-3.1.7)
release rather than fetched from `downloads.iedb.org`. That host became
unreachable from GitHub runners on 2026-09-28 (`curl: (28) Connection timeout`,
four attempts) with the Actions cache evicted, which blocked merges on PRs that
had nothing to do with SMM. The release unpacks to 1031 MB across 38,236
members, dominated by bundled DTU executables under `method/` (netmhc-4.0 is
210 MB, netmhc-3.4 192 MB, netmhcpan-4.1 114 MB) plus 192 MB of per-method
training data under `data/`. mhctools does not run those executables from this
bundle; the netMHC family is wrapped through its own licensed distribution.
The paths SMM needs come to 9.4 MB.

The subset is a verbatim copy of `LIAI_license.txt`, `Copenhagen_license.txt`,
the upstream `README`, `src/`, `method/allele-info/`,
`method/iedbtools-utilities/` and the `smm/`, `smmpmbec/` and `consensus/`
training data. `consensus/` is required despite
the consensus method being unused: percentile ranks read
`distribution_consensus_bin.cpickle` from it. `scripts/build_iedb_smm_subset.py`
derives it from the official archive, verifying that archive's own SHA-256 and
normalizing entry order and metadata. Rebuilding with the same compression
runtime reproduces the checksum; different zlib versions may produce different
compressed bytes even when the extracted files match upstream exactly:

```sh
curl -o env/test-backends/IEDB_MHC_I-3.1.7.tar.gz \
    https://downloads.iedb.org/tools/mhci/3.1.7/IEDB_MHC_I-3.1.7.tar.gz
python scripts/build_iedb_smm_subset.py \
    --archive env/test-backends/IEDB_MHC_I-3.1.7.tar.gz \
    --output dist/IEDB_MHC_I-3.1.7-smm-subset.tar.gz
```

Both paths are already git-ignored. The build exits with an error if the output
digest differs from `EXPECTED_SUBSET_SHA256`; do not publish that output. It
also fails if any allowlisted prefix contains no regular files or the archive
grows a member type it does not copy.
