# Getting models: `fetch`, `ls`, and `predictors`

Most predictors need something downloaded before they will run: model weights,
reference files, or a snapshot of the upstream tool. One command covers all of it:

```sh
mhctools fetch <name>
```

Predictors that ship their own download manager keep using it. mhctools reports
which manager owns the files and where they landed, rather than copying them
into a second cache.

## See what is installed

`mhctools ls` lists packaged and optional artifacts, where they live, and who
manages them. Status and manager fields go in a compact table, with installed
paths and missing-artifact notes below it.

```sh
mhctools ls                               # every artifact
mhctools ls mixtcrpred                    # one artifact
mhctools ls mixtcrpred --models --downloaded
mhctools ls --json                        # full structured records
mhctools ls --data-dir /shared/models     # inventory rooted somewhere else
```

## Download an artifact

```sh
mhctools fetch mhcflurry                    # the upstream package's default compatible release
mhctools fetch mhcflurry --version 2.2.0    # an explicit release, for a reproducibility run
mhctools fetch eramer                       # pinned snapshot: code plus the wrapper's model files
mhctools fetch phagescout                   # optional CC BY 4.0 peptide tables (data only)
```

Academic licenses must be reviewed and accepted explicitly:

```sh
mhctools fetch nettcr --accept-license
mhctools fetch mixtcrpred --accept-license
mhctools fetch mixtcrpred --model A0201_NLVPMVATV   # MixTCRpred ships two checkpoints; name others
```

NetCleave and TLimmuno2 publish no license at all. mhctools can still fetch a
pinned snapshot, but only when you confirm that your own use is authorized.
Here `--accept-license` records that acknowledgement; it does not grant rights
mhctools does not have.

```sh
mhctools fetch netcleave --accept-license
mhctools fetch tlimmuno2 --accept-license
```

## Check what can run

`mhctools predictors` reports launchability, which is not the same as
reproduced inference:

```sh
mhctools predictors

# Run the registered reference-inference probes and fail if either backend does
# not reproduce its reference (unknown and not-checked also fail here).
mhctools predictors calis netchop --check reproduced --strict --json
```

## How `fetch` behaves

`fetch` works the same way for every artifact, whichever tier it belongs to.

Already-available counts as success, whether mhctools, a native
downloader, the package itself, or you installed it. Re-running is a no-op, so
a provisioning script can call `fetch` over a whole list without
special-casing manual tools. The `manager` and `fetchable` fields say who owns
each one. Naming a different destination or revision with `--data-dir` or
`--version` is a request to install that exact thing, so it is not satisfied by
an install found elsewhere; `--version` still errors if it disagrees with what
a foreign manager already has.

An artifact that is missing and cannot be fetched fails with one shape of
message: what to install, then the environment variable or `PATH` entry the
wrapper actually reads.

With `--json`, only JSON goes to stdout. Downloader and git progress go to
stderr, so `mhctools fetch <name> --json | jq` is always safe.

## The same thing from Python

```python
from mhctools import ERAMER, MHCflurry, fetch, list_artifacts

MHCflurry.fetch()
fetch("mhcflurry-affinity")
ERAMER.fetch()
for artifact in list_artifacts():
    print(artifact.name, artifact.manager, artifact.version, artifact.path)
```

`fetch()` obtains every safely and legally downloadable artifact the named
wrapper needs. It stops there: it does not install Python packages, execute
upstream setup scripts, or duplicate a cache owned by another package.

## Where snapshots live

mhctools-managed snapshots default to the platform's user data directory:
`~/Library/Application Support/mhctools` on macOS, `~/.local/share/mhctools` on
Linux. To put them on shared or scratch storage, set `MHCTOOLS_DATA_DIR`, pass
`--data-dir`, or use the Python `data_dir=` argument.

Every snapshot lives under `artifacts/<tool>/<git-commit>/` and includes a
`.mhctools-artifact.json` recording its source repository, exact commit, sparse
paths, and license provenance.

[PhageScout's optional full peptide tables](cleavage/phagescout.md#optional-full-deseq2-peptide-profiles)
are data-only assets under `artifacts/phagescout/21387981/`. Their manifest
records fixed source URLs, byte counts and SHA-256 hashes instead of a git
commit; both tables, attribution and license are published atomically.

## Who owns what: the MANAGER column

| Manager | Meaning |
|---|---|
| `mhctools package` / `<package> package` | Weights shipped inside an installed Python package |
| `mhcflurry` | [MHCflurry](predictors/binding.md#mhcflurry)'s own native download cache |
| `mhctools` | A pinned snapshot fetched into the data directory above |
| `user` / `manual` | An existing checkout or licensed executable you own |

Manual artifacts are listed but `fetch` will not redistribute them.

Small published models such as [Calis](predictors/immunogenicity.md#calis) are fully embedded
in the mhctools package, appear as `mhctools package`, and never need a fetch.

## Inventory is not capability

`mhctools ls` is an artifact inventory: `ready` means the required path was
located, not that an executable works.

For a capability report, use `mhctools predictors`. Its `LOCATED`, `RUNNABLE`,
and `REPRODUCED` columns are independent observations, and `not checked` is
never promoted to success. `--check` sets the highest level it will attempt,
while `--strict` returns a nonzero exit unless every selected integration
reaches that level. Reproduction is available only for registered,
reference-backed probes. A successful help command establishes `runnable`, not
`reproduced`.

The earlier `mhctools integrations` spelling still works as an alias.

## Licensing

Some upstream tools are academic or non-commercial, and a few publish no
license at all. `--accept-license` records that you reviewed the terms and
confirmed your own use is authorized. It does not grant rights mhctools does
not have, and it cannot stand in for a license you must request yourself.

The DTU [NetMHC](predictors/binding.md#netmhc)-family downloads are the clearest example: they are
identity-bound licenses. DTU requires a name, position, academic email,
affiliation, and acceptance, then sends a private download link. So
`--accept-license` cannot substitute for the official DTU request form, and
those installations stay `manual` in the inventory.
