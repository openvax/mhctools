# Optional backend conformance

Optional peptide pharmacokinetics and uptake adapters use a deliberately small
contract. It records what a wrapper was developed against separately from the
files it actually resolved, and it does not fetch models during prediction.

`BackendSpec` declares the endpoint, upstream compatibility target, license,
serialization format, prediction-only entry point, and reviewed platform and
interpreter support. `ArtifactIdentity` records the role, size, and SHA-256 of
each exact code, weight, preprocessing, or auxiliary file. `BackendInventory`
combines those records with output-affecting settings into a path-independent
identity suitable for prediction metadata and cache keys.

Artifact state and runtime capability are intentionally separate:

| Artifact status | Meaning |
|---|---|
| `missing` | At least one required file was not found. |
| `mismatch` | A file differs from the reviewed checksum. |
| `unverified` | Files were hashed, but no reviewed checksum was declared. |
| `verified` | Every file matches its reviewed checksum. |

| Capability | Meaning |
|---|---|
| `blocked` | Required files are missing or mismatched. |
| `artifacts_located` | Files exist, but their provenance is unverified. |
| `artifacts_verified` | Files match the reviewed inventory; inference has not run. |
| `inference_reproduced` | The verified backend completed an offline prediction. |

These labels prevent file discovery from being presented as reproduced
inference. A dataset, training script, interactive design program, web-only
service, or checkpoint for another endpoint cannot satisfy this contract.

`run_python_sidecar()` supplies the shared runtime boundary. It runs a reviewed,
noninteractive Python prediction entry point under the selected isolated
interpreter, disables user-site imports, activates standard offline modes,
blocks Python socket connections, captures failures, and enforces a timeout.
This is defense in depth for reviewed sidecars, not a sandbox for hostile code:
pickle, joblib, and unrestricted torch checkpoints can execute arbitrary code.
They must match reviewed checksums by default, and opting into an unverified
artifact is an explicit trust decision.

## PeptiVerse inventory

`PeptiVerse` verifies the pinned upstream `inference.py`, the exact
`transformer_wt_log` checkpoint/configuration/calibration files, and the ESM2
weights, configuration, and tokenizer files. Its manifest names
`Transformer_WT_Log` directly; it never accepts upstream's fallback to
`transformer_wt`. The ESM2 snapshot must be local before construction, and the
sidecar replaces the two unused SMILES embedders so prediction cannot trigger
PeptideCLM or ChemBERTa downloads.

`predictor_version` contains both the revisions mhctools was developed against
and a path-independent SHA-256 identity of the files actually supplied. After
a successful call, `artifact_inventory.capability` changes from
`artifacts_verified` to `inference_reproduced`; merely finding the files does
not make that claim.

## PlifePred2 inventory

`PlifePred2` verifies the natural-peptide forest from the official 1.0 wheel,
Pfeature's QSO implementation at revision
`93636eb95bed9df2893b7a0c56b1215e648ecdbf`, both QSO distance matrices, and
the three additional data files that Pfeature reads at process startup. Only
those reviewed resources are copied into the isolated feature workspace.

The official wheel pins scikit-learn 1.4.2. The sidecar checks that version
before loading the joblib forest, so an interpreter that merely happens to
unpickle it with warnings is not reported as exact reproduction. Both the
outer predictor and nested Pfeature process run with socket access disabled and
timeouts. The output-setting choice (`assume_log10_seconds`) is part of the
identity because it changes whether a duration is emitted.

## Capability matrix

This release makes no blanket platform claim for real-model inference. The
conformance suite exercises artifact checks, input/output contracts, offline
execution, and isolated stub sidecars on Linux (Python 3.10-3.12 in CI) and
macOS (Python 3.12 locally); actual third-party inference remains an opt-in
smoke test. A verified predictor instance advances to `inference_reproduced`
only after that local smoke succeeds.

| Backend/candidate | Endpoint-specific release | Static artifact gate | Real inference |
|---|---|---|---|
| PeptiVerse | Human-serum half-life, exact `transformer_wt_log` | Exact source/model/calibration/ESM2 inventory; automated conformance | Opt-in smoke; no published platform combination yet |
| PlifePred2 | Undocumented blood-half-life native regression | Exact forest/QSO/resources inventory; automated conformance | Reproduced on macOS arm64 / Python 3.12 / CPU; otherwise opt-in |
| POSEIDON | None established | Blocked: research/training repository is not an inference-complete endpoint release | Not run |
| PERSEU | None established | Blocked: interactive design path and serialized models do not provide a reviewed prediction-only entry point | Not run |

Finding a repository, dataset, or checkpoint is therefore never shown as an
inference capability. A platform claim can be added only alongside a recorded
real-model smoke result for that exact inventory.

The PlifePred2 entry records a local smoke on 2026-09-14 with scikit-learn
1.4.2 and asset identity
`d1fb877035f93b1302a975c47c08a084d8757964692c64f89f108792da0915c0`.
The verified backend returned native scores `3.157266404183623` and
`3.7934723194625297` for the two opt-in reference peptides and advanced to
`inference_reproduced`. This is evidence only for the listed runtime; it does
not establish the endpoint's undocumented biological semantics.
