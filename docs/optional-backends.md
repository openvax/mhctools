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
