# CleaveNet MMP substrate scores

Use `CleaveNet` to score short peptide substrates against 18 matrix
metalloproteinases (MMPs). It returns whole-substrate evidence in dedicated
`CleaveNetResult` records. It does not return bond tracks.

## Install and run

Fetch the pinned prediction assets, then install their requirements into an
isolated Python 3.11 or 3.12 environment. TensorFlow stays outside the mhctools
interpreter. CPU inference is exercised on Linux and macOS arm64.

```bash
mhctools fetch cleavenet
python scripts/setup_test_backends.py cleavenet --python python3.11
source env/test-backends/activate.sh
```

The provisioning script is in the source checkout. For a package-only install,
run `mhctools fetch cleavenet --json`, create a separate virtual environment,
install `-r <returned-path>/requirements.txt`, and set `CLEAVENET_PYTHON` to its
`bin/python`. Set `CLEAVENET_HOME` to that returned path if you use a custom
artifact directory. No GPU is required. TensorFlow is pinned to 2.18.0.

```python
from mhctools import CleaveNet, PeptideInput

model = CleaveNet()
results = model.predict([
    "PRVFQLRVFL",
    PeptideInput("LRVFL", source_sequence_name="candidate", source_start=12),
])
for score in results[0].scores:
    print(score.enzyme, score.z_score, score.ensemble_sd)

record = results[1].to_dict()  # exact input, padding, assay and asset inventory
frame = model.predict_dataframe(["PRVFQLRVFL", "LRVFL"])
```

Inputs must be canonical uppercase sequences of 1–10 residues. Modified
termini, attachments, ambiguous residues and longer substrates raise before
inference. Short sequences receive centered padding to ten positions;
`LRVFL` becomes `--LRVFL---` and `AHA` becomes `---AHA----`. The extra gap goes
on the right. This follows the [upstream README](https://github.com/microsoft/cleavenet/tree/4dac67defc99ca35d967ddc76eca0fe8b74afdad#cleavenet-predictor),
rather than its CLI's batch-dependent right padding. Short inputs carry an
explicit limited-validation annotation.

For longer sequences, score overlapping ten-residue windows:

```python
windows = model.predict_windows(PeptideInput(
    "ACDEFGHIKLMNPQRSTVWY", source_sequence_name="candidate", source_start=20))
assert windows[1].peptide_input.source_start == 21
```

Coordinates are zero-based, and DataFrame `source_end` is exclusive. A window
score does not assert that the window is released from the parent or identify
a central bond. Duplicate inputs remain separate occurrences.

## Read the evidence

The native regression target is a dimensionless cleavage Z-score from the
Kukreja mRNA-display dataset. Five transformer models supply the mean and
population standard deviation. Higher scores indicate stronger relative
substrate cleavage in that assay; the standard deviation describes ensemble
spread. Neither quantity is a cleavage probability or a calibrated confidence
interval. The [primary publication](https://doi.org/10.1038/s41467-025-67226-1)
describes the training endpoint and its dependence on assay/library context.

The output includes all 18 native heads in upstream order, including MMP24,
which the pinned README's enzyme enumeration omits. No cross-enzyme survival
score is calculated. The model has no enzyme concentration, incubation time,
pH or experimental-matrix inputs. Its display-assay target does not establish
cleavage of free peptides in serum, cathepsin/AEP activity, APC uptake,
antigen presentation or peptide half-life.

Each result preserves the exact `PeptideInput`, padded model input, assay
endpoint, source reference, runtime versions, and SHA-256 inventory of code,
weights and vocabulary resources. Files must match the reviewed revision
`4dac67defc99ca35d967ddc76eca0fe8b74afdad`. Prediction runs offline and never
trains or generates sequences.

The real-model regression reproduces the official source on a ten-mer and a
padded five-mer. This establishes adapter conformance. Independent benchmark
curation, sequence-overlap audits and downstream batch-report integration
remain tracked in [#476](https://github.com/openvax/mhctools/issues/476).

## Licenses

The fetch retains the [MIT code license](https://github.com/microsoft/cleavenet/blob/4dac67defc99ca35d967ddc76eca0fe8b74afdad/LICENSE)
and the separate [CDLA-Permissive-2.0 data agreement](https://github.com/microsoft/cleavenet/blob/4dac67defc99ca35d967ddc76eca0fe8b74afdad/data/LICENSE).
Preserve that data agreement when sharing upstream data. mhctools vendors no
upstream code, weights or training data.
