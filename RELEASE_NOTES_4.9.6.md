# MHCTools 4.9.6

The optional `nettcr` extra uses TensorFlow 2.21's supported TFLite interpreter
on Python 3.11–3.13, avoiding LiteRT's unconditional Python<3.11 backport
dependency (#651). Both pip and uv dependency checks are included in its
real-model CI job. Current TensorFlow wheels do not support Python 3.14;
that optional runtime remains dependent on upstream wheel/metadata releases.
The core package still supports Python 3.11–3.14.

Only the known TensorFlow interpreter migration notice is filtered during
construction; other warnings and inference errors propagate. Model weights,
encoding and inference remain unchanged and are checked against NetTCR's
published ensemble predictions.
