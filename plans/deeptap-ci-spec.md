# Required DeepTAP inference in CI (#562)

Add a `deeptap` setup group that fetches the existing pinned Apache-2.0 snapshot
into the managed test root and provisions a separate Python 3.11 runtime.
Pin the tested torch, Lightning, NumPy and pandas versions; use the official
CPU torch wheel on Linux and run `pip check` in the isolated environment.
Export `DEEPTAP_HOME` and `DEEPTAP_PYTHON` through the existing activation file.

Add DeepTAP to the existing optional-backend CI matrix, where the complete
test module runs with `--require-all`. Keep offline parser/validation tests in
the public Python matrix. Document one setup-and-test command sequence.

Verify fresh provisioning and both actual classification/regression models,
lint, the offline suite, all CI jobs, and the full release gate. Bump the patch
version and release after the preceding report fix.
