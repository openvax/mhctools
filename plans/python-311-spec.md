# Python 3.11 support-policy migration

Raise the installed mhctools host runtime minimum to Python 3.11 while retaining
the existing isolated runtimes needed by external predictors. Keep the public
API and predictor behavior unchanged.

1. Update package metadata, version, lint target, support docs, and the primary
   CI matrix to Python 3.11–3.14. Keep pandas 2 and 3 coverage on Python 3.11.
2. Resolve the package and optional CI dependencies without forcing versions.
   Run lint and the offline suite on local Python 3.11 and 3.14, then build and
   inspect both distributions.
3. Open a PR closing #601, require green CI, merge, and deploy from clean master.

The older Python interpreters documented for third-party predictor sidecars are
separate from mhctools' installed host runtime and remain pinned as required by
those upstream tools.
