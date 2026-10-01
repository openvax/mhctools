# Releasing mhctools

Every PR bumps the [semantic version](https://semver.org/) in
`mhctools/__init__.py` — at minimum a patch bump, even for docs.

Once the PR is merged, a core developer runs `./deploy.sh` from a clean
`master`. It re-runs lint and tests, builds, uploads to PyPI, tags, and pushes.
Passing a version (`./deploy.sh 3.44.56`) makes the bump and commits it for you.

`deploy.sh` refuses to run off `main`/`master`, refuses a dirty tree, stops on
an existing tag or an already-published version, and gates on lint and tests.
If it stops, fix the cause rather than working around the gate.

The documentation site (<https://openvax.github.io/mhctools/>) is built from `docs/` and
deployed to GitHub Pages by the `Docs` workflow on every merge to `master`, so
it needs no separate release step. A new predictor needs a row in
`scripts/predictor_matrix.py` (then `python scripts/predictor_matrix.py`) before
its PR can pass.
