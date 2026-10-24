# Deployment coordination (#610)

## Design

Use one temporary `refs/heads/release-lock` ref on `origin` for the repository.
Create it atomically with an explicit empty `--force-with-lease` expectation.
Each invocation creates a unique commit containing the version and a UUID;
the commit author and timestamp identify the owner for inspection. This
coordinates separate clones, worktrees and machines without a new service.

Hold the lock through lint, tests, build, upload and tagging. On normal exit,
failure or handled interruption, delete it only with a lease matching this
invocation's owner commit. Never automatically steal or expire a lock. Document
inspection and conditional removal after the operator confirms the owner has
stopped; SIGKILL and lost network acknowledgements can require that recovery.

Keep the early tag/PyPI checks, then repeat them after building, immediately
before upload. Verify ownership at that final boundary. An existing release or
an unverifiable lookup stops publication. Lock operations do not push local tags.

## Verification and release

Exercise the real deploy script against temporary bare Git remotes, replacing
only lint/test/build/upload and PyPI responses with deterministic test doubles.
Cover competing clones, failure/interruption cleanup, changed ownership, and
publication or lookup failure during the long gates. Preserve branch, dirty-tree,
existing-release and version-bump behavior with executable tests.

Release 4.4.2 after lint, full backend tests and green PR CI. Deploy from clean
master, then verify published hashes, the release tag and removal of the lock.
