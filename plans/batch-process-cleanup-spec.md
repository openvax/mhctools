# Batch process cleanup (#618)

On command failure, startup failure, or interruption, the batch runner must stop
all commands it started, reap direct children, and preserve the original error.
Use owned process groups on POSIX so shell launchers cannot leave descendants
running. Stop submitting new work after a failure is observed. Successful output
and existing concurrency limits remain unchanged.

Add real-process regression coverage for failure while filling a bounded batch,
failure while draining it, failed startup, interruption, and a child that ignores
SIGTERM. Check that cancellation reaches an actual Docker container through the
repository's legacy test launcher; do not assume killing its host client is
enough. Keep any launcher correction local to its ownership boundary.

The Docker regression confirmed a surviving container after the host group was
killed. Give the launcher a unique container name, create before starting, and
remove that exact container on completion or handled interruption. Keep Docker
knowledge out of the generic batch runner and exercise this path in CI.

An immediate logging failure also exposed macOS EPERM when the group leader
exits between polling and signaling. A 100-child probe and a deterministic
already-exited-leader test reproduced it. Reap the exited leader before retrying
the group operation, while preserving real permission errors.

Use the existing process helper rather than a new execution framework. Bump the
patch version, run lint and the full required-backend suite, obtain green CI,
merge, deploy from clean master, and verify PyPI artifacts and the release tag.
