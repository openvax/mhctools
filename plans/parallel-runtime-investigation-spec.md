# Parallel runtime investigation (#599, #547)

Record the host, Docker VM and backend runtime identities before reproducing
the existing tests. Compare isolated execution with controlled parallel loads,
capturing elapsed time, process-tree CPU/RSS and diagnostic output. Run only one
load experiment at a time; do not change model inputs, prediction semantics or
the native reference's 180-second timeout to obtain a pass.

Check process cleanup independently with real launcher/child fixtures, including
a child that ignores SIGTERM after its parent exits. Repair verified lifecycle
or diagnostic defects with small shared changes and behavioral regression tests.
Use measured evidence before changing worker defaults or runtime scheduling.

Record observations and their limits on the existing issues. A passing rerun or
lower worker count alone does not resolve an intermittent crash or timeout.
Any PR bumps the version, passes lint/full tests/CI, then merges and deploys from
clean master with published artifacts, tag and release-lock cleanup verified.
