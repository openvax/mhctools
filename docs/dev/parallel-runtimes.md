# Parallel runtime investigation

The October 2026 investigation of [#599](https://github.com/openvax/mhctools/issues/599)
and [#547](https://github.com/openvax/mhctools/issues/547) found a reproducible
timeout-cleanup bug. It did **not** reproduce the historical QEMU segmentation
fault, native alignment timeout, or alignment `IndexError`. Those causes remain
unresolved. Passing under a particular worker count does not establish a fix.

## Verified cleanup defect

`AsyncProcess` previously waited only for the group leader after sending
`SIGTERM`. If that leader exited but a child ignored the signal, cleanup returned
without killing the child. A real-process regression failed before the fix and
passes afterward. Cleanup now checks the owned process group through the grace
period and sends `SIGKILL` if members remain, even after the leader exits.

The native MixMHCpred reference runner now uses this group cleanup instead of
`subprocess.run`, which terminates only its direct child on timeout. Its limit
remains **180 seconds**. Failures report the exit status or timeout limit, elapsed
wall time including cleanup, and the final 2,000 bytes of combined stdout/stderr.
Offline tests verify both diagnostics and child termination. Children that
create a different session and Docker containers owned by the daemon are outside
this process-group guarantee.

A separate deterministic failure leaves already-started sibling commands alive
when a command batch raises. [#618](https://github.com/openvax/mhctools/issues/618)
tracks that issue, including Docker cancellation; this change does not fix it.

## Environment and measurements

The baseline was mhctools 4.4.3, commit
`dde001ffb8c49936377422913933448c728e0d1d`, on macOS 26.5 arm64 with 10 logical CPUs.
The repository-local test environment used Python 3.12.6, pandas 3.1.0 and
pytest 9.1.1. The shared development environment was not modified.

- Docker: aarch64 VM, 10 CPUs, 8,218,034,176 bytes of memory, Docker 28.1.1.
- Legacy NetMHC image: amd64 under emulation, image ID
  `sha256:03da179ae2bdabf795f5f06009f266a2c2217f8cc4371d237b53c9250abdd34a`;
  netmhc-bundle commit `7cbfaf34d8f023f8cce2384eb43af4942e7e1ae4`.
- MixMHCpred commit `0a7f9b9e20d1cf02236f4a0a90d16735be879b38`, isolated
  Python 3.11.12, numpy 2.4.6 with Apple Accelerate, pandas 2.3.3, MAFFT 7.526.

Experiments ran sequentially with all required backends provisioned and
`--require-all`. The focused selection contained the four NetMHCcons tests and
three native-reference integration tests. Durations below are pytest wall times;
alignment durations cover the entire regression, including native and adapter
calls, rather than just the native command subject to the 180-second limit.

| Selection | Workers | Result | Suite wall time | Alignment regression |
|---|---:|---|---:|---:|
| Focused baseline | 1 | 7 passed | 101.88 s | 28.33 s |
| Focused baseline | 4 | 7 passed | 55.55 s | — |
| Full baseline, new cleanup regression added | 7 | 2,168 passed; cleanup regression failed | 571.69 s | 133.16 s |
| Same full suite, numerical-library thread limits of 1 | 7 | Same sole cleanup failure | 566.92 s | 87.72 s |
| Two-allele NetMHCcons regression, 32 parametrized trials | 8 | 32 passed | 98.18 s | — |

The failing baseline regression specifically had a parent that exited on
`SIGTERM` and a child that ignored it. Its own teardown stopped that child;
neither full run reproduced either original issue. The seven-worker choice came
from `test.sh`'s memory calculation; the comparison explicitly held it at seven.
The NetMHCcons stress run used eight concurrent tests with independent temporary
directories. It observed up to 41 host processes, with no observed survivors at
completion; all 32 trials passed without prediction retries.

The thread-limited comparison set `OMP_NUM_THREADS`, `OPENBLAS_NUM_THREADS`,
`MKL_NUM_THREADS`, `VECLIB_MAXIMUM_THREADS`, `NUMEXPR_NUM_THREADS`,
`TF_NUM_INTRAOP_THREADS` and `TF_NUM_INTEROP_THREADS` to `1`. Alignment improved,
but TLimmuno2 slowed from 301.85 to 335.49 seconds. Total time changed by less
than 1%. This single comparison does not justify changing global thread defaults.

One-second host process-tree samples reached 7,597.6 MiB summed RSS and 28
processes in the baseline, versus 8,000.5 MiB and 24 processes with thread limits.
One baseline pytest worker and its descendants reached 5,069.3 MiB summed RSS;
the host's one-minute load average reached 30.08. These measurements show that
a worker can own several expensive subprocesses. RSS sums include shared pages,
host load includes unrelated activity, and host samples exclude Docker guest
memory. They do not prove the cause of a historical crash or timeout.

## Reproduction and next evidence

After activating the repository-local test environment and provisioning backends:

```sh
source env/test-backends/activate.sh
python -m pytest -n 1 tests/test_netmhc_cons.py \
    tests/test_mixmhcpred_runtime_integration.py --require-all -vv --durations=15
# Repeat the focused selection with -n 4 for the parallel comparison.
TEST_SH_MIN=8 TEST_SH_MAX=8 ./test.sh --require-all -vv --durations=25
```

The current focused selection also includes the two offline diagnostic tests.
Run only one load experiment at a time. Preserve failure output, runtime
identities, actual collected test counts, process-tree observations, and Docker
resource observations. Repeating the same node ID on the pytest command line
does not necessarily collect repeated tests, even with `--keep-duplicates`.

For another alignment failure, retain `final_alignment.fasta` and inspect the
MAFFT exit status and diagnostics. The pinned upstream sequence script suppresses
MAFFT stderr and does not check its exit status; this is a diagnostic limitation,
not evidence that MAFFT caused the reported `IndexError`.

Worker defaults and the documented two-worker release configuration remain
unchanged. No prediction retries, skips, or increased timeout were introduced.
The next scheduling change needs repeated measurements or a captured failure
that identifies the exhausted resource or faulty runtime.
