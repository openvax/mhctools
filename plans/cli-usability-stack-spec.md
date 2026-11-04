# CLI usability stack (#633, #636)

## Scope and release order

1. **#633 — version reporting, 4.7.1:** add top-level `--version` / `-V`,
   printing `mhctools <version>` without requiring a predictor. Preserve
   `fetch --version` as an upstream artifact revision. Document and test both.
2. **#636 — completion and startup, 4.8.0:** make package exports lazy without
   changing the public API; keep parser construction free of scientific/runtime
   imports; enable argcomplete for prediction flags, all subcommands and aliases,
   and each chosen command's options. Apply the existing OpenVax cyan help style
   only when printing to a capable terminal; keep machine output plain.

Stack the second branch on the first. Review each final diff, pass lint, the full
real-backend suite with `--require-all` and `TEST_SH_MAX=2`, docs checks and GitHub
CI, then merge/deploy each release from clean master in order. Use repository
release environments, preserve the shared environment, verify `pip check` and
published PyPI artifact hashes. Measure fresh-process completion startup against
#636's 100 ms budget and enforce import boundaries in subprocess tests.

## Follow-up plan

After implementing usability, update the existing #631 canonical allele naming
spec with completed prerequisites, exact compatibility decisions, adapter/cache
impacts, release staging and validation gates. This task plans that change; it
does not change allele semantics.

## Release coordination

The concurrent #628, #630 and #645 reserve 4.6.2, 4.6.3 and 4.7.0.
This stack therefore uses 4.7.1 then 4.8.0 and follows those releases.
Rebase and revalidate against the resulting main before merging.
