# Canonical allele identity migration (#631)

Status: plan only, refreshed 2026-11-04 after the CLI usability implementation.
Do not close #631 until the implementation, downstream checks and minor release
ship. The usability stack changes no allele identities.

## Objective and completed prerequisites

Use mhcgnomes' parsed `Allele`/`Pair` objects and `to_string()` as the source of
canonical prediction identities, supported-allele lookup keys and cache inputs.
Stop selecting an arbitrary allele for an ambiguous serotype or supertype.

The formatter prerequisites shipped in #637 and #638: `allele_formats.py`
provides public NetMHC, NetMHC 4, NetMHCIIpan, DeepImmuno and CapHLA formatters,
plus `prefer_allele_spelling`. These operate on mhcgnomes objects independently
of the legacy parsed tuple. Keep that predictor-specific knowledge in mhctools.
A canonical identity and a native executable's accepted input spelling serve
different purposes.

The remaining compatibility layer is `allele_normalization.py`. It reconstructs
names, truncates fields, removes meaningful leading zeros and silently chooses
representatives for groups. CapHLA and DeepImmuno still normalize then parse a
second time. The MHCflurry adapter still has special handling for original allele
identities; audit the exact #544 path before removing it.

## Proposed identity contract

1. Parse with mhcgnomes, retaining the existing documented input aliases (`;`
   between chains and `_` in place of `*`) through one shared parsing entry point.
   Only a specific `Allele`, or a `Pair` of specific alleles, is a prediction
   identity. Use parsed types and properties; never infer species from a string
   prefix or select a group's first/representative member.
2. Return `parsed.to_string()` without local re-spelling. Preserve all supplied
   allele fields, species-specific leading zeros and expression annotations.
   A predictor may explicitly select a supported resolution at its adapter
   boundary; canonicalization must not silently discard that information.
3. Preserve the current default class-II pairing policy
   (`infer_class2_pair=True`) during this migration. A DR beta chain therefore
   retains its inferred alpha chain, now spelled by mhcgnomes. Keep the explicit
   `infer_class2_pair=False` behavior. Audit `omit_dra1` callers and retain its
   documented projection using parsed chain identity during a deprecation
   period; native DR spelling stays in the NetMHCIIpan/CapHLA formatters.
4. Reject `Serotype`, `Supertype`, gene-only, unresolved and grouped results
   before predictor execution. The error must explain the ambiguity and request
   a specific allele, without suggesting a representative as if equivalent.
5. Distinguish a parser failure from a parsed-but-ambiguous name. The permissive
   helper must never catch an ambiguity error and return the original group name
   as a successful raw fallback. Keep a dedicated ambiguity error or confine
   the fallback's `except ParseError` to mhcgnomes' parsing call itself.
6. Retain the existing raw-name path for names mhcgnomes genuinely cannot parse
   but a native tool lists. Preserve case and the existing whitespace/asterisk
   matching rule; retain validation against the tool's actual supported list.
   Do not grow guessed aliases or a second nomenclature implementation.

## Concrete before/after fixtures

The following were probed against baseline commit `921b791` and mhcgnomes 3.67.0.
The implementation must commit fixtures with the tested parser version and
verify them against primary nomenclature sources and each native tool's lists.

| Input | Current canonical output | Planned canonical output / result |
|---|---|---|
| `HLA-A*02:01` | `HLA-A*02:01` | unchanged |
| `H2-Kb` | `H-2-Kb` | `H2-K*b` |
| `H2-IAb` | `H-2-IAb` | `H2-AA*b/AB*b` |
| `HLA-DQA1*01:02-DQB1*06:02` | hyphen-separated pair | `HLA-DQA1*01:02/DQB1*06:02` |
| `DLA-88*001:01` | `DLA-88*01:01` | `DLA-88*001:01` |
| `HLA-DRB1*01:01` (default pairing) | `HLA-DRA1*01:01-DRB1*01:01` | `HLA-DRA*01:01/DRB1*01:01` |
| `HLA-A*02:01:01:01` | `HLA-A*02:01` | preserve all four fields |
| `HLA-A*30:14L` | strict error; permissive raw path | canonical parsed allele with `L` retained |
| `A2` / `A02` | `HLA-A*02:01` | ambiguity error; no prediction launched |

HLA suffixes and higher-resolution fields carry defined information; preserving
these in identity does not claim every predictor supports them. Unsupported
model inputs must still raise the existing capability/allele error.

## Implementation order

### 1. Contract and dependency audit

- Inspect all uses of `normalize_allele_name`, `normalize_allele_name_or_raw`,
  `parse_classi_or_classii_allele_name`, `omit_dra1`, supported-allele caches and
  serialized allele keys in mhctools and topiary. Check MHCflurry's two-field
  model normalization separately from the lossless canonical identity contract.
- Capture full current native allele lists and classify parser failures,
  ambiguous results, annotations and high-resolution names. Report counts and
  examples; file parser defects in pirl-unc/mhcgnomes rather than compensating
  for them locally. Existing parser issues remain separate work.
- Test the minimum supported mhcgnomes against this corpus and raise the declared
  minimum to the earliest stable version that passes. The current audited
  reference is 3.67.0; do not rely on a newer local checkout while declaring an
  untested old minimum.
- Inventory old cache formats and the original-input information retained in
  saved tables before selecting the migration procedure.

### 2. One coordinated mhctools minor-release PR

- Replace the compatibility tuple and spelling helpers with the contract above.
  Use bounded caching only after parsing/error behavior is settled. Remove the
  tuple-based API once internal consumers are gone; document its removal and
  public replacements rather than silently changing its return shape.
- Route input names, supported-allele maps and parsed predictor output through
  the same canonical identity function. Deduplicate aliases by parsed identity;
  retain `prefer_allele_spelling` and original native accepted spellings for
  command construction. Preserve request/output identity when tools emit a
  lower-resolution or alternative spelling.
- Remove the redundant normalize/parse cycle in CapHLA and DeepImmuno and the
  DLA workaround only after a real MHCflurry DLA regression demonstrates the
  generic path preserves the original identity and scores.
- Update stored-schema/cache namespaces where canonicalization affects keys.
  Include canonicalization-policy version and mhcgnomes version in relevant
  cache provenance. Explicitly invalidate/rebuild legacy supported-allele and
  prediction caches instead of mixing old/new keys or editing cache files in
  place. Confirm unaffected cache formats need no change.
- Add migration notes with the complete before/after table, ambiguity errors,
  high-resolution handling, parameter deprecations and downstream join guidance.
  Original inputs or trustworthy upstream identities are required to repair
  lossy old rows; do not invent missing DLA zeros, fields or serotype members.

Ship as the next available minor version after the usability stack (tentatively
4.9.0). Do not expose a half-migrated default identity across separate releases.
If a behavior-preserving prerequisite is discovered, ship it first as a patch
with its own full verification and deployment.

### 3. Downstream compatibility and release gate

Before merging the minor release, exercise topiary's supported-allele selection,
ranking joins, cached predictions and serialization against the candidate wheel.
Verify class-I joins with MHCflurry using explicit model-resolution mapping,
without pretending a four-field input and a two-field model key retain identical
resolution. File and link any necessary downstream PR, then land in dependency
order. Record required cache invalidation/recomputation in release instructions.

## Acceptance tests

- Idempotence and round trips for human class I, DR/DQ/DP, mouse class I and II,
  DLA and other supported non-human alleles; preserve zeros, full resolution,
  suffixes and chain order. Independently asserted expected identities, not
  expectations computed using the same helper under test.
- `A2`, `A02`, other groups and invalid pair components fail through both strict
  and permissive APIs, CLI and adapter entry points before any subprocess runs.
- Raw fallbacks still accept tool-listed unparseable names and remain
  case-sensitive. Parser failures, unsupported annotations and ambiguity have
  distinct tested behavior.
- Full native-list parse/format/lookup comparisons; real NetMHCpan and
  NetMHCIIpan predictions, class-II pairs, DLA/MHCflurry, CapHLA and DeepImmuno.
  Compare exact peptides/identities and numerical outputs with existing native
  tolerances. Do not replace scientific expected scores with new snapshots just
  because canonical strings changed.
- Alias collisions, duplicate requests, high-resolution model projection,
  mixed old/new cache inputs, cache-version boundaries, serialized round trips
  and downstream joins have explicit regressions.
- Lint, complete real-backend suite (`--require-all`, two workers), docs checks,
  tested dependency minimum/current versions, supported Python CI and relevant
  downstream tests pass. Merge only with green CI; deploy clean master and
  verify the PyPI version and artifact hashes before closing #631.

## Sources and decision boundaries

- [Issue #631 and its scope clarification](https://github.com/openvax/mhctools/issues/631): canonical identity changes; predictor-specific formats stay.
- [mhcgnomes implementation](https://github.com/pirl-unc/mhcgnomes): parsed result types, pair inference and `to_string()`; the fixture probe used 3.67.0.
- [Official HLA allele naming](https://hla.alleles.org/pages/nomenclature/naming_alleles/): field meanings and expression suffixes.
- [Official HLA suffix assignment](https://hla.alleles.org/pages/nomenclature/assigning_suffixes/): expression annotations must not be discarded as punctuation.
- [IPD-MHC DLA group](https://www.ebi.ac.uk/ipd/mhc/group/DLA/): species-specific nomenclature and authoritative allele records.

This plan preserves the existing pairing default and proposes lossless canonical
identity. The audit must surface any downstream requirement that conflicts with
those choices before implementation; it must not silently collapse identities
to make a native tool or an old fixture pass.
