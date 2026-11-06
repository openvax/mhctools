# Canonical allele identity migration (#631)

Status: implementation in progress, 2026-11-05; release gates below remain required.
Do not close #631 until the implementation, downstream checks and minor release
ship. The implementation covers accepted types, pair validation, cross-repository
pairing, G/P groups and historical aliases; release 4.9 changes runtime identities
as described in `docs/alleles.md`.

## Objective and completed prerequisites

Use validated mhcgnomes `Allele`, `AlleleWithoutGene` and `Pair` objects and
`to_string()` as the source of canonical prediction identities, supported-allele
lookup keys and cache inputs.
Stop selecting an arbitrary allele for an ambiguous serotype or supertype.

The formatter prerequisites shipped in #637 and #638: `allele_formats.py`
provides public NetMHC, NetMHC 4, NetMHCIIpan, DeepImmuno and CapHLA formatters,
plus `prefer_allele_spelling`. These operate on mhcgnomes objects independently
of the legacy parsed tuple. Keep that predictor-specific knowledge in mhctools.
A canonical identity and a native executable's accepted input spelling serve
different purposes.

The previous `allele_normalization.py` reconstructed names, truncated fields,
removed meaningful leading zeros and silently chose representatives for groups.
The implementation removes that layer and the redundant CapHLA/DeepImmuno parse
cycle. MHCflurry keeps canonical requests separate from its native model keys,
including the exact supported-key behavior required by #544.

## Proposed identity contract

### Parsing and accepted identities

Use one shared parsing entry point, retaining the documented `;` pair separator
and `_` spelling of `*`. Set the identity-affecting options explicitly:

```python
parse(name, use_allele_aliases=False, infer_class2_pairing=False,
      collapse_singleton_haplotypes=True, collapse_singleton_serotypes=False,
      max_allele_fields=None, raise_on_error=True)
```

Keep the existing species-context policy. Validate the result before returning
`parsed.to_string()`; preserve all fields, leading zeros and expression suffixes.
Accept these three cases:

- A specific `Allele` which passes the group checks below.
- A standalone `AlleleWithoutGene`, including `BoLA-T2b` and `BoLA-T2c` from the
  committed native lists. Preserve its own canonical spelling without inventing
  a gene. Absence of an assigned gene does not make this an ambiguous allele.
  Model support still requires the normal native-list/capability validation.
- A `Pair` whose two components are specific `Allele` objects of the same parsed
  species, with `alpha.is_class2_alpha` and `beta.is_class2_beta` both true, in
  alpha/beta order. Validate both chains for groups and annotations. Reject
  class-I/class-I, beta/alpha, same-role, gene-only and unresolved components;
  never silently reorder them. `Pair` alone is not proof of these conditions.
  Actual support for a particular alpha/beta combination remains the predictor's
  responsibility; do not add guessed locus-pairing rules to the normalizer.

For each `Allele`, reject `is_allele_group` and either `G` or `P` in the parsed
`annotations`. These are different tests: mhcgnomes 3.68.1 deliberately reports
`is_allele_group=False` for G/P typings, and `annotation_group` covers only G.
G/P typings are rejected by this single-allele contract, even when their binding
domains agree. Do not remove the suffix, choose a member or treat the group as a
specific allele. A future group-aware prediction API would need its own identity
and aggregation contract. Retain expression annotations such as `N`, `L` and `Q`;
canonical identity does not override a model's or downstream application's
expression/capability restrictions.

Reject remaining `Serotype`, `Supertype`, `Gene`, unresolved `Haplotype` and other
result types with an error requesting a specific allele or explicit valid pair.
Use parsed properties, never species-prefix tests or representative selection.

### Class-II pairing and model inputs

Canonical identity does **not** infer a missing alpha chain. This deliberately
replaces mhctools' old default and agrees with Topiary's current parsing default
and its proposed stricter normalizer in topiary#578. `HLA-DRB1*01:01` and
`HLA-DQB1*06:02` stay single-chain identities; an explicitly supplied pair stays
a pair. A beta-only input and an explicit pair have different cache/join keys.

The default singleton-haplotype collapse remains enabled: mhcgnomes already
resolves the curated `H2-IAb` shorthand to `H2-AA*b/AB*b` without enabling
`infer_class2_pairing`. Preserve this case and test it separately from adding a
missing human alpha chain. Disabling collapse unexpectedly turns `H2-IAb` into
`Gene H2-AB` in 3.68.1; that separate parser defect is tracked in
[mhcgnomes#253](https://github.com/pirl-unc/mhcgnomes/issues/253).

A predictor may accept a beta chain directly or require an explicit pair. If it
requires a pair, report the missing chain rather than guessing one. A caller may
explicitly resolve a pair before prediction; that pair is then a new, explicit
input identity, with the resolution choice recorded in provenance. Native
formatters may omit DRA or project model resolution without changing the
canonical request/output identity. When different requests share a native
spelling, retain a request-to-model mapping and restore each requested identity
in output; include effective model inputs in prediction-cache provenance/keys.

Change `normalize_allele_name`'s pairing default to false in the coordinated minor
release. Retain explicitly false `infer_class2_pair`/`omit_dra1` keywords for one
documented deprecation cycle. True values must raise a migration error directing
callers to resolve an explicit pair or use a native formatter; they must not
produce alternate canonical keys or be silently ignored. Migrate all internal
callers before removing the old tuple and projection code.

### Historical aliases and native fallback

Keep `use_allele_aliases=False` in canonicalization, in both mhctools and Topiary.
This preserves `HLA-C*04:22` instead of renaming it to `HLA-C*04:15:02` and avoids
an unannounced historical-name migration. `Mamu-A01` remains a group-looking
input and is rejected; an error may suggest the explicit current name resolved
by mhcgnomes, without substituting it. Generic inputs and imported canonical
tables must not enable aliases merely because the first parse was rejected.

At a native-tool boundary, a catalog entry that parses as an allele group may
have a verified historical spelling, such as NetMHC's `Mamu-A01`. An adapter may
resolve that entry with mhcgnomes' curated alias lookup only if it yields one
identity passing the full contract. Map the explicit current identity to the
tool's accepted spelling, and map native output back through that same catalog.
Keep the native name, resolved identity, tool version and parser/alias provenance.
This scoped adapter mapping must not rewrite already-specific catalog identities
such as `HLA-C*04:22` or accept ambiguous user input. Never maintain a second alias
table in mhctools. Gate this behavior with the cases in
[topiary#582](https://github.com/openvax/topiary/issues/582).

Restrict permissive raw fallback to an actual parsing failure for a name validated
against the native tool's list. Preserve case and the existing whitespace/asterisk
matching rule. A parsed `Gene` such as `H-2-Qa1` or `Mamu-B17` is not a parser
failure: it remains an explicit unsupported-identity error, tracked upstream in
[mhcgnomes#190](https://github.com/pirl-unc/mhcgnomes/issues/190) and topiary#582.
Document this deliberate loss of the old permissive path; do not claim those
names are preserved by raw fallback. Keep `BoLA-T2b`/`BoLA-T2c` working through
their accepted `AlleleWithoutGene` type.

Use a dedicated identity-validation error, and restrict `except ParseError` to
the parsing call itself. Strict and permissive user-input paths must propagate
group/pair/identity errors. Catalog ingestion instead classifies and reports
unsupported entries while retaining valid entries: an unrelated `Gene` or group
in a tool's full list must not prevent a valid HLA or BoLA request from running.
Audit all existing `AlleleParseError` handlers for this distinction.

## Concrete before/after fixtures

The baseline is mhctools 4.8.1 (`8a35062`); parser probes use mhcgnomes 3.68.1
with the explicit options above. Planned rejections are policy decisions, not
claims that mhcgnomes itself refuses those strings.

The implementation must commit fixtures with the tested parser version and
verify them against primary nomenclature sources and each native tool's lists.

| Input | Current behavior | Planned canonical output / result |
|---|---|---|
| `HLA-A*02:01` | `HLA-A*02:01` | unchanged |
| `H2-Kb` | `H-2-Kb` | `H2-K*b` |
| `H2-IAb` | `H-2-IAb` | `H2-AA*b/AB*b` |
| `HLA-DQA1*01:02-DQB1*06:02` | hyphen-separated pair | `HLA-DQA1*01:02/DQB1*06:02` |
| `DLA-88*001:01` | `DLA-88*01:01` | `DLA-88*001:01` |
| `BoLA-T2b` / `BoLA-T2c` | verbatim native names | unchanged, accepted `AlleleWithoutGene` |
| `HLA-DRB1*01:01` | `HLA-DRA1*01:01-DRB1*01:01` | `HLA-DRB1*01:01` |
| `HLA-DQB1*06:02` | `HLA-DQA1*01:02-DQB1*06:02` | `HLA-DQB1*06:02` |
| explicit `HLA-DRA*01:01/DRB1*01:01` | inferred/explicit DR pair share one key | `HLA-DRA*01:01/DRB1*01:01`, distinct from beta-only |
| `HLA-A*02:01:01:01` | `HLA-A*02:01` | preserve all four fields |
| `HLA-A*30:14L` | strict error; permissive raw path | canonical parsed allele with `L` retained |
| `A2` / `A02` | `HLA-A*02:01` | ambiguity error; no prediction launched |
| `HLA-A*02:01/HLA-A*03:01` | accepted by the formatter parser as a `Pair` | invalid-pair error |
| `HLA-DRB1*01:01/HLA-DRA*01:01` | accepted by the formatter parser as a `Pair` | invalid-order error; no reordering |
| `HLA-DQA1*01/DQB1*06:02` | accepted by the formatter parser as a `Pair` | group-in-chain error |
| `HLA-A*02:01:01G` / `HLA-A*02:01P` | strict error; permissive raw path | group error in both APIs |
| `Mamu-A01` as user input | `Mamu-A01` | group error; explicit `Mamu-A1*001:01` can use a verified native alias mapping |
| `HLA-C*04:22` | `HLA-C*04:22` | unchanged; historical aliases disabled |
| `H-2-Qa1` / `Mamu-B17` | permissive raw path | gene-only identity error; no raw escape |

HLA suffixes and higher-resolution fields carry defined information; preserving
these in identity does not claim every predictor supports them. Unsupported
model inputs must still raise the existing capability/allele error.

## Implementation order

### 1. Contract and dependency audit

- Inspect all uses of `normalize_allele_name`, `normalize_allele_name_or_raw`,
  `parse_classi_or_classii_allele_name`, `infer_class2_pair`, `omit_dra1`,
  supported-allele caches and serialized allele keys in mhctools and topiary.
  Check MHCflurry's two-field model normalization separately from the lossless
  canonical identity contract.
- Capture full current native allele lists and classify parser failures,
  ambiguous results, annotations and high-resolution names. Report counts and
  examples, including `AlleleWithoutGene`, G/P annotations and invalid pairs;
  file parser defects in pirl-unc/mhcgnomes rather than compensating for them
  locally. Separate verified native aliases, genuine parse failures and parsed
  unsupported identities. Existing parser issues remain separate work.
- Test the minimum supported mhcgnomes against this corpus and raise the declared
  minimum to the earliest stable version that passes. The current audited
  reference is 3.68.1; do not rely on a newer local checkout while declaring an
  untested old minimum.
- Inventory old cache formats and the original-input information retained in
  saved tables before selecting the migration procedure.

### 2. One coordinated mhctools minor-release PR

- Replace the compatibility tuple and spelling helpers with the contract above.
  Use bounded caching only after parsing/error behavior is settled. Remove the
  tuple-based API once internal consumers are gone; document its removal and
  public replacements rather than silently changing its return shape.
- Route input names, supported-allele maps and parsed predictor output through
  the same canonical identity function, after any explicitly scoped native alias
  resolution described above. Deduplicate names by validated parsed identity;
  retain `prefer_allele_spelling` and original native accepted spellings for
  command construction. Preserve request/output identity when tools emit a
  lower-resolution or alternative spelling.
- Remove the redundant normalize/parse cycle in CapHLA and DeepImmuno and the
  DLA workaround only after a real MHCflurry DLA regression demonstrates the
  generic path preserves the original identity and scores.
- Update stored-schema/cache namespaces where canonicalization affects keys.
  Include canonicalization-policy version and mhcgnomes version in relevant
  cache provenance, including pairing and alias options and native alias maps.
  Explicitly invalidate/rebuild legacy supported-allele and prediction caches
  instead of mixing old/new keys or editing cache files in
  place. Confirm unaffected cache formats need no change.
- Add migration notes with the complete before/after table, ambiguity errors,
  high-resolution handling, parameter deprecations and downstream join guidance.
  Original inputs or trustworthy upstream identities are required to repair
  lossy old rows; do not invent missing DLA zeros, fields or serotype members.
  Old inferred pairs cannot reveal whether the user supplied an alpha chain;
  retain the recorded pair identity with unknown original-input provenance or
  recompute from trustworthy original inputs. Do not relabel an old pair as a
  beta-only request by guessing.

Ship as the next available minor version after the usability stack (tentatively
4.9.0). Do not expose a half-migrated default identity across separate releases.
If a behavior-preserving prerequisite is discovered, ship it first as a patch
with its own full verification and deployment.

### 3. Downstream compatibility and release gate

Before merging the minor release, exercise topiary's supported-allele selection,
ranking joins, cached predictions and serialization against the candidate wheel.
Use one shared fixture corpus in both repositories. Coordinate with
[topiary#578](https://github.com/openvax/topiary/pull/578), which already accepts
`AlleleWithoutGene`, validates pair roles and leaves beta-only inputs unpaired.
Its reviewed candidate does not reject G/P groups; require that downstream change
and matching alias-policy tests before shipping this identity migration.

Assert exact canonical-key equality across mhctools and Topiary for DR and DQ
beta-only inputs and explicit pairs, DP pairs, mouse class II and BoLA identities.
Assert that beta-only and explicit-pair keys stay distinct through live/native
output, cached predictions, occurrence joins and serialized reloads. Test a batch
containing both a DR beta-only request and its explicit pair, so restoring native
output cannot collapse one request into the other. A backend that requires a
DQ pair must reject a beta-only request clearly, without creating a hidden alpha.
Keep Topiary's expression/eligibility rules separate: a preserved null suffix
does not require downstream code to accept that allele for presentation.

Verify class-I joins with MHCflurry using explicit model-resolution mapping,
without pretending a four-field input and a two-field model key retain identical
resolution. File and link any necessary downstream PR, then land in dependency
order. Record required cache invalidation/recomputation in release instructions.

## Acceptance tests

- Idempotence and round trips for human class I, DR/DQ/DP, mouse class I and II,
  DLA and other supported non-human alleles; preserve zeros, full resolution,
  suffixes and chain order. Independently asserted expected identities, not
  expectations computed using the same helper under test.
- `BoLA-T2b` and `BoLA-T2c` stay accepted as `AlleleWithoutGene`; compare their
  canonical keys, native list spellings and parsed predictions. Inference-disabled
  `H2-IAb` still yields its explicit pair with singleton-haplotype collapse on.
- `A2`, `A02`, one-field groups, G/P groups (including either chain of a pair),
  class-I/class-I pairs, reversed pairs, same-role pairs and invalid components
  fail through strict/permissive APIs, CLI and adapters before prediction starts.
  Keep `HLA-MICA*001` accepted: a single field alone is not a group test.
- DR/DQ beta-only identities remain unpaired with matching Topiary keys;
  explicit pairs stay distinct. Test the legacy keyword migration errors and
  prevent adapter inference/projection from overwriting request identities.
- `HLA-C*04:22` stays unchanged; user `Mamu-A01` fails while explicit
  `Mamu-A1*001:01` can round-trip through NetMHC's verified `Mamu-A01` catalog
  entry. Alias lookup must not rescue an unverified group, gene, G/P typing or
  arbitrary imported table value. Record the alias provenance in cache fixtures.
- Raw fallbacks still accept tool-listed unparseable names and remain
  case-sensitive. Parser failures, unsupported annotations and ambiguity have
  distinct tested behavior. A mixed native list containing rejected identities
  still supports unrelated valid HLA and BoLA requests; requesting the rejected
  identity itself raises, including `H-2-Qa1` and `Mamu-B17`.
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
- [mhcgnomes implementation](https://github.com/pirl-unc/mhcgnomes): parsed result types, pair inference and `to_string()`; this review's probes used 3.68.1.
- [Committed NetMHCpan list](../tests/data/allele_lists/netmhcpan42_listmhc.txt): native `BoLA-T2b` and `BoLA-T2c` spellings.
- [Topiary normalizer reviewed at 00a9baa](https://github.com/openvax/topiary/blob/00a9baaff797c8e90f56bd7c13596ce231659a05/topiary/alleles.py): accepted result types, ordered class-II roles, unpaired beta identity and historical-alias defaults; recheck the final downstream implementation before release.
- [Official HLA allele naming](https://hla.alleles.org/pages/nomenclature/naming_alleles/): field meanings and expression suffixes.
- [Official HLA suffix assignment](https://hla.alleles.org/pages/nomenclature/assigning_suffixes/): expression annotations must not be discarded as punctuation.
- [Official G groups](https://hla.alleles.org/pages/p%26g_groups/g_groups/) and [P groups](https://hla.alleles.org/pages/p%26g_groups/p_groups/): typing groups describe shared binding-domain sequences, not one fully specified allele.
- [IPD-MHC DLA group](https://www.ebi.ac.uk/ipd/mhc/group/DLA/): species-specific nomenclature and authoritative allele records.

This plan changes mhctools' canonical pairing default deliberately and requires
the same identity policy in both repositories. Changes to that policy, including
globally resolving historical aliases or accepting G/P group identities, require
a separate documented migration decision and cache-policy version. The audit
must surface conflicting downstream requirements before implementation; it must
not silently collapse identities to make a native tool or an old fixture pass.


## Implementation audit (4.9.0 candidate)

The shared fixture corpus is `tests/data/canonical_allele_identities.json`.
The 41 core contract regressions pass with both the minimum mhcgnomes 3.68.1
and current 3.70.1. The full captured native catalogs (NetMHCpan 4.2: 12,997;
NetMHCIIpan 4.3: 11,048; NetMHC 4.0: 132) were classified with both parsers.
The recorded counts, excluded entries, curated aliases and catalog hashes are
in `tests/data/allele_lists/canonical-audit.json`.

There is no persistent MHC prediction-cache format in mhctools. Its supported
allele map is process-local and now keys on canonicalization policy, parser
version and native catalog content. Model-weight caches and allele-independent
cleavage caches are unaffected. External prediction caches must use the new
policy namespace and preserve the request-to-model provenance exposed by the
adapters. Migration notes in `docs/alleles.md` require invalidation/recomputation
of legacy lossy rows and prohibit guessing missing original inputs.

Topiary's existing PR #578 is being extended to share this corpus and reject
G/P typings, with composed cache/serialization/occurrence coverage. Its null
expression eligibility remains separate. It must pass review, full gates and
release before the mhctools minor release ships.
