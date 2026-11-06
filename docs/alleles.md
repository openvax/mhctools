# Allele names

mhctools parses allele names with [mhcgnomes](https://github.com/pirl-unc/mhcgnomes)
and converts them to each tool's own spelling, so you can write alleles the way
you have them. Do not match alleles with string tricks such as
`startswith("HLA-")`: alleles are not always human.

## What is accepted

Since 4.9, canonical identity is mhcgnomes' lossless spelling. Leading zeros,
all allele fields, and expression suffixes survive. No missing alpha chain is
inferred and no ambiguous group is replaced with one member.

| You write | Canonical prediction identity |
|---|---|
| `HLA-A*02:01`, `A0201`, `hla-a*02:01` | `HLA-A*02:01` |
| `HLA-A*02:01:01:01` | `HLA-A*02:01:01:01` |
| `DLA-88*001:01` | `DLA-88*001:01` |
| `HLA-DRB1*15:01`, `DRB1_1501` | `HLA-DRB1*15:01` |
| `HLA-DQB1*06:02` | `HLA-DQB1*06:02`; a tool requiring a pair needs an explicit alpha chain |
| `HLA-DRA1*01:01-DRB1*15:01` | `HLA-DRA*01:01/DRB1*15:01`, distinct from beta-only |
| `DPA1*01:03/DPB1*04:01`, `DPA1*01:03;DPB1*04:01` | `HLA-DPA1*01:03/DPB1*04:01` |
| `H-2-Kb`, `H2-Kb` | `H2-K*b` |
| `H2-IAb` | `H2-AA*b/AB*b` (curated singleton haplotype) |
| `BoLA-T2b`, `BoLA-T2c` | unchanged; a specific allele need not have an assigned gene |
| `HLA-A30:14L` | `HLA-A*30:14L` |
| `HLA-C*04:22` | unchanged; historical aliases are disabled |
| `A2`, `A02`, `HLA-A*02`, `Mamu-A01` | ambiguity error |
| `HLA-A*02:01:01G`, `HLA-A*02:01P` | group error |
| class-I/class-I pairs, beta/alpha pairs, genes such as `H-2-Qa1` | identity error |

The [official HLA nomenclature](https://hla.alleles.org/pages/nomenclature/naming_alleles/)
defines the fields and expression suffixes. [G groups](https://hla.alleles.org/pages/p%26g_groups/g_groups/)
and [P groups](https://hla.alleles.org/pages/p%26g_groups/p_groups/)
are typing groups, not one full allele identity. Preserving an expression suffix
does not imply that every model or downstream application supports that allele.

MHCflurry preserves the requested canonical identity in its outputs, while
passing the loaded model's actual key to inference. A four-field request can use
an explicitly mapped two-field model without losing its request resolution.
Two distinct requests that use the same model key remain distinct outputs.
For canine `DLA-88*001:01`, mhcgnomes also recognizes `DLA-88*01:01`; mhctools
preserves the canonical three-digit family. See [IPD-MHC's canine section](https://www.ebi.ac.uk/ipd/mhc/group/DLA/).
Supported alleles and percentile calibration still depend on the loaded bundle.

Constructors take a list (`alleles=["HLA-A*02:01", "HLA-B*07:02"]`), a single
string, or a comma-separated string. The command line takes `--mhc-alleles`
(comma- or space-separated) or `--mhc-alleles-file` (one per line); there is no
default allele.

## In each tool's spelling

Command-line and sidecar tools each want their own spelling. mhctools prefers
the spelling in the tool's own allele list (`netMHCpan -listMHC`); only when an
allele is not found there does it format the name itself. Those formatters are
public. Each takes the result of
[`parse_allele`][mhctools.parse_allele] (or call
[`format_allele`][mhctools.format_allele] with a name) and raises
`mhcgnomes.ParseError` when the tool has no spelling for the allele, such as a
class I allele for NetMHCIIpan.

| Formatter | Tool | `HLA-A*02:01` | Class II / mouse |
|---|---|---|---|
| [`netmhc_allele_name`][mhctools.netmhc_allele_name] | NetMHCpan family, most wrappers | `HLA-A02:01` | pairs joined by `-` |
| [`netmhc4_allele_name`][mhctools.netmhc4_allele_name] | NetMHC 4.0 | `HLA-A0201` | class I only |
| [`netmhciipan_allele_name`][mhctools.netmhciipan_allele_name] | NetMHCIIpan | not accepted | `DRB1_0101`, `HLA-DQA10501-DQB10636`, `H-2-IAb` |
| [`deepimmuno_allele_name`][mhctools.deepimmuno_allele_name] | DeepImmuno | `HLA-A*0201` | class I only |
| [`caphla_allele_key`][mhctools.caphla_allele_key] | CapHLA | `HLA-A*02:01` | `HLA-DRB1*15:01`, `HLA-DPA1*01:03/DPB1*04:01`, `H2-IAb` |

mhcgnomes pads non-human families to three digits (`Mamu-B*010:01`) and
renames some species (`Chi-B0401` is `Gaga-BF2*04:01`), while the tools keep
older spellings. For non-human species the tool's own list therefore decides
the spelling and the formatters are only a fallback. Names mhcgnomes cannot
parse follow the tool's list verbatim only after a catalog match. Parsed
genes and groups do not enter this fallback. A verified native catalog alias
such as `Mamu-A01` can serve an explicit request for `Mamu-A1*001:01`; the
ambiguous user spelling itself is still rejected. Already-specific catalog
identities such as `HLA-C*04:22` are never renamed.
[`prefer_allele_spelling`][mhctools.prefer_allele_spelling] picks between
several list spellings of one allele.

## What happens with a bad allele

- A string that mhcgnomes cannot parse (`foo`) raises `mhcgnomes.ParseError`
  from predictors using the shared normalizer. [MHCflurry](predictors/binding.md#mhcflurry) and command-line predictors check an
  unparseable name against the tool's own allele list and raise
  `UnsupportedAllele` when it is not there.
- An allele that parses but the predictor does not support (`HLA-A*99:99`, or a
  class II allele passed to a class I predictor) raises `UnsupportedAllele`
  naming the predictor and, for the [NetMHC](predictors/binding.md#netmhc) family, the command that lists what
  it does support (for example `netMHCpan-4.2 -listMHC`).
- Predictors never silently drop an unsupported allele or return a shorter
  result list.

`AlleleIdentityError` distinguishes a parsed but ambiguous or invalid identity
from `mhcgnomes.ParseError`. Neither strict nor raw-fallback APIs accept a
serotype, G/P group, gene, or invalid pair as an allele. Native catalog ingestion
reports excluded identities and retains unrelated valid entries.

[TLimmuno2](predictors/immunogenicity.md#tlimmuno2) also accepts its native [NetMHCIIpan](predictors/binding.md#netmhciipan)-style keys (`DRB1_0803`,
`HLA-DPA10103-DPB10101`), and [MixMHC2pred](predictors/binding.md#mixmhc2pred) accepts its own spelling
(`DRB1_15_01`, `DQA1_01_02__DQB1_06_02`).

## Which class does an allele need?

Class I predictors ([NetMHCpan](predictors/binding.md#netmhcpan), [MHCflurry](predictors/binding.md#mhcflurry), ...) need class I alleles and class II
predictors ([NetMHCIIpan](predictors/binding.md#netmhciipan), [MixMHC2pred](predictors/binding.md#mixmhc2pred), [TLimmuno2](predictors/immunogenicity.md#tlimmuno2)) need class II. See the
**MHC class** column of the [predictor matrix](predictor-matrix.md). Predictors
that report `mhc_dependence` of `none` ignore alleles entirely; see
[MHC dependence and class](kinds.md#mhc-dependence-and-class).

## Migrating from 4.8 and earlier

Prediction allele columns, joins and cache keys now use the identities in the
table above. Class-II pairs use `/`, mouse alleles use mhcgnomes names, and
beta-only requests stay unpaired. The obsolete tuple-returning
`parse_classi_or_classii_allele_name` API is removed; use `parse_allele` and the
public native formatters. `normalize_allele_name` accepts the old
`infer_class2_pair=False` and `omit_dra1=False` keywords during the 4.9 cycle;
true values raise a migration error. Supply explicit pairs or use a native
formatter when preparing a tool command.

Rebuild supported-allele maps and invalidate old prediction caches before
switching a pipeline to 4.9. The in-process native lookup cache has a new key
including the identity policy, parser version, and catalog contents. mhctools
has no persistent MHC prediction cache to rewrite. Its model-weight caches and
allele-independent cleavage cache formats are unaffected.

External caches should include `allele_identity_provenance()` plus predictor
version, settings, and effective model inputs in their namespace. Command-line
and MHCflurry adapters expose `allele_provenance()` for the exact
request-to-model mapping. This includes historical alias policy at the native
boundary. Keep the provenance alongside serialized tables; do not combine
legacy and canonical-policy rows in one cache namespace.

For existing tables, re-run from trustworthy original inputs when the old
normalizer lost resolution or selected a serotype member. Normalizing an old
output cannot recover discarded fields, the original ambiguity, or whether an
alpha chain was supplied or inferred. If original inputs are unavailable, retain
the recorded pair as a pair with unknown original-input provenance; do not guess
that it represents a beta-only request. Apply the same identity policy to every
side of a downstream join. Topiary's expression eligibility rules remain
separate from this lossless naming contract.


The MixTCRpred catalog also contains four targets specified only as a haplotype,
allele group or serotype. These are excluded with a diagnostic while valid
models remain available; their precise identities require upstream evidence
([#649](https://github.com/openvax/mhctools/issues/649)).
