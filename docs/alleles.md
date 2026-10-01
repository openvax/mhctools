# Allele names

mhctools parses allele names with [mhcgnomes](https://github.com/pirl-unc/mhcgnomes)
and converts them to each tool's own spelling, so you can write alleles the way
you have them. Do not match alleles with string tricks such as
`startswith("HLA-")`: alleles are not always human.

## What is accepted

| You write | mhctools uses |
|---|---|
| `HLA-A*02:01`, `A*02:01`, `HLA-A0201`, `A0201`, `hla-a*02:01` | `HLA-A*02:01` |
| `A2` | `HLA-A*02:01` (a serotype resolves to its representative allele) |
| `HLA-DRB1*15:01`, `DRB1*15:01`, `DRB1_1501` | `HLA-DRA1*01:01-DRB1*15:01` (the invariant DRA1 chain is added) |
| `HLA-DPA1*01:03-DPB1*04:01`, `DPA1*01:03/DPB1*04:01`, `DPA1*01:03;DPB1*04:01` | `HLA-DPA1*01:03-DPB1*04:01` |
| `H-2-Kb`, `H2-Kb` | `H-2-Kb` |

Write class II DP and DQ molecules as an alpha-beta pair. A single-chain DR
allele is expanded to the pair automatically.

Constructors take a list (`alleles=["HLA-A*02:01", "HLA-B*07:02"]`), a single
string, or a comma-separated string. The command line takes `--mhc-alleles`
(comma- or space-separated) or `--mhc-alleles-file` (one per line); there is no
default allele.

## What happens with a bad allele

- A string that mhcgnomes cannot parse (`foo`) raises `mhcgnomes.ParseError`
  from predictors such as `MHCflurry`. Command-line predictors instead check an
  unparseable name against the tool's own allele list and raise
  `UnsupportedAllele` when it is not there.
- An allele that parses but the predictor does not support (`HLA-A*99:99`, or a
  class II allele passed to a class I predictor) raises `UnsupportedAllele`
  naming the predictor and, for the NetMHC family, the command that lists what
  it does support (for example `netMHCpan-4.2 -listMHC`).
- Predictors never silently drop an unsupported allele or return a shorter
  result list.

Some tools accept names mhcgnomes cannot parse, such as `H-2-Qa1` or
`BoLA-amani.1`. Command-line predictors validate these against the tool's own
list, so they can be requested; `MHCflurry` rejects `H-2-Kb` with
`UnsupportedAllele`.

`TLimmuno2` also accepts its native NetMHCIIpan-style keys (`DRB1_0803`,
`HLA-DPA10103-DPB10101`), and `MixMHC2pred` accepts its own spelling
(`DRB1_15_01`, `DQA1_01_02__DQB1_06_02`).

## Which class does an allele need?

Class I predictors (NetMHCpan, MHCflurry, ...) need class I alleles and class II
predictors (NetMHCIIpan, MixMHC2pred, TLimmuno2) need class II. See the
**MHC class** column of the [predictor matrix](predictor-matrix.md). Predictors
that report `mhc_dependence` of `none` ignore alleles entirely; see
[MHC dependence and class](kinds.md#mhc-dependence-and-class).
