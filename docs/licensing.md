# Licensing

## The mhctools license

mhctools adapter code uses the [Apache License 2.0](https://github.com/openvax/mhctools/blob/master/LICENSE).
It permits use, modification, and redistribution, including commercial use.
You can include mhctools in an application without publishing that application's
source code.

When redistributing mhctools or a modified version:

- Include a copy of the license.
- Retain applicable copyright and attribution notices, including notices from
  a supplied NOTICE file.
- Mark files you changed.

The license includes a limited contributor patent grant. It does not grant
trademark rights, and the software comes without warranties. The
[full terms](https://www.apache.org/licenses/LICENSE-2.0) and
[Apache licensing FAQ](https://www.apache.org/foundation/license-faq.html)
explain these conditions.

Bundled predictor data retain their own licenses: the
[ITCell profiles](cleavage/itcell.md) are LGPL-2.1-only, and the
[PhageScout matrices and peptide profiles](cleavage/phagescout.md) are
CC BY 4.0. Their license texts and attribution notices are included in the
distribution. Package metadata records the combined license expression;
the independent PhageScout adapter does not include the unlicensed author
notebook code.

<a id="the-five-tiers"></a>

## Predictor licenses

Each upstream program, model, or dataset has its own terms. The mhctools
license does not grant permission to use those materials. Check the
predictor's linked upstream license before installing it or redistributing
its code or weights.

The [predictor matrix](predictor-matrix.md) gives every predictor one of these:

| Tier | Meaning | What you do |
|---|---|---|
| `builtin` | Built into mhctools; nothing to download | Nothing. [Calis](predictors/immunogenicity.md#calis) and `RandomBindingPredictor` |
| `open` | Open source; code and weights can be fetched or installed freely | `mhctools fetch <name>` or `pip install` |
| `academic` | Academic / non-commercial terms, often with no redistribution | Read the upstream license, then `mhctools fetch <name> --accept-license` or install it yourself |
| `dtu` | DTU academic license, bound to an identity | Request it from DTU; mhctools calls your installation |
| `unlicensed` | Upstream publishes no license | `--accept-license` records that you have confirmed your own use is authorized |

## What `--accept-license` means

It records that you reviewed the terms and confirmed that your own use is
authorized. It does not grant rights mhctools does not have, and it cannot stand
in for a license you must request yourself.

The DTU downloads ([NetMHCpan](predictors/binding.md#netmhcpan), [NetMHC](predictors/binding.md#netmhc), [NetMHCcons](predictors/binding.md#netmhccons), [NetMHCIIpan](predictors/binding.md#netmhciipan), [NetMHCstabpan](predictors/binding.md#netmhcstabpan),
[NetChop](predictors/processing.md#netchop)) are the clearest example. DTU requires a name, position, academic
email, affiliation and acceptance, then sends a private link, so those
installations stay `manual` in `mhctools ls` and `fetch` will not install them.

[NetCleave](predictors/processing.md#netcleave) and [TLimmuno2](predictors/immunogenicity.md#tlimmuno2) publish no license at all. mhctools can fetch a
pinned snapshot, but only after `--accept-license`; the recorded manifest says
`"license": "none published"`.

## Copyleft tools

[ERAMER](predictors/processing.md#eramer), [Tulip](predictors/tcr.md#tulip) and [PlifePred2](predictors/peptide-pk.md#plifepred2) are GPLv3. mhctools vendors none of them and
does not import them into its own interpreter: ERAMER's workbook is read at
runtime, TULIP runs out of process through its own `predict.py`, and
PlifePred2 and Pfeature run in a separate interpreter.

## Ambiguous upstream licenses

[PeptiVerse](predictors/peptide-pk.md#peptiverse) declares Apache-2.0 on its model card and MIT in its README;
mhctools lists it as `open` but you should confirm which applies to your use.

## Per-predictor details

| Topic | Where |
|---|---|
| [MixMHCpred](predictors/binding.md#mixmhcpred), [MixMHC2pred](predictors/binding.md#mixmhc2pred), [PRIME](predictors/immunogenicity.md#prime), [MixTCRpred](predictors/tcr.md#mixtcrpred) terms | their sections in [binding](predictors/binding.md#mixmhcpred), [immunogenicity](predictors/immunogenicity.md#prime) and [TCR](predictors/tcr.md#mixtcrpred) |
| [SMM](predictors/binding.md#smm-and-smm-pmbec) bundle (Non-Profit Open Software License 3.0) | [SMM setup](backends.md#local-smm-and-smm-pmbec) |
| How `fetch` records provenance | [getting models](artifacts.md#licensing) |
