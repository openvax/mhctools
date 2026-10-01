# Licensing

mhctools is Apache-2.0. The predictors it wraps are not, and mhctools never
bundles another project's code or weights: it fetches a pinned snapshot, calls
an installation you provide, or reimplements a published model. Check the
upstream license before you use any of them, particularly commercially.

## The five tiers

The [predictor matrix](predictor-matrix.md) gives every predictor one of these:

| Tier | Meaning | What you do |
|---|---|---|
| `builtin` | Built into mhctools; nothing to download | Nothing. `Calis` and `RandomBindingPredictor` |
| `open` | Open source; code and weights can be fetched or installed freely | `mhctools fetch <name>` or `pip install` |
| `academic` | Academic / non-commercial terms, often with no redistribution | Read the upstream license, then `mhctools fetch <name> --accept-license` or install it yourself |
| `dtu` | DTU academic license, bound to an identity | Request it from DTU; mhctools calls your installation |
| `unlicensed` | Upstream publishes no license | `--accept-license` records that you have confirmed your own use is authorized |

## What `--accept-license` means

It records that you reviewed the terms and confirmed that your own use is
authorized. It does not grant rights mhctools does not have, and it cannot stand
in for a license you must request yourself.

The DTU downloads (NetMHCpan, NetMHC, NetMHCcons, NetMHCIIpan, NetMHCstabpan,
NetChop) are the clearest example. DTU requires a name, position, academic
email, affiliation and acceptance, then sends a private link, so those
installations stay `manual` in `mhctools ls` and `fetch` will not install them.

`NetCleave` and `TLimmuno2` publish no license at all. mhctools can fetch a
pinned snapshot, but only after `--accept-license`; the recorded manifest says
`"license": "none published"`.

## Copyleft tools

`ERAMER`, `Tulip` and `PlifePred2` are GPLv3. mhctools vendors none of them and
does not import them into its own interpreter: ERAMER's workbook is read at
runtime, TULIP runs out of process through its own `predict.py`, and
PlifePred2 and Pfeature run in a separate interpreter.

## Ambiguous upstream licenses

`PeptiVerse` declares Apache-2.0 on its model card and MIT in its README;
mhctools lists it as `open` but you should confirm which applies to your use.

## Per-predictor details

| Topic | Where |
|---|---|
| MixMHCpred, MixMHC2pred, PRIME, MixTCRpred terms | their sections in [binding](predictors/binding.md#mixmhcpred), [immunogenicity](predictors/immunogenicity.md#prime) and [TCR](predictors/tcr.md#mixtcrpred) |
| SMM bundle (Non-Profit Open Software License 3.0) | [SMM setup](backends.md#local-smm-and-smm-pmbec) |
| How `fetch` records provenance | [getting models](artifacts.md#licensing) |
