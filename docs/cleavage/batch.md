# Batch cleavage assessments

`predict_cleavage_batch(inputs, scenarios)` assesses explicit sequence
occurrences and overlays evidence on intended epitopes. The CLI uses the
same implementation:

```sh
mhctools cleavage --input docs/cleavage-batch-example.json \
  --out cleavage.json --html cleavage.html
```

The example uses only built-in models. Add exact optional model names from
`mhctools cleavage --list-models` to a scenario to run proteasome or ERAP1
models. Missing assets and failed runtimes produce visible failed records.
Python callers can set `raise_on_error=True` to fail immediately.

```python
from mhctools import predict_cleavage_batch, write_cleavage_batch, load_cleavage_batch

report = predict_cleavage_batch(
    inputs=[{
        "id": "native-occurrence-1",
        "peptide": "RPPGFSPFR", "n_flank": "AA", "c_flank": "GG",
        "source_id": "protein-1", "source_start": 12,
        "evidence": [{"source": "input-table", "affinity_nM": 40.0}],
    }],
    scenarios=[{
        "id": "tumor-cytosol", "context": "tumor",
        "compartments": ["cytosol"],
        "models": ["pepsickle-in-vitro-2-human-only-constitutive"],
    }],
)
write_cleavage_batch(report, "cleavage.json", html_path="cleavage.html")
assert load_cleavage_batch("cleavage.json") == report
```

## Inputs and coordinates

Each input needs a unique `id` and either a complete `sequence`, or a
`peptide` plus explicit `n_flank` and `c_flank` strings (empty strings are
allowed). Peptide mode creates one epitope interval automatically. Sequence
mode accepts `epitopes`, each with an `id`, `start`, `end`, optional exact
`sequence`, and any additional JSON annotations.

**All batch intervals are zero-based and half-open.** Bond `b` separates
`sequence[:b]` from `sequence[b:]`. `source_start` is the offset of the whole
supplied sequence, including any N flank; it is not the epitope's offset.
Source bond coordinates add `source_start`. The older vaccine-report schema
uses one-based inclusive epitope windows; convert explicitly when using it.

`scope` is `native_window` (default), `protein`, or `construct`. A native
reconstruction window must retain unknown molecular termini. A construct
can declare `n_term` and `c_term`; both default to `unknown`. Do not turn
native-window edges into exposed substrate ends to obtain peptidase scores.
Pepsickle can assess native, protein, construct or conditional-fragment sequence
context with free or unknown termini (including mixed free/unknown states);
the result explicitly records that sequence-only assumption. Known chemical
modifications still cause abstention in sequence-only models.

For a final vaccine, supply the entire manufactured/translated sequence,
including linkers and junctions, and its intended epitope intervals. Retain
delivery/routing/formulation facts in input annotations and choose applicable
scenarios explicitly. The batch runner does not infer intracellular routing
from a sequence or assume an injected peptide reaches circulation.

Imported `evidence` and other JSON fields survive normalization, inference,
save/reload and reporting unchanged. They do not automatically become scores
or independent corroboration. Generalized source-file reconciliation and
ranking are Topiary responsibilities, and construct selection is Vaxrank's; see
[known gaps](../known-gaps.md#downstream-integration).

## Scenarios and trimming

Each scenario has a unique `id`, `context` (`tumor`, `apc`, `extracellular`),
explicit `models`, optional `compartments`, and optional `enzyme_states`
(currently activation-gated CPB2). Compartments filter known enzyme locations;
the scenario label does not change predictor weights or simulate a cell.
Additional scenario annotations, including pH, are descriptive; they do not
silently tune a model. A future condition-dependent predictor must declare
and validate its supported settings explicitly.

An input may contain explicit hypothetical `fragments`, each with `id`,
`start`, `end`, `n_term`, `c_term` and a nonempty production `assumption`.
Retained parent termini must retain their chemistry. Results for these
fragments retain their condition on every epitope overlay. No cascade order,
production probability, equilibrium or competition is inferred.

Each overlay preserves N and C boundaries separately, all internal bonds,
native scores/statuses and whole-substrate observations. Sequence endpoints
are `sequence_endpoint`; unreturned bonds are `unassessed`. A motif non-match
does not become a zero probability or proof of resistance. Failed models
have explicit error records rather than empty successful tracks.

## Existing vaccine reports

`CleavageTrack.from_result(result, context, full_sequence, start=0)` creates a
track without inventing numbers for categorical evidence. A fragment track
also requires `conditional_on="..."`. Full canonical evidence survives the
normalized manifest and placement assessments; the PDF uses circles for
motif matches, crosses for non-matches and squares for source-reported sites.
Native numerical scores retain their own optional thresholds.

## Experimental source panels

The Python API and CLI batch request accept `reference_panels`: a list of
`{"model": ..., "cases": [...]}` catalogs using the existing
`PeptidaseSubstrateReference` schema. Name each assay/condition panel separately
and include its model name in the requested scenario. The model's evidence
must be `substrate_reference`; every case requires exact sequence/termini,
source URL, source measurement ID, assay conditions, supported bonds,
interpretation and substrate observation. The complete panel is retained in
the output. Catalog names cannot override built-in models.

This provides an input mode for cathepsin/AEP experiments while transferable
prediction remains unavailable. A source lookup abstains on every unobserved
chemical form. Never fill omitted cleavage sites with experimental negatives.
See [the biological coverage review](validation.md).

## Guarantees

- Preserve input occurrence IDs, source offsets, chemistry, sequence scope,
  epitope intervals, imported predictions and arbitrary source annotations.
  Native windows do not establish exposed molecular termini.
- Use zero-based, half-open epitope/fragment intervals and existing cleavage
  bond coordinates. Reject inconsistent peptide/sequence/interval combinations.
- Scenarios explicitly name models and compartments. Default panels report
  their limitations and missing biological coverage; optional models never
  silently disappear or substitute for requested models.
- Hypothetical fragments require explicit terminal chemistry and an assumption
  explaining their production. Their assessments remain conditional.
- Preserve categorical motif decisions, quantitative scores, substrate-only
  observations, unsupported inputs and runtime failures separately. Overlay
  internal bonds and the N/C boundaries without inventing terminal bonds.
- JSON round-trips retain the original inputs and evidence. Human reports
  consume these same records and do not assign aggregate protection scores.
- Expose upstream Pepsickle epitope and in-vitro model families with explicit
  constitutive/immunoproteasome selection, exact artifact identity and native
  score semantics. Reject settings the upstream model ignores.
- Review cathepsin/AEP evidence with pH, activation and assay scope. Include
  source-linked observations only where the original sequence and bond are
  verified. Missing transferable models remain explicit coverage gaps.
