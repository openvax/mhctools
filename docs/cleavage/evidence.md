# Reading cleavage evidence

How to interpret what the cleavage models return. None of it is a probability
that a peptide survives or is presented.

## Interpreting rule-based evidence

The shared result contract also supports `motif_rule` evidence, reported as
`matched` or `not_matched` with no numerical score. A rule match indicates a
known recognition pattern; a non-match does not establish resistance.
`unsupported_reason` means the input could not be assessed and carries no
score. These states must remain separate in downstream displays and ranking.

## Building a position track for a full sequence

A downstream tool (topiary, vaxrank) that wants to overlay cleavage evidence
from several models onto one parent sequence, as a track indexed by
position, reads `source_bond` from `to_dict()["sites"]` and keys everything
on it. That coordinate is `source_start + bond`, so it is the same absolute
number no matter which model or which fragment produced the site:

```python
track = {}
for result in predict_cleavage(peptide, compartment="cytosol"):
    for site in result.to_dict()["sites"]:
        track.setdefault(site["source_bond"], []).append(
            {"model": result.model.name, "status": site["status"], "score": site["score"]})
```

Two model shapes contribute to that track differently, and a caller building
one has to model both:

- **Internal topology** (`mme-hydrophobic`, `fap-endo-gp`, `prep-pro`)
  assesses every internal bond of whatever peptide it is given in one
  `predict()` call. Calling it once on the full input already produces a
  multi-position run of track entries.
- **Terminal topology** (the aminopeptidases, carboxypeptidases, DPP-family,
  `dpp4-qpisa` and `eramer-step`) only ever assesses the *currently exposed*
  end of its input: always exactly one bond, regardless of peptide length.
  It cannot tell you whether a bond in the middle of a long precursor is a
  plausible trimming stop; it can only assess a candidate fragment you
  supply. To extend a track with these models, model the hypothesized
  trimming step explicitly with `parent.fragment(start, end, n_term=...,
  c_term=...)` and predict on that fragment. Its `source_bond` still lands
  on the parent's absolute coordinates, so it merges into the same track.
  Producing a track over a whole precursor this way means enumerating
  candidate fragments yourself; the model does not search for them.

The same rule applies to `substrate_reference` evidence (THOP1, neurolysin,
IRAP): it never invents a bond, so a source case that reports degradation
without pinning one contributes no track entry, not a guessed position.

Because motif and source-reference models never emit numerical scores, and
because `not_matched`/`no_cleavage_detected` are not probabilities of
resistance, a track built this way is evidence to overlay and inspect, not a
single per-position score to rank or threshold.

## How strict is each motif?

A recognition pattern is only as informative as the evidence behind it, so
every motif rule carries a `motif_strictness` grade and a `strictness_basis`
naming the source observation that supports the grade. Read a decision through
the grade rather than treating all matches alike:

| Grade | What a match means | What a non-match means |
| --- | --- | --- |
| `required` | The pattern is necessary for this activity, so a match clears a real gate. Rates, exposure and competition still decide the outcome. | Meaningful evidence against cleavage by this enzyme through this route. Still not proof of resistance. |
| `preferred` | A favoured context, not a gate. | Weak evidence. Non-matching bonds are cleaved, usually more slowly. |
| `permissive` | Little information; the rule mostly describes topology or a broad enzyme. | Almost no information. |

Grades apply to motif rules only. Scored models report native units instead,
and source references report what an experiment observed. Neither carries a
grade, and both leave `motif_strictness` null.
