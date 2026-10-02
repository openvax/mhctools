# Documentation style

Write for someone who wants to run a predictor or interpret its output.
Start with the task, show the smallest useful example, and explain what the
result means. Put complete lists of parameters and model metadata in reference
pages.

## Structure

The [scikit-learn getting started guide](https://scikit-learn.org/stable/getting_started.html)
introduces its interface through short examples and points readers to the user
guide and API reference for detail. Follow that progression here:

- The homepage describes the library and links to the main guides.
- Getting started takes a reader from installation to a first result.
- User guides explain a task or predictor family, including its limitations.
- Reference pages provide complete fields, classes, commands, and metadata.

Use descriptive headings such as “Predict for peptides” or “Read the results”.
The [pandas tutorials](https://pandas.pydata.org/docs/getting_started/intro_tutorials/index.html)
are a useful model for organizing pages around reader questions.

Avoid repeating the sidebar or table of contents in the page body. A short
list of related pages is useful when it explains where to go next.

## Prose

Use ordinary paragraphs for explanations, numbered lists for ordered steps,
and bullets for parallel choices. Use tables when readers need to compare
short entries. Move long explanations below a table instead of squeezing
paragraphs into cells.

Describe behavior directly. Replace “Every predictor answers predict” with
“Each predictor provides a prediction method”. Remove introductions that
repeat the heading, promotional claims, and commentary about implementation
unless it affects how a reader uses the API.

Keep scientific qualifications close to the claim they constrain. Preserve
units, input restrictions, assay context, and source links when shortening a
page. Distinguish author-reported performance from independent validation.

## Inline code

Use inline code for literal syntax a reader can type or look up:

| Use code | Use ordinary text |
|---|---|
| Methods: `predict()`, `predict_dataframe()` | Model names: MHCflurry, NetMHCpan, Calis |
| Parameters: `alleles=`, `peptide_lengths=` | Concepts: binding affinity, presentation, peptide length |
| Fields: `result.affinity`, `percentile_rank` | Units: nM, hours, residues |
| Commands: `mhctools fetch mhcflurry` | Tool names in prose and navigation links |
| Literal values: `None`, `"pMHC_affinity"` | Descriptions of those values |

A Python class belongs in code when discussing the interface itself, such as
“`PeptideResult` contains a tuple of predictions” or listing the exported
classes in the predictor matrix. When referring to the upstream model, use
its ordinary name. Write “NetMHCpan 4.2”, rather than the wrapper spelling
`NetMHCpan42`, unless the wrapper is the subject.

Do not combine bold and code styling just to emphasize a field name. Explain
the field in a sentence or give it a row in a reference table.

## Examples and presentation

Prefer complete, copyable examples with imports and real inputs. If an example
needs an existing predictor or downloaded model, state that requirement and
link to its setup. Separate code from its explanation; do not turn comments
into paragraphs. Label illustrative fragments that are not runnable.

Use the site's existing typography and navigation. Keep prose to a readable
line length and let wide reference tables scroll. Add a callout only when the
reader needs a distinct warning or note to use the example correctly.

## Checks

Before submitting documentation changes, run:

```sh
python scripts/predictor_matrix.py --check
python scripts/check_docs_links.py
mkdocs build --strict
./lint.sh
./test.sh
```

Edit the predictor matrix in
[scripts/predictor_matrix.py](https://github.com/openvax/mhctools/blob/master/scripts/predictor_matrix.py)
and regenerate it. Review the built pages in a browser, including a narrow
window. Keep existing page URLs and linked heading anchors when practical;
update internal links when a heading changes.
