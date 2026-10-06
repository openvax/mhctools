# Cleavage mechanism scorecards

Render the verified, frozen Sid vaccine cleavage results as a site-by-model
PDF and a long-form CSV, preserving every native score and assessment state.
Each model declares its biological route, enzyme attribution, training/assay
endpoint, bond topology, units, exposure assumptions, and primary references.
NetChop Cterm is ligand-trained processing evidence; NetChop 20S is an
in-vitro proteasome model. Neither assigns a catalytic subunit. DPP4 scores
only the exposed N-terminal dipeptide bond of the intact input. Motif matches
remain categorical; unavailable/unassessed sites never become zero.

For the human Sid dataset, prefer human-only Pepsickle. An explicit organism
option selects all-mammal for nonhuman, mixed, or unknown context, while both
recorded models remain available without combining their scores or votes.
Source-listed minimal epitopes are highlighted only where their recorded
offset matches the input sequence; no missing epitope is inferred.

The renderer performs no inference or network requests. It verifies source
checksums, writes a new immutable output directory, and records input/model
metadata and output checksums. Validate coordinate joins, scores and missing
states with focused offline tests; visually inspect the PDF. Run repository
lint/tests and docs checks, then ship the versioned PR through CI, merge and
PyPI deployment as required by AGENTS.md.
