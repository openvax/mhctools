# Simple input validation (#535)

Replace mhctools' undeclared `typechecks` dependency with two small private
`isinstance` helpers. Keep runtime `TypeError` validation and existing string
and integer acceptance, including Python's existing bool/int relationship.
Errors identify the argument and item index without dumping input objects.

Materialize iterable inputs before validation wherever they are reused. Check
peptides in the existing validation loop. Annotate the affected input boundaries;
repository-wide typing and a `py.typed` declaration require a separate migration.

Verify generator alleles, lengths, flags, peptides and flanks, invalid elements,
and importing/using mhctools with `typechecks` unavailable. Run lint and the full
required-backend suite, bump the patch version, and ship through green CI and
the normal release script.
