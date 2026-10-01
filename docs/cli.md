# Command line

`mhctools` is one executable with subcommands. Run `mhctools <command> --help`
for every option.

| Command | Does | More |
|---|---|---|
| `mhctools` | Predict for peptides or FASTA sequences | below |
| `mhctools ls` | List model artifacts and where they live | [getting models](artifacts.md) |
| `mhctools fetch` | Download model weights and tool snapshots | [getting models](artifacts.md) |
| `mhctools predictors` | Report which predictors can actually run (`integrations` is an alias) | [getting models](artifacts.md#inventory-is-not-capability) |
| `mhctools predict-table` | Append predictor scores to a CSV | [below](#annotate-a-table-predict-table) |
| `mhctools cleavage` | Per-bond peptidase evidence | [cleavage](cleavage/index.md) |
| `mhctools vaccine-report` | Route-aware vaccine construct report | [vaccine reports](vaccine-reports.md) |
| `mhctools benchmark` | Assay-aware model evaluation | [benchmarks](benchmarks.md) |
| `mhctools mixtcrpred` | Score paired TCRs against a fixed target | [TCR predictors](predictors/tcr.md#mixtcrpred) |

The default command takes `--mhc-predictor` (one or more names, space or comma
separated), `--mhc-alleles` or `--mhc-alleles-file`, and one input. The
accepted names are listed by `mhctools --help` and in the **CLI names** column
of the [predictor matrix](predictor-matrix.md). TCR predictors have no name
here because their input is a peptide plus a TCR.

## Predict for peptides you supply

```sh
mhctools --sequence SIINFEKL SIINFEKLQ --mhc-predictor netmhc --mhc-alleles A0201
```

`--sequence` may be repeated and all occurrences accumulate. Or use
`--input-peptides-file` for one peptide per line (blank lines ignored), or
`--input-fasta-file` for protein sequences. Pick exactly one of the three.

## Extract subsequences automatically

```sh
mhctools --sequence AAAQQQSIINFEKL --extract-subsequences \
    --mhc-peptide-lengths 8-10 --mhc-predictor mhcflurry --mhc-alleles A0201
```

## Annotate a table (`predict-table`)

Reads a CSV, runs each requested predictor once, and appends one score column
per predictor — choosing the best allele per row — while preserving every input
column:

```sh
mhctools predict-table \
    --input benchmark.csv.bz2 \
    --peptide-column peptide \
    --alleles-column hla \
    --predictor netmhcpan42-ba:netmhcpan4.2.ba:affinity \
    --predictor netmhcpan42-el:netmhcpan4.2.el:score \
    --out benchmark.with_scores.csv.bz2
```

Each `--predictor` spec is `NAME[:OUTPUT_COLUMN[:FIELD]]`, where `FIELD` is
`affinity`, `score`, or `percentile_rank`. Lower is better for `affinity` and
`percentile_rank`, higher for `score`.

A row may hold several alleles per cell (whitespace-, comma-, or
semicolon-separated); the best one per peptide is chosen and recorded in a
`<OUTPUT_COLUMN>_best_allele` provenance column. Missing or blank
peptide/allele cells stay unscored — they are never coerced into a literal
sequence or allele string and sent to a predictor.

Pass `--predictor-info info.csv` to also write a sidecar describing each
column's `score_field`, `units`, and `higher_is_better`. Empty `units` means
the field is dimensionless or predictor-specific.

## Output conventions

CLI prediction tables follow one convention across every predictor. Plain
peptide inputs get an empty `source_sequence_name` and offset `0`, while FASTA
and subsequence inputs keep their source and zero-based offset and are ordered
by those coordinates. `prediction_method_name` is the exact CLI predictor name
you selected, including version and mode. `affinity` is IC50 in nM,
`percentile_rank` is a 0–100 percentile, and `score` stays
predictor-specific. CSV floats are serialized with six significant digits.

Default stdout is streamed as tab-separated values, so an empty source name
survives as an empty field and large tables don't need a second formatted copy
in memory. A downstream closed pipe (`| head`, say) exits cleanly.
