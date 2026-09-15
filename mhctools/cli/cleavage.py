"""JSON reporting for per-bond cleavage evidence."""

import argparse
from dataclasses import asdict
import json
from pathlib import Path

from mhctools.cleavage import CleavageInput
from mhctools.peptidases import cleavage_models, predict_cleavage


def _format_models(models):
    """Render the model catalog as a compact, aligned table."""
    headers = ("NAME", "ENZYME", "COMPARTMENTS", "EVIDENCE", "SCORE UNITS")
    rows = [
        (
            model.name,
            model.enzyme,
            ",".join(model.compartments),
            model.evidence,
            model.score_units or "-",
        )
        for model in models
    ]
    widths = [
        max(len(headers[i]), *(len(row[i]) for row in rows))
        for i in range(len(headers))
    ]
    lines = [
        "  ".join(value.ljust(widths[i]) for i, value in enumerate(headers))]
    lines.extend(
        "  ".join(value.ljust(widths[i]) for i, value in enumerate(row))
        for row in rows)
    return "\n".join(lines)


def _compact_report(predictions):
    """Build schema v2 with model metadata stored once per model name."""
    models = {}
    results = []
    for prediction in predictions:
        row = prediction.to_dict()
        model = row.pop("model")
        name = model["name"]
        if name in models and models[name] != model:
            raise ValueError("Conflicting metadata for cleavage model %r" % name)
        models[name] = model
        peptide = row.pop("peptide")
        results.append({"peptide": peptide, "model": name, **row})
    return {"schema_version": 2, "models": models, "results": results}


def _round_floats(value):
    """Round JSON float leaves to six significant digits without stringifying."""
    if isinstance(value, float):
        return float("%.6g" % value)
    if isinstance(value, dict):
        return {key: _round_floats(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [_round_floats(item) for item in value]
    return value


def main(argv=None):
    parser = argparse.ArgumentParser(
        prog="mhctools cleavage",
        description="Report peptidase motif evidence and native model scores "
                    "per bond, as JSON.")
    parser.add_argument("--list-models", action="store_true",
                        help="List the curated cleavage models with their "
                             "enzymes, compartments and evidence, then exit")
    parser.add_argument("--json", action="store_true",
                        help="With --list-models, emit the full JSON catalog "
                             "instead of the human-readable table")
    parser.add_argument("--sequence", action="append", help="Canonical linear L-peptide; repeat for multiple inputs")
    parser.add_argument("--model", action="append",
                        help="Exact model name (see --list-models); repeat to "
                             "select multiple models. Defaults to every "
                             "always-available model for --compartment; "
                             "models needing a fetched asset (eramer-step) "
                             "must be named explicitly.")
    parser.add_argument("--compartment", help="Filter enzyme locations: serum, plasma, extracellular, cytosol, er, endosome")
    parser.add_argument("--n-term", choices=("free", "acetylated", "unknown"), default="free",
                        help="N-terminal chemistry of the input peptide "
                             "(default: %(default)s)")
    parser.add_argument("--c-term", choices=("free", "amidated", "unknown"), default="free",
                        help="C-terminal chemistry of the input peptide "
                             "(default: %(default)s)")
    parser.add_argument("--source-id",
                        help="Identifier of the protein the peptide came from, "
                             "echoed back in the report")
    parser.add_argument("--source-start", type=int, default=0,
                        help="Offset of the peptide within --source-id, used "
                             "to report source coordinates "
                             "(default: %(default)s)")
    parser.add_argument("--enzyme-state", action="append", default=[],
                        metavar="ENZYME=STATE",
                        help="Declare one enzyme's state, e.g. "
                             "'CPB2=active'; STATE is active, zymogen, "
                             "inactive or unknown. Repeat per enzyme.")
    parser.add_argument("--out", help="Write JSON to this path instead of stdout")
    args = parser.parse_args(argv)
    if args.list_models:
        if args.sequence or args.model or args.compartment or args.enzyme_state:
            parser.error("--list-models cannot be combined with prediction options")
        models = cleavage_models(include_optional=True)
        if not args.json and not args.out:
            print(_format_models(models))
            return
        # Preserve the original machine-readable catalog shape. ``--out``
        # implies JSON, matching its longstanding "Write JSON" contract.
        result = {
            "models": [asdict(model) for model in models],
            "schema_version": 1,
        }
    else:
        if not args.sequence:
            parser.error("Provide --sequence or --list-models")
        try:
            states = {}
            for value in args.enzyme_state:
                enzyme, state = value.split("=", 1)
                if enzyme in states:
                    raise ValueError("Duplicate enzyme state for %s" % enzyme)
                states[enzyme] = state
            predictions = [
                prediction
                for seq in args.sequence
                for prediction in predict_cleavage(CleavageInput(
                    seq, args.n_term, args.c_term, args.source_id, args.source_start),
                    models=args.model, compartment=args.compartment,
                    enzyme_states=states)]
            result = _compact_report(predictions)
        except Exception as error:
            # Deliberately broad: this is the CLI's user-facing error
            # boundary. Resolving or running a named model can load an
            # external asset (e.g. eramer-step's PWM workbook), and any
            # failure there should exit cleanly with parser.error() rather
            # than a raw traceback, regardless of the exception type raised
            # by that asset's own loader.
            parser.error(str(error))
    output = json.dumps(_round_floats(result), indent=2, allow_nan=False) + "\n"
    if args.out:
        Path(args.out).write_text(output, encoding="utf-8")
    else:
        print(output, end="")
