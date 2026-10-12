"""Assess each candidate vaccine window and every supplied scored core."""

import argparse
import csv
import hashlib
import json
from pathlib import Path

from ..vaccine_windows import assess_vaccine_windows, format_vaccine_window_report, vaccine_window_rows
from .errors import CLI_ERROR_TYPES


def main(argv=None):
    parser = argparse.ArgumentParser(prog="mhctools vaccine-windows", description=__doc__)
    parser.add_argument("--input", required=True, help="Window schema JSON with explicit kinetic profiles")
    parser.add_argument("--output", required=True, help="New JSON with full per-core trajectory audits")
    parser.add_argument("--table", required=True, help="New per-core CSV")
    parser.add_argument("--report", required=True, help="New readable Markdown report")
    args = parser.parse_args(argv)
    try:
        paths = [Path(args.output), Path(args.table), Path(args.report)]
        if len({p.resolve() for p in paths}) != len(paths):
            raise ValueError("Output paths must differ")
        for path in paths:
            if path.exists():
                raise FileExistsError("Output already exists: %s" % path)
            if not path.parent.is_dir():
                raise FileNotFoundError("Output parent directory does not exist: %s" % path.parent)
        data = Path(args.input).read_bytes()
        result = assess_vaccine_windows(json.loads(data))
        result["input_file"] = {"path": args.input, "sha256": hashlib.sha256(data).hexdigest()}
        rows = vaccine_window_rows(result)
        with paths[0].open("x") as stream:
            json.dump(result, stream, indent=2, allow_nan=False)
            stream.write("\n")
        with paths[1].open("x", newline="") as stream:
            writer = csv.DictWriter(stream, fieldnames=list(rows[0]) if rows else ["candidate_id", "status"])
            writer.writeheader()
            writer.writerows(rows)
        with paths[2].open("x") as stream:
            stream.write(format_vaccine_window_report(result))
    except CLI_ERROR_TYPES as error:
        parser.error(str(error))
    print(args.report)


if __name__ == "__main__":
    main()
