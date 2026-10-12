"""Compare existing vaccine trajectory results with explicit endpoint matching."""

import argparse
import hashlib
import json
from pathlib import Path

from ..vaccine_trajectory_report import compare_vaccine_trajectories, format_vaccine_comparison_report
from .errors import CLI_ERROR_TYPES


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__, prog="mhctools vaccine-compare")
    parser.add_argument("--input", required=True, nargs="+", help="Two or more trajectory RESULT JSON files, in display order")
    parser.add_argument("--output", required=True, help="New comparison JSON")
    parser.add_argument("--report", required=True, help="New readable Markdown comparison")
    args = parser.parse_args(argv)
    try:
        outputs = [Path(args.output), Path(args.report)]
        if outputs[0].resolve() == outputs[1].resolve():
            raise ValueError("JSON and report paths must differ")
        for path in outputs:
            if path.exists():
                raise FileExistsError("Output already exists: %s" % path)
        inputs = [Path(path).read_bytes() for path in args.input]
        comparison = compare_vaccine_trajectories(json.loads(data) for data in inputs)
        comparison["input_files"] = [{"path": path, "sha256": hashlib.sha256(data).hexdigest()}
                                     for path, data in zip(args.input, inputs)]
        with outputs[0].open("x") as stream:
            json.dump(comparison, stream, indent=2, allow_nan=False)
            stream.write("\n")
        with outputs[1].open("x") as stream:
            stream.write(format_vaccine_comparison_report(comparison))
    except CLI_ERROR_TYPES as error:
        parser.error(str(error))
    print(args.report)


if __name__ == "__main__":
    main()
