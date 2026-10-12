"""JSON interface to explicitly parameterized SLP/mRNA pMHC scenarios."""

import argparse
import json
from pathlib import Path

from ..vaccine_trajectory import (
    VaccineTrajectoryInput, simulate_vaccine_trajectory,
)
from ..vaccine_windows import fragment_cleavage_table
from ..vaccine_trajectory_report import format_vaccine_trajectory_report
from .errors import CLI_ERROR_TYPES


def main(argv=None):
    parser = argparse.ArgumentParser(
        prog="mhctools vaccine-trajectory",
        description="Conditional local SLP/secreted-mRNA trajectories to APC pMHC; requires explicit kinetics.")
    parser.add_argument("--input", required=True, help="Trajectory schema JSON, including observation times")
    parser.add_argument("--output", required=True, help="Output JSON with curves and full mechanism/assumption audit")
    parser.add_argument("--report", help="Optional new Markdown report with dose, presentation exposure and missing inputs")
    args = parser.parse_args(argv)
    try:
        paths = [Path(args.output)] + ([Path(args.report)] if args.report else [])
        if len({path.resolve() for path in paths}) != len(paths):
            raise ValueError("JSON and report paths must differ")
        for path in paths:
            if path.exists():
                raise FileExistsError("Output already exists: %s" % path)
        value = json.loads(Path(args.input).read_text())
        model = VaccineTrajectoryInput.from_dict(value)
        sequence = model.validate().sequence
        output = simulate_vaccine_trajectory(
            model, value["times_hours"], fragment_cleavage_table(sequence, value.get("cleavage_assessments", [])))
        # A path's parent must exist; do not silently overwrite a prior analysis.
        with Path(args.output).open("x") as stream:
            json.dump(output, stream, indent=2, allow_nan=False)
            stream.write("\n")
        if args.report:
            with Path(args.report).open("x") as stream:
                stream.write(format_vaccine_trajectory_report(output))
    except CLI_ERROR_TYPES as error:
        parser.error(str(error))
    print(args.output)
    if args.report:
        print(args.report)


if __name__ == "__main__":
    main()
