"""JSON interface to explicitly parameterized SLP/mRNA pMHC scenarios."""

import argparse
import json
from pathlib import Path

from ..serum_contributions import EnzymeCutRate
from ..vaccine_trajectory import (
    VaccineCleavageRates, VaccineTrajectoryInput, simulate_vaccine_trajectory,
)
from .errors import CLI_ERROR_TYPES


def main(argv=None):
    parser = argparse.ArgumentParser(
        prog="mhctools vaccine-trajectory",
        description="Conditional local SLP/secreted-mRNA trajectories to APC pMHC; requires explicit kinetics.")
    parser.add_argument("--input", required=True, help="Trajectory schema JSON, including observation times")
    parser.add_argument("--output", required=True, help="Output JSON with curves and full mechanism/assumption audit")
    args = parser.parse_args(argv)
    try:
        value = json.loads(Path(args.input).read_text())
        model = VaccineTrajectoryInput.from_dict(value)
        sequence = model.validate().sequence
        table = {}
        for item in value.get("cleavage_assessments", []):
            start, end = item["start"], item["end"]
            if (type(start) is not int or type(end) is not int or
                    not 0 <= start < end <= len(sequence) or
                    item["sequence"] != sequence[start:end]):
                raise ValueError("Cleavage assessment must match its original-parent interval")
            key = item["compartment"], start, end
            if key in table:
                raise ValueError("Duplicate fragment cleavage assessment")
            table[key] = VaccineCleavageRates(
                tuple(EnzymeCutRate(**channel) for channel in item["channels"]), item["basis"], item["source"])
        output = simulate_vaccine_trajectory(
            model, value["times_hours"], lambda compartment, fragment: table.get((compartment, fragment.start, fragment.end)))
        # A path's parent must exist; do not silently overwrite a prior analysis.
        with Path(args.output).open("x") as stream:
            json.dump(output, stream, indent=2, allow_nan=False)
            stream.write("\n")
    except CLI_ERROR_TYPES as error:
        parser.error(str(error))
    print(args.output)


if __name__ == "__main__":
    main()
