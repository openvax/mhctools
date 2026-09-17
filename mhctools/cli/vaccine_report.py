"""Command-line entry point for route-aware vaccine processing reports."""

import argparse
from pathlib import Path

from ..vaccine_report import VaccineReportInput, generate_vaccine_report
from .errors import CLI_ERROR_TYPES


def main(argv=None):
    parser = argparse.ArgumentParser(
        prog="mhctools vaccine-report",
        description=(
            "Generate a timestamped, route-aware peptide cleavage and MHC-window "
            "report from structured JSON input."
        ),
    )
    parser.add_argument("--input", required=True, help="Vaccine report schema JSON")
    parser.add_argument(
        "--output-dir", required=True,
        help="Base directory; a date/time-stamped child directory is created",
    )
    parser.add_argument(
        "--max-mhc-windows", type=int, default=10,
        help="Maximum displayed windows per MHC class and construct (default: 10)",
    )
    args = parser.parse_args(argv)
    try:
        report = VaccineReportInput.from_json(args.input)
        output = generate_vaccine_report(
            report, Path(args.output_dir), maximum_mhc_windows=args.max_mhc_windows
        )
    except CLI_ERROR_TYPES as error:
        parser.error(str(error))
    print(output)


if __name__ == "__main__":
    main()
