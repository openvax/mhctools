"""
Capture a tool's own allele list as a test fixture.

    python scripts/capture_allele_lists.py OUTPUT EXECUTABLE --list-flag=FLAG

for example::

    python scripts/capture_allele_lists.py \\
        tests/data/allele_lists/netmhc4_listmhc.txt netMHC-4.0 --list-flag=-listMHC

Run the tool, drop comment/header lines (which can contain local paths and
user names), and write a sample: every unusual name (non-human species,
expression variants, old-style spellings) plus a deterministic one-in-ten of
the routine human names. Pass ``--all`` to keep every name.
"""

import argparse
import re
import subprocess
import zlib

ROUTINE = re.compile(r"^(HLA-.*|DRB\d_\d+)\d$")
SAMPLE_ONE_IN = 10


def allele_names(text):
    names = []
    for line in text.splitlines():
        line = line.strip()
        if not line or line.startswith("#") or line.startswith("PLATFORM="):
            continue
        names.append(line)
    return names


def sample(names):
    kept = []
    for name in names:
        routine = ROUTINE.match(name) and not re.search(r"\d[A-Z]$", name)
        if not routine or zlib.crc32(name.encode()) % SAMPLE_ONE_IN == 0:
            kept.append(name)
    return kept


def main():
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("output")
    parser.add_argument("executable")
    parser.add_argument("--list-flag", required=True)
    parser.add_argument("--all", action="store_true")
    args = parser.parse_args()
    result = subprocess.run(
        [args.executable, args.list_flag], capture_output=True, text=True,
        check=True)
    names = allele_names(result.stdout)
    if not args.all:
        names = sample(names)
    with open(args.output, "w") as f:
        f.write("\n".join(names) + "\n")
    print("wrote %d names to %s" % (len(names), args.output))


if __name__ == "__main__":
    main()
