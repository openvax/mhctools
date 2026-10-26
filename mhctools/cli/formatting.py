# Licensed under the Apache License, Version 2.0 (the "License");
# https://www.apache.org/licenses/LICENSE-2.0

"""Plain-text formatting shared by CLI catalogs."""


def format_table(headers, rows):
    """Align string columns, retaining headers when there are no rows."""
    table = [headers, *rows]
    widths = [max(map(len, column)) for column in zip(*table)]
    return "\n".join(
        "  ".join(value.ljust(width) for value, width in zip(row, widths)).rstrip()
        for row in table
    )
