# Licensed under the Apache License, Version 2.0 (the "License");
# https://www.apache.org/licenses/LICENSE-2.0

"""Aligned CLI catalogs, with optional styling at the output boundary."""

from .help import color_enabled, style_status, style_table_header


def format_table(headers, rows, *, stream=None):
    """Align columns; style only an explicitly supplied capable terminal."""
    table = [headers, *rows]
    widths = [max(map(len, column)) for column in zip(*table)]
    colored = stream is not None and color_enabled(stream)
    lines = []
    for index, row in enumerate(table):
        cells = []
        for header, value, width in zip(headers, row, widths):
            cell = value.ljust(width)
            if colored and index:
                if value in ("-", "not checked", "not applicable"):
                    cell = style_status(cell, level="dim", use_color=True)
                elif header == "STATUS" and value in ("missing", "unavailable", "failed", "error"):
                    cell = style_status(cell, level="warning", use_color=True)
            cells.append(cell)
        line = "  ".join(cells).rstrip()
        if colored and index == 0:
            line = style_table_header(line, use_color=True)
        lines.append(line)
    return "\n".join(lines)
