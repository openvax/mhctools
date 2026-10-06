"""Shared display labels; native prediction values remain unchanged."""

import math


def dpp4_loss_label(score: float) -> str:
    """Format native log2 control/treated signal as a simple assay-loss label."""
    if not math.isfinite(score):
        raise ValueError("DPP4 loss display requires a finite native score")
    if score < 0:
        return "No predicted loss"
    percent = -100 * math.expm1(-score * math.log(2))
    if 0 < percent < 1:
        return "<1% predicted loss"
    return f"{percent:.0f}% predicted loss"
