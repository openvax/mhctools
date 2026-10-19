"""Small runtime checks for predictor inputs, independent of optional backends."""

from collections.abc import Sequence


def require_type(value: object, expected_type: type, name: str) -> None:
    """Reject a wrong type with a short, argument-specific error."""
    if not isinstance(value, expected_type):
        raise TypeError(
            f"{name} must be {expected_type.__name__}, "
            f"got {type(value).__name__}")


def require_items(values: Sequence, expected_type: type, name: str) -> None:
    """Validate materialized inputs in one pass, retaining the failing index."""
    for index, value in enumerate(values):
        require_type(value, expected_type, f"{name}: element at index {index}")
