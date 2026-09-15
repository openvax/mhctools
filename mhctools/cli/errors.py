"""Shared formatting for expected command-line failures."""

from subprocess import SubprocessError


CLI_ERROR_TYPES = (
    ImportError,
    KeyError,
    OSError,
    RuntimeError,
    SubprocessError,
    SystemError,
    TypeError,
    ValueError,
)


def cli_error_message(error):
    """Return a user-facing message without Python API-only wording."""
    if isinstance(error, KeyError) and error.args:
        message = str(error.args[0])
    else:
        message = str(error)
    return message.replace(
        "pass overwrite=True", "pass --overwrite").replace(
        " (or accept_license=True)", "")


def predictor_error_message(error, predictor_specs):
    """Add artifact inventory guidance for a missing optional dependency."""
    if isinstance(error, ModuleNotFoundError):
        from ..artifacts import artifact_status

        for raw_spec in predictor_specs:
            name = raw_spec.split(":", 1)[0].strip().lower()
            try:
                status = artifact_status(name)
            except ValueError:
                continue
            if status.status != "ready" and status.detail:
                return status.detail
    return cli_error_message(error)

