import numpy as np


def valid_affinity(x):
    """
    Check that an IC50 affinity value is valid.

    Parameters
    ----------
    x : float

    Returns
    -------
    bool
    """
    if x is None:
        return False
    if np.isnan(x) or np.isinf(x):
        return False
    return x >= 0


def valid_percentile_rank(x):
    """
    Check whether a percentile rank is valid.

    Parameters
    ----------
    x : float

    Returns
    -------
    bool
    """
    if x is None:
        return False
    if np.isnan(x) or np.isinf(x):
        return False
    return 0 <= x <= 100
