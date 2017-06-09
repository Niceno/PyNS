"""
An interface to the standard NumPy's function "concatenate".  With this, one
transforms a tuple with a number of arrays and/or matrices into a contiguous
memory space.
"""

from numpy import concatenate

from scrins.constants.coordinates import X, Y, Z
# =============================================================================
def cat(d, tup):
    """
    Args:
      d:   Direction for catenation (X, Y or Z)
      tup: Tuple containing whatever has to be sent to NumPy's function.

    Returns:
      Concatenad arrays.

    Note:
      Actually, it would make perfect sense to derive separate functions for
      catenating in "x", "y" and "z" directions.  Less arguments passing and
      shorter syntax.
    """

    return concatenate(tup, d)  # end of function

def cat_x(tup):
    # -----------------------------------------------------------------------------
    """
    Args:
      tup: Tuple containing whatever has to be sent to NumPy's function.

    Returns:
      Concatenad arrays or matrices.
    """

    return concatenate(tup, X)  # end of function

def cat_y(tup):
    # -----------------------------------------------------------------------------
    """
    Args:
      tup: Tuple containing whatever has to be sent to NumPy's function.

    Returns:
      Concatenad arrays or matrices.
    """

    return concatenate(tup, Y)  # end of function


def cat_z(tup):
    # -----------------------------------------------------------------------------
    """
    Args:
      tup: Tuple containing whatever has to be sent to NumPy's function.

    Returns:
      Concatenad arrays or matrices.
    """

    return concatenate(tup, Z)  # end of function

