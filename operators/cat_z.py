"""
An interface to the standard NumPy's function "concatenate" by third index.
"""

# Standard Python modules
from standard import *

# PyNS modules
from scrins.constants.coordinates import Z

# =============================================================================


def cat_z(tup):
    # -----------------------------------------------------------------------------
    """
    Args:
      tup: Tuple containing whatever has to be sent to NumPy's function.

    Returns:
      Concatenad arrays or matrices.
    """

    return concatenate(tup, Z)  # end of function
