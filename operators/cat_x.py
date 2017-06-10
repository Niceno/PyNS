"""
An interface to the standard NumPy's function "concatenate" by first index.
"""

# Standard Python modules
from pyns.standard import *

# PyNS modules
from pyns.constants.all import X

# =============================================================================
def cat_x(tup):
# -----------------------------------------------------------------------------
    """
    Args:
      tup: Tuple containing whatever has to be sent to NumPy's function.

    Returns:
      Concatenad arrays or matrices.
    """

    return concatenate(tup, X)  # end of function
