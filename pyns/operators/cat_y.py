"""
An interface to the standard NumPy's function "concatenate" by second index.
"""

# Standard Python modules
from pyns.standard import concatenate

# PyNS modules
from pyns.constants import Y

# =============================================================================
def cat_y(tup):
# -----------------------------------------------------------------------------
    """
    Args:
      tup: Tuple containing whatever has to be sent to NumPy's function.

    Returns:
      Concatenad arrays or matrices.
    """

    return concatenate(tup, Y)  # end of function
