"""
An interface to the standard NumPy's function "concatenate".  With this, one
transforms a tuple with a number of arrays and/or matrices into a contiguous
memory space.
"""

# Standard Python modules
from pyns.standard import *

# PyNS modules
from pyns.constants.all import *

# =============================================================================
def cat(dir, tup):
# -----------------------------------------------------------------------------
    """
    Args:
      dir: Direction for catenation (X, Y or Z)
      tup: Tuple containing whatever has to be sent to NumPy's function.

    Returns:
      Concatenad arrays.

    Note:
      Actually, it would make perfect sense to derive separate functions for
      catenating in "x", "y" and "z" directions.  Less arguments passing and
      shorter syntax.
    """

    return concatenate(tup, dir)  # end of function
