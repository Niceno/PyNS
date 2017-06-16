"""
An interface to the standard NumPy's function "concatenate".  With this, one
transforms a tuple with a number of arrays and/or matrices into a contiguous
memory space.
"""

# Standard Python modules
from pyns.standard import concatenate

# PyNS modules
from pyns.constants import X, Y, Z

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
      Separate functions for catenating in "x", "y" and "z" directions
      are also defined.  They require fewer arguments being passed and
      shorter syntax.  This function, however, is still usefull in
      instances in which the catenating direction is not predefined.
    """

    return concatenate(tup, dir)  # end of function

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
