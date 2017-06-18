"""
Vector-dot product of two vectors stored as three-dimensional arrays.
"""

# Standard Python modules
from pyns.standard import *

# =============================================================================
def vec_vec(x, y):
# -----------------------------------------------------------------------------
    """
    Args:
      x: Three-dimensional array holding vector for multiplication.
      y: Three-dimensional array holding vector for multiplication.

    Returns:
      Result of the vector-dot product.

    Note:
      Try to find a better way to summ the elements of a matrix than
      sum(sum(sum()))
    """

    return sum( sum( sum( multiply(x, y) ) ) )  # end of function
