"""
Compute the norm of a vector.
"""

from pyns.standard import *

from pyns.solvers.vec_vec import vec_vec

# =============================================================================
def norm(x):
# -----------------------------------------------------------------------------
    """
    Args:
      x: Vector.

    Returns:
      Norm of the vector.
    """

    return sqrt( vec_vec(x, x) )  # end of function
