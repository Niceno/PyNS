"""
Matrix-vector product for PyNS matrix format.

Note:
  This function doesn't take care of periodic boundary conditions, and is
  therefore superseeded by the "mat_vec_bnd".  Actually, it is a candidate
  for deletion.
"""

# Standard Python modules
from pyns.standard import *

# =============================================================================
def mat_vec(a, x):
# -----------------------------------------------------------------------------
    """
    Args:
      a: Matrix sent for multiplication in PyNS format (essentially that
         is storing a bundle of non-zero diagonals in compas directions)
      x: Three-dimensional array holding a vector for multiplication.

    Returns:
      r: Result of the matrix-vector product, which is a vector stored
         in a three-dimensional array.
    """

    r = zeros(x.shape)

    r[  :,  :,  :]  = a.P[  :,  :,  :] * x[  :,  :,  :]
    r[ 1:,  :,  :] -= a.W[ 1:,  :,  :] * x[:-1,  :,  :]
    r[:-1,  :,  :] -= a.E[:-1,  :,  :] * x[ 1:,  :,  :]
    r[  :, 1:,  :] -= a.S[  :, 1:,  :] * x[  :,:-1,  :]
    r[  :,:-1,  :] -= a.N[  :,:-1,  :] * x[  :, 1:,  :]
    r[  :,  :, 1:] -= a.B[  :,  :, 1:] * x[  :,  :,:-1]
    r[  :,  :,:-1] -= a.T[  :,  :,:-1] * x[  :,  :, 1:]

    return r  # end of function
