"""
Matrix-vector product, including booundary values, for PyNS matrix format.

Notes:
  This function superseedes the plain "mat_vec".
"""

# Standard Python modules
from pyns.standard import *

# PyNS modules
from pyns.constants import W, E, S, N, B, T
from pyns.operators import cat_x, cat_y, cat_z

# =============================================================================
def mat_vec_bnd(a, phi):
# -----------------------------------------------------------------------------
    """
    Args:
      a: Matrix sent for multiplication in PyNS format (essentially that
         is storing a bundle of non-zero diagonals in compas directions)
      x: Object of the type "Unknown" to be solved.

    Returns:
      r: Result of the matrix-vector product, which is a vector stored
         in a three-dimensional array.
    """

    r = zeros(phi.val.shape)
    
    phi.exchange()

    r[:]  = a.P[:] * phi.val[:]
    
    r[:] -= a.W[:] * cat_x( (phi.bnd[W].val[ :1,:,:], 
                             phi.val       [:-1,:,:]) )
    r[:] -= a.E[:] * cat_x( (phi.val       [ 1:,:,:], 
                             phi.bnd[E].val[ :1,:,:]) )
    
    r[:] -= a.S[:] * cat_y( (phi.bnd[S].val[:, :1,:], 
                             phi.val       [:,:-1,:]) )
    r[:] -= a.N[:] * cat_y( (phi.val       [:, 1:,:], 
                             phi.bnd[N].val[:, :1,:]) )
    
    r[:] -= a.B[:] * cat_z( (phi.bnd[B].val[:,:, :1], 
                             phi.val       [:,:,:-1]) )
    r[:] -= a.T[:] * cat_z( (phi.val       [:,:, 1:], 
                             phi.bnd[T].val[:,:, :1]) )

    return r  # end of function
