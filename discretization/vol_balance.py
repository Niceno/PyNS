"""
Computes volume balance inside the domain, which is, in essence, the right
hand side in the pressure-Poisson equation.
"""

# Standard Python modules
from standard import *

# ScriNS modules
from scrins.constants.coordinates import X, Y, Z
from scrins.constants.compass import W, E, S, N, B, T, C
from scrins.discretization.obst_zero_val import obst_zero_val
from scrins.operators.cat import cat
from scrins.operators.dif import dif

# =============================================================================


def vol_balance(uvwf, dxyz, obst):
    """
    Args:
      uvwf: Tuple with three staggered velocity components (where each
            component is created with "create_unknown" function.
      dxyz: Tuple holding cell dimensions in "x", "y" and "z" directions.
            Each cell dimension is a three-dimensional array.
      obst: Obstacle, three-dimensional matrix with zeros and ones.  It is
            zero in fluid, one in solid.

    Note:
      "obst" is an optional parameter.  If it is not sent, the source won't be
      set to zero inside the obstacle.  That is important for calculation of
      pressure, see function "pyns.calc_p" as well.
    """
    # Unpack tuples
    dx, dy, dz = dxyz
    uf, vf, wf = uvwf

    # Compute it throughout the domain
    src = - dif(X, cat(X, (uf.bnd[W].val, uf.val, uf.bnd[E].val))) * dy * dz  \
          - dif(Y, cat(Y, (vf.bnd[S].val, vf.val, vf.bnd[N].val))) * dx * dz  \
          - dif(Z, cat(Z, (wf.bnd[B].val, wf.val, wf.bnd[T].val))) * dx * dy

    # Zero it inside obstacles, if obstacle is sent as parameter
    if obst.any() != 0:
        src = obst_zero_val(C, src, obst)

    return src  # end of function
