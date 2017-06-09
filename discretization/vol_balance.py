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
from scrins.operators.cat import cat, cat_x, cat_y, cat_z
#from scrins.operators.cat_x import cat_x
#from scrins.operators.cat_y import cat_y
#from scrins.operators.cat_z import cat_z
from scrins.operators.dif import dif, dif_x, dif_y, dif_z
#from scrins.operators.dif_x import dif_x
#from scrins.operators.dif_y import dif_y
#from scrins.operators.dif_z import dif_z

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
    src = - dif_x(cat_x((uf.bnd[W].val, uf.val, uf.bnd[E].val))) * dy * dz  \
          - dif_y(cat_y((vf.bnd[S].val, vf.val, vf.bnd[N].val))) * dx * dz  \
          - dif_z(cat_z((wf.bnd[B].val, wf.val, wf.bnd[T].val))) * dx * dy

    # Zero it inside obstacles, if obstacle is sent as parameter
    if obst.any() != 0:
        src = obst_zero_val(C, src, obst)

    return src  # end of function
