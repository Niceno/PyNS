"""
Computes volume balance inside the domain, which is, in essence, the right
hand side in the pressure-Poisson equation.
"""

# Standard Python modules
from pyns.standard import *

# PyNS modules
from pyns.constants.compass import *
from pyns.operators         import *

from pyns.discretization.obst_zero_val  import obst_zero_val

# =============================================================================
def vol_balance(uvwf, dxyz, 
                obstacle = None):
# -----------------------------------------------------------------------------
    """
    Args:
      uvwf: ... Tuple with three staggered velocity components (where each
                component is object of type "Unknown").
      dxyz: ... Tuple holding cell dimensions in "x", "y" and "z" directions.
                Each cell dimension is a three-dimensional array.
      obstacle: Obstacle, three-dimensional array with zeros and ones.
                It is zero in fluid, one in solid.

    Note:
      "obstacle" is an optional parameter.  If it is not sent, the source
      won't be set to zero inside the obstacle.  That is important for
      calculation of pressure, see function "calc_p" as well.
    """

    # Unpack tuples
    dx, dy, dz = dxyz
    uf, vf, wf = uvwf
    
    # Refresh buffers here
    uf.exchange()
    vf.exchange()
    wf.exchange()

    # Compute it throughout the domain
    src = - dif_x(cat_x((uf.bnd[W].val, uf.val, uf.bnd[E].val)))*dy*dz  \
          - dif_y(cat_y((vf.bnd[S].val, vf.val, vf.bnd[N].val)))*dx*dz  \
          - dif_z(cat_z((wf.bnd[B].val, wf.val, wf.bnd[T].val)))*dx*dy

    # Zero it inside obstacles, if obstacle is sent as parameter
    if obstacle is not None:
        src = obst_zero_val(C, src, obstacle)

    return src  # end of function
