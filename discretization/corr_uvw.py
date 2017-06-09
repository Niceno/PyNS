"""
Projects velocity into a divergence-free velocity field with the pressure
correction gradient (provided pressure correction is accuratelly solved).
"""

# Standard Python modules
from standard import *

# ScriNS modules
from scrins.constants.coordinates import X, Y, Z
from scrins.constants.compass import W, E, S, N, B, T, C
from scrins.discretization.obst_zero_val import obst_zero_val
from scrins.operators.avg import avg, avg_x, avg_y, avg_z
#from scrins.operators.avg_x import avg_x
#from scrins.operators.avg_y import avg_y
#from scrins.operators.avg_z import avg_z
from scrins.operators.cat import cat, cat_x, cat_y, cat_z
#from scrins.operators.cat_x import cat_x
#from scrins.operators.cat_y import cat_y
#from scrins.operators.cat_z import cat_z
from scrins.operators.dif import dif, dif_x, dif_y, dif_z
#from scrins.operators.dif_x import dif_x
#from scrins.operators.dif_y import dif_y
#from scrins.operators.dif_z import dif_z

# =============================================================================


def corr_uvw(uvw, p, rho, dt, dxyz, obst):
    # -----------------------------------------------------------------------------
    """
    Args:
      uvw:  Tuple with three staggered or centered velocity components
            (each component is created with "create_unknown" function.
      p:    Unknown holding the pressure correction.
      rho:  Three-dimensional matrix holding density for all cells.
      dt:   Time step
      dxyz: Tuple holding cell dimensions in "x", "y" and "z" directions.
            Each cell dimension is a three-dimensional array.
      obst: Obstacle, three-dimensional matrix with zeros and ones.  It is
            zero in fluid, one in solid.

    Returns:
      none, but input argument uvw is modified.
    """

    # Unpack received tuples
    dx, dy, dz = dxyz

    # Compute pressure correction gradients
    p_x = dif_x(p.val) / avg_x(dx)
    p_y = dif_y(p.val) / avg_y(dy)
    p_z = dif_z(p.val) / avg_z(dz)

    # Set to zero in obst
    if obst.any() != 0:
        p_x = obst_zero_val(X, p_x, obst)
        p_y = obst_zero_val(Y, p_y, obst)
        p_z = obst_zero_val(Z, p_z, obst)

    # Pad with boundary values by expanding from interior
    # (This is done only for collocated formulation)
    if uvw[X].pos == C:
        p_x = avg_x(cat_x((p_x[:1, :, :], p_x, p_x[-1:, :, :])))
        p_y = avg_y(cat_y((p_y[:, :1, :], p_y, p_y[:, -1:, :])))
        p_z = avg_z(cat_z((p_z[:, :, :1], p_z, p_z[:, :, -1:])))

    # Correct the velocities
    uvw[X].val[:] = uvw[X].val[:] - dt / avg(uvw[X].pos, rho) * p_x
    uvw[Y].val[:] = uvw[Y].val[:] - dt / avg(uvw[Y].pos, rho) * p_y
    uvw[Z].val[:] = uvw[Z].val[:] - dt / avg(uvw[Z].pos, rho) * p_z

    return  # end of function
