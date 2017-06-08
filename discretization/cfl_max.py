"""
Computes maximum Courant-Friedrich-Levy (CFD) number for the given velocity.
"""

# ScriNS modules
from scrins.constants.coordinates import X, Y, Z
from scrins.constants.compass import W, E, S, N, B, T, C
from scrins.operators.avg import avg
from scrins.operators.avg_x import avg_x 
from scrins.operators.avg_y import avg_y
from scrins.operators.avg_z import avg_z

# =============================================================================


def cfl_max(uvw, dt, dxyz):
    # -----------------------------------------------------------------------------
    """
    Args:
      uvw:  Tuple with three staggered or centered velocity components
            (each component is created with "create_unknown" function.
      dt:   Time step
      dxyz: Tuple holding cell dimensions in "x", "y" and "z" directions.
            Each cell dimension is a three-dimensional array.

    Returns:
      cfl: Floating point number holding the maximum value of CFL.

    Note:
      It could be written in a more compact way, without unpacking the tuples,
      but that would only lead to poorer readability of the code.
    """

    # Unpack received tuples
    u, v, w = uvw
    dx, dy, dz = dxyz

    # Take velocity's position
    d = u.pos

    # Mesh is cell-centered
    if d == C:
        cfl = dt * max(abs(u.val / dx).max(),
                       abs(v.val / dy).max(),
                       abs(w.val / dz).max())

    # Mesh is staggered
    else:
        cfl = dt * max(abs(u.val / avg_x(dx)).max(),
                       abs(v.val / avg_y(dy)).max(),
                       abs(w.val / avg_z(dz)).max())

    return cfl
