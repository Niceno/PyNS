"""
Computes maximum Courant-Friedrich-Levy (CFD) number for the given velocity.
"""

# PyNS modules
from pyns.constants import *
from pyns.operators import *
from pyns.display   import write

# =============================================================================
def cfl_max(uvw, dt, dxyz,
            verbose = True):
# -----------------------------------------------------------------------------
    """
    Args:
      uvw:  Tuple with three staggered or centered velocity components.
            (Each component is an object of type "Unknown".
      dt:   Time step
      dxyz: Tuple holding cell dimensions in "x", "y" and "z" directions.
            Each cell dimension is a three-dimensional array.

    Returns:
      cfl: Floating point number holding the maximum value of CFL.

    Note:
      It could be written in a more compact way, without unpacking the
      tuples, but that would only lead to poorer readability of the code.
    """

    # Unpack received tuples
    u,  v,  w  = uvw
    dx, dy, dz = dxyz

    # Take velocity's position
    d = u.pos

    # Mesh is cell-centered
    if d == C:
        cfl = dt * max( abs(u.val/dx).max(),   \
                        abs(v.val/dy).max(),   \
                        abs(w.val/dz).max() )

    # Mesh is staggered
    else:
        cfl = dt * max( abs(u.val/avg_x(dx)).max(),   \
                        abs(v.val/avg_y(dy)).max(),   \
                        abs(w.val/avg_z(dz)).max() )

    if verbose is True:
        write.at(__name__)
        print("  Maximum CFL number: %12.5e" % cfl)

    return cfl
