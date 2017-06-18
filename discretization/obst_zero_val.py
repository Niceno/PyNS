"""
Set values of an unknown to zero inside the obstacle.
"""

# Standard Python modules
from pyns.standard import *

# PyNS modules
from pyns.constants import *
from pyns.operators import *

# =============================================================================
def obst_zero_val(pos, val, obst):
# -----------------------------------------------------------------------------
    """
    Args:
      pos:  Position of the variable (C, X, Y or Z).
      val:  Value to be set in the obstacle.
      obst: Obstacle, three-dimensional array with zeros and ones.
            It is zero in fluid, one in solid.
    """

    if pos == C:
        val = val * lnot(obst)

    elif pos == X:
        obst_x = mx(obst[:-1,:,:], obst[1:,:,:])
        val = val * lnot(obst_x)

    elif pos == Y:
        obst_y = mx(obst[:,:-1,:], obst[:,1:,:])
        val = val * lnot(obst_y)

    elif pos == Z:
        obst_z = mx(obst[:,:,:-1], obst[:,:,1:])
        val = val * lnot(obst_z)

    return val  # end of function
