"""
Set values of an unknown to zero inside the obstacle.
"""

# Standard Python modules
from pyns.standard import *

# PyNS modules
from pyns.constants.all import *
from pyns.operators.all import *

# =============================================================================
def obst_zero_val(d, val, obst):
# -----------------------------------------------------------------------------
  """
  Args:
    d:    Position of the variable (C, X, Y or Z).
    val:  Value to be set in the obstacle.
    obst: Obstacle, three-dimensional matrix with zeros and ones.  It is
          zero in fluid, one in solid.
  """

  if d == C:  
    val = val * lnot(obst)
    
  elif d == X:
    obst_x = mx(obst[:-1,:,:], obst[1:,:,:])
    val = val * lnot(obst_x)
    
  elif d == Y: 
    obst_y = mx(obst[:,:-1,:], obst[:,1:,:])
    val = val * lnot(obst_y)
    
  elif d == Z: 
    obst_z = mx(obst[:,:,:-1], obst[:,:,1:])
    val = val * lnot(obst_z)

  return val  # end of function
