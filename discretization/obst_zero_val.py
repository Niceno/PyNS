"""
Set values of an unknown to zero inside the obstacle.
"""

# Standard Python modules
from standard import *

# ScriNS modules
from constants.all import *
from operators.all import *

# =============================================================================
def obst_zero_val(d, val, obst):
# -----------------------------------------------------------------------------
  """
  Args:
    d:    position of the variable (C, X, Y or Z)
    val:  value to be set in the obstacle
    obst: obstacle
  """

  if d == C:  
    val = val * lnot(obst)
    
  elif d==X:
    obst_x = mx(obst[:-1,:,:], obst[1:,:,:])
    val = val * lnot(obst_x)
    
  elif d==Y: 
    obst_y = mx(obst[:,:-1,:], obst[:,1:,:])
    val = val * lnot(obst_y)
    
  elif d==Z: 
    obst_z = mx(obst[:,:,:-1], obst[:,:,1:])
    val = val * lnot(obst_z)

  return val  # end of function