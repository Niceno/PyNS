"""
Compute the norm of a vector.
"""

from standard import *

from solvers.vec_vec import vec_vec

# =============================================================================
def norm(x):
# -----------------------------------------------------------------------------
  """
  Args:
    x: Vector.
    
  Returns:
    Norm of the vector. 
  """
  
  return sqrt( vec_vec(x, x) )  # end of function
