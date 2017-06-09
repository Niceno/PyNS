"""
Matrix-vector product for PyNS matrix format.
"""

from standard import *

# =============================================================================
def mat_vec(a, x):
# -----------------------------------------------------------------------------
  """
  Args:
    a: Matrix sent for multiplication in PyNS format (essentially that 
       is storing a bundle of non-zero diagonals in compas directions)
    x: Three-dimensional matrix with a vector for multiplication.
    
  Returns:
    r: Result of the matrix-vector product as a three-dimensional matrix
  """
 
  r = zeros(x.shape)
   
  r[  :,  :,  :]  = a.P[  :,  :,  :] * x[  :,  :,  :]
  r[ 1:,  :,  :] -= a.W[ 1:,  :,  :] * x[:-1,  :,  :]
  r[:-1,  :,  :] -= a.E[:-1,  :,  :] * x[ 1:,  :,  :]
  r[  :, 1:,  :] -= a.S[  :, 1:,  :] * x[  :,:-1,  :]
  r[  :,:-1,  :] -= a.N[  :,:-1,  :] * x[  :, 1:,  :]
  r[  :,  :, 1:] -= a.B[  :,  :, 1:] * x[  :,  :,:-1]
  r[  :,  :,:-1] -= a.T[  :,  :,:-1] * x[  :,  :, 1:]
  
  return r  # end of function
