"""
Vector-vector product for PyNS matrix format
"""

from standard import *

# =============================================================================
def vec_vec(x, y):
# -----------------------------------------------------------------------------
  """
  Args:
    x: Three-dimensional matrix with a vector for multiplication.
    y: Three-dimensional matrix with a vector for multiplication.
    
  Returns:
    Result of the matrix-vector product as a three-dimensional matrix. 
    
  Note:
    Try to find a better way to summ the elements of a matrix than 
    sum(sum(sum()))
  """
  
  return sum( sum( sum( multiply(x, y) ) ) )  # end of function
