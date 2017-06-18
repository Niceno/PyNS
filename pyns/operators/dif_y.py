"""
Returns runnig difference in "y" direction of the matrix sent as parameter.

Returning matrix will have one element less in "y" direction.

Note:
  Difference is NOT derivative.  To find a derivative, you should divide
  difference with a proper array with distances between elements!
"""

# =============================================================================
def dif_y(a):
# -----------------------------------------------------------------------------
    """
    Args:
      a: matrix for differencing.

    Returns:
      Matrix with differenced values.
    """

    return (a[:,1:,:] - a[:,:-1,:])  # end of function
