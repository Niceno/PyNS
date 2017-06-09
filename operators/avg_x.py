"""
Returns runnig average in "x" direction of the matrix sent as parameter.

Returning matrix will have one element less in "x" direction.

Example:

  Array sent:    |-------|-------|-------|-------|-------| => six elements
                 0       1       2       3       4       5

  Returnig array:    o-------o-------o-------o-------o     => five elements
                     0       1       2       3       4

  Element 0 in the returning array will be the average between the values
  in elements 1 and 0 in the sending array.
"""

# =============================================================================


def avg_x(a):
    # -----------------------------------------------------------------------------
    """
    Args:
      a: matrix for averaging.

    Returns:
      Matrix with averaged values.
    """

    return (a[1:, :, :] + a[:-1, :, :]) * 0.5  # end of function
