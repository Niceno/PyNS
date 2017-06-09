"""
Returns runnig difference of the array or matrix sent as parameter.

If array is sent for differencing, only one argument is needed, array itself.

If matrix is sent for differencing, two arguments are needed, matrix itself,
and the desired direction for differencing.

It is clear that the returning difference will have one element less in the
direction in which it was differenced:

Example:

  Array sent:    |-------|-------|-------|-------|-------| => six elements
                 0       1       2       3       4       5

  Returnig array:    o-------o-------o-------o-------o     => five elements
                     0       1       2       3       4

  Element 0 in the returning array will be the difference between the values
  in elements 1 and 0 in the sending array.

Note:
  Difference is NOT derivative.  To find a derivative, you should divide
  difference with a proper array with distances between elements!
"""

from scrins.constants.coordinates import X, Y, Z

# =============================================================================


def dif(*args):
    # -----------------------------------------------------------------------------
    """
    Args: Number of input arguments can be one or two, depending if one wants to
          run a difference on array or matrix.

      One argument provided (for arrays)
        args[0]: Array for differencing.

      Two arguments provided (for matrices)
        args[0]: Matrix for differencing.
        args[1]: Direction for differencing (X, Y or Z)

    Returns:
      Array or matrix with differenced values.

    Note:
      Actually, it would make perfect sense to derive separate functions for
      differencing in "x", "y" and "z" directions.  Less arguments passing, less
      arguments parsing, and fewer "if" statements.
    """

    # Only one argument is sent - perform
    # differencing of one-dimensional array
    if len((args)) == 1:
        a = args[0]
        return (a[1:] - a[:-1])

    # Two arguments are sent - perform differencing of a
    # three-dimensional array, with a given direction
    elif len((args)) == 2:
        d = args[0]  # direction
        a = args[1]  # array

        if d == X:
            return (a[1:, :, :] - a[:-1, :, :])
        elif d == Y:
            return (a[:, 1:, :] - a[:, :-1, :])
        elif d == Z:
            return (a[:, :, 1:] - a[:, :, :-1])

    # Some error message might
    # be printed if number of
    # arguments is wrong

    return  # end of function

def dif_x(a):
    # -----------------------------------------------------------------------------
    """
    Args:
      a: matrix for differencing.

    Returns:
      Matrix with differenced values.
    """

    return (a[1:, :, :] - a[:-1, :, :])  # end of function

def dif_y(a):
    # -----------------------------------------------------------------------------
    """
    Args:
      a: matrix for differencing.

    Returns:
      Matrix with differenced values.
    """

    return (a[:, 1:, :] - a[:, :-1, :])  # end of function

def dif_z(a):
    # -----------------------------------------------------------------------------
    """
    Args:
      a: matrix for differencing.

    Returns:
      Matrix with differenced values.
    """

    return (a[:, :, 1:] - a[:, :, :-1])  # end of function

