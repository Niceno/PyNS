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

# PyNS modules
from pyns.constants.coordinates import X, Y, Z

# =============================================================================
def dif(*args):
# -----------------------------------------------------------------------------
    """
    Args: Number of input arguments can be one or two, depending if one
          wants to run a difference on array or matrix.

      One argument provided (for arrays)
        args[0]: Array for differencing.

      Two arguments provided (for matrices)
        args[0]: Matrix for differencing.
        args[1]: Direction for differencing (X, Y or Z)

    Returns:
      Array or matrix with differenced values.

    Note:
      Separate functions for differencing in "x", "y" and "z" directions
      are also defined.  They require fewer arguments being passed and
      shorter syntax.  This function, however, is still usefull in
      instances in which the differencing direction is not predefined.
    """

    # Only one argument is sent - perform
    # differencing of one-dimensional array
    if len((args)) == 1:
        a = args[0]
        return (a[1:] - a[:-1])

    # Two arguments are sent - perform differencing of a
    # three-dimensional array, with a given direction
    elif len((args)) == 2:
        dir = args[0]  # direction
        a   = args[1]  # array

        if dir == X:
            return (a[1:,:,:] - a[:-1,:,:])
        elif dir == Y:
            return (a[:,1:,:] - a[:,:-1,:])
        elif dir == Z:
            return (a[:,:,1:] - a[:,:,:-1])

    # Some error message might
    # be printed if number of
    # arguments is wrong

    return  # end of function
