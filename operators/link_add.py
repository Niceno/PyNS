"""
Links the values of array in "x" direction.

Example:

  Arrays sent:   |---o---|---o---|---o---|---o---|---o---|---o---| 
                     0       1       2       3       4       5

  Returnig:      |---o---|---o---|---o---|---o---|---o---|---o---| 
                    0+5      1       2       3       4      0+5   

First and last value in the returning array will be the average of themselves.                   
"""

# =============================================================================
def link_add_x(a):
# -----------------------------------------------------------------------------
    """
    Args:
      a: matrices for linking.

    Returns:
      None, but input arguments are changed.
    """
 
    a[ :1,:,:] += a[-1:,:,:]
    a[-1:,:,:]  = a[ :1,:,:]

    return  # end of function

# =============================================================================
def link_add_y(a):
# -----------------------------------------------------------------------------
    """
    Args:
      a: matrices for linking.

    Returns:
      None, but input arguments are changed.
    """
 
    a[:, :1,:] += a[:,-1:,:]
    a[:,-1:,:]  = a[:, :1,:]

    return  # end of function

# =============================================================================
def link_add_z(a):
# -----------------------------------------------------------------------------
    """
    Args:
      a: matrices for linking.

    Returns:
      None, but input arguments are changed.
    """
 
    a[:,:, :1] += a[:,:,-1:]
    a[:,:,-1:]  = a[:,:, :1]

    return  # end of function
