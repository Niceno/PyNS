"""
Links the values of array in "x", "y" or "z" direction.

Example:

  Arrays sent:   |---o---|---o---|---o---|---o---|---o---|---o---| 
                     0       1       2       3       4       5

  Returnig:      |---o---|---o---|---o---|---o---|---o---|---o---| 
                  (0+5)/2    1       2       3       4    (0+5)/2  

First and last value in the returning array will be the average of themselves.
"""

# =============================================================================
def link_avg_x(a):
# -----------------------------------------------------------------------------
    """
    Args:
      a: matrix for linking.

    Returns:
      None, but input arguments are changed.
    """
 
    a[ :1,:,:] = (a[ :1,:,:] + a[-1:,:,:]) * 0.5
    a[-1:,:,:] =  a[ :1,:,:]
      
    return  # end of function

# =============================================================================
def link_avg_y(a):
# -----------------------------------------------------------------------------
    """
    Args:
      a: matrix for linking.

    Returns:
      None, but input arguments are changed.
    """
 
    a[:, :1,:] = (a[:, :1,:] + a[:,-1:,:]) * 0.5
    a[:,-1:,:] =  a[:, :1,:]
      
    return  # end of function

# =============================================================================
def link_avg_z(a):
# -----------------------------------------------------------------------------
    """
    Args:
      a: matrix for linking.

    Returns:
      None, but input arguments are changed.
    """
 
    a[:,:, :1] = (a[:,:, :1] + a[:,:,-1:]) * 0.5
    a[:,:,-1:] =  a[:,:, :1]
      
    return  # end of function
