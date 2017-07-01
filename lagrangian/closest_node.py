"""
This script calculates the closest node to the particles position. 
It is important when interpolating the particles velocity. 
"""

# =============================================================================
def closest_node(xn, xp):
# -----------------------------------------------------------------------------  
    """
    Args:
      xn: Nodes coordinates in the "x", "y" or "z" direction.
      xp: Particle location in the "x", "y" or "z" direction.
        
    Returns:
        The index of the cell which the particle is closest to, in the x, y or 
        z direction.
    """
    
    # Returns the index of the closest node.
    idx = (abs(xn-xp)).argmin()
    
    # For the interpolater we need to find the index of the node on the 
    # other side to the closest node.
    #
    # o = particle
    # | = nodes
    #
    #
    # ---|------|------|------|------|----o-|------|------|------|------|
    #                               idx1   idx

    
    if xn[idx] > xp:
        idx1 = idx - 1 
    else:
        idx1 = idx + 1 
        
    return idx, idx1  # end of function

