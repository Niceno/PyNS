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
        The index of the node which the particle is closest to and the nodes 
        either side of it, in the x, y or z direction.
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
    #                               idx1   idx    idx2
    #
    #
    #
    # ---|------|------|------|------|-o----|------|------|------|------|
    #                       idx2    idx    idx1 
    #
    #
    # ---|------|------|------|------|-o----|------|------|------|----o-|
    #                                                           idx1   idx
    # idx2 = None in last case. 

    
    if xn[idx] > xp:
        idx1 = idx - 1
        idx2 = idx + 1
    else:
        idx1 = idx + 1
        idx2 = idx - 1
    
    # Have to deal with boundaries 
    if idx2 > xn.argmax() or idx2 < 0:
        
        return idx1, idx , None
    
    else:
        
       return idx1, idx, idx2  # end of function

