"""
Calculates the closest cell centre to the particles position. 
It is important when checking whether the particle is inside the obstacle or 
not.
"""

# =============================================================================
def closest_cell(xc, xp):
# -----------------------------------------------------------------------------  
    """
    Args:
      xc: Cell coordinates in the "x", "y" or "z" direction.
      xp: Particle location in the "x", "y" or "z" direction.
        
    Returns:
      The index of the cell which the particle is closest to.
    """
    idxc = (abs(xc-xp)).argmin()
        
    return idxc  # end of function