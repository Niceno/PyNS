"""
Generates the particles with certain initial position and velocities. 
"""

# Standard Python modules
from pyns.standard import *

# Sisters from this module
from pyns.lagrangian import *

import random

# =============================================================================
def initialiser(n, verbose = False):
# -----------------------------------------------------------------------------  
    """
    Args:
        n: ..... The number particle's that one wishes to generate. 
        verbose: Logical variable setting if particle's initial properties 
                 are printed.     
    Returns:
      pt: Particles.
    """
    pt = [Particles(random.uniform(0.4,  0.6), 
                    random.uniform(0.1,  0.9), 
                    random.uniform(0.01, 0.24), 
                    0, 0, 0) for p in range(0, n)]
    
    if verbose is True:
        for p in range(0,x):
            print(pt[p].x, pt[p].y, pt[p].z, pt[p].u, pt[p].v, pt[p].w)
        
    return pt  # end of function
