"""
Generates the particles with certain initial position and velocities. 
"""

# Standard Python modules
from pyns.standard import *

# Sisters from this module
from pyns.lagrangian import *

import random

# =============================================================================
def initialiser(n, rho_p=[], d=[], verbose = False):
# -----------------------------------------------------------------------------  
    """
    Args:
        n:        The number of particles that one wishes to generate.
        rho_p:    Density of the particles. 
        r:        Radius of the particles.
        verbose: Logical variable setting, if particle's initial properties 
                  are to be printed.     
    Returns:
      Particles
    """

    pt = [Particles(random.uniform(0.3875, 0.4125), 0.09,
                    random.uniform(0.0001,0.0249), 
                    0, 0, 0, rho_p, d) for p in range(0, n)]
    
    if verbose:
        for p in range(0,n):
            print(pt[p].x, pt[p].y, pt[p].z, pt[p].u, pt[p].v, pt[p].w)
    
    # Writing initial conditions to a text file. To be inputted into 
    # Fluent.
    #
    # Format of the File:
    # ((x, y, z, u, v, w, d, temp, mass-flow))
    
    f = open("initial_particles.txt", "w")

    for p in range(0,n):
        f.write("((%12.5e %12.5e %12.5e %12.5e %12.5e %12.5e %12.5e %12.5e %12.5e)) \n" 
                % (pt[p].x, pt[p].y, pt[p].z, pt[p].u, pt[p].v, pt[p].w, d, 300, 1e-20))

    f.close()        
    
    return pt  # end of function
