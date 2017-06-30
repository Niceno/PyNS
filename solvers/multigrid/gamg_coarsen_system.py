"""
Base for the geometrically agglomerated multigrid method.
"""

# Standard Python modules
from pyns.standard import *

# PyNS modules
from pyns.constants      import *
from pyns.operators      import *
from pyns.display        import write
from pyns.discretization import Unknown

# Sisters from this module
from pyns.solvers.Matrix import Matrix

# =============================================================================
def gamg_coarsen_system(a, phi, b,
                        verbatim = True):
# -----------------------------------------------------------------------------
    """
    Args:

    Returns:
 
    """

    # ---------------------------------------------------------
    # Copy Matrix shape into array for easier (shorter) syntax    
    # ---------------------------------------------------------
    shape = array(a.C.shape)    
        
    # -------------------------------------------------------
    # Check if any coarsening is possible
    #
    # Essentially, you do it by ensuring that all dimensions 
    # are greater than four, and all divideable by two, 
    #
    # If that is the case, also store the zeroth level.
    # -------------------------------------------------------
    if (shape     >= 4).all() and \
       (shape % 2 == 0).all():
            shape_ = (shape,)                # resolution
            a_     = (a,)                    # system matrix
            phi_   = (phi,)                  # unknown
            r_     = (Unknown("res_%d" % 0,  # residual at level "0"
                             phi.pos, 
                             shape_[0], -1, 
                             per = phi.per, 
                             verbatim = False),)
            b_     = (b,)                    # source
            d_     = (zeros(shape),)         # to store diffusion term 
            i_     = (zeros(shape),)         # to store innertial term 
    else:
        if verbatim:
            print("  Coarsening not possible!")
        return 1, None, None, None, None, None, None, None

    print(shape_[0])
        
    # ------------------------------
    # Count the levels and allocate 
    # memory for coarser levels
    # ------------------------------
    grid = 1  # level counter
    while (array(shape_[grid-1])      > 2).all() and  \
          (array(shape_[grid-1]) % 2 == 0).all() :
        print("  coarsening!")

        shape_ += ((shape_[grid-1][0]//2,       # resolution at level "grid"
                    shape_[grid-1][1]//2, 
                    shape_[grid-1][2]//2),)
        a_     += (Matrix(shape_[grid]),)
        phi_   += (Unknown("phi_%d" % grid,     # unknown at level "grid"
                           phi.pos, 
                           shape_[grid], -1, 
                           per = phi.per, 
                           verbatim = False),)
        r_     += (Unknown("res_%d" % grid,     # residual at level "grid"
                           phi.pos, 
                           shape_[grid], -1, 
                           per = phi.per, 
                           verbatim = False),)
        b_     += (zeros(shape_[grid]),)        # source at level "grid"
        d_     += (zeros(shape_[grid]),)        # diffusion at level "grid"
        i_     += (zeros(shape_[grid]),)        # innertia at level "grid"

        grid += 1                               # increase the level counter 
                              
    n_ = grid
    
    # -------------------------------------------------------------------------
    #
    # Browse through all grid levels
    #
    # -------------------------------------------------------------------------
    for grid in range(1, n_):
        print("  Resolution at level %d = " % grid, shape_[grid])

        # Compute diffusion term at previous (finer)
        d_[grid-1][:] = a_[grid-1].W[:] + a_[grid-1].E[:]  \
                      + a_[grid-1].S[:] + a_[grid-1].N[:]  \
                      + a_[grid-1].B[:] + a_[grid-1].T[:]
        
        # Store innertial term at previous (finer) level 
        # (Essentially than means: remove diffusion)
        i_[grid-1][:] = a_[grid-1].C[:] - d_[grid-1][:]

        #------------------------------------------------------------------
        # Think of the order of magnitude for neighboring coefficients 
        # derived from finite volume method in, a geometric progression
        #
        # Grid 0: delta = 1, coeff = delta^2 / delta ~ 1
        # Grid 1: delta = 2, coeff = delta^2 / delta ~ 2
        # Grid 3: delta = 4, coeff = delta^2 / delta ~ 4
        # Grid 4: delta = 8, coeff = delta^2 / delta ~ 8
        #
        # That means that neighboring coefficients scale with grid size 
        # on coarser levels.  That explains "0.5" in summation expressions 
        # for neighbouring coefficients below.
        #------------------------------------------------------------------
        
        # West coefficient
        a_[grid].W[:,:,:] =  a_[grid-1].W[0::2, 0::2, 0::2]
        a_[grid].W[:,:,:] += a_[grid-1].W[0::2, 1::2, 0::2]
        a_[grid].W[:,:,:] += a_[grid-1].W[0::2, 0::2, 1::2]
        a_[grid].W[:,:,:] += a_[grid-1].W[0::2, 1::2, 1::2]
        a_[grid].W[:,:,:] = a_[grid].W[:,:,:] * 0.5

        # West coefficient
        a_[grid].E[:,:,:] =  a_[grid-1].E[1::2, 0::2, 0::2]
        a_[grid].E[:,:,:] += a_[grid-1].E[1::2, 1::2, 0::2]
        a_[grid].E[:,:,:] += a_[grid-1].E[1::2, 0::2, 1::2]
        a_[grid].E[:,:,:] += a_[grid-1].E[1::2, 1::2, 1::2]
        a_[grid].E[:,:,:] = a_[grid].E[:,:,:] * 0.5

        # South coefficient
        a_[grid].S[:,:,:] =  a_[grid-1].S[0::2, 0::2, 0::2]
        a_[grid].S[:,:,:] += a_[grid-1].S[1::2, 0::2, 0::2]
        a_[grid].S[:,:,:] += a_[grid-1].S[0::2, 0::2, 1::2]
        a_[grid].S[:,:,:] += a_[grid-1].S[1::2, 0::2, 1::2]
        a_[grid].S[:,:,:] = a_[grid].S[:,:,:] * 0.5

        # North coefficient
        a_[grid].N[:,:,:] =  a_[grid-1].N[0::2, 1::2, 0::2]
        a_[grid].N[:,:,:] += a_[grid-1].N[1::2, 1::2, 0::2]
        a_[grid].N[:,:,:] += a_[grid-1].N[0::2, 1::2, 1::2]
        a_[grid].N[:,:,:] += a_[grid-1].N[1::2, 1::2, 1::2]
        a_[grid].N[:,:,:] = a_[grid].N[:,:,:] * 0.5

        # Bottom coefficient
        a_[grid].B[:,:,:] =  a_[grid-1].B[0::2, 0::2, 0::2]
        a_[grid].B[:,:,:] += a_[grid-1].B[1::2, 0::2, 0::2]
        a_[grid].B[:,:,:] += a_[grid-1].B[0::2, 1::2, 0::2]
        a_[grid].B[:,:,:] += a_[grid-1].B[1::2, 1::2, 0::2]
        a_[grid].B[:,:,:] = a_[grid].B[:,:,:] * 0.5

        # Top coefficient
        a_[grid].T[:,:,:] =  a_[grid-1].T[0::2, 0::2, 1::2]
        a_[grid].T[:,:,:] += a_[grid-1].T[1::2, 0::2, 1::2]
        a_[grid].T[:,:,:] += a_[grid-1].T[0::2, 1::2, 1::2]
        a_[grid].T[:,:,:] += a_[grid-1].T[1::2, 1::2, 1::2]
        a_[grid].T[:,:,:] = a_[grid].T[:,:,:] * 0.5

        # Central coefficient goes in two stages:
        # 1. add the innertial term like a simple summation
        a_[grid].C[:,:,:] =  i_[grid-1][0::2, 0::2, 0::2]
        a_[grid].C[:,:,:] += i_[grid-1][0::2, 1::2, 0::2]
        a_[grid].C[:,:,:] += i_[grid-1][0::2, 0::2, 1::2]
        a_[grid].C[:,:,:] += i_[grid-1][0::2, 1::2, 1::2]
        a_[grid].C[:,:,:] += i_[grid-1][1::2, 0::2, 0::2]
        a_[grid].C[:,:,:] += i_[grid-1][1::2, 1::2, 0::2]
        a_[grid].C[:,:,:] += i_[grid-1][1::2, 0::2, 1::2]
        a_[grid].C[:,:,:] += i_[grid-1][1::2, 1::2, 1::2]

        # 2. add newly formed diffusion terms 
        a_[grid].C[:,:,:] += (  a_[grid].W[:,:,:] + a_[grid].E[:,:,:]
                              + a_[grid].S[:,:,:] + a_[grid].N[:,:,:]
                              + a_[grid].B[:,:,:] + a_[grid].T[:,:,:]) 

    return n_, shape_, a_, i_, d_, phi_, b_, r_  # end of function
