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
from pyns.solvers.Matrix              import Matrix
from pyns.solvers.multigrid.coarsable import coarsable
from pyns.solvers.multigrid.coarsen   import coarsen

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
    if coarsable(shape):
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
    while coarsable(shape_[grid-1]):
      
        print("  coarsening!")

        shape_ += (coarsen(shape_[grid-1]),)    # resolution at level "grid"
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

        # Compute ratio between grid levels 
        rx = shape_[grid-1][X] // shape_[grid][X]
        ry = shape_[grid-1][Y] // shape_[grid][Y]
        rz = shape_[grid-1][Z] // shape_[grid][Z]
        
        # Lower bounds for browsing through grid levels
        wl = 0
        el = 1
        sl = 0
        nl = 1
        bl = 0
        tl = 1

        if rx < 2:
            el = 0
        if ry < 2:
            nl = 0
        if rz < 2:
            tl = 0

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
        a_[grid].W[:] =  a_[grid-1].W[wl::rx, sl::ry, bl::rz]
        a_[grid].W[:] += a_[grid-1].W[wl::rx, nl::ry, bl::rz]
        a_[grid].W[:] += a_[grid-1].W[wl::rx, sl::ry, tl::rz]
        a_[grid].W[:] += a_[grid-1].W[wl::rx, nl::ry, tl::rz]
        a_[grid].W[:] =  a_[grid].W[:] / (rx * (3-ry) * (3-rz))

        # East coefficient
        a_[grid].E[:] =  a_[grid-1].E[el::rx, sl::ry, bl::rz]
        a_[grid].E[:] += a_[grid-1].E[el::rx, nl::ry, bl::rz]
        a_[grid].E[:] += a_[grid-1].E[el::rx, sl::ry, tl::rz]
        a_[grid].E[:] += a_[grid-1].E[el::rx, nl::ry, tl::rz]
        a_[grid].E[:] =  a_[grid].E[:] / (rx * (3-ry) * (3-rz))

        # South coefficient
        a_[grid].S[:] =  a_[grid-1].S[wl::rx, sl::ry, bl::rz]
        a_[grid].S[:] += a_[grid-1].S[el::rx, sl::ry, bl::rz]
        a_[grid].S[:] += a_[grid-1].S[wl::rx, sl::ry, tl::rz]
        a_[grid].S[:] += a_[grid-1].S[el::rx, sl::ry, tl::rz]
        a_[grid].S[:] =  a_[grid].S[:] / ((3-rx) * ry * (3-rz))

        # North coefficient
        a_[grid].N[:] =  a_[grid-1].N[wl::rx, nl::ry, bl::rz]
        a_[grid].N[:] += a_[grid-1].N[el::rx, nl::ry, bl::rz]
        a_[grid].N[:] += a_[grid-1].N[wl::rx, nl::ry, tl::rz]
        a_[grid].N[:] += a_[grid-1].N[el::rx, nl::ry, tl::rz]
        a_[grid].N[:] =  a_[grid].N[:] / ((3-rx) * ry * (3-rz))

        # Bottom coefficient
        a_[grid].B[:] =  a_[grid-1].B[wl::rx, sl::ry, bl::rz]
        a_[grid].B[:] += a_[grid-1].B[el::rx, sl::ry, bl::rz]
        a_[grid].B[:] += a_[grid-1].B[wl::rx, nl::ry, bl::rz]
        a_[grid].B[:] += a_[grid-1].B[el::rx, nl::ry, bl::rz]
        a_[grid].B[:] =  a_[grid].B[:] / ((3-rx) * (3-ry) * rz)

        # Top coefficient
        a_[grid].T[:] =  a_[grid-1].T[wl::rx, sl::ry, tl::rz]
        a_[grid].T[:] += a_[grid-1].T[el::rx, sl::ry, tl::rz]
        a_[grid].T[:] += a_[grid-1].T[wl::rx, nl::ry, tl::rz]
        a_[grid].T[:] += a_[grid-1].T[el::rx, nl::ry, tl::rz]
        a_[grid].T[:] =  a_[grid].T[:] / ((3-rx) * (3-ry) * rz)

        # Central coefficient goes in two stages:
        # 1. add the innertial term like a simple summation
        a_[grid].C[:] =  i_[grid-1][wl::rx, sl::ry, bl::rz]
        a_[grid].C[:] += i_[grid-1][wl::rx, nl::ry, bl::rz]
        a_[grid].C[:] += i_[grid-1][wl::rx, sl::ry, tl::rz]
        a_[grid].C[:] += i_[grid-1][wl::rx, nl::ry, tl::rz]
        a_[grid].C[:] += i_[grid-1][el::rx, sl::ry, bl::rz]
        a_[grid].C[:] += i_[grid-1][el::rx, nl::ry, bl::rz]
        a_[grid].C[:] += i_[grid-1][el::rx, sl::ry, tl::rz]
        a_[grid].C[:] += i_[grid-1][el::rx, nl::ry, tl::rz]
        a_[grid].C[:] =  a_[grid].C[:] / ((3-rx) * (3-ry) * (3-rz))

        # 2. add newly formed diffusion terms 
        a_[grid].C[:] += (  a_[grid].W[:] + a_[grid].E[:]
                          + a_[grid].S[:] + a_[grid].N[:]
                          + a_[grid].B[:] + a_[grid].T[:]) 

    return n_, shape_, a_, i_, d_, phi_, b_, r_  # end of function
