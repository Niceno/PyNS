"""
Creates system matrix and right hand side for both cell-centered and
staggered variables.
"""

# Standard Python modules
from pyns.standard import *

# PyNS modules
from pyns.constants import *
from pyns.operators import *

from pyns.discretization.obst_mod_matrix import obst_mod_matrix
from pyns.solvers                        import Matrix

# =============================================================================
def diffusion(phi, inn, mu, dxyz, obst, obc):
# -----------------------------------------------------------------------------
    """
    Args:
      phi:  Object of the type "Unknown" for which we need the system matrix.
            It is needed for its resolution, essentially.
      inn:  Three-dimensional array holding the innertial term (innertial
            term is whatever multiplies the time derivative)
      mu:   Three-dimensional array holding diffusion coefficient.
      dxyz: Tuple holding cell dimensions in "x", "y" and "z" directions.
            Each cell dimension is a three-dimensional array.
      obst: Obstacle, three-dimensional array with zeros and ones.
            It is zero in fluid, one in solid.
      obc:  Obstacle's boundary condition (NEUMANN or DIRICHLET)

    Returns:
      c: Object of the type "Matrix", holding the system matrix.
    """

    # Unpack tuples
    dx, dy, dz = dxyz

    res = phi.val.shape

    # -----------------------------------
    # Create default matrix coefficients
    # -----------------------------------
    matrix = Matrix(res)

    d = phi.pos

    # Pre-compute geometrical quantities
    sx = dy * dz
    sy = dx * dz
    sz = dx * dy

    # Allocate memory for specification of boundary conditions
    nx, ny, nz = res
    c_bc_x = zeros(( 1, ny, nz))
    c_bc_y = zeros((nx,  1, nz))
    c_bc_z = zeros((nx, ny,  1))

    # In the following lines, coefficients are computed simply by multiplying
    # diffusivity (mu) with cell-face areas (sx, sy and sz) and dividing by
    # distance between cells (dx, dy and dz).  Near the boundaries, only
    # half-distance between cells is taken into account.
    #
    #    DIRICHLET:
    #
    #      ///+---------+---------+-- -   - --+---------+---------+///
    #      ///|         |         |           |         |         |///
    #      ///|    o    |    o    |           |    o    |    o    |///  
    #      ///|    0    |    1    |           |  nx-2   |  nx-1   |///
    #      ///+---------+---------+-- -   - --+---------+---------+///
    #         :    :    __   :                     :   __    :    :
    #         |dx/2|    dx   |                     |   dx    |dx/2|
    #         |<-->|<------->|                     |<------->|<-->|
    #
    #    PERIODIC:
    #
    #    - First cell (0) is the same as last cell (nx-1)
    #    - Thaty makes domain a bit shorter
    #
    #         +---------+---------+-- -   - --+---------+---------+
    #   [W]   |         |         |           |         |         |   [E]
    #    o    |    o    |    o    |           |    o    |    o    |    o
    #  nx-2   |    0    |    1    |           |  nx-2   |  nx-1   |    1
    #         +---------+---------+-- -   - --+---------+---------+
    #    :   __    :   __    :                     :   __    :   __    :
    #    |   dx    |   dx    |                     |   dx    |   dx    |
    #    |<------->|<------->|                     |<------->|<------->|
    #       (nx-2)     (0)                           (nx-2)      (0) 
    #
    #
    if d != X:
        # West
        c_bc_x[:] = avg(d, mu[ :1,:,:]) * avg(d,  sx[ :1,:,:])         \
                                        / avg(d, (dx[ :1,:,:]) / 2.0)
        matrix.W[:] = cat_x((c_bc_x, avg(d, avg_x(mu)) * avg(d, avg_x(sx))
                                                       / avg(d, avg_x(dx))))
        # East
        c_bc_x[:] = avg(d, mu[-1:,:,:]) * avg(d,  sx[-1:,:,:])         \
                                        / avg(d, (dx[-1:,:,:]) / 2.0)      
        matrix.E[:] = cat_x((avg(d, avg_x(mu)) * avg(d, avg_x(sx))
                                               / avg(d, avg_x(dx)), c_bc_x))
        # Correct for periodicity
        if phi.per[X] == True:
            matrix.W[ :1,:,:] = matrix.W[-1:,:,:]
            matrix.E[-1:,:,:] = matrix.E[ :1,:,:]

    if d != Y:
        # South
        c_bc_y[:] = avg(d, mu[:, :1,:]) * avg(d,  sy[:, :1,:])         \
                                        / avg(d, (dy[:, :1,:]) / 2.0)
        matrix.S[:] = cat_y((c_bc_y, avg(d, avg_y(mu)) * avg(d, avg_y(sy))
                                                  / avg(d, avg_y(dy))))
        # North
        c_bc_y[:] = avg(d, mu[:,-1:,:]) * avg(d,  sy[:,-1:,:])         \
                                        / avg(d, (dy[:,-1:,:]) / 2.0)      
        matrix.N[:] = cat_y((avg(d, avg_y(mu)) * avg(d, avg_y(sy))
                                          / avg(d, avg_y(dy)), c_bc_y))
        # Correct for periodicity
        if phi.per[Y] == True:
            matrix.S[:, :1,:] = matrix.S[:,-1:,:]
            matrix.N[:,-1:,:] = matrix.N[:, :1,:]

    if d != Z:
        # Bottom
        c_bc_z[:] = avg(d, mu[:,:, :1]) * avg(d,  sz[:,:, :1])         \
                                        / avg(d, (dz[:,:, :1]) / 2.0)            
        matrix.B[:] = cat_z((c_bc_z, avg(d, avg_z(mu)) * avg(d, avg_z(sz))
                                                  / avg(d, avg_z(dz))))
        # Top
        c_bc_z[:] = avg(d, mu[:,:,-1:]) * avg(d,  sz[:,:,-1:])         \
                                        / avg(d, (dz[:,:,-1:]) / 2.0)            
        matrix.T[:] = cat_z((avg(d, avg_z(mu)) * avg(d, avg_z(sz))
                                          / avg(d, avg_z(dz)), c_bc_z))
        # Correct for periodicity
        if phi.per[Z] == True:
            matrix.B[:,:, :1] = matrix.B[:,:,-1:]
            matrix.T[:,:,-1:] = matrix.T[:,:, :1]

    # Correct for staggered variables.  For staggered variables, near wall
    # cells are not at a half-distance from the wall, but at full distance.
    #
    #    DIRICHLET:
    #
    #      ///+---------+---------+-- -   - --+---------+---------+///
    #      ///|         |         |           |         |         |///
    #  [W] ///|        -->       -->         -->       -->        |/// [E]
    #      ///|         |0        |1          |nx-2     |nx-1     |///
    #      ///+---------+---------+-- -   - --+---------+---------+///
    #         :         :         :           :         :         :
    #         |   dx    |   dx    |           |   dx    |   dx    |
    #         |<------->|<------->|           |<------->|<------->|
    #
    #    PERIODIC:
    #
    #         +---------+---------+-- -   - --+---------+---------+
    #         |         |         |           |         |         |
    #  [W]   -->       -->       -->         -->       -->       -->   [E]
    #         |nx-2     |0        |1          |         |nx-2     |0
    #         +---------+---------+-- -   - --+---------+---------+   
    #         :         :         :           :         :         :
    #         |   dx    |   dx    |           |   dx    |   dx    |
    #         |<------->|<------->|           |<------->|<------->|
    #
    if d == X:
        matrix.W[:] = mu[0:-1,:,:] * sx[0:-1,:,:] / dx[0:-1,:,:]
        matrix.E[:] = mu[1:,  :,:] * sx[1:,  :,:] / dx[1:,  :,:]
    elif d == Y:
        matrix.S[:] = mu[:,0:-1,:] * sy[:,0:-1,:] / dy[:,0:-1,:]
        matrix.N[:] = mu[:,1:,  :] * sy[:,1:,  :] / dy[:,1:,  :]
    elif d == Z:
        matrix.B[:] = mu[:,:,0:-1] * sz[:,:,0:-1] / dz[:,:,0:-1]
        matrix.T[:] = mu[:,:,1:  ] * sz[:,:,1:  ] / dz[:,:,1:  ]

    # ----------------------------------------------------------------------
    # Zero them (correct them) for vanishing derivative boundary condition.
    # ----------------------------------------------------------------------

    # The values defined here will be false (numerical value 0)
    # wherever there is Neumann boundary condition.
    matrix.W[ :1,  :,  :] *= ( phi.bnd[W].typ[:] != NEUMANN )
    matrix.E[-1:,  :,  :] *= ( phi.bnd[E].typ[:] != NEUMANN )
    matrix.S[  :, :1,  :] *= ( phi.bnd[S].typ[:] != NEUMANN )
    matrix.N[  :,-1:,  :] *= ( phi.bnd[N].typ[:] != NEUMANN )
    matrix.B[  :,  :, :1] *= ( phi.bnd[B].typ[:] != NEUMANN )
    matrix.T[  :,  :,-1:] *= ( phi.bnd[T].typ[:] != NEUMANN )

    # --------------------------------------
    # Correct system matrices for obstacles
    # --------------------------------------
    if obst is not None:
        c = obst_mod_matrix(phi, c, obst, obc)

    # ----------------------------------------------
    # Add all neighbours to the central matrix,
    # and zero the coefficients towards boundaries
    # ----------------------------------------------
    matrix.C[:] += matrix.W[:] + matrix.E[:]  \
                +  matrix.S[:] + matrix.N[:]  \
                +  matrix.B[:] + matrix.T[:]

    return matrix  # end of function
