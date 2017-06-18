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

# =============================================================================
def create_matrix(phi, inn, mu, dxyz, obst, obc):
# -----------------------------------------------------------------------------
    """
    Args:
      phi:  Unknown (created by "pyns.create_unknown") for which we need the
            system matrix.  It is needed for its resolution, essentially.
      inn:  Three-dimensional matrix holding the innertial term (innertial
            term is whatever multiplies the time derivative)
      mu:   Three-dimensional array holding diffusion coefficient.
      dxyz: Tuple holding cell dimensions in "x", "y" and "z" directions.
            Each cell dimension is a three-dimensional matrix.
      obst: Obstacle, three-dimensional matrix with zeros and ones.
            It is zero in fluid, one in solid.
      obc:  Obstacle's boundary condition (NEUMANN or DIRICHLET)

    Returns:
      A: Matrix in sparse diagonal format.
      b: Three-dimensional matrix for the right hand side vector.

    Note:
      It is a bit of missnomer, it is called "create_matrix", but it also
      creates right hand side vector.
    """

    # Unpack tuples
    dx, dy, dz = dxyz

    res = phi.val.shape

    # -----------------------------------
    # Create default matrix coefficients
    # -----------------------------------
    coefficients = namedtuple('matrix_diagonal', 'W E S N B T P')
    c = coefficients(zeros(res), zeros(res), zeros(res),   \
                     zeros(res), zeros(res), zeros(res),   \
                     zeros(res))

    d = phi.pos

    # Handle central coefficient due to innertia
    c.P[:] = avg(d, inn) * avg(d, dx*dy*dz)

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
    #     - --+---------+---------+-- -   - --+---------+---------+-- -
    #         |         |         |           |         |         |   
    #    o    |    o    |    o    |           |    o    |    o    |    o
    #  nx-1   |    0    |    1    |           |  nx-2   |  nx-1   |    0
    #     - --+---------+---------+-- -   - --+---------+---------+-- -
    #    :         :         :                     :         :         :
    #    |   dx    |   dx    |                     |   dx    |   dx    |
    #    |<------->|<------->|                     |<------->|<------->|
    #
    if d != X:
        # West
        c_bc_x[:] = avg(d, mu[ :1,:,:]) * avg(d,  sx[ :1,:,:])         \
                                        / avg(d, (dx[ :1,:,:]) / 2.0)
        c.W[:] = cat_x((c_bc_x, avg(d, avg_x(mu)) * avg(d, avg_x(sx))
                                                  / avg(d, avg_x(dx))))
        # East
        c_bc_x[:] = avg(d, mu[-1:,:,:]) * avg(d,  sx[-1:,:,:])         \
                                        / avg(d, (dx[-1:,:,:]) / 2.0)      
        c.E[:] = cat_x((avg(d, avg_x(mu)) * avg(d, avg_x(sx))
                                          / avg(d, avg_x(dx)), c_bc_x))

    if d != Y:
        # South
        c_bc_y[:] = avg(d, mu[:, :1,:]) * avg(d,  sy[:, :1,:])         \
                                        / avg(d, (dy[:, :1,:]) / 2.0)
        c.S[:] = cat_y((c_bc_y, avg(d, avg_y(mu)) * avg(d, avg_y(sy))
                                                  / avg(d, avg_y(dy))))
        # North
        c_bc_y[:] = avg(d, mu[:,-1:,:]) * avg(d,  sy[:,-1:,:])         \
                                        / avg(d, (dy[:,-1:,:]) / 2.0)      
        c.N[:] = cat_y((avg(d, avg_y(mu)) * avg(d, avg_y(sy))
                                          / avg(d, avg_y(dy)), c_bc_y))

    if d != Z:
        # Bottom
        c_bc_z[:] = avg(d, mu[:,:, :1]) * avg(d,  sz[:,:, :1])         \
                                        / avg(d, (dz[:,:, :1]) / 2.0)            
        c.B[:] = cat_z((c_bc_z, avg(d, avg_z(mu)) * avg(d, avg_z(sz))
                                                  / avg(d, avg_z(dz))))
        # Top
        c_bc_z[:] = avg(d, mu[:,:,-1:]) * avg(d,  sz[:,:,-1:])         \
                                        / avg(d, (dz[:,:,-1:]) / 2.0)            
        c.T[:] = cat_z((avg(d, avg_z(mu)) * avg(d, avg_z(sz))
                                          / avg(d, avg_z(dz)), c_bc_z))

    # Correct for staggered variables.  For staggered variables, near wall
    # cells are not at a half-distance from the wall, but at full distance.
    #
    #    DIRICHLET:
    #
    #      ///+---------+---------+-- -   - --+---------+---------+///
    #      ///|         |         |           |         |         |///
    #      ///|        -->       -->         -->       -->        |/// 
    #      ///|         |0        |1          |nx-2     |nx-1     |///
    #      ///+---------+---------+-- -   - --+---------+---------+///
    #         :         :         :           :         :         :
    #         |   dx    |   dx    |           |   dx    |   dx    |
    #         |<------->|<------->|           |<------->|<------->|
    #
    #    PERIODIC:
    #
    #         +---------+---------+-- -   - --+---------+---------+-- -   
    #         |         |         |           |         |         |  
    #        -->       -->       -->         -->       -->       -->       -->   
    #         |nx-1     |0        |1          |         |nx-2     |nx-1      0
    #         +---------+---------+-- -   - --+---------+---------+-- -   
    #         :         :         :           :         :         :         :
    #         |   dx    |   dx    |           |   dx    |   dx    |   dx    |
    #         |<------->|<------->|           |<------->|<------->|<------->|
    #
    if d == X:
        c.W[:] = mu[0:-1,:,:] * sx[0:-1,:,:] / dx[0:-1,:,:]
        c.E[:] = mu[1:,  :,:] * sx[1:,  :,:] / dx[1:,  :,:]
    elif d == Y:
        c.S[:] = mu[:,0:-1,:] * sy[:,0:-1,:] / dy[:,0:-1,:]
        c.N[:] = mu[:,1:,  :] * sy[:,1:,  :] / dy[:,1:,  :]
    elif d == Z:
        c.B[:] = mu[:,:,0:-1] * sz[:,:,0:-1] / dz[:,:,0:-1]
        c.T[:] = mu[:,:,1:  ] * sz[:,:,1:  ] / dz[:,:,1:  ]

    # ----------------------------------------------------------------------
    # Zero them (correct them) for vanishing derivative boundary condition.
    # ----------------------------------------------------------------------

    # The values defined here will be false (numerical value 0)
    # wherever there is Neumann boundary condition.
    c.W[ :1,  :,  :] *= ( phi.bnd[W].typ[:] != NEUMANN )
    c.E[-1:,  :,  :] *= ( phi.bnd[E].typ[:] != NEUMANN )
    c.S[  :, :1,  :] *= ( phi.bnd[S].typ[:] != NEUMANN )
    c.N[  :,-1:,  :] *= ( phi.bnd[N].typ[:] != NEUMANN )
    c.B[  :,  :, :1] *= ( phi.bnd[B].typ[:] != NEUMANN )
    c.T[  :,  :,-1:] *= ( phi.bnd[T].typ[:] != NEUMANN )

    # --------------------------------------
    # Correct system matrices for obstacles
    # --------------------------------------
    if obst.any() != 0:
        c = obst_mod_matrix(phi, c, obst, obc)

    # ----------------------------------------------
    # Add all neighbours to the central matrix,
    # and zero the coefficients towards boundaries
    # ----------------------------------------------
    c.P[:] += c.W[:] + c.E[:] + c.S[:] + c.N[:] + c.B[:] + c.T[:]

    return c  # end of function
