"""
Adjusts the system matrix for obstacles.
"""

# Standard Python modules
from pyns.standard import *

# PyNS modules
from pyns.constants import *
from pyns.operators import *

# =============================================================================
def obst_mod_matrix(phi, matrix, obstacle, obstacle_bc):
# -----------------------------------------------------------------------------
    """
    Args:
      phi: ....... Object of the type "Unknown".
      matrix: .... Object of the type "Matrix", holding the system matrix.
      obstacle: .. Obstacle, three-dimensional array with zeros and ones.  
                   It is zero in fluid, one in solid.
      obstacle_bc: Obstacle's boundary condition (NEUMANN or DIRICHLET).

    Returns:
      matrix: Modified system matrix.
    """

    pos = phi.pos

    # -------------------------
    #
    # For collocated variables
    #
    # -------------------------
    if pos == C:

        # -----------------------------------
        # Neumann's boundary on the obstacle
        # -----------------------------------
        if obstacle_bc == NEUMANN:

            # Correct west and east
            sol_x = dif_x(obstacle)  # +1 east of obstacle, -1 west of it
            corr = 1 - (sol_x < 0)
            matrix.W[1:,:,:] = matrix.W[1:,:,:] * corr
            corr = 1 - (sol_x > 0)
            matrix.E[:-1,:,:] = matrix.E[:-1,:,:] * corr

            # Correct south and north
            sol_y = dif_y(obstacle)  # +1 north of obstacle, -1 south of it
            corr = 1 - (sol_y < 0)
            matrix.S[:,1:,:] = matrix.S[:,1:,:] * corr
            corr = 1 - (sol_y > 0)
            matrix.N[:,:-1,:] = matrix.N[:,:-1,:] * corr

            # Correct bottom and top
            sol_z = dif_z(obstacle)  # +1 north of obstacle, -1 south of it
            corr = 1 - (sol_z < 0)
            matrix.B[:,:,1:] = matrix.B[:,:,1:] * corr
            corr = 1 - (sol_z > 0)
            matrix.T[:,:,:-1] = matrix.T[:,:,:-1] * corr

        # -------------------------------------
        # Dirichlet's boundary on the obstacle
        # -------------------------------------
        elif obstacle_bc == DIRICHLET:

            # Set central coefficient to 1 in obst, unchanged elsewhere
            matrix.C[:] = matrix.C[:] * lnot(obstacle) + obstacle

            # Set neighbour coefficients to zero in obst
            matrix.W[:] = matrix.W[:] * lnot(obstacle)
            matrix.E[:] = matrix.E[:] * lnot(obstacle)
            matrix.S[:] = matrix.S[:] * lnot(obstacle)
            matrix.N[:] = matrix.N[:] * lnot(obstacle)
            matrix.B[:] = matrix.B[:] * lnot(obstacle)
            matrix.T[:] = matrix.T[:] * lnot(obstacle)

            # Increase coefficients close to obst (makes sense for momentum)
            sol_x = dif_x(obstacle)  # +1 east of obstacle, -1 west of it
            corr = 1 + (sol_x > 0)
            matrix.E[:-1,:,:] = matrix.E[:-1,:,:] * corr
            corr = 1 + (sol_x < 0)
            matrix.W[1:,:,:] = matrix.W[1:,:,:] * corr

            sol_y = dif_y(obstacle)  # +1 north of obstacle, -1 south of it
            corr = 1 + (sol_y > 0)
            matrix.N[:,:-1,:] = matrix.N[:,:-1,:] * corr
            corr = 1 + (sol_y < 0)
            matrix.S[:,1:,:] = matrix.S[:,1:,:] * corr

            sol_z = dif_z(obstacle)  # +1 top of obstacle, -1 bottom of it
            corr = 1 + (sol_z > 0)
            matrix.T[:,:,:-1] = matrix.T[:,:,:-1] * corr
            corr = 1 + (sol_z < 0)
            matrix.B[:,:,1:] = matrix.B[:,:,1:] * corr

    # ------------------------
    #
    # For staggered variables
    #
    # ------------------------
    elif pos == X:

        # Set central coefficient to 1 in obst, unchanged elsewhere
        obst_x = mx(obstacle[:-1,:,:], obstacle[1:,:,:])
        matrix.C[:] = matrix.C[:] * lnot(obst_x) + obst_x

        # Set neighbour coefficients to zero in obst
        matrix.W[:] = matrix.W[:] * lnot(obst_x)
        matrix.E[:] = matrix.E[:] * lnot(obst_x)
        matrix.S[:] = matrix.S[:] * lnot(obst_x)
        matrix.N[:] = matrix.N[:] * lnot(obst_x)
        matrix.B[:] = matrix.B[:] * lnot(obst_x)
        matrix.T[:] = matrix.T[:] * lnot(obst_x)

        # Increase coefficients close to obst (makes sense for momentum)
        sol_y = dif_y(obst_x)  # will be +1 north of obst, -1 south of obst
        corr = 1 + (sol_y > 0)
        matrix.N[:,:-1,:] = matrix.N[:,:-1,:] * corr
        corr = 1 + (sol_y < 0)
        matrix.S[:,1:,:] = matrix.S[:,1:,:] * corr

        sol_z = dif_z(obst_x)  # will be +1 top of obst, -1 bottom of obst
        corr = 1 + (sol_z > 0)
        matrix.T[:,:,:-1] = matrix.T[:,:,:-1] * corr
        corr = 1 + (sol_z < 0)
        matrix.B[:,:,1:] = matrix.B[:,:,1:] * corr

    elif pos == Y:

        # Set central coefficient to 1 in obst, unchanged elsewhere
        obst_y = mx(obstacle[:,:-1,:], obstacle[:,1:,:])
        matrix.C[:] = matrix.C[:] * lnot(obst_y) + obst_y

        # Set neighbour coefficients to zero in obst
        matrix.W[:] = matrix.W[:] * lnot(obst_y)
        matrix.E[:] = matrix.E[:] * lnot(obst_y)
        matrix.S[:] = matrix.S[:] * lnot(obst_y)
        matrix.N[:] = matrix.N[:] * lnot(obst_y)
        matrix.B[:] = matrix.B[:] * lnot(obst_y)
        matrix.T[:] = matrix.T[:] * lnot(obst_y)

        # Increase coefficients close to obst (makes sense for momentum)
        sol_x = dif_x(obst_y)  # will be +1 north of obst, -1 south of obst
        corr = 1 + (sol_x > 0)
        matrix.E[:-1,:,:] = matrix.E[:-1,:,:] * corr
        corr = 1 + (sol_x < 0)
        matrix.W[1:,:,:] = matrix.W[1:,:,:] * corr

        sol_z = dif_z(obst_y)  # will be +1 north of obst, -1 south of obst
        corr = 1 + (sol_z > 0)
        matrix.T[:,:,:-1] = matrix.T[:,:,:-1] * corr
        corr = 1 + (sol_z < 0)
        matrix.B[:,:,1:] = matrix.B[:,:,1:] * corr

    elif pos == Z:

        # Set central coefficient to 1 in obst, unchanged elsewhere
        obst_z = mx(obstacle[:,:,:-1], obstacle[:,:,1:])
        matrix.C[:] = matrix.C[:] * lnot(obst_z) + obst_z

        # Set neighbour coefficients to zero in obst
        matrix.W[:] = matrix.W[:] * lnot(obst_z)
        matrix.E[:] = matrix.E[:] * lnot(obst_z)
        matrix.S[:] = matrix.S[:] * lnot(obst_z)
        matrix.N[:] = matrix.N[:] * lnot(obst_z)
        matrix.B[:] = matrix.B[:] * lnot(obst_z)
        matrix.T[:] = matrix.T[:] * lnot(obst_z)

        # Increase coefficients close to obst (makes sense for momentum)
        sol_x = dif_x(obst_z)  # will be +1 north of obst, -1 south of obst
        corr = 1 + (sol_x > 0)
        matrix.E[:-1,:,:] = matrix.E[:-1,:,:] * corr
        corr = 1 + (sol_x < 0)
        matrix.W[1:,:,:] = matrix.W[1:,:,:] * corr

        sol_y = dif_y(obst_z)  # will be +1 north of obst, -1 south of obst
        corr = 1 + (sol_y > 0)
        matrix.N[:,:-1,:] = matrix.N[:,:-1,:] * corr
        corr = 1 + (sol_y < 0)
        matrix.S[:,1:,:] = matrix.S[:,1:,:] * corr

    return matrix  # end of function
