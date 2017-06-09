"""
Creates system matrix and right hand side for both cell-centered and
staggered variables.
"""

# Standard Python modules
from standard import *

# ScriNS modules
from scrins.constants.boundary_conditions import DIRICHLET, NEUMANN
from scrins.constants.coordinates import X, Y, Z
from scrins.constants.compass import W, E, S, N, B, T, C
from scrins.discretization.obst_mod_matrix import obst_mod_matrix
from scrins.operators.avg import avg, avg_x, avg_y, avg_z
#from scrins.operators.avg_x import avg_x
#from scrins.operators.avg_y import avg_y
#from scrins.operators.avg_z import avg_z
from scrins.operators.cat import cat, cat_x, cat_y, cat_z
#from scrins.operators.cat_x import cat_x
#from scrins.operators.cat_y import cat_y
#from scrins.operators.cat_z import cat_z

# =============================================================================
def create_matrix(phi, inn, mu, dxyz, obst, obc):
    """
    Args:
      phi:  Unknown (created by "pyns.create_unknown") for which we need the
            system matrix.  It is needed for its resolution, essentially.
      inn:  Three-dimensional matrix holding the innertial term (innertial term
            is whatever multiplies the time derivative)
      mu:   Three-dimensional array holding diffusion coefficient.
      dxyz: Tuple holding cell dimensions in "x", "y" and "z" directions.
            Each cell dimension is a three-dimensional matrix.
      obst: Obstacle, three-dimensional matrix with zeros and ones.  It is
            zero in fluid, one in solid.
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

    # ------------------------------
    # Create right hand side vector
    # ------------------------------
    b = zeros(res)

    # -----------------------------------
    # Create default matrix coefficients
    # -----------------------------------
    coefficients = namedtuple('matrix_diagonal', 'W E S N B T P')
    c = coefficients(zeros(res), zeros(res), zeros(res),
                     zeros(res), zeros(res), zeros(res),
                     zeros(res))

    d = phi.pos

    # Handle central coefficient due to innertia
    c.P[:] = avg(d, inn) * avg(d, dx * dy * dz)

    # Pre-compute geometrical quantities
    sx = dy * dz
    sy = dx * dz
    sz = dx * dy

    # In the following lines, coefficients are computed simply by multiplying
    # diffusivity (mu) with cell-face areas (sx, sy and sz) and dividing by
    # distance between cells (dx, dy and dz).  Near the boundaries, only
    # half-distance between cells is taken into account.
    #
    #    +-----------+-----------+//
    #    |           |           |//
    #    |     o     |     o     |//  wall
    #    |           |           |//
    #    +-----------+-----------+//
    #          |           |     |
    #          |    dx     | dx/2|
    #          |<--------->|<--->|
    #
    if d != X:
        c.W[:] = cat_x((
            avg(d, mu[:1, :, :]) * avg(d, sx[:1, :, :])
                                 / avg(d, (dx[:1, :, :])/2.0),
            avg(d, avg_x(mu)) * avg(d, avg_x(sx))
                              / avg(d, avg_x(dx))))

        c.E[:] = cat_x((
            avg(d, avg_x(mu)) * avg(d, avg_x(sx))
                              / avg(d, avg_x(dx)),
            avg(d, mu[-1:, :, :]) * avg(d, sx[-1:, :, :])
                                  / avg(d, (dx[-1:, :, :]) / 2.0)))

    if d != Y:
        c.S[:] = cat_y((
            avg(d, mu[:, :1, :]) * avg(d, sy[:, :1, :])
                                 / avg(d, (dy[:, :1, :]) / 2.0),
            avg(d, avg_y(mu)) * avg(d, avg_y(sy))
                              / avg(d, avg_y(dy))))

        c.N[:] = cat_y((
            avg(d, avg_y(mu)) * avg(d, avg_y(sy))
                              / avg(d, avg_y(dy)),
            avg(d, mu[:, -1:, :]) * avg(d, sy[:, -1:, :])
                                  / avg(d, (dy[:, -1:, :]) / 2.0)))

    if d != Z:
        c.B[:] = cat_z((
            avg(d, mu[:, :, :1]) * avg(d, sz[:, :, :1])
                                 / avg(d, (dz[:, :, :1]) / 2.0),
            avg(d, avg_z(mu)) * avg(d, avg_z(sz))
                              / avg(d, avg_z(dz))))

        c.T[:] = cat_z((
            avg(d, avg_z(mu)) * avg(d, avg_z(sz))
                              / avg(d, avg_z(dz)),
            avg(d, mu[:, :, -1:]) * avg(d, sz[:, :, -1:])
                                  / avg(d, (dz[:, :, -1:]) / 2.0)))

    # Correct for staggered variables.  For staggered variables, near wall
    # cells are not at a half-distance from the wall, but at full distance.
    #
    #    +-----------+-----------+//
    #    |           |           |//
    #   -->         -->          |//  wall
    #    |           |           |//
    #    +-----------+-----------+//
    #    |           |           |
    #    |    dx     |     dx    |
    #    |<--------->|<--------->|
    #
    if d == X:
        c.W[:] = mu[0:-1, :, :] * sx[0:-1, :, :] / dx[0:-1, :, :]
        c.E[:] = mu[1:, :, :] * sx[1:, :, :] / dx[1:, :, :]
    elif d == Y:
        c.S[:] = mu[:, 0:-1, :] * sy[:, 0:-1, :] / dy[:, 0:-1, :]
        c.N[:] = mu[:, 1:, :] * sy[:, 1:, :] / dy[:, 1:, :]
    elif d == Z:
        c.B[:] = mu[:, :, 0:-1] * sz[:, :, 0:-1] / dz[:, :, 0:-1]
        c.T[:] = mu[:, :, 1:] * sz[:, :, 1:] / dz[:, :, 1:]

    # ----------------------------------------------------------------------
    # Zero them (correct them) for vanishing derivative boundary condition.
    # ----------------------------------------------------------------------

    # The values defined here will be false (numerical value 0)
    # wherever there is or Neumann boundary condition.
    c.W[:1, :, :] *= (phi.bnd[W].typ[:] == DIRICHLET)
    c.E[-1:, :, :] *= (phi.bnd[E].typ[:] == DIRICHLET)
    c.S[:, :1, :] *= (phi.bnd[S].typ[:] == DIRICHLET)
    c.N[:, -1:, :] *= (phi.bnd[N].typ[:] == DIRICHLET)
    c.B[:, :, :1] *= (phi.bnd[B].typ[:] == DIRICHLET)
    c.T[:, :, -1:] *= (phi.bnd[T].typ[:] == DIRICHLET)

    # -------------------------------------------
    # Fill the source terms with boundary values
    # -------------------------------------------
    b[:1, :, :] += c.W[:1, :, :] * phi.bnd[W].val[:1, :, :]
    b[-1:, :, :] += c.E[-1:, :, :] * phi.bnd[E].val[:1, :, :]
    b[:, :1, :] += c.S[:, :1, :] * phi.bnd[S].val[:, :1, :]
    b[:, -1:, :] += c.N[:, -1:, :] * phi.bnd[N].val[:, :1, :]
    b[:, :, :1] += c.B[:, :, :1] * phi.bnd[B].val[:, :, :1]
    b[:, :, -1:] += c.T[:, :, -1:] * phi.bnd[T].val[:, :, :1]

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

    c.W[:1, :, :] = 0.0
    c.E[-1:, :, :] = 0.0
    c.S[:, :1, :] = 0.0
    c.N[:, -1:, :] = 0.0
    c.B[:, :, :1] = 0.0
    c.T[:, :, -1:] = 0.0

    # ---------------------
    # Create sparse matrix
    # ---------------------
    nx, ny, nz = res
    n = nx * ny * nz

    data = array([reshape(c.P, n),
                  reshape(-c.W, n), reshape(-c.E, n),
                  reshape(-c.S, n), reshape(-c.N, n),
                  reshape(-c.B, n), reshape(-c.T, n)])

    diag = array([0, +ny * nz, -ny * nz, +nz, -nz, +1, -1])

    return spdiags(data, diag, n, n), b  # end of function
