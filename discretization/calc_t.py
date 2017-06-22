"""
Discretizes and solves equation for temperature.

Note:
  It should, however, be usable for any scalar.
"""

# Standard Python modules
from pyns.standard import *

# PyNS modules
from pyns.constants import *
from pyns.operators import *

from pyns.discretization.adj_n_bnds     import adj_n_bnds
from pyns.discretization.advection      import advection
from pyns.discretization.create_matrix  import create_matrix
from pyns.solvers                       import cg, cgs, bicgstab

# =============================================================================
def calc_t(t, uvwf, rho_cap, kappa, dt, dxyz, obst, 
           source = None):
# -----------------------------------------------------------------------------
    """
    Args:
      t:       Temperature unknown (from "create_unknown" function)
      uvwf:    a tuple with three staggered velocity components (where each
               component is created with "create_unknown" function.
      rho_cap: Three-dimensional array holding density times thermal
               capactity for all cells.
      dt:      Time step.
      dxyz:    Tuple holding cell dimensions in "x", "y" and "z" directions.
               Each cell dimension is a three-dimensional array.
      obst:    Obstacle, three-dimensional array with zeros and ones.
               It is zero in fluid, one in solid.

    Returns:
      none, but input argument t is modified!

    Note:
      Source (or sink) term is missing.
    """

    # Unpack tuple(s)
    dx, dy, dz = dxyz

    # Fetch the resolution
    rc = t.val.shape

    # Discretize the diffusive part
    A_t = create_matrix(t, rho_cap/dt, kappa, dxyz, obst, NEUMANN)
    b_t = zeros(rc)

    # The advective fluxes
    c_t = advection(rho_cap, t, uvwf, dxyz, dt, 'minmod')

    # Innertial term for enthalpy
    i_t = t.old * avg(t.pos, rho_cap) * avg(t.pos, dx*dy*dz) / dt

    # Handle external source
    if source is None:
        s_t = zeros(r_phi)
    else:
        s_t = source * avg(t.pos, dx*dy*dz)

    # The entire source term
    f_t = b_t - c_t + i_t + s_t

    # Solve for temperature
    t.val[:] = cgs(A_t, t, f_t, TOL, False)

    adj_n_bnds(t)

    return  # end of function
