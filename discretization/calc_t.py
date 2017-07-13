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

from pyns.discretization.adj_n_bnds import adj_n_bnds
from pyns.discretization.advection  import advection
from pyns.discretization.diffusion  import diffusion
from pyns.solvers.nonstationary     import cg, cgs, bicgstab

# =============================================================================
def calc_t(t, uvwf, rho_cap, kappa, dt, dxyz, 
           obstacle = None, 
           source = None,
           under_relaxation = 1.0,
           advection_scheme = "minmod"):
# -----------------------------------------------------------------------------
    """
    Args:
      t: .............. Object of type "Unknown" holding the temperature.
      uvwf: ........... Tuple with three staggered velocity components (where 
                        each component is object of the type "Unknown").
      rho_cap: ........ Three-dimensional array holding density times thermal 
                        capacity for cells.
      kappa: .......... Three-dimensional array holding conductivity for cells.
      dt: ............. Time step.
      dxyz: ........... Tuple holding cell dimensions in "x", "y" and "z" 
                        directions.  Each component (each cell dimension) is 
                        a three-dimensional array.
      obstacle: ....... Obstacle, three-dimensional array with zeros and ones.
                        It is zero in fluid, one in solid.
      source: ......... External source.
      under_relaxation: Under relaxation factor.
      advection_scheme: Advection scheme.

    Returns:
      None, but input argument t is modified!
    """

    # Unpack tuple(s)
    dx, dy, dz = dxyz

    # Fetch the resolution
    rt = t.val.shape

    # Discretize the diffusive part
    A_t = diffusion(t, rho_cap/dt, kappa, dxyz, obstacle, NEUMANN)
    b_t = zeros(rt)

    # The advective fluxes
    c_t = advection(rho_cap, t, uvwf, dxyz, dt, advection_scheme,
                    matrix = A_t)

    # Innertial term for enthalpy
    A_t.C       += avg(t.pos, rho_cap) * avg(t.pos, dx*dy*dz) / dt
    i_t = t.old  * avg(t.pos, rho_cap) * avg(t.pos, dx*dy*dz) / dt

    # Handle external source
    if source is None:
        s_t = zeros(rt)
    else:
        s_t = source * avg(t.pos, dx*dy*dz)

    # The entire source term
    f_t = b_t - c_t + i_t + s_t

    # Solve for temperature
    ur = under_relaxation
    t.val[:] = (1-ur)*t.val[:] + ur * bicgstab(A_t, t, f_t, TOL, False)

    adj_n_bnds(t)

    return  # end of function
