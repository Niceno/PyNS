"""
Discretizes and solves equation for temperature.  

Note:
  It should, however, be usable for any scalar.
"""

# Standard Python modules
from standard import *

# ScriNS modules
from constants.all      import *
from operators.all      import *

from discretization.adj_n_bnds     import adj_n_bnds
from discretization.advection      import advection
from discretization.create_matrix  import create_matrix

# =============================================================================
def calc_t(t, uvwf, rho_cap, kappa, dt, dxyz, obst):
# -----------------------------------------------------------------------------
  """
  Args:
    t:       temperature unknown (from "pyns.create_unknown" function)
    uvwf:    a tuple with three staggered velocity components (where each 
             component is created with "scrins.create_unknown" function.
    rho_cap: three-dimensional matrix holding density times thermal capactity 
             for all cells.
    dt:      time step
    dxyz:    a tuple holding cell dimensions in "x", "y" and "z" directions.
             Each cell dimension is a three-dimensional matrix.
    obst:    obstacle

  Returns:
    none, but input argument t is modified!      
  """

  # Unpack tuple(s)
  dx, dy, dz = dxyz

  # Fetch the resolution  
  rc = t.val.shape 

  # Discretize the diffusive part
  A_t, b_t = create_matrix(t, rho_cap/dt, kappa, dxyz, obst, NEUMANN)

  # The advective fluxes 
  c_t = advection(rho_cap, t, uvwf, dxyz, dt, 'minmod')
  
  # Innertial term for enthalpy
  i_t = t.old * rho_cap * dx*dy*dz / dt

  # The entire source term
  f_t = b_t - c_t + i_t

  # Solve for temperature 
  res0 = bicgstab( A_t, reshape(f_t, prod(rc)), tol=TOL )
  t.val[:] = reshape(res0[0], rc) 

  adj_n_bnds(t)

  return  # end of function
