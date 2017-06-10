"""
Discretizes and solves equation for temperature.  

Note:
  It should, however, be usable for any scalar.
"""

# Standard Python modules
from pyns.standard import *

# PyNS modules
from pyns.constants.all      import *
from pyns.operators.all      import *

from pyns.discretization.adj_n_bnds     import adj_n_bnds
from pyns.discretization.advection      import advection
from pyns.discretization.create_matrix  import create_matrix
from pyns.solvers.all                   import cg, cgs, bicgstab

# =============================================================================
def calc_t(t, uvwf, rho_cap, kappa, dt, dxyz, obst):
# -----------------------------------------------------------------------------
  """
  Args:
    t:       Temperature unknown (from "create_unknown" function)
    uvwf:    a tuple with three staggered velocity components (where each 
             component is created with "create_unknown" function.
    rho_cap: Three-dimensional matrix holding density times thermal capactity 
             for all cells.
    dt:      Time step.
    dxyz:    Tuple holding cell dimensions in "x", "y" and "z" directions.
             Each cell dimension is a three-dimensional matrix.
    obst:    Obstacle, three-dimensional matrix with zeros and ones.  It is
             zero in fluid, one in solid.

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
  A_t, b_t = create_matrix(t, rho_cap/dt, kappa, dxyz, obst, NEUMANN)

  # The advective fluxes 
  c_t = advection(rho_cap, t, uvwf, dxyz, dt, 'minmod')
  
  # Innertial term for enthalpy
  i_t = t.old * rho_cap * dx*dy*dz / dt

  # The entire source term
  f_t = b_t - c_t + i_t

  # Solve for temperature 
  t.val[:] = cgs(A_t, t, f_t, TOL, False)

  adj_n_bnds(t)

  return  # end of function
