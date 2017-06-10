"""
Projects velocity into a divergence-free velocity field with the pressure 
correction gradient (provided pressure correction is accuratelly solved).
"""

# Standard Python modules
from pyns.standard import *

# PyNS modules
from pyns.constants.all      import *
from pyns.operators.all      import *

from pyns.discretization.obst_zero_val import obst_zero_val

# =============================================================================
def corr_uvw(uvw, p, rho, dt, dxyz, obst):
# -----------------------------------------------------------------------------
  """
  Args:
    uvw:  Tuple with three staggered or centered velocity components 
          (each component is created with "create_unknown" function.
    p:    Unknown holding the pressure correction.     
    rho:  Three-dimensional matrix holding density for all cells.
    dt:   Time step
    dxyz: Tuple holding cell dimensions in "x", "y" and "z" directions.
          Each cell dimension is a three-dimensional array.
    obst: Obstacle, three-dimensional matrix with zeros and ones.  It is
          zero in fluid, one in solid.
    
  Returns:
    none, but input argument uvw is modified.
  """

  # Unpack received tuples
  dx, dy, dz = dxyz

  # Compute pressure correction gradients
  p_x = dif_x(p.val) / avg_x(dx)  
  p_y = dif_y(p.val) / avg_y(dy)  
  p_z = dif_z(p.val) / avg_z(dz)  
  
  # Set to zero in obst
  if obst.any() != 0:
    p_x = obst_zero_val(X, p_x, obst)
    p_y = obst_zero_val(Y, p_y, obst)
    p_z = obst_zero_val(Z, p_z, obst)
  
  # Pad with boundary values by expanding from interior 
  # (This is done only for collocated formulation)
  if uvw[X].pos == C:
    p_x = avg_x(cat_x((p_x[:1,:,:], p_x, p_x[-1:,:,:])))
    p_y = avg_y(cat_y((p_y[:,:1,:], p_y, p_y[:,-1:,:])))
    p_z = avg_z(cat_z((p_z[:,:,:1], p_z, p_z[:,:,-1:])))
  
  # Correct the velocities
  uvw[X].val[:] = uvw[X].val[:] - dt / avg(uvw[X].pos, rho) * p_x
  uvw[Y].val[:] = uvw[Y].val[:] - dt / avg(uvw[Y].pos, rho) * p_y
  uvw[Z].val[:] = uvw[Z].val[:] - dt / avg(uvw[Z].pos, rho) * p_z
  
  return  # end of function
