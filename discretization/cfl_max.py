"""
Computes maximum Courant-Friedrich-Levy (CFD) number for the given velocity. 
"""

# ScriNS modules
from constants.all      import *
from operators.all      import *

# =============================================================================
def cfl_max(uvw, dt, dxyz):
# -----------------------------------------------------------------------------
  """
  Args:
    uvw: a tuple with three staggered or centered velocity components 
         (each component is created with "scrins.create_unknown" function.
    dt:  time step
    dxyz: a tuple holding cell dimensions in "x", "y" and "z" directions.
          Each cell dimension is a three-dimensional array.

  Returns:
    cfl: a floating point number holding the maximum value of CFL.

  Note:
    It could be written in a more compact way, without unpacking the tuples,
    but that would only lead to poorer readability of the code.    
  """

  # Unpack received tuples
  u,  v,  w  = uvw
  dx, dy, dz = dxyz

  # Take velocity's position
  d = u.pos
  
  # Mesh is cell-centered
  if d == C:
    cfl = dt * max( abs(u.val/dx).max(),   \
                    abs(v.val/dy).max(),   \
                    abs(w.val/dz).max() )
  
  # Mesh is staggered
  else:
    cfl = dt * max( abs(u.val/avg(X,dx)).max(),   \
                    abs(v.val/avg(Y,dy)).max(),   \
                    abs(w.val/avg(Z,dz)).max() )
  
  return cfl