"""
Spawns a 3D Cartesian grid from three arrays with node coordinates.
"""

# Standard Python modules
from standard import *

# PyNS modules
from constants.all      import *
from operators.all      import *

# =============================================================================
def cartesian_grid(xn, yn, zn):
# -----------------------------------------------------------------------------
  """
  Args:
    xn: One-dimensional array with node coordinates in "x" direction.
    yn: One-dimensional array with node coordinates in "y" direction.
    zn: One-dimensional array with node coordinates in "z" direction.
    
  Returns:
    nx: Number of cells (one less than nodes) in "x" direction.
    ny: Number of cells (one less than nodes) in "y" direction.
    nz: Number of cells (one less than nodes) in "z" direction.
    dx: Three-dimensional matrix holding "dx" for all cells.
    dy: Three-dimensional matrix holding "dy" for all cells.
    dz: Three-dimensional matrix holding "dz" for all cells.
    rc: Resolution tuple (nx, ny and nz) for centered variable.
    rx: Resolution tuple (nx, ny and nz) for variable staggered in "x".
    ry: Resolution tuple (nx, ny and nz) for variable staggered in "y".
    rz: Resolution tuple (nx, ny and nz) for variable staggered in "z".
    
  Note:
    Cell resolutions are for one smaller than node resolutions for each 
    direction.  Staggered cell resolutions are smaller than centered.
  """

  # Compute cell resolutions
  nx = len(xn)-1
  ny = len(yn)-1
  nz = len(zn)-1
  
  # Create matrices for cell dimensions ...
  dx = empty((nx,ny,nz))
  dy = empty((nx,ny,nz))
  dz = empty((nx,ny,nz))
  
  # ... and fill them up!
  dx[:] = dif(xn).reshape(nx,1,1)
  dy[:] = dif(yn).reshape(1,ny,1)
  dz[:] = dif(zn).reshape(1,1,nz)
  
  # Compute resolutions for cell-centered and all collocated variables
  rc = nx,   ny,   nz
  ru = nx-1, ny,   nz
  rv = nx,   ny-1, nz
  rw = nx,   ny,   nz-1
  
  return nx,ny,nz, dx,dy,dz, rc,ru,rv,rw  # end of function
