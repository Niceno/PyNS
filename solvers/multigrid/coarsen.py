# PyNS modules
from pyns.constants      import *

# =============================================================================
def coarsen(shape):
# -----------------------------------------------------------------------------

    nx, ny, nz = shape

    if (shape[X] > 2 and shape[X] % 2 == 0):
        nx = nx // 2
        
    if (shape[Y] > 2 and shape[Y] % 2 == 0):
        ny = ny // 2
        
    if (shape[Z] > 2 and shape[Z] % 2 == 0):
        nz = nz // 2
      
    return (nx, ny, nz)  # end of function