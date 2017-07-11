# PyNS modules
from pyns.constants      import *

# =============================================================================
def coarsable(shape):
# -----------------------------------------------------------------------------

    if ( (shape[X] > 2 and shape[X] % 2 == 0) or
         (shape[Y] > 2 and shape[Y] % 2 == 0) or
         (shape[Z] > 2 and shape[Z] % 2 == 0) ):
        return True
    
    return False  # end of function