"""
Function to generate a parabolic profile over a set of cell centers,
useful for specifying parabolic inlet velocity profiles.  

It expects nodal coordinates as input, but sends values at cell centers back:

Input coordinates at:    |-----|-----|-----|-----|-----|-----|-----|  

Output values at:           o-----o-----o-----o-----o-----o-----o
"""

# Standard Python modules
from pyns.standard import *

# PyNS modules
from pyns.constants.all import *
from pyns.display.all   import *

from pyns.operators.all import avg

# =============================================================================
def par(mean_val, x_nodes):
# -----------------------------------------------------------------------------
  """
  Args:
    mean_val: Mean value the parabola will have.
    x_nodes:  Nodal coordinatels.
    
  Returns:
    An array with values of parabolic function.
    
  Note:
    It would really be useful to have a two-dimensional variant of this.
  """

  # It is known that maximum of a parabola is 3/2 of its mean value
  max_val = mean_val * 3/2

  # Normalized x coordinates (from -1 to +1)
  xn = copy(x_nodes)
  
  xn -= xn.min()                  
  xn /= (xn.max()-xn.min())
  xn *= 2                               
  xn -= 1  
  
  xc = avg(xn)

  yc = (1.0-xc*xc) * max_val

  return yc  # end of function
