"""
Returns physical properties of air for a given domain (resolution).

Properties are take for 60 deg from:
  http://www.engineeringtoolbox.com/air-properties-d_156.html
"""

# Standard Python modules
from pyns.standard import *

# PyNS modules
from pyns.constants.all      import *
from pyns.operators.all      import *

# =============================================================================
def properties_for_water(rc):
# -----------------------------------------------------------------------------
    """
    Args:
      rc: Tuple holding the resolution of the computational domain
          (number of cells in "x", "y" and "z" direction).

    Returns:
      rho:   Three-dimensional matrix holding density for cells.
      mu:    Three-dimensional matrix holding dymanic viscosity for cells.
      cp:    Three-dimensional matrix holding thermal capacity for cells.
      kappa: Three-dimensional matrix holding thermal conductivity for cells.

    Note:
      It is a no-brainer, but comes in handy.
    """

    # Create and fill matrices for all properties
    rho   = ones(rc) *  983.0       # density              [kg/m^3]
    mu    = ones(rc) *    0.466E-3  # viscosity            [Pa s]
    cp    = ones(rc) * 4185         # thermal capacity     [J/kg/K]
    kappa = ones(rc) *    0.654     # thermal conductivity [W/m/K]

    return rho, mu, cp, kappa  # end of function
