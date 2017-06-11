"""
Returns physical properties of air for a given domain (resolution).

Properties are take for 60 deg from:
  http://www.engineeringtoolbox.com/air-properties-d_156.html
"""

# Standard Python modules
from pyns.standard import *

# =============================================================================
def air(rc):
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

    # Create and fill matrice for all properties
    rho   = ones(rc) *    1.067     # density              [kg/m^3]
    mu    = ones(rc) *   20.17E-06  # viscosity            [Pa s]
    cp    = ones(rc) * 1009         # thermal capacity     [J/kg/K]
    kappa = ones(rc) *    0.0285    # thermal conductivity [W/m/K]

    return rho, mu, cp, kappa  # end of function
