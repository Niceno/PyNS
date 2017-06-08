"""
Returns physical properties of air for a given domain (resolution).

Properties are take for 60 deg from:
  http://www.engineeringtoolbox.com/air-properties-d_156.html
"""

from numpy import ones

# =============================================================================
def properties_for_air(rc):
# -----------------------------------------------------------------------------
  """
  Args:
    rc: a tuple holding the resolution of the computational domain (number of
        cells in "x", "y" and "z" direction.

  Returns:
    rho:   three-dimensional matrix holding density for all cells
    mu:    three-dimensional matrix holding dymanic viscosity for all cells
    cp:    three-dimensional matrix holding thermal capacity for all cells
    kappa: three-dimensional matrix holding thermal conductivity for all cells

  Note:
    It is a no-brainer, but comes in handy.
  """

  # Create and fill matrice for all properties
  rho   = ones(rc) *    1.067     # density              [kg/m^3]
  mu    = ones(rc) *   20.17E-06  # viscosity            [Pa s]
  cp    = ones(rc) * 1009         # thermal capacity     [J/kg/K]
  kappa = ones(rc) *    0.0285    # thermal conductivity [W/m/K]

  return rho, mu, cp, kappa  # end of function
