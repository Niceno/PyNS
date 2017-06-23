"""
This is ruthless violation of Python's ethical code of conduct (best practice).
It includes all the functionality defined in the module "discretization".
Very bad, but it makes code development faster for non-IT oriendted minds.
"""

# Classes in "discretization":
from .Unknown import Unknown
  
# Functions in "discretization":
from .adj_n_bnds      import adj_n_bnds
from .calc_p          import calc_p
from .calc_phi        import calc_phi
from .calc_t          import calc_t
from .calc_uvw        import calc_uvw
from .cartesian_grid  import cartesian_grid
from .cfl_max         import cfl_max
from .corr_uvw        import corr_uvw
from .diffusion       import diffusion
from .nodes           import nodes
from .nodal_uvw       import nodal_uvw
from .vol_balance     import vol_balance
