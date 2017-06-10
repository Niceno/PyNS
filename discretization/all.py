"""
This is ruthless violation of Python's ethical code of conduct (best practice).  
It includes all the functionality defined in the module "discretization".  
Very bad, but it makes code development faster for non-IT oriendted minds.
"""

# Discretization
from pyns.discretization.adj_n_bnds      import adj_n_bnds
from pyns.discretization.adj_o_bnds      import adj_o_bnds
from pyns.discretization.advection       import advection
from pyns.discretization.calc_p          import calc_p
from pyns.discretization.calc_t          import calc_t
from pyns.discretization.calc_uvw        import calc_uvw
from pyns.discretization.cartesian_grid  import cartesian_grid
from pyns.discretization.cfl_max         import cfl_max
from pyns.discretization.corr_uvw        import corr_uvw
from pyns.discretization.create_matrix   import create_matrix
from pyns.discretization.create_unknown  import create_unknown
from pyns.discretization.nodes           import nodes
from pyns.discretization.obst_mod_matrix import obst_mod_matrix
from pyns.discretization.obst_zero_val   import obst_zero_val
from pyns.discretization.vol_balance     import vol_balance
