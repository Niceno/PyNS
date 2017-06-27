"""
This is ruthless violation of Python's ethical code of conduct (best practice).
It includes all the solvers defined in the code.  Very bad, but it makes
code development faster for non-IT oriendted minds.
"""

# Classes in "discretization":
from .Matrix import Matrix

# Functions in "discretization":
from .bicgstab     import bicgstab
from .cg           import cg
from .cgs          import cgs
from .jacobi       import jacobi
from .coarsen      import coarsen