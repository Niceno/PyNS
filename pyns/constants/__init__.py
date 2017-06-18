"""
This is ruthless violation of Python's ethical code of conduct (best practice).
It includes all the constants defined in the code.  Very bad, but it makes
code development faster for non-IT oriendted minds.
"""

# Constants
from .boundary_conditions    import DIRICHLET, NEUMANN, OUTLET
from .compass                import W, E, S, N, B, T, C
from .coordinates            import X, Y, Z
from .solver                 import TOL
from .tiny_and_huge          import TINY, HUGE
