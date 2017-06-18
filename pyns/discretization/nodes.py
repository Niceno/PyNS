"""
Creates uniform or non-uniform node coordinates in one dimension.

Uniform distribution is obtained if only three parameters are sent
(for starting and ending coordinate and the number of cells).

Non-uniform distribution is obtained if, in addition to the above,
two more parameters are sent (starting and ending cell size).
"""

# Standard Python modules
from pyns.standard import *

# PyNS modules
from pyns.constants import *
from pyns.operators import *

# =============================================================================
def nodes(*args):
# -----------------------------------------------------------------------------
    """
    Args: Number of input arguments can be three or five, depending if
          uniform or non-uniform node distribution is desired.

      Three arguments provided (for uniform meshes)
        args[0]: Starting coordinate.
        args[1]: Ending coordinate.
        args[2]: Number of cells (number of nodes minus one).

      Five arguments provided (for non-uniform meshes)
        args[0]: Same as above.
        args[1]: Same as above.
        args[2]: Same as above.
        args[3]: Starting cell size.
        args[4]: Ending cell size.

    Returns:
      An array with node coordinates (size is number of cells plus one).

    Examples:

      An example of the node distribution for a mesh with five cells
      is given here:

       starting                                                    ending
      coordinate                                                coordinate
            |                                                         |
      ----->|                                                     --->|

            0     1         2           3             4               5   nodes
            |--o--|----o----|-----o-----|------o------|-------o-------|
               0       1          2            3              4           cells
            |     |                                   |               |
            |<--->|                                   |<------------->|
        starting cell size                            ending cell size
      """

    # --------------------------------------------
    # Just a simple constant-spacing distribution
    # --------------------------------------------
    if len((args)) == 3:
        s = args[0]
        e = args[1]
        n = args[2]
        return linspace(s, e, n+1)

    # -----------------------------------------------------------------------
    # Grid distribution is determined from the following function:
    #   x = a + b y + c y^2 + d y^3
    #
    # One should imagine this function as taking integer arguments.
    #
    # Boundary conditions:
    #
    #   x(0)   = s         => a                                   = s
    #   x(1)   = s + del_s => a + b       + c          + d        = s + del_s
    #   x(n-1) = e - del_e => a + b*(n-1) + c*(n-1)^2 + d*(n-1)^3 = e - del_e
    #   x(n)   = e         => a + b*n     + c*n^2     + d*n^3     = e
    #
    # It follows that:
    #
    #   |1   0      0        0     |  |a|  |s      |
    #   |1   1      1        1     |  |b|  |s+del_s|
    #   |1   n-1   (n-1)^2  (n-1)^3|  |c|  |e-del_e|
    #   |1   n      n^2      n^3   |  |d|  |e      |
    # -----------------------------------------------------------------------
    elif len((args)) == 5:
        s     = args[0]
        e     = args[1]
        n     = args[2]
        del_s = args[3]
        del_e = args[4]

        # Form the system matrix
        A = matrix( [ [1,  0,    0,           0         ],     \
                      [1,  1,    1,           1         ],     \
                      [1,  n-1,  pow(n-1,2),  pow(n-1,3)],     \
                      [1,  n,    pow(n,  2),  pow(n,  3)] ] )

        # Form the right hand side
        f = array( [ s        ,     \
                     s + del_s,     \
                     e - del_e,     \
                     e         ] )

        # Solve for coefficients
        abcd = solve(A, f)

        # Pre-allocate the array for coordinates ...
        x = zeros(n+1);

        # ... and fill it up.
        for i in range(0, n+1):
            x[i] =   abcd[0]          \
                   + abcd[1] * i      \
                   + abcd[2] * i*i    \
                   + abcd[3] * i*i*i

        return x

    # -------------------------
    # Some error message might
    # be printed if number of
    # arguments is wrong
    # -------------------------

    return  # end of function
