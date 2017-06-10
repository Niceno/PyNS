"""
Preconditioned Conjugate Gradient (CG) solver.

Source:
  http://www.netlib.org/templates/templates.pdf
"""

# Standard Python modules
from pyns.standard import *

from pyns.solvers.mat_vec import mat_vec
from pyns.solvers.vec_vec import vec_vec
from pyns.solvers.norm    import norm

# =============================================================================
def cg(a, phi, b, tol, ver):
# -----------------------------------------------------------------------------
    """
    Args:
      a:   System matrix in PyNS format (which ssentially means storing a
           bundle of non-zero diagonals in compas directions)
      phi: Unknown to be solved (from "create_unknown" function)
      b:   Three-dimensional matrix holding the source term.
      tol: Absolute solver tolerance
      ver: Logical variable setting if solver will be verbatim (print
           info on Python console) or not.

    Returns:
      x: Three-dimensional matrix with solution.

    Note:
      One should also consider implementing periodic boundary conditions
      in this version of the solver.
    """

    if ver:
        print("Solver: CG")

    # Helping variables
    x = phi.val
    n = prod(x.shape)

    # Intitilize arrays
    p = zeros(x.shape)
    q = zeros(x.shape)
    r = zeros(x.shape)
    z = zeros(x.shape)

    # r = b - A * x
    r[:,:,:] = b[:,:,:] - mat_vec(a, x)

    # ---------------
    # Iteration loop
    # ---------------
    for i in range(1,n):

        if ver:
            print("  iteration: %3d:" % (i), end = "" )

        # Solve M z = r
        z[:,:,:] = r[:,:,:] / a.P[:,:,:]

        # rho = r * z
        rho = vec_vec(r, z)

        if i == 1:
            # p = z
            p[:,:,:] = z[:,:,:]

        else:
            # beta = rho / rho_old
            beta = rho / rho_old

            # p = z + beta p
            p[:,:,:] = z[:,:,:] + beta * p[:,:,:]

        # q = A * p
        q[  :,  :,  :]  = mat_vec(a, p)

        # alfa = rho / (p * q)
        alfa = rho / vec_vec(p, q)

        # x = x + alfa p
        x[:,:,:] += alfa * p[:,:,:]

        # r = r - alfa q
        r[:,:,:] -= alfa * q[:,:,:]

        # Compute residual
        res = norm(r)

        if ver:
            print("%12.5e" %res)

        # If tolerance has been reached, get out of here
        if res < tol:
            return x

        # Prepare for next iteration
        rho_old = rho

    return x  # end of function
