"""
Preconditioned Conjugate Gradient Squared (CGS) solver.

Source:
  http://www.netlib.org/templates/templates.pdf
"""

# Standard Python modules
from pyns.standard import *

from pyns.constants       import TINY
from pyns.solvers.mat_vec import mat_vec
from pyns.solvers.vec_vec import vec_vec
from pyns.solvers.norm    import norm

# =============================================================================
def cgs(a, phi, b, tol, ver):
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
        print("Solver: CGS")

    # Helping variables
    x = phi.val
    n = prod(x.shape)

    # Intitilize arrays
    p       = zeros(x.shape)
    p_hat   = zeros(x.shape)
    q       = zeros(x.shape)
    r       = zeros(x.shape)
    r_tilda = zeros(x.shape)
    u       = zeros(x.shape)
    u_hat   = zeros(x.shape)
    v_hat   = zeros(x.shape)

    # r = b - A * x
    r[:,:,:] = b[:,:,:] - mat_vec(a, x)

    # Chose r~
    r_tilda[:,:,:] = r[:,:,:]

    # ---------------
    # Iteration loop
    # ---------------
    for i in range(1,n):

        if ver:
            print("  iteration: %3d:" % (i), end = "" )

        # rho = r~ * r
        rho = vec_vec(r_tilda, r)

        # If rho == 0 method fails
        if abs(rho) < TINY * TINY:
            print("cgs fails becuase rho = %12.5e" % rho)
            return x

        if i == 1:
            # u = r
            u[:,:,:] = r[:,:,:]

            # p = u
            p[:,:,:] = u[:,:,:]

        else:
            # beta = rho / rho_old
            beta = rho / rho_old

            # u = r + beta q
            u[:,:,:] = r[:,:,:] + beta * q[:,:,:]

            # p = u + beta (q + beta p)
            p[:,:,:] = u[:,:,:] + beta * (q[:,:,:] + beta * p[:,:,:])

        # Solve M p_hat = p
        p_hat[:,:,:] = p[:,:,:] / a.P[:,:,:]

        # v^ = A * p^
        v_hat[:,:,:] = mat_vec(a, p_hat)

        # alfa = rho / (r~ * v^)
        alfa = rho / vec_vec(r_tilda, v_hat)

        # q = u - alfa v^
        q[:,:,:] = u[:,:,:] - alfa * v_hat[:,:,:]

        # Solve M u^ = u + q
        u_hat[:,:,:] = (u[:,:,:] + q[:,:,:]) / a.P[:,:,:]

        # x = x + alfa u^
        x[:,:,:] += alfa * u_hat[:,:,:]

        # q^ = A u^
        q_hat = mat_vec(a, u_hat)

        # r = r - alfa q^
        r[:,:,:] -= alfa * q_hat[:,:,:]

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
