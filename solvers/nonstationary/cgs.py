"""
Preconditioned Conjugate Gradient Squared (CGS) solver.

Source:
  http://www.netlib.org/templates/templates.pdf
"""

# Standard Python modules
from pyns.standard import *

# PyNS modules
from pyns.constants      import *
from pyns.display        import write
from pyns.discretization import Unknown

# Sisters from this module
from pyns.solvers.mat_vec_bnd import mat_vec_bnd
from pyns.solvers.vec_vec     import vec_vec
from pyns.solvers.norm        import norm

# =============================================================================
def cgs(a, phi, b, tol, 
        verbatim = False,
        max_iter = -1):
# -----------------------------------------------------------------------------
    """
    Args:
      a:        Object of the type "Matrix", holding the system matrix.
      phi:      Object of the type "Unknown" to be solved.
      b:        Three-dimensional array holding the source term.
      tol:      Absolute solver tolerance
      verbatim: Logical variable setting if solver will be verbatim (print
                info on Python console) or not.
      max_iter: Maxiumum number of iterations.
      
    Returns:
      x: Three-dimensional array with solution.
    """

    if verbatim:
        write.at(__name__)

    # Helping variable
    x = phi.val

    # Intitilize arrays
    p       = zeros(x.shape)
    p_hat   = Unknown("vec_p_hat", phi.pos, x.shape, -1, per=phi.per, 
                      verbatim=False)
    q       = zeros(x.shape)
    r       = zeros(x.shape)
    r_tilda = zeros(x.shape)
    u       = zeros(x.shape)
    u_hat   = Unknown("vec_u_hat", phi.pos, x.shape, -1, per=phi.per, 
                      verbatim=False)
    v_hat   = zeros(x.shape)

    # r = b - A * x
    r[:,:,:] = b[:,:,:] - mat_vec_bnd(a, phi)

    # Chose r~
    r_tilda[:,:,:] = r[:,:,:]

    # ---------------
    # Iteration loop
    # ---------------
    if max_iter == -1:
        max_iter = prod(phi.val.shape)
        
    for i in range(1, max_iter):

        if verbatim:
            print("  iteration: %3d:" % (i), end = "" )

        # rho = r~ * r
        rho = vec_vec(r_tilda, r)

        # If rho == 0 method fails
        if abs(rho) < TINY * TINY:
            if verbatim == True:  
                write.at(__name__)
                print("  Fails becuase rho = %12.5e" % rho)
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
        p_hat.val[:,:,:] = p[:,:,:] / a.C[:,:,:]

        # v^ = A * p^
        v_hat[:,:,:] = mat_vec_bnd(a, p_hat)

        # alfa = rho / (r~ * v^)
        alfa = rho / vec_vec(r_tilda, v_hat)

        # q = u - alfa v^
        q[:,:,:] = u[:,:,:] - alfa * v_hat[:,:,:]

        # Solve M u^ = u + q
        u_hat.val[:,:,:] = (u[:,:,:] + q[:,:,:]) / a.C[:,:,:]

        # x = x + alfa u^
        x[:,:,:] += alfa * u_hat.val[:,:,:]

        # q^ = A u^
        q_hat = mat_vec_bnd(a, u_hat)

        # r = r - alfa q^
        r[:,:,:] -= alfa * q_hat[:,:,:]

        # Compute residual
        res = norm(r)

        if verbatim:
            print("%12.5e" %res)

        # If tolerance has been reached, get out of here
        if res < tol:
            return x

        # Prepare for next iteration
        rho_old = rho

    return x  # end of function
