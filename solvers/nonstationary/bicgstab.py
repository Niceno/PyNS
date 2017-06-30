"""
Preconditioned Bi-Conjugate Gradient Stabilized (BiCGStab) solver.

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
def bicgstab(a, phi, b, tol, 
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
    r       = zeros(x.shape)
    r_tilda = zeros(x.shape)
    s       = zeros(x.shape)
    s_hat   = Unknown("vec_s_hat", phi.pos, x.shape, -1, per=phi.per, 
                      verbatim=False)
    v       = zeros(x.shape)

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
            write.at(__name__)
            print("  Fails becuase rho = %12.5e" % rho)
            return x

        if i == 1:
            # p = r
            p[:,:,:] = r[:,:,:]

        else:
            # beta = (rho / rho_old)(alfa/omega)
            beta = rho / rho_old * alfa / omega

            # p = r + beta (p - omega v)
            p[:,:,:] = r[:,:,:] + beta * (p[:,:,:] - omega * v[:,:,:])

        # Solve M p_hat = p
        p_hat.val[:,:,:] = p[:,:,:] / a.C[:,:,:]

        # v = A * p^
        v[:,:,:] = mat_vec_bnd(a, p_hat)

        # alfa = rho / (r~ * v)
        alfa = rho / vec_vec(r_tilda, v)

        # s = r - alfa v
        s[:,:,:] = r[:,:,:] - alfa * v[:,:,:]

        # Check norm of s, if small enough set x = x + alfa p_hat and stop
        res = norm(s)
        if res < tol:
            if verbatim == True:  
                write.at(__name__)
                print("  Fails becuase rho = %12.5e" % rho)
            x[:,:,:] += alfa * p_hat.val[:,:,:]
            return x

        # Solve M s^ = s
        s_hat.val[:,:,:] = s[:,:,:] / a.C[:,:,:]

        # t = A s^
        t = mat_vec_bnd(a, s_hat)

        # omega = (t * s) / (t * t)
        omega = vec_vec(t, s) / vec_vec(t, t)

        # x = x + alfa p^ + omega * s^
        x[:,:,:] += alfa * p_hat.val[:,:,:] + omega * s_hat.val[:,:,:]

        # r = s - omega q^
        r[:,:,:] = s[:,:,:] - omega * t[:,:,:]

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
