"""
Weighted Jacobi Method for solution of linear systems of equations.

Source:
  https://en.wikipedia.org/wiki/Jacobi_method
  
Note:
  It is not very fast but will serve as the basis for multigrid solver.
"""

from __future__ import print_function

# Standard Python modules
from pyns.standard import *

# PyNS modules
from pyns.constants      import *
from pyns.display        import write
from pyns.discretization import Unknown

# Modules from the parent's directory
from pyns.solvers.mat_vec_bnd import mat_vec_bnd
from pyns.solvers.norm        import norm

# =============================================================================
def jacobi(a, phi, b, tol, 
           verbose = False,
           max_iter = -1):
# -----------------------------------------------------------------------------
    """
    Args:
      a: ...... Object of the type "Matrix", holding the system matrix.
      phi: .... Object of the type "Unknown" to be solved.
      b: ...... Three-dimensional array holding the source term.
      tol: .... Absolute solver tolerance
      verbose:  Logical variable setting if solver will be verbose (print
                info on Python console) or not.
      max_iter: Maxiumum number of iterations.

    Returns:
      phi.val: Three-dimensional array with solution.
    """

    if verbose is True:
        write.at(__name__)

    sum = zeros(phi.val.shape)
    r   = zeros(phi.val.shape)

    if max_iter == -1:
        max_iter = prod(phi.val.shape)

    # Under-relaxation factor
    alfa = 0.9

    # Main iteration loop
    for iter in range(0, max_iter):

        if verbose is True:
            print("  iteration: %3d:" % (iter), end = "" )
      
        # Add source term 
        sum[:] = b[:]
        
        # Add contribution from west and east
        sum[ :1,:,:] += a.W[ :1,:,:] * phi.bnd[W].val[:1,:,:]
        sum[-1:,:,:] += a.E[-1:,:,:] * phi.bnd[E].val[:1,:,:]

        # Add up west and east neighbors from the inside
        sum[ 1:,:,:] += a.W[ 1:,:,:] * phi.val[:-1,:,:]
        sum[:-1,:,:] += a.E[:-1,:,:] * phi.val[ 1:,:,:]
        
        # Add contribution from south and north
        sum[:, :1,:] += a.S[:, :1,:] * phi.bnd[S].val[:,:1,:]
        sum[:,-1:,:] += a.N[:,-1:,:] * phi.bnd[N].val[:,:1,:]

        # Add up south and north neighbors from the inside
        sum[:, 1:,:] += a.S[:, 1:,:] * phi.val[:,:-1,:]
        sum[:,:-1,:] += a.N[:,:-1,:] * phi.val[:, 1:,:]
        
        # Add contribution from bottom and top
        sum[:,:, :1] += a.B[:,:, :1] * phi.bnd[B].val[:,:,:1]
        sum[:,:,-1:] += a.T[:,:,-1:] * phi.bnd[T].val[:,:,:1]

        # Add up south and north neighbors from the inside
        sum[:,:, 1:] += a.B[:,:, 1:] * phi.val[:,:,:-1]
        sum[:,:,:-1] += a.T[:,:,:-1] * phi.val[:,:, 1:]
        
        phi.val[:,:,:] = (1.0 - alfa) * phi.val[:,:,:]   \
                       +        alfa  * sum[:,:,:] / a.C[:,:,:]
        phi.exchange()

        # r = b - A * x
        r[:,:,:] = b[:,:,:] - mat_vec_bnd(a, phi)

        # Compute residual
        res = norm(r)

        if verbose is True:
            print("%12.5e" %res)

        # If tolerance has been reached, get out of here
        if res < tol:
            return phi.val

    return phi.val  # end of function
