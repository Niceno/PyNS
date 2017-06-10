"""
Preconditioned Bi-Conjugate Gradient Stabilized (BiCGStab) solver.

Source:
  http://www.netlib.org/templates/templates.pdf
"""

# Standard Python modules
from pyns.standard import *

# PyNS modules
from pyns.constants.all   import TINY
from pyns.solvers.mat_vec import mat_vec
from pyns.solvers.vec_vec import vec_vec
from pyns.solvers.norm    import norm

# =============================================================================
def bicgstab(a, phi, b, tol, ver):
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
    print("Solver: BiCGStab")

  # Helping variables
  x = phi.val
  n = prod(x.shape)
  
  # Intitilize arrays
  p       = zeros(x.shape)
  p_hat   = zeros(x.shape)
  r       = zeros(x.shape)
  r_tilda = zeros(x.shape)
  s       = zeros(x.shape)
  s_hat   = zeros(x.shape)
  v       = zeros(x.shape)
  
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
      print("bicgstab fails becuase rho = %12.5e" % rho)
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
    p_hat[:,:,:] = p[:,:,:] / a.P[:,:,:]         
    
    # v = A * p^
    v[:,:,:] = mat_vec(a, p_hat)

    # alfa = rho / (r~ * v)
    alfa = rho / vec_vec(r_tilda, v)

    # s = r - alfa v
    s[:,:,:] = r[:,:,:] - alfa * v[:,:,:]

    # Check norm of s, if small enough set x = x + alfa p_hat and stop
    res = norm(s)
    if res < tol:
      if ver:  
        print("%12.5e" %res)
      x[:,:,:] += alfa * p_hat[:,:,:]
      return x

    # Solve M s^ = s
    s_hat[:,:,:] = s[:,:,:] / a.P[:,:,:]         

    # t = A s^
    t = mat_vec(a, s_hat)
    
    # omega = (t * s) / (t * t)
    omega = vec_vec(t, s) / vec_vec(t, t)
    
    # x = x + alfa p^ + omega * s^
    x[:,:,:] += alfa * p_hat[:,:,:] + omega * s_hat[:,:,:]
      
    # r = s - omega q^
    r[:,:,:] = s[:,:,:] - omega * t[:,:,:]
      
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
