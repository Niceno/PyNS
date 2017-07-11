"""
Base for the geometrically agglomerated multigrid method.
"""

# Standard Python modules
from pyns.standard import *

# PyNS modules
from pyns.constants      import *
from pyns.operators      import *
from pyns.display        import write
from pyns.discretization import Unknown

# Modules from the parent's directory
from pyns.solvers.mat_vec_bnd import mat_vec_bnd
from pyns.solvers.norm        import norm

# Sisters from this module
from pyns.solvers.stationary import jacobi
from pyns.solvers.multigrid.gamg_coarsen_system import gamg_coarsen_system

# =============================================================================
def gamg_v_cycle(a, phi, b, tol, 
                 verbatim = False):
# -----------------------------------------------------------------------------
    """
    Args:
      a:        Object of the type "Matrix", holding the system matrix.
      phi:      Object of the type "Unknown" to be solved.
      b:        Three-dimensional array holding the source term.
      tol:      Absolute solver tolerance
      verbatim: Logical variable setting if solver will be verbatim (print
                info on Python console) or not.

    Returns:
      x: Three-dimensional array with solution.
    """

    if verbatim:
        write.at(__name__)

    # ------------------
    #    
    # Get system levels
    #    
    # ------------------
    n_, shape_, a_, i_, d_, phi_, b_, r_ = gamg_coarsen_system(a, phi, b)

    # Set a couple of parameters for the solver
    n_steps  = n_-1  # number of steps down the "V" cycle      
    n_smooth = 4     # number of smoothing sweeps   

    # -----------------------------------------------------
    #
    # Solve a bit on the finest level and compute residual
    #
    # -----------------------------------------------------
    grid = 0  # finest level
    phi_[grid].val = jacobi(a_[grid], phi_[grid], b_[grid], TOL, 
                         verbatim = True, 
                         max_iter = 4)
    # r = b - A * x
    r_[grid].val[:,:,:] = b_[grid][:,:,:] - mat_vec_bnd(a_[grid], phi_[grid])
    if verbatim:
        print("  residual at level %d" % grid, norm(r_[grid].val))

    # =========================================================================
    #
    # Start the V-cycle
    #
    # =========================================================================
    for cycle in range(0, 9):    

        # --------------------
        #
        # Go down a few steps
        #
        # --------------------        
        for level in range(1, n_steps+1):
            grid = grid + 1    

            # Compute ratio between grid levels 
            rx = shape_[grid-1][X] // shape_[grid][X]
            ry = shape_[grid-1][Y] // shape_[grid][Y]
            rz = shape_[grid-1][Z] // shape_[grid][Z]
            
            # Lower bounds for browsing through grid levels
            wl = 0
            el = 1
            sl = 0
            nl = 1
            bl = 0
            tl = 1
            
            # -------------------------------------
            # Restrict 
            #
            # Computes r.h.s. on the coarser level 
            # from residual on the finer level.
            # -------------------------------------
            b_[grid][:,:,:] =  r_[grid-1].val[wl::rx, sl::ry, bl::rz]
            b_[grid][:,:,:] += r_[grid-1].val[wl::rx, nl::ry, bl::rz]
            b_[grid][:,:,:] += r_[grid-1].val[wl::rx, sl::ry, tl::rz]
            b_[grid][:,:,:] += r_[grid-1].val[wl::rx, nl::ry, tl::rz]
            b_[grid][:,:,:] += r_[grid-1].val[el::rx, sl::ry, bl::rz]
            b_[grid][:,:,:] += r_[grid-1].val[el::rx, nl::ry, bl::rz]
            b_[grid][:,:,:] += r_[grid-1].val[el::rx, sl::ry, tl::rz]
            b_[grid][:,:,:] += r_[grid-1].val[el::rx, nl::ry, tl::rz]
            
            # ------------------------------------------------
            # Solve on the coarser level and compute residual
            # ------------------------------------------------
            phi_[grid].val[:] = 0  # nulify to forget previous corrections
            phi_[grid].val = jacobi(a_[grid], phi_[grid], b_[grid], TOL, 
                                 verbatim = False, 
                                 max_iter = level * n_smooth)
            # r = b - A * x
            r_[grid].val[:,:,:] = b_[grid][:,:,:]  \
                                - mat_vec_bnd(a_[grid], phi_[grid])
            if verbatim:
                print("  residual at level %d" % grid, norm(r_[grid].val))
    
        # ==================================
        #
        # This is the bottom of the V-cycle
        #
        # ==================================
    
        # ------------------
        #
        # Go up a few steps
        #
        # ------------------
        for level in range(1, n_steps+1):
            grid = grid - 1

            # Compute ratio between grid levels 
            rx = shape_[grid][X] // shape_[grid+1][X]
            ry = shape_[grid][Y] // shape_[grid+1][Y]
            rz = shape_[grid][Z] // shape_[grid+1][Z]
            
            # Lower bounds for browsing through grid levels
            wl = 0
            el = 1
            sl = 0
            nl = 1
            bl = 0
            tl = 1
                
            # -------------------------------------------
            # Prolongation 
            #
            # Interpolates the solution from the coarser 
            # level as the correction to the current. 
            # It also smoooths it out a little bit.
            # -------------------------------------------
            r_[grid].val[:] = 0
              
            # First copy in each available cell on fine level  
            r_[grid].val[wl::rx, sl::ry, bl::rz] = phi_[grid+1].val[:,:,:]
            
            # Then spread arond
            r_[grid].val[el::rx, sl::ry, bl::rz] = r_[grid].val[wl::rx, sl::ry, bl::rz]
            r_[grid].val[wl::rx, nl::ry, bl::rz] = r_[grid].val[wl::rx, sl::ry, bl::rz]
            r_[grid].val[el::rx, nl::ry, bl::rz] = r_[grid].val[wl::rx, sl::ry, bl::rz]

            r_[grid].val[wl::rx, sl::ry, tl::rz] = r_[grid].val[wl::rx, sl::ry, bl::rz]
            r_[grid].val[el::rx, sl::ry, tl::rz] = r_[grid].val[wl::rx, sl::ry, bl::rz]
            r_[grid].val[wl::rx, nl::ry, tl::rz] = r_[grid].val[wl::rx, sl::ry, bl::rz]
            r_[grid].val[el::rx, nl::ry, tl::rz] = r_[grid].val[wl::rx, sl::ry, bl::rz]

            # Then smooth them out a little bit
            for smooth in range(0, n_smooth):
                r_[grid].exchange()
                summ = zeros((shape_[grid]))
                summ[:,:,:] += cat_x((r_[grid].bnd[W].val, 
                                      r_[grid].val[:-1,:,:])) * a_[grid].W
                summ[:,:,:] += cat_x((r_[grid].val[ 1:,:,:], 
                                      r_[grid].bnd[E].val  )) * a_[grid].E
                summ[:,:,:] += cat_y((r_[grid].bnd[S].val, 
                                      r_[grid].val[:,:-1,:])) * a_[grid].S
                summ[:,:,:] += cat_y((r_[grid].val[:, 1:,:], 
                                      r_[grid].bnd[N].val  )) * a_[grid].N
                summ[:,:,:] += cat_z((r_[grid].bnd[B].val, 
                                      r_[grid].val[:,:,:-1])) * a_[grid].B
                summ[:,:,:] += cat_z((r_[grid].val[:,:, 1:], 
                                      r_[grid].bnd[T].val  )) * a_[grid].T
                r_[grid].val[:] = summ[:] / d_[grid][:]
        
            # -----------------------------------------------
            # Correction on the finer level, followed by a 
            # bit of smoothing and computation of residuals.
            # -----------------------------------------------
            phi_[grid].val[:] += r_[grid].val[:]
            phi_[grid].val = jacobi(a_[grid], phi_[grid], b_[grid], TOL, 
                                 verbatim = False, 
                                 max_iter = n_smooth)
            # r = b - A * x
            r_[grid].val[:,:,:] = b_[grid][:,:,:]  \
                                - mat_vec_bnd(a_[grid], phi_[grid])
            print("  residual at level %d" % grid, norm(r_[grid].val))

    return phi_[0].val  # end of function
