"""
Interpolates velocities to nodes.  Important for Lagrangian particle tracking.
"""

# Standard Python modules
from pyns.standard import *

# PyNS modules
from pyns.constants      import *
from pyns.operators      import *

# =============================================================================
def nodal_uvw(xyzn, uvwf, obst):
# -----------------------------------------------------------------------------

    # Unpack tuples
    xn, yn, zn = xyzn
    uf, vf, wf = uvwf

    nx = xn.shape[0] 
    ny = yn.shape[0] 
    nz = zn.shape[0] 

    # Reserve space for nodal values of velocities
    un = zeros( (nx, ny, nz) )
    vn = zeros( (nx, ny, nz) )
    wn = zeros( (nx, ny, nz) )

    # ---------------------------------------------
    # Nodal values of velocities inside the domain
    # ---------------------------------------------
    
    # "u" velocity component
    un[:, 1:-1, 1:-1] = avg_z(avg_y(cat_x((uf.bnd[W].val[:], 
                                           uf.val[:], 
                                           uf.bnd[E].val[:])))) 

    # "v" velocity component
    vn[1:-1, :, 1:-1] = avg_z(avg_x(cat_y((vf.bnd[S].val[:], 
                                           vf.val[:], 
                                           vf.bnd[N].val[:])))) 

    # "w" velocity component
    wn[1:-1, 1:-1, :] = avg_y(avg_x(cat_z((wf.bnd[B].val[:], 
                                           wf.val[:], 
                                           wf.bnd[T].val[:])))) 

    # ---------------------------------------------
    # Nodal values of velocities on the boundaries
    # ---------------------------------------------
    
    # Nodal values for "u" velocity components on S, N, B and T
    un[1:-1,  :1, 1:-1] = avg_z(uf.bnd[S].val)
    un[1:-1, -1:, 1:-1] = avg_z(uf.bnd[N].val)
    un[1:-1, 1:-1,  :1] = avg_y(uf.bnd[B].val)
    un[1:-1, 1:-1, -1:] = avg_y(uf.bnd[T].val)

    # Nodal values for "v" velocity components on W, E, B and T
    vn[ :1, 1:-1, 1:-1] = avg_z(vf.bnd[W].val)
    vn[-1:, 1:-1, 1:-1] = avg_z(vf.bnd[E].val)
    vn[1:-1, 1:-1,  :1] = avg_x(vf.bnd[B].val)
    vn[1:-1, 1:-1, -1:] = avg_x(vf.bnd[T].val)

    # Nodal values for "v" velocity components on W, E, S and N
    wn[ :1, 1:-1, 1:-1] = avg_y(wf.bnd[W].val)
    wn[-1:, 1:-1, 1:-1] = avg_y(wf.bnd[E].val)
    wn[1:-1,  :1, 1:-1] = avg_x(wf.bnd[S].val)
    wn[1:-1, -1:, 1:-1] = avg_x(wf.bnd[N].val)

    # -------------------------------------------
    # Take care of edges 
    # This is complicated because I don't have  
    # separate boundary condition type for wall.
    # -------------------------------------------
    for comp in (un, vn, wn):

        # Edges in "x" direction 
        m_0 = abs(comp[ :, 1: 2, :1]) > TINY
        m_1 = abs(comp[ :, :1, 1: 2]) > TINY
        comp[ :, :1, :1] = (   comp[ :, 1: 2, :1] 
                             + comp[ :, :1, 1: 2]) * 0.5 * m_0 * m_1
            
        m_0 = abs(comp[ :, 1: 2,-1:]) > TINY
        m_1 = abs(comp[ :, :1,-2:-1]) > TINY              
        comp[ :, :1,-1:] = (   comp[ :, 1: 2,-1:] 
                             + comp[ :, :1,-2:-1]) * 0.5 * m_0 * m_1

        m_0 = abs(comp[ :,-2:-1, :1]) > TINY
        m_1 = abs(comp[ :,-1:, 1: 2]) > TINY 
        comp[ :,-1:, :1] = (   comp[ :,-2:-1, :1] 
                             + comp[ :,-1:, 1: 2]) * 0.5 * m_0 * m_1
            
        m_0 = abs(comp[ :,-2:-1,-1:]) > TINY
        m_1 = abs(comp[ :,-1:,-2:-1]) > TINY
        comp[ :,-1:,-1:] = (   comp[ :,-2:-1,-1:] 
                             + comp[ :,-1:,-2:-1]) * 0.5 * m_0 * m_1

        # Edges in "y" direction 
        m_0 = abs(comp[ :1, :, 1: 2]) > TINY
        m_1 = abs(comp[ 1: 2, :, :1]) > TINY
        comp[ :1, :, :1] = (   comp[ :1, :, 1: 2] 
                             + comp[ 1: 2, :, :1]) * 0.5 * m_0 * m_1
        
        m_0 = abs(comp[-1:, :, 1: 2]) > TINY
        m_1 = abs(comp[-2:-1, :, :1]) > TINY
        comp[-1:, :, :1] = (   comp[-1:, :, 1: 2] 
                             + comp[-2:-1, :, :1]) * 0.5 * m_0 * m_1
        
        m_0 = abs(comp[ :1, :,-2:-1]) > TINY
        m_1 = abs(comp[ 1 :2, :,-1:]) > TINY
        comp[ :1, :,-1:] = (   comp[ :1, :,-2:-1] 
                             + comp[ 1 :2, :,-1:]) * 0.5 * m_0 * m_1
        
        m_0 = abs(comp[-1:, :,-2:-1]) > TINY
        m_1 = abs(comp[-2:-1, :,-1:]) > TINY
        comp[-1:, :,-1:] = (   comp[-1:, :,-2:-1] 
                             + comp[-2:-1, :,-1:]) * 0.5 * m_0 * m_1

        # Edges in "z" direction 
        m_0 = abs(comp[ 1: 2, :1, :]) > TINY 
        m_1 = abs(comp[ :1, 1: 2, :]) > TINY
        comp[ :1, :1, :] = (   comp[ 1: 2, :1, :] 
                             + comp[ :1, 1: 2, :]) * 0.5 * m_0 * m_1
        
        m_0 = abs(comp[ 1: 2,-1:, :]) > TINY 
        m_1 = abs(comp[ :1,-2:-1, :]) > TINY
        comp[ :1,-1:, :] = (   comp[ 1: 2,-1:, :] 
                             + comp[ :1,-2:-1, :]) * 0.5 * m_0 * m_1
        
        m_0 = abs(comp[-2:-1, :1, :]) > TINY 
        m_1 = abs(comp[-1:, 1 :2, :]) > TINY
        comp[-1:, :1, :] = (   comp[-2:-1, :1, :] 
                             + comp[-1:, 1 :2, :]) * 0.5 * m_0 * m_1
        
        m_0 = abs(comp[-2:-1,-1:, :]) > TINY
        m_1 = abs(comp[-1:,-2:-1, :]) > TINY
        comp[-1:,-1:, :] = (   comp[-2:-1,-1:, :] 
                             + comp[-1:,-2:-1, :]) * 0.5 * m_0 * m_1
        
    # -----------------------------------------------------
    # I am not sure what to do about the corners yet, but 
    # when I see coding for edges, I get a bit discouraged
    # -----------------------------------------------------

    # ------------------------
    # Filter the obstacle out
    # ------------------------
    if obst is not None:
        on = zeros( (nx, ny, nz) )
 
        # Nodify the obstacle
        on[0:-1, 0:-1, 0:-1] = obst[:, :, :]     
        on[1:  , 0:-1, 0:-1] = mx(on[1:   , 0:-1, 0:-1], on[0:-1, 0:-1, 0:-1])
        on[0:-1, 1:  , 0:-1] = mx(on[0:-1:, 1:  , 0:-1], on[0:-1, 0:-1, 0:-1])
        on[0:-1, 0:-1, 1:  ] = mx(on[0:-1:, 0:-1, 1:  ], on[0:-1, 0:-1, 0:-1])

        # Obstacle is non-zero in fluid, therefore "lnot"
        un[:,:,:] *= lnot(on[:,:,:])
        vn[:,:,:] *= lnot(on[:,:,:])
        wn[:,:,:] *= lnot(on[:,:,:])

    return un, vn, wn
