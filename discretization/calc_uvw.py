"""
Discretizes and solves momentum equation (all three components).
"""

# Standard Python modules
from pyns.standard import *

# PyNS modules
from pyns.constants import *
from pyns.operators import *

from pyns.discretization.adj_n_bnds     import adj_n_bnds
from pyns.discretization.adj_o_bnds     import adj_o_bnds
from pyns.discretization.advection      import advection
from pyns.discretization.create_matrix  import create_matrix
from pyns.discretization.obst_zero_val  import obst_zero_val
from pyns.solvers                       import cg, cgs, bicgstab

# =============================================================================
def calc_uvw(uvw, uvwf, rho, mu, dt, dxyz, obst,
             pressure = None,
             force    = None):
# -----------------------------------------------------------------------------
    """
    Args:
      uvw:      Tuple with three velocity components, staggered or collocated.
                (Each component is created with "create_unknown" function.)
      uvwf:     Tuple with three staggered velocity components (where each
                component is created with "pyns.create_unknown" function.
      rho:      Three-dimensional array holding density for all cells.
      mu:       Three-dimensional array holding dynamic viscosity.
      pressure: Object "Unknown" holding total pressure.
      force:    Tuple containing three-dimensional matrices holding external
                forces in each direction.
      dt:       Time step.
      dxyz:     Tuple holding cell dimensions in "x", "y" and "z" directions.
                Each cell dimension is a three-dimensional array.
      obst:     Obstacle, three-dimensional array with zeros and ones.
                It is zero in fluid, one in solid.

    Returns:
      none, but input argument uvw is modified!
    """

    # Unpack tuples
    u,   v,   w   = uvw
    uf,  vf,  wf  = uvwf
    dx,  dy,  dz  = dxyz

    # Fetch resolutions
    ru = u.val.shape
    rv = v.val.shape
    rw = w.val.shape

    # Pre-compute geometrical quantities
    dv = dx * dy * dz

    d = u.pos

    # Create system matrices and right hand sides
    A_u = create_matrix(u, rho/dt, mu, dxyz, obst, DIRICHLET)
    A_v = create_matrix(v, rho/dt, mu, dxyz, obst, DIRICHLET)
    A_w = create_matrix(w, rho/dt, mu, dxyz, obst, DIRICHLET)
    b_u = zeros(ru)
    b_v = zeros(rv)
    b_w = zeros(rw)

    # Advection terms for momentum
    c_u = advection(rho, u, uvwf, dxyz, dt, 'superbee');
    c_v = advection(rho, v, uvwf, dxyz, dt, 'superbee');
    c_w = advection(rho, w, uvwf, dxyz, dt, 'superbee');

    # Innertial term for momentum (this works for collocated and staggered)
    i_u = u.old * avg(u.pos, rho) * avg(u.pos, dv) / dt
    i_v = v.old * avg(v.pos, rho) * avg(v.pos, dv) / dt
    i_w = w.old * avg(w.pos, rho) * avg(w.pos, dv) / dt

    # Full force terms for momentum equations (collocated and staggered)
    f_u = b_u - c_u + i_u 
    f_v = b_v - c_v + i_v 
    f_w = b_w - c_w + i_w 

    # Compute staggered pressure gradients
    if pressure is not None:
        pressure.exchange()
        p_tot_x = dif_x(pressure.val) / avg_x(dx)
        p_tot_y = dif_y(pressure.val) / avg_y(dy)
        p_tot_z = dif_z(pressure.val) / avg_z(dz)

        # Make pressure gradients cell-centered
        if d == C:
            p_tot_x = avg_x(cat_x((p_tot_x[:1,:,:], p_tot_x, p_tot_x[-1:,:,:])))
            p_tot_y = avg_y(cat_y((p_tot_y[:,:1,:], p_tot_y, p_tot_y[:,-1:,:])))
            p_tot_z = avg_z(cat_z((p_tot_z[:,:,:1], p_tot_z, p_tot_z[:,:,-1:])))
    
        # Total pressure gradients (this works for collocated and staggered)
        p_st_u = p_tot_x * avg(u.pos, dv)
        p_st_v = p_tot_y * avg(v.pos, dv)
        p_st_w = p_tot_z * avg(w.pos, dv)
    
        f_u -= p_st_u 
        f_v -= p_st_v 
        f_w -= p_st_w 

    # If external force is given
    if force is not None:
        
        # Unpack the tuple 
        e_u, e_v, e_w = force
        
        # Add to right hand side
        f_u += e_u * avg(u.pos, dv)
        f_v += e_v * avg(v.pos, dv)
        f_w += e_w * avg(w.pos, dv)

    # Take care of obsts in the domian
    if obst is not None:
        f_u = obst_zero_val(u.pos, f_u, obst)
        f_v = obst_zero_val(v.pos, f_v, obst)
        f_w = obst_zero_val(w.pos, f_w, obst)

    # Solve for velocities
    u.val[:] = cgs(A_u, u, f_u, TOL, False)
    v.val[:] = cgs(A_v, v, f_v, TOL, False)
    w.val[:] = cgs(A_w, w, f_w, TOL, False)

    # Update velocities in boundary cells
    adj_o_bnds((u,v,w), (dx,dy,dz), dt)
    adj_n_bnds(u)
    adj_n_bnds(v)
    adj_n_bnds(w)

    # Update face velocities
    # (For collocated arrangement also substract cell-centered
    # pressure gradients and add staggered pressure gradients)
    if d == C:
        uf.val[:]  = avg_x(u.val)
        vf.val[:]  = avg_y(v.val)
        wf.val[:]  = avg_z(w.val)

        if pressure is not None:        
            uf.val[:] += avg_x(dt / rho  * p_tot_x)
            uf.val[:] -= dt / avg_x(rho) * (dif_x(pressure.val) / avg_x(dx))
    
            vf.val[:] += avg_y(dt / rho  * p_tot_y)
            vf.val[:] -= dt / avg_y(rho) * (dif_y(pressure.val) / avg_y(dy))
    
            wf.val[:] += avg_z(dt / rho  * p_tot_z)
            wf.val[:] -= dt / avg_z(rho) * (dif_z(pressure.val) / avg_z(dz))

        for j in (W,E):
            uf.bnd[j].val[:] = u.bnd[j].val[:]
            vf.bnd[j].val[:] = avg_y(v.bnd[j].val[:])
            wf.bnd[j].val[:] = avg_z(w.bnd[j].val[:])
        for j in (S,N):
            uf.bnd[j].val[:] = avg_x(u.bnd[j].val[:])
            vf.bnd[j].val[:] = v.bnd[j].val[:]
            wf.bnd[j].val[:] = avg_z(w.bnd[j].val[:])
        for j in (B,T):
            uf.bnd[j].val[:] = avg_x(u.bnd[j].val[:])
            vf.bnd[j].val[:] = avg_y(v.bnd[j].val[:])
            wf.bnd[j].val[:] = w.bnd[j].val[:]

    else:
        uf.val[:] = u.val[:]
        vf.val[:] = v.val[:]
        wf.val[:] = w.val[:]
        for j in (W,E,S,N,B,T):
            uf.bnd[j].val[:] = u.bnd[j].val[:]
            vf.bnd[j].val[:] = v.bnd[j].val[:]
            wf.bnd[j].val[:] = w.bnd[j].val[:]

    if obst is not None:
        uf.val[:] = obst_zero_val(X, uf.val, obst)
        vf.val[:] = obst_zero_val(Y, vf.val, obst)
        wf.val[:] = obst_zero_val(Z, wf.val, obst)

    

    return  # end of function
