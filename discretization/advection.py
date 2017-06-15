"""
Computes advection fluxes for a given variable, with various advection schemes.

Following schemes are supported: "upwind", "minmod", "koren" and "superbee"
"""

# Standard Python modules
from pyns.standard import *

# PyNS modules
from pyns.constants import *
from pyns.operators import *

# =============================================================================
def advection(rho, phi, uvwf, dxyz, dt, lim_name):
# -----------------------------------------------------------------------------
    """
    Args:
      rho:      Three-dimensional matrix holding physical property in advection
                term (density or density times capacity ...) for cells.
      phi:      Unknown transported by advection (from "create_unknown").
      uvwf:     Tuple with three staggered velocity components (where each
                component is created with "create_unknown" function.
      dxyz:     Tuple holding cell dimensions in "x", "y" and "z" directions.
                Each cell dimension is a three-dimensional matrix.
      dt:       Time step.
      lim_name: Limiter name.
                Can be: "upwind", "minmod", "koren" and "superbee"

    Returns:
      Three-dimensional matrix with advection term.
    """

    res = phi.val.shape
    nx, ny, nz = res

    # Unpack tuples
    uf, vf, wf = uvwf
    dx, dy, dz = dxyz

    pos = phi.pos

    # Pre-compute geometrical quantities
    sx = dy * dz
    sy = dx * dz
    sz = dx * dy

    # ------------------------------------------------
    # Specific for cell-centered transported variable
    # ------------------------------------------------
    if pos == C:

        # Facial values of physical properties including boundary cells
        rho_x_fac = cat_x((rho[:1,:,:], avg_x(rho), rho[-1:,:,:]))  
        rho_y_fac = cat_y((rho[:,:1,:], avg_y(rho), rho[:,-1:,:]))  
        rho_z_fac = cat_z((rho[:,:,:1], avg_z(rho), rho[:,:,-1:]))  

        # Facial values of areas including boundary cells
        a_x_fac = cat_x((sx[:1,:,:], avg_x(sx), sx[-1:,:,:]))
        a_y_fac = cat_y((sy[:,:1,:], avg_y(sy), sy[:,-1:,:]))
        a_z_fac = cat_z((sz[:,:,:1], avg_z(sz), sz[:,:,-1:]))

        # Distance between cell centers, defined at faces including boundaries
        # MODIFY FOR PERIODIC!!!
        del_x = cat_x((dx[:1,:,:]*0.5, avg_x(dx), dx[-1:,:,:]*0.5))
        del_y = cat_y((dy[:,:1,:]*0.5, avg_y(dy), dy[:,-1:,:]*0.5))
        del_z = cat_z((dz[:,:,:1]*0.5, avg_z(dz), dz[:,:,-1:]*0.5))

        # Velocities defined at faces including boundaries
        # TAKE CARE YOU HAVE FRESH VALUES FOR PERIODIC
        u_fac = cat_x((uf.bnd[W].val, uf.val, uf.bnd[E].val))
        v_fac = cat_y((vf.bnd[S].val, vf.val, vf.bnd[N].val))
        w_fac = cat_z((wf.bnd[B].val, wf.val, wf.bnd[T].val))
        
    # -----------------------------------------------------------
    # Specific for transported variable staggered in x direction
    # -----------------------------------------------------------
    if pos == X:

        # Facial values of physical properties including boundary cells
        rho_x_fac = rho                             
        rho_nod_y = avg_x(avg_y(rho))               
        rho_y_fac = cat_y((rho_nod_y[:, :1,:], rho_nod_y, rho_nod_y[:,-1:,:]))     
        rho_nod_z = avg_x(avg_z(rho))               
        rho_z_fac = cat_z((rho_nod_z[:,:, :1], rho_nod_z, rho_nod_z[:,:,-1:]))     

        # Facial values of areas including boundary cells
        a_x_fac = sx
        a_y_fac = cat_y((avg_x(sy[:,:1,:]),
                         avg_x(avg_y(sy)),
                         avg_x(sy[:,-1:,:])))
        a_z_fac = cat_z((avg_x(sz[:,:,:1]),
                         avg_x(avg_z(sz)),
                         avg_x(sz[:,:,-1:])))

        # Distance between cell centers, defined at faces including boundaries
        # MODIFY FOR PERIODIC
        del_x = dx[:,:,:]
        del_y = avg_x(cat_y((dy[:,:1,:]*0.5, avg_y(dy), dy[:,-1:,:]*0.5)))
        del_z = avg_x(cat_z((dz[:,:,:1]*0.5, avg_z(dz), dz[:,:,-1:]*0.5)))

        # Velocities defined at faces including boundaries
        # TAKE CARE YOU HAVE FRESH VALUES FOR PERIODIC
        u_fac = cat_x((uf.bnd[W].val, avg_x(uf.val), uf.bnd[E].val))  
        v_fac = avg_x(cat_y((vf.bnd[S].val, vf.val, vf.bnd[N].val)))
        w_fac = avg_x(cat_z((wf.bnd[B].val, wf.val, wf.bnd[T].val)))  

    # -----------------------------------------------------------
    # Specific for transported variable staggered in y direction
    # -----------------------------------------------------------
    if pos == Y:

        # Facial values of physical properties including boundary cells
        rho_nod_x = avg_y(avg_x(rho) )
        rho_x_fac = cat_x((rho_nod_x[ :1,:,:], rho_nod_x, rho_nod_x[-1:,:,:]))
        rho_y_fac = rho
        rho_nod_z = avg_y(avg_z(rho) )
        rho_z_fac = cat_z((rho_nod_z[:,:, :1], rho_nod_z, rho_nod_z[:,:,-1:]))

        # Facial values of areas including boundary cells
        a_x_fac = cat_x((avg_y(sx[:1,:,:]),
                         avg_y(avg_x(sx)),
                         avg_y(sx[-1:,:,:])))
        a_y_fac = sy
        a_z_fac = cat_z((avg_y(sz[:,:,:1]),
                         avg_y(avg_z(sz)),
                         avg_y(sz[:,:,-1:])))

        # Distance between cell centers, defined at faces including boundaries
        # MODIFY FOR PERIODIC
        del_x = avg_y(cat_x((dx[:1,:,:]*0.5, avg_x(dx), dx[-1:,:,:]*0.5)))
        del_y = dy[:,:,:]
        del_z = avg_y(cat_z((dz[:,:,:1]*0.5, avg_z(dz), dz[:,:,-1:]*0.5)))

        # Velocities defined at faces including boundaries
        # TAKE CARE YOU HAVE FRESH VALUES FOR PERIODIC
        u_fac = avg_y(cat_x((uf.bnd[W].val, uf.val, uf.bnd[E].val)))
        v_fac = cat_y((vf.bnd[S].val, avg_y(vf.val), vf.bnd[N].val)) 
        w_fac = avg_y(cat_z((wf.bnd[B].val, wf.val, wf.bnd[T].val)))
        
    # -----------------------------------------------------------
    # Specific for transported variable staggered in z direction
    # -----------------------------------------------------------
    if pos == Z:

        # Facial values of physical properties including boundary cells
        rho_nod_x = avg_z(avg_x(rho) )
        rho_x_fac = cat_x((rho_nod_x[ :1,:,:], rho_nod_x, rho_nod_x[-1:,:,:]))
        rho_nod_y = avg_z(avg_y(rho) )
        rho_y_fac = cat_y((rho_nod_y[:, :1,:], rho_nod_y, rho_nod_y[:,-1:,:]))
        rho_z_fac = rho                             # nx, ny, nz

        # Facial values of areas including boundary cells
        a_x_fac = cat_x((avg_z(sx[:1,:,:]),
                         avg_z(avg_x(sx)),
                         avg_z(sx[-1:,:,:])))
        a_y_fac = cat_y((avg_z(sy[:,:1,:]),
                         avg_z(avg_y(sy)),
                         avg_z(sy[:,-1:,:])))
        a_z_fac = sz

        # Facial values of distance between cell centers
        # MODIFY FOR PERIODIC
        del_x = avg_z(cat_x((dx[:1,:,:]*0.5, avg_x(dx), dx[-1:,:,:]*0.5)))
        del_y = avg_z(cat_y((dy[:,:1,:]*0.5, avg_y(dy), dy[:,-1:,:]*0.5)))
        del_z = dz[:,:,:]

        # Facial values of velocities without boundary values
        # TAKE CARE YOU HAVE FRESH VALUES FOR PERIODIC
        u_fac = avg_z(cat_x((uf.bnd[W].val, uf.val, uf.bnd[E].val)))  
        v_fac = avg_z(cat_y((vf.bnd[S].val, vf.val, vf.bnd[N].val)))
        w_fac = cat_z((wf.bnd[B].val, avg_z(wf.val), wf.bnd[T].val))  

    # -----------------------------
    # Common part of the algorithm
    # -----------------------------

    # ------------------------------------------------------------------
    #
    #    |-W-|-W-|-o-|-o-|-o-|-o-|-o-|-o-|-o-|-o-|-o-|-o-|-E-|-E-|
    #      0   1   2   3   4   5   6   7   8   9  10  11  12  13    phi
    #        x---x---x---x---x---x---x---x---x---x---x---x---x
    #        0   1   2   3   4   5   6   7   8   9  10  11  12      d_x
    #
    #            x---x---x---x---x---x---x---x---x---x---x
    #            0   1   2   3   4   5   6   7   8   9  10          r_x
    #
    # ------------------------------------------------------------------

    # Compute consecutive differences (and avoid division by zero)
    # FILL THE BUFFERS PROPERLY FOR PERIODIC
    d_x = dif_x(cat_x((phi.bnd[W].val, 
                       phi.bnd[W].val, 
                       phi.val, 
                       phi.bnd[E].val,
                       phi.bnd[E].val)))  
    d_x[(d_x >  -TINY) & (d_x <=   0.0)] = -TINY
    d_x[(d_x >=   0.0) & (d_x <  +TINY)] = +TINY

    d_y = dif_y(phi.val)  
    d_y = dif_y(cat_y((phi.bnd[S].val, 
                       phi.bnd[S].val, 
                       phi.val, 
                       phi.bnd[N].val,
                       phi.bnd[N].val)))  
    d_y[(d_y >  -TINY) & (d_y <=   0.0)] = -TINY
    d_y[(d_y >=   0.0) & (d_y <  +TINY)] = +TINY

    d_z = dif_z(phi.val)  
    d_z = dif_z(cat_z((phi.bnd[B].val, 
                       phi.bnd[B].val, 
                       phi.val, 
                       phi.bnd[T].val,
                       phi.bnd[T].val)))  
    d_z[(d_z >  -TINY) & (d_z <=   0.0)] = -TINY
    d_z[(d_z >=   0.0) & (d_z <  +TINY)] = +TINY

    # Ratio of consecutive gradients for positive and negative flow
    r_x_we = d_x[1:-1,:,:] / d_x[0:-2,:,:]  
    r_x_ew = d_x[2:,  :,:] / d_x[1:-1,:,:]  
    r_y_sn = d_y[:,1:-1,:] / d_y[:,0:-2,:]  
    r_y_ns = d_y[:,2:,  :] / d_y[:,1:-1,:]  
    r_z_bt = d_z[:,:,1:-1] / d_z[:,:,0:-2]  
    r_z_tb = d_z[:,:,2:  ] / d_z[:,:,1:-1] 

    flow_we = u_fac >= 0
    flow_ew = lnot(flow_we)
    flow_sn = v_fac >= 0
    flow_ns = lnot(flow_sn)
    flow_bt = w_fac >= 0
    flow_tb = lnot(flow_bt)

    r_x = r_x_we * flow_we + r_x_ew * flow_ew
    r_y = r_y_sn * flow_sn + r_y_ns * flow_ns
    r_z = r_z_bt * flow_bt + r_z_tb * flow_tb

    # Apply a limiter
    if lim_name == 'upwind':
        psi_x = r_x * 0.0
        psi_y = r_y * 0.0
        psi_z = r_z * 0.0
    elif lim_name == 'minmod':
        psi_x = mx(zeros(r_x.shape),mn(r_x,ones(r_x.shape)))
        psi_y = mx(zeros(r_y.shape),mn(r_y,ones(r_y.shape)))
        psi_z = mx(zeros(r_z.shape),mn(r_z,ones(r_z.shape)))
    elif lim_name == 'superbee':
        psi_x = mx(zeros(r_x.shape),mn(2.*r_x, ones(r_x.shape)),mn(r_x, 2.))
        psi_y = mx(zeros(r_y.shape),mn(2.*r_y, ones(r_y.shape)),mn(r_y, 2.))
        psi_z = mx(zeros(r_z.shape),mn(2.*r_z, ones(r_z.shape)),mn(r_z, 2.))
    elif lim_name == 'koren':
        psi_x = mx(zeros(r_x.shape),mn(2.*r_x,(2.+r_x)/3.,2.*ones(r_x.shape)))
        psi_y = mx(zeros(r_y.shape),mn(2.*r_y,(2.+r_y)/3.,2.*ones(r_y.shape)))
        psi_z = mx(zeros(r_z.shape),mn(2.*r_z,(2.+r_z)/3.,2.*ones(r_z.shape)))

    flux_fac_lim_x =   cat_x((phi.bnd[W].val, phi.val)) * u_fac * flow_we  \
                   +   cat_x((phi.val, phi.bnd[E].val)) * u_fac * flow_ew  \
                   +   0.5 * abs(u_fac) * (1 - abs(u_fac) * dt / del_x)    \
                   *  (   psi_x[:,:,:] * d_x[0:-2,:,:] * flow_we           \
                        + psi_x[:,:,:] * d_x[1:-1,:,:] * flow_ew )
    flux_fac_lim_y =   cat_y((phi.bnd[S].val, phi.val)) * v_fac * flow_sn  \
                   +   cat_y((phi.val, phi.bnd[N].val)) * v_fac * flow_ns  \
                   +   0.5 * abs(v_fac) * (1 - abs(v_fac) * dt / del_y)    \
                   *  (   psi_y[:,:,:] * d_y[:,0:-2,:] * flow_sn           \
                        + psi_y[:,:,:] * d_y[:,1:-1,:] * flow_ns )
    flux_fac_lim_z =   cat_z((phi.bnd[B].val, phi.val)) * w_fac * flow_bt  \
                   +   cat_z((phi.val, phi.bnd[T].val)) * w_fac * flow_tb  \
                   +   0.5 * abs(w_fac) * (1 - abs(w_fac) * dt / del_z)    \
                   *  (   psi_z[:,:,:] * d_z[:,:,0:-2] * flow_bt           \
                        + psi_z[:,:,:] * d_z[:,:,1:-1] * flow_tb )

    # Multiply with face areas
    flux_fac_lim_x = rho_x_fac * flux_fac_lim_x * a_x_fac
    flux_fac_lim_y = rho_y_fac * flux_fac_lim_y * a_y_fac
    flux_fac_lim_z = rho_z_fac * flux_fac_lim_z * a_z_fac

    # Sum contributions from all directions up
    c = dif_x(flux_fac_lim_x) + \
        dif_y(flux_fac_lim_y) + \
        dif_z(flux_fac_lim_z)

    return c  # end of function
