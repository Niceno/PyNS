"""
Computes advection fluxes for a given variable, with various advection schemes.

Following schemes are supported: "upwind", "minmod", "koren" and "superbee"
"""

# Standard Python modules
from pyns.standard import *

# PyNS modules
from pyns.constants.all import *
from pyns.operators.all import *

# =============================================================================
def advection(rho, phi, uvwf, dxyz, dt, lim_name):
# -----------------------------------------------------------------------------
  """
  Args:
    rho:      Three-dimensional matrix holding physical property in advection 
              term (density or density times capacity ...) for all cells.
    phi:      Unknown transported by advection (from "create_unknown").
    uvwf:     Tuple with three staggered velocity components (where each 
              component is created with "create_unknown" function.
    dxyz:     Tuple holding cell dimensions in "x", "y" and "z" directions.
              Each cell dimension is a three-dimensional matrix.
    dt:       Time step.
    lim_name: Limiter name.  Can be: "upwind", "minmod", "koren" and "superbee"
    
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
    rho_x_fac = cat_x((rho[:1,:,:], avg_x(rho), rho[-1:,:,:]))  # nxp,ny, nz
    rho_y_fac = cat_y((rho[:,:1,:], avg_y(rho), rho[:,-1:,:]))  # nx, nyp,nz 
    rho_z_fac = cat_z((rho[:,:,:1], avg_z(rho), rho[:,:,-1:]))  # nx, ny, nzp 

    # Facial values of areas including boundary cells
    a_x_fac = cat_x((sx[:1,:,:], avg_x(sx), sx[-1:,:,:]))
    a_y_fac = cat_y((sy[:,:1,:], avg_y(sy), sy[:,-1:,:]))
    a_z_fac = cat_z((sz[:,:,:1], avg_z(sz), sz[:,:,-1:]))
  
    del_x = avg_x(dx)
    del_y = avg_y(dy)
    del_z = avg_z(dz)

    # Facial values of velocities without boundary values
    u_fac = uf.val  # nxm,ny, nz
    v_fac = vf.val  # nx, nym,nz
    w_fac = wf.val  # nx, ny, nzm

    # Boundary velocity values
    u_bnd_W = uf.bnd[W].val
    u_bnd_E = uf.bnd[E].val
    v_bnd_S = vf.bnd[S].val  
    v_bnd_N = vf.bnd[N].val
    w_bnd_B = wf.bnd[B].val  
    w_bnd_T = wf.bnd[T].val

  # -----------------------------------------------------------
  # Specific for transported variable staggered in x direction
  # -----------------------------------------------------------
  if pos == X:  
    
    # Facial values of physical properties including boundary cells
    rho_x_fac = rho                             # nx, ny, nz
    rho_nod_y = avg_x(avg_y(rho))               # nxm,nym,nz
    rho_y_fac = cat_y((rho_nod_y[:, :1,:],      \
                       rho_nod_y[:,  :,:],      \
                       rho_nod_y[:,-1:,:]))     # nxm,nyp,nz                      
    rho_nod_z = avg_x(avg_z(rho))               # nxm,ny,nzm
    rho_z_fac = cat_z((rho_nod_z[:,:, :1],      \
                       rho_nod_z[:,:,  :],      \
                       rho_nod_z[:,:,-1:]))     # nxm,ny,nzp    

    # Facial values of areas including boundary cells
    a_x_fac = sx
    a_y_fac = cat_y((avg_x(sy[:,:1,:]),  \
                     avg_x(avg_y(sy)),   \
                     avg_x(sy[:,-1:,:])))
    a_z_fac = cat_z((avg_x(sz[:,:,:1]),  \
                     avg_x(avg_z(sz)),   \
                     avg_x(sz[:,:,-1:])))

    del_x = dx[1:-1,:,:]
    del_y = avg_x(avg_y(dy))
    del_z = avg_x(avg_z(dz))
  
    # Facial values of velocities without boundary values
    u_fac = avg_x(uf.val)  # nxmm,ny, nz
    v_fac = avg_x(vf.val)  # nxm, nym,nz
    w_fac = avg_x(wf.val)  # nxm, ny, nzm

    # Boundary velocity values
    u_bnd_W = uf.bnd[W].val         
    u_bnd_E = uf.bnd[E].val
    v_bnd_S = avg_x(vf.bnd[S].val)  
    v_bnd_N = avg_x(vf.bnd[N].val)
    w_bnd_B = avg_x(wf.bnd[B].val) 
    w_bnd_T = avg_x(wf.bnd[T].val)

  # -----------------------------------------------------------
  # Specific for transported variable staggered in y direction
  # -----------------------------------------------------------
  if pos == Y:  
    
    # Facial values of physical properties including boundary cells
    rho_nod_x = avg_y(avg_x(rho) )              # nxm,nym,nz
    rho_x_fac = cat_x((rho_nod_x[ :1,:,:],      \
                       rho_nod_x[  :,:,:],      \
                       rho_nod_x[-1:,:,:]))     # nxp,nym,nz
    rho_y_fac = rho                             # nx, ny, nz
    rho_nod_z = avg_y(avg_z(rho) )              # nx, nym,nzm
    rho_z_fac = cat_z((rho_nod_z[:,:, :1],      \
                       rho_nod_z[:,:,  :],      \
                       rho_nod_z[:,:,-1:]))     # nx, nym,nzp

    # Facial values of areas including boundary cells
    a_x_fac = cat_x((avg_y(sx[:1,:,:]),  \
                     avg_y(avg_x(sx)),   \
                     avg_y(sx[-1:,:,:])))
    a_y_fac = sy
    a_z_fac = cat_z((avg_y(sz[:,:,:1]),  \
                     avg_y(avg_z(sz)),   \
                     avg_y(sz[:,:,-1:])))

    del_x = avg_y(avg_x(dx))
    del_y = dy[:,1:-1,:]
    del_z = avg_y(avg_z(dz))
    
    # Facial values of velocities without boundary values
    u_fac = avg_y(uf.val)  # nxm,nym, nz
    v_fac = avg_y(vf.val)  # nx, nymm,nz
    w_fac = avg_y(wf.val)  # nx, nym, nzm

    # Facial values of velocities with boundary values
    u_bnd_W = avg_y(uf.bnd[W].val)  
    u_bnd_E = avg_y(uf.bnd[E].val)
    v_bnd_S = vf.bnd[S].val          
    v_bnd_N = vf.bnd[N].val
    w_bnd_B = avg_y(wf.bnd[B].val)  
    w_bnd_T = avg_y(wf.bnd[T].val)

  # -----------------------------------------------------------
  # Specific for transported variable staggered in z direction
  # -----------------------------------------------------------
  if pos == Z:
    
    # Facial values of physical properties including boundary cells
    rho_nod_x = avg_z(avg_x(rho) )              # nxm,ny, nzm
    rho_x_fac = cat_x((rho_nod_x[ :1,:,:],      \
                       rho_nod_x[  :,:,:],      \
                       rho_nod_x[-1:,:,:]))     # nxp,ny, nzm
    rho_nod_y = avg_z(avg_y(rho) )              # nx, nym,nzm
    rho_y_fac = cat_y((rho_nod_y[:, :1,:],      \
                       rho_nod_y[:,  :,:],      \
                       rho_nod_y[:,-1:,:]))     # nx, nyp,nzm
    rho_z_fac = rho                             # nx, ny, nz

    # Facial values of areas including boundary cells
    a_x_fac = cat_x((             \
              avg_z(sx[:1,:,:]),  \
              avg_z(avg_x(sx)),   \
              avg_z(sx[-1:,:,:])))
    a_y_fac = cat_y((             \
              avg_z(sy[:,:1,:]),  \
              avg_z(avg_y(sy)),   \
              avg_z(sy[:,-1:,:])))
    a_z_fac = sz

    del_x = avg_z(avg_x(dx))
    del_y = avg_z(avg_y(dy))
    del_z = dz[:,:,1:-1]
    
    # Facial values of velocities without boundary values
    u_fac = avg_z(uf.val)  # nxm,ny,  nzm
    v_fac = avg_z(vf.val)  # nx, nym, nzm
    w_fac = avg_z(wf.val)  # nx, ny,  nzmm

    # Facial values of velocities with boundary values
    u_bnd_W = avg_z(uf.bnd[W].val)  
    u_bnd_E = avg_z(uf.bnd[E].val)
    v_bnd_S = avg_z(vf.bnd[S].val)  
    v_bnd_N = avg_z(vf.bnd[N].val)
    w_bnd_B = wf.bnd[B].val          
    w_bnd_T = wf.bnd[T].val

  # -----------------------------
  # Common part of the algorithm
  # -----------------------------
    
  # ------------------------------------------------------------
  #
  #    |-o-|-o-|-o-|-o-|-o-|-o-|-o-|-o-|-o-|-o-|
  #      1   2   3   4   5   6   7   8   9   10     phi
  #        x---x---x---x---x---x---x---x---x      
  #        1   2   3   4   5   6   7   8   9        d_x initial
  #    0---x---x---x---x---x---x---x---x---x---0      
  #    1   2   3   4   5   6   7   8   9  10  11    d_x padded
  #
  # ------------------------------------------------------------
 
  # Compute consecutive differences (and avoid division by zero)
  d_x = dif_x(phi.val)  # nxm, ny, nz  
  d_x[(d_x >  -TINY) & (d_x <=   0.0)] = -TINY 
  d_x[(d_x >=   0.0) & (d_x <  +TINY)] = +TINY 
  d_x = cat_x((d_x[:1,:,:], d_x, d_x[-1:,:,:]))  # nxp, ny, nz
    
  d_y = dif_y(phi.val)  # nx, nym, nz  
  d_y[(d_y >  -TINY) & (d_y <=   0.0)] = -TINY 
  d_y[(d_y >=   0.0) & (d_y <  +TINY)] = +TINY 
  d_y = cat_y((d_y[:,:1,:], d_y, d_y[:,-1:,:]))  # nx, nyp, nz
   
  d_z = dif_z(phi.val)  # nx, ny, nzm  
  d_z[(d_z >  -TINY) & (d_z <=   0.0)] = -TINY 
  d_z[(d_z >=   0.0) & (d_z <  +TINY)] = +TINY 
  d_z = cat_z((d_z[:,:,:1], d_z, d_z[:,:,-1:]))  # nx, ny, nzp
      
  # Ratio of consecutive gradients for positive and negative flow
  r_x_we = d_x[1:-1,:,:] / d_x[0:-2,:,:]  # nxm,ny, nz
  r_x_ew = d_x[2:,  :,:] / d_x[1:-1,:,:]  # nxm,ny, nz 
  r_y_sn = d_y[:,1:-1,:] / d_y[:,0:-2,:]  # nx, nym,nz
  r_y_ns = d_y[:,2:,  :] / d_y[:,1:-1,:]  # nx, nym,nz
  r_z_bt = d_z[:,:,1:-1] / d_z[:,:,0:-2]  # nx, ny, nzm
  r_z_tb = d_z[:,:,2:  ] / d_z[:,:,1:-1]  # nx, ny, nzm
    
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
    psi_x = mx(zeros(r_x.shape), mn(r_x, ones(r_x.shape)))
    psi_y = mx(zeros(r_y.shape), mn(r_y, ones(r_y.shape)))
    psi_z = mx(zeros(r_z.shape), mn(r_z, ones(r_z.shape)))
  elif lim_name == 'superbee':
    psi_x = mx(zeros(r_x.shape), mn(2.*r_x, ones(r_x.shape)), mn(r_x, 2.))
    psi_y = mx(zeros(r_y.shape), mn(2.*r_y, ones(r_y.shape)), mn(r_y, 2.))
    psi_z = mx(zeros(r_z.shape), mn(2.*r_z, ones(r_z.shape)), mn(r_z, 2.))
  elif lim_name == 'koren':
    psi_x = mx(zeros(r_x.shape), mn(2.*r_x, (2.+r_x)/3., 2.*ones(r_x.shape)))
    psi_y = mx(zeros(r_y.shape), mn(2.*r_y, (2.+r_y)/3., 2.*ones(r_y.shape)))
    psi_z = mx(zeros(r_z.shape), mn(2.*r_z, (2.+r_z)/3., 2.*ones(r_z.shape)))

  flux_fac_lim_x =   phi.val[0:-1,:,:] * u_fac * flow_we               \
                 +   phi.val[1:,  :,:] * u_fac * flow_ew               \
                 +   0.5 * abs(u_fac) * (1 - abs(u_fac) * dt / del_x)  \
                 *  (   psi_x[:,:,:] * d_x[0:nx-1,:,:] * flow_we       \
                      + psi_x[:,:,:] * d_x[1:nx,  :,:] * flow_ew )   
  flux_fac_lim_y =   phi.val[:,0:-1,:] * v_fac * flow_sn               \
                 +   phi.val[:,1:  ,:] * v_fac * flow_ns               \
                 +   0.5 * abs(v_fac) * (1 - abs(v_fac) * dt / del_y)  \
                 *  (   psi_y[:,:,:] * d_y[:,0:ny-1,:] * flow_sn       \
                      + psi_y[:,:,:] * d_y[:,1:ny,  :] * flow_ns ) 
  flux_fac_lim_z =   phi.val[:,:,0:-1] * w_fac * flow_bt               \
                 +   phi.val[:,:,1:  ] * w_fac * flow_tb               \
                 +   0.5 * abs(w_fac) * (1 - abs(w_fac) * dt / del_z)  \
                 *  (   psi_z[:,:,:] * d_z[:,:,0:nz-1] * flow_bt       \
                      + psi_z[:,:,:] * d_z[:,:,1:nz  ] * flow_tb ) 

  # Pad with boundary values
  flux_fac_lim_x = cat_x((phi.bnd[W].val * u_bnd_W,      \
                          flux_fac_lim_x,                \
                          phi.bnd[E].val * u_bnd_E))
  flux_fac_lim_y = cat_y((phi.bnd[S].val * v_bnd_S,      \
                          flux_fac_lim_y,                \
                          phi.bnd[N].val * v_bnd_N))
  flux_fac_lim_z = cat_z((phi.bnd[B].val * w_bnd_B,      \
                          flux_fac_lim_z,                \
                          phi.bnd[T].val * w_bnd_T))
                       
  # Multiply with face areas                     
  flux_fac_lim_x = rho_x_fac * flux_fac_lim_x * a_x_fac
  flux_fac_lim_y = rho_y_fac * flux_fac_lim_y * a_y_fac
  flux_fac_lim_z = rho_z_fac * flux_fac_lim_z * a_z_fac
      
  # Sum contributions from all directions up
  c = dif_x(flux_fac_lim_x) + \
      dif_y(flux_fac_lim_y) + \
      dif_z(flux_fac_lim_z)
    
  return c  # end of function
